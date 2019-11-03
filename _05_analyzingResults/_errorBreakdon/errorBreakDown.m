%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the
aliasing phase

INPUT VARS
 r0,L0          :: Respectively the Fried's parameter and outer in meter
  d                :: The sub-aperture size
nActu          :: the \# DM actuators
nT               :: number of multiple of the Dm cut-off frequency for
defining the maximal reconstructed frequency
 
OUTPUT VARS
 sfFit_2D             :: bi-dimensional phase structure function map
(Toeplitz)
 psd               :: bi-dimensional fitting PSD
sfFit_4D              : bi-dimensional phase structure matrix
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function out = errorBreakDown(bestFit,varargin)
inputs = inputParser;
inputs.addRequired('bestFit',@(x) isa(x,'bestFitting'));
inputs.addParameter('display', true,@islogical);
inputs.parse(bestFit);
display= inputs.Results.display;

wvl = bestFit.psfr_.trs.cam.wavelength;
otf   = bestFit.psfr_.otf;
sf     = bestFit.psfr_.sf;
S      = sum(otf.otfDL(:));
r053 = bestFit.r0_fit^(-5/3);
r053Init = bestFit.psfr_.res.seeing.r0^(-5/3);

%1\ Get the static error from the static map (NCPA)
srStatic = sum(real(otf.otfStat(:)))/S;
wfeStatic = tools.sr2wfe(srStatic,wvl);

%2\ Get the fitting error by integrating the Fitting PSD
srFit = sum(sum(otf.otfDL.*exp(-0.5*sf.Dfit*r053/r053Init)))/S;
wfeFit = tools.sr2wfe(srFit,wvl);

%3\ Get the aliasing error from the covariance matrix
if ~isempty(bestFit.idxDal)
    gAl  = bestFit.xao_fit(bestFit.idxDal - length(bestFit.idxR0));
else
    gAl = 1;
end

srAl = sum(sum(otf.otfDL.*exp(-0.5*sf.Dal*r053/r053Init*gAl)))/S;
wfeAl = tools.sr2wfe(srAl,wvl);

%4\ Get the noise errors
%4.1 HO WFS Noise
if ~isempty(bestFit.idxR0) && ~isempty(bestFit.idxDho)
    gHO = bestFit.xao_fit(bestFit.idxDho- length(bestFit.idxR0));
else
    gHO = 1;
end

wfeNoise = 0;
cte = 1e9*sqrt(bestFit.psfr_.trs.holoop.tf.pn/size(bestFit.psfr_.res.noise.Cn_ho,1));
cn_tmp = bestFit.psfr_.res.noise.Cn_ho;

for i=1:bestFit.nGainsHO-1
    tmp = bestFit.psfr_.trs.mat.Hz{i}*bestFit.psfr_.res.noise.Cn_ho*bestFit.psfr_.trs.mat.Hz{i}';
    wfeNoise  = wfeNoise + gHO(i)*trace(tmp);
    cn_tmp = cn_tmp - tmp;
end
wfeNoise = sqrt(wfeNoise + gHO(end)*trace(cn_tmp))*cte;

%4.2 TT WFS Noise
if ~isempty(bestFit.idxR0) && ~isempty(bestFit.idxDtt)
    gTT = bestFit.xao_fit(bestFit.idxDtt - length(bestFit.idxR0));
else
    gTT = 1;
end
wfeNoiseTT = sqrt(trace(gTT*bestFit.psfr_.trs.ttloop.tf.pn*bestFit.psfr_.res.noise.Cn_tt))*1e9;

%5\. AO Bandwidth errors
srLag    = sum(sum(otf.otfDL.*exp(-0.5*sum(bsxfun(@times,sf.Dho_z,reshape(gHO,1,1,[])), 3))))/S;
wfeLag = sqrt(tools.sr2wfe(srLag,wvl)^2 - wfeNoise^2) ;

%6.  Residual tip-tilt
srTT    = sum(sum(otf.otfDL.*exp(-0.5*sf.Dtt*gTT)))/S;
wfeTT = sqrt(tools.sr2wfe(srTT,wvl)^2 - wfeNoiseTT^2);

%7\  Anisoplanatism
srAni = real(sum(sum(otf.Kani.*otf.otfDL))/S);
wfeAni= tools.sr2wfe(srAni, wvl);

%8\  Total residual
wfeTot = sqrt(sum(wfeStatic^2 + wfeFit^2 + wfeLag^2 ...
    + wfeNoise^2 + wfeAni^2 + wfeTT^2 +  wfeNoiseTT^2 + wfeAl^2));

srMar = tools.wfe2sr(wfeTot,wvl);
SRho  = srStatic*srFit.*srAni*srLag*srAl;
srPar = 1e2*SRho/(1+tools.sr2wfe(srTT,wvl)^2*(2*pi*1e-9/wvl)^2) + (1-SRho)/(1+(bestFit.psfr_.trs.tel.Dcircle/ bestFit.r0_fit)^2);

% -----------------  DISPLAY
if display
    fprintf('-------------------------------\n');
    fprintf('Wavelength\t\t%.3g micron\n', wvl*1e6)
    fprintf('Strehl Image\t\t%.4g%s\t\n',1e2*bestFit.psfr_.trs.cam.SR,'%');
    fprintf('Strehl Marechal \t%.4g%s\t\n',srMar,'%');
    fprintf('Strehl Parenti\t\t%.4g%s\t\n',srPar,'%');
    fprintf('Wavefront error (nm)\t%.4g\t\n',wfeTot);
    fprintf('-------------------------------\n');
    fprintf('NCPA calibration\t%.4g\n',wfeStatic);
    fprintf('Atmospheric Fitting\t%.4g\n',wfeFit);
    fprintf('-------------------------------\n');
    fprintf('Servo-lag\t\t%.4g\n',wfeLag);
    fprintf('WFS noise\t\t%.4g\n',wfeNoise);
    fprintf('WFS aliasing\t\t%.4g\n',wfeAl);
    fprintf('Anisoplanatism\t%.4g\n',max(wfeAni));
    %fprintf('Focal anisoplanatism.\t%.4g\n',wfeAnisoLGS);
    fprintf('-------------------------------\n');
    fprintf('Tip-tilt bandwidth\t%.4g\n',wfeTT);
    fprintf('Tip-tilt noise\t\t%.4g\n',wfeNoiseTT);
    %fprintf('Tip-tilt anisoplanatism\t%g\n',wfeAnisoTT);
    fprintf('-------------------------------\n');
end

out = [srMar,srPar,wfeTot,wfeStatic,wfeFit,wfeLag,wfeNoise,wfeAl,wfeAni,...
    wfeTT,wfeNoiseTT];

psfr.wfe.sr_tot      = srMar;
psfr.wfe.wfe_tot     = wfeTot;
psfr.wfe.wfe_ncpa    = wfeStatic;
psfr.wfe.wfe_fit     = wfeFit;
psfr.wfe.wfe_lag      = wfeLag;
psfr.wfe.wfe_noise   = wfeNoise;
psfr.wfe.wfe_alias   = wfeAl;
psfr.wfe.wfe_aniso   = wfeAni;
psfr.wfe.wfe_tt  = wfeTT;
psfr.wfe.wfe_noiseTT = wfeNoiseTT;
