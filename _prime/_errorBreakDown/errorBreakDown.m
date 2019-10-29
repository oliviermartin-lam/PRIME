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

function out = errorBreakDown(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'prime'));
inputs.addParameter('zS', 0,@isnumeric);
inputs.addParameter('aS', 0,@isnumeric);
inputs.addParameter('display', true,@islogical);
inputs.parse(psfr);
display= inputs.Results.display;
zS     = inputs.Results.zS;
aS     = inputs.Results.aS;

S      = sum(psfr.otfDL(:));
r053 = psfr.r0_fit^(-5/3);

%1\ Get the static error from the static map (NCPA)
srStatic = sum(real(psfr.otfStat(:)))/S;
wfeStatic = tools.sr2wfe(srStatic,psfr.wvl);

%2\ Get the fitting error by integrating the Fitting PSD
srFit = sum(sum(psfr.otfDL.*exp(-0.5*psfr.Dfit*r053)))/S;
wfeFit = tools.sr2wfe(srFit,psfr.wvl);

%3\ Get the aliasing error from the covariance matrix
if ~isempty(psfr.idxDal)
    gAl  = psfr.xao_fit(psfr.idxDal - length(psfr.idxR0));
else
    gAl = 1;
end

srAl = sum(sum(psfr.otfDL.*exp(-0.5*psfr.Dal*r053*gAl)))/S;
wfeAl = tools.sr2wfe(srAl,psfr.wvl);

%4\ Get the noise errors
%4.1 HO WFS Noise
if ~isempty(psfr.idxR0) && ~isempty(psfr.idxDho)
    gHO = psfr.xao_fit(psfr.idxDho- length(psfr.idxR0));
else
    gHO = 1;
end

wfeNoise = 0;
cte = 1e9*sqrt(psfr.trs.holoop_pn/size(psfr.Cn_ho,1));
cn_tmp = psfr.Cn_ho;

for i=1:psfr.nGainsHO-1
    tmp = psfr.Hz{i}*psfr.Cn_ho*psfr.Hz{i}';
    wfeNoise  = wfeNoise + gHO(i)*trace(tmp);
    cn_tmp = cn_tmp - tmp;
end
wfeNoise = sqrt(wfeNoise + gHO(end)*trace(cn_tmp))*cte;

%4.2 TT WFS Noise
if ~isempty(psfr.idxR0) && ~isempty(psfr.idxDtt)
    gTT = psfr.xao_fit(psfr.idxDtt - length(psfr.idxR0));
else
    gTT = 1;
end
wfeNoiseTT = sqrt(trace(gTT*psfr.trs.ttloop_pn*psfr.Cn_tt))*1e9;

%5\. AO Bandwidth errors
srLag    = sum(sum(psfr.otfDL.*exp(-0.5*sum(bsxfun(@times,psfr.Dho_z,reshape(gHO,1,1,[])), 3))))/S;
wfeLag = sqrt(tools.sr2wfe(srLag,psfr.wvl)^2 - wfeNoise^2) ;

%6.  Residual tip-tilt
srTT    = sum(sum(psfr.otfDL.*exp(-0.5*psfr.Dtt*gTT)))/S;
wfeTT = sqrt(tools.sr2wfe(srTT,psfr.wvl)^2 - wfeNoiseTT^2);

%7\  Anisoplanatism
if zS~=psfr.src(1).zenith*constants.radian2arcsec
    sci   = source('zenith',zS*constants.arcsec2radian,'azimuth',aS*pi/180);
    psfr.atm_fit.wavelength = psfr.src(1).photometry;
    SFani = psfr.anisoplanaticStructureFunctionModel(psfr.gs);
    SFani = sum(bsxfun(@times, SFani, reshape( psfr.r0_fit^(-5/3)*[psfr.atm_fit.layer.fractionnalR0],1,1,[])), 3);
    srAni = real(sum(sum(exp(-0.5*SFani).*psfr.otfDL))/S);
    wfeAni= tools.sr2wfe(srAni,psfr.wvl);
elseif zS==psfr.src(1).zenith*constants.radian2arcsec && psfr.flagAniso
    srAni = real(sum(sum(psfr.Kani.*psfr.otfDL))/S);
    wfeAni= tools.sr2wfe(srAni, psfr.wvl);
else
    wfeAni = 0;
    srAni  = 1;
end

%8\  Total residual
wfeTot = sqrt(sum(wfeStatic^2 + wfeFit^2 + wfeLag^2 ...
    + wfeNoise^2 + wfeAni^2 + wfeTT^2 +  wfeNoiseTT^2 + wfeAl^2));

srMar = tools.wfe2sr(wfeTot,psfr.wvl);
SRho  = srStatic*srFit.*srAni*srLag*srAl;
srPar = 1e2*SRho/(1+tools.sr2wfe(srTT,psfr.wvl)^2*(2*pi*1e-9/psfr.wvl)^2) + (1-SRho)/(1+(psfr.tel.D/ psfr.r0_fit)^2);

% -----------------  DISPLAY
if display
    fprintf('-------------------------------\n');
    fprintf('Wavelength\t\t%.3g micron\n', psfr.wvl*1e6)
    fprintf('Strehl Image\t\t%.4g%s\t\n',1e2*psfr.SR_sky,'%');
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
