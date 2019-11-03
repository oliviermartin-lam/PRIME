%{
------------HEADER-----------------
Objective          ::  Estimate the noise covariance matrix
INPUT VARS
psfr              :: The PRIME top-level class
nMin, nMax  :: (optional) 'Min and Max Zernike modes order the wavefront is decomposed onto (defaut: 4, 100)
fitL0             :: (optional) fit the L0 as well if set to true (default: true)
best              :: (optional) Repeat the fitting process for different values of nMin from nMin to nMax/4 and pick up the most precise one
wvl_ref         :: (optional) Wavelength the r0 is given at.
mskPhase     :: (optional) Phase mask to reject non-valid actuators.
OUTPUT VARS
r0               :: The estimated line-of-sight r0 in meters at wvl
L0               :: The estimated L0 in meters
fwhm          :: The estimated PSF FWHM in arcsec
dr0             :: The 3-sigma precision on the r0 estimation
dL0             :: The 3-sigma precision on the L0 estimation
dfwhm        :: The 3-sigma precision on the fwhm estimation
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}


function res = estimateSeeingFromTelemetry(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.addParameter('nMin',4,@isnumeric);
inputs.addParameter('nMax',120,@isnumeric);
inputs.addParameter('fitL0',true,@islogical);
inputs.addParameter('best',true,@islogical);
inputs.addParameter('wvl',0.5e-6,@isnumeric);
inputs.addParameter('aoMode','NGS',@ischar);
inputs.addParameter('D1',9,@isnumeric);
inputs.addParameter('D2',2.65,@isnumeric);
inputs.parse(psfr,varargin{:});

%1\ Parsing inputs 
D1 = inputs.Results.D1; % Equivalent outer diameter of a circular pupil within which actuators commands are not sensitive to edges effects
D2 = inputs.Results.D2; % Same for the inner diameter
fitL0 = inputs.Results.fitL0;
wvl = inputs.Results.wvl;
nMin = inputs.Results.nMin;
nMax= inputs.Results.nMax;


%2\ Getting the reconstructed wavefront and noise covariance at 500 nm
uout = (2*pi/wvl)*psfr.trs.dm.com;%selectActuators(psfr.trs.hodm_pos,D1,D2);
D1 = 11.25;D2 = 0;
validActu = find(uout(:,1));
nValid = numel(validActu);

%2\ Get the noise covariance matrix and variance
Cn   = 0*psfr.res.noise.Cn_ho;
Cn(validActu,validActu) = 0*psfr.res.noise.Cn_ho(validActu,validActu)*(2*pi/wvl)^2;
varN = trace(Cn(validActu,validActu))/nValid;

%3\ Preliminary r0 estimation sig^2 = 0.111*(D(1-o)/r0)^5/3
varPh = sum(std(uout(validActu,:),[],2).^2)/nValid - varN;
r0_ = (0.111/varPh)^(3/5)*D1*(1-D2/D1);

%4\ Estimating the r0/outer scale
[res.r0,res.L0,res.dr0,res.dL0,res.zern.jindex,res.zern.model,res.zern.meas] =  getr0L0FromDMcommands(psfr.trs.mat.dmIF_hr,uout,psfr.trs.tel,Cn,...
    'best',true,'fitL0',fitL0,'initR0',r0_,'D1',D1,'D2',D2,'aoMode',psfr.trs.aoMode,'nMin',nMin,'nMax',nMax);

fprintf('r_0 estimated from the telemetry: %.3g cm\n',1e2*res.r0);
fprintf('L_0 estimated from the telemetry: %.3g m\n',res.L0);

%5\ Estimating the atmosphere-limited PSF FWHM
k1 = 1.03*wvl*constants.radian2arcsec;
res.fwhm = k1/res.r0;
res.dfwhm= k1*res.dr0/res.r0^2;

if res.L0 ~=0
    k2     = 2.813;
    k3     = 0.356;
    sqfact = sqrt(1-k2*(res.r0/res.L0)^k3);
    res.fwhm = res.fwhm*sqfact;
    
    dfdr01 = res.dr0/res.r0*res.fwhm;
    dfdr02 = 0.5*k1*res.dr0/res.r0/sqfact*k2*k3*res.L0^(-k3)*res.r0^(k3-1);
    dfdr0 = hypot(dfdr01,dfdr02);        
    dfdL0 = 0.5*k1*res.dL0/res.r0/sqfact*k2*k3*res.r0^(k3)*res.L0^(-k3-1);    
    res.dfwhm=hypot(dfdr0,dfdL0);
end

