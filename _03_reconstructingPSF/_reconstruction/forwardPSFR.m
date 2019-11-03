function psfr = forwardPSFR(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.addParameter('flagToeplitz',false,@islogical);
inputs.addParameter('flagAnisoMethod','oomao',@ischar);
inputs.addParameter('flagNoiseMethod','autocorrelation',@ischar);
inputs.addParameter('flagAoPattern','circle',@ischar);
inputs.parse(psfr,varargin{:});

psfr.flags.toeplitz         = inputs.Results.flagToeplitz;
psfr.flags.anisoMethod = inputs.Results.flagAnisoMethod;
psfr.flags.noiseMethod = inputs.Results.flagNoiseMethod;
psfr.flags.aoPattern = inputs.Results.flagAoPattern;

%1\ Noise covariance matrices estimation
fprintf('Estimating noise covariance matrices\n');
psfr.res.noise = estimateNoiseCovarianceFromTelemetry(psfr.trs,'method',psfr.flags.noiseMethod);

%\2 Seeing estimation
fprintf('Estimating the seeing');
psfr.res.seeing = estimateSeeingFromTelemetry(psfr);

%3\ Diffraction-limit OTF - Nyquist sampling
psfr.otf.otfStat= computeStaticOpticalTransferFunction(psfr);

%4\ Normalized Fitting SF
[psfr.sf.Dfit,psfr.cov.psdFit] = computeFittingPhaseStructureFunction(psfr,'aoPattern',psfr.flags.aoPattern);

%5\ Normalized Alasing SF
[psfr.sf.Dal,psfr.cov.Cal] = computeAliasingPhaseStructureFunction(psfr,'aoPattern',psfr.flags.aoPattern);

%6\ AO residual SF
[psfr.sf.Dho_z,psfr.cov.Cho_z] = computeResidualPhaseStructureFunction(psfr);
psfr.sf.Dho = psfr.sf.Dho_z;
psfr.cov.Cho = psfr.cov.Cho_z;

%7\ Tip-tilt SF
[psfr.sf.Dtt,psfr.cov.Ctt] = computeTipTiltPhaseStructureFunction(psfr);

%8\ Anisoplanatism
[psfr.sf.Dani_l,psfr.sf.Dani,psfr.otf.Kani] = computeAnisoplanatismPhaseStructureFunction(psfr);

%9\ Reconstruct the PSF
psfr.psf.image = zeros(psfr.psf.fov,psfr.psf.fov,psfr.trs.src.nSrc);
for iSrc = 1:psfr.trs.src.nSrc
    psfr.otf.otfShannon = psfr.otf.otfStat.*exp(-0.5*(psfr.sf.Dfit+ psfr.sf.Dal+ psfr.sf.Dho + psfr.sf.Dtt + psfr.sf.Dani(:,:,iSrc)));  
    psf_ij = tools.otfShannon2psf(psfr.otf.otfShannon,psfr.trs.cam.samp(iSrc),psfr.psf.fov);
    psfr.psf.rec(:,:,iSrc) = psf_ij/sum(psf_ij(:));    
end

