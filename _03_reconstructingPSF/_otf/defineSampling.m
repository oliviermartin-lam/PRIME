function [dk,nT,nPup,nOtf] = defineSampling(nActuators,Samp,fov)

% Number of pixels in the PSD to describe the AO-corrected area
dk = 2*nActuators - 1;
% Number of pixels in the PSF to describe the AO-corrected area
aoBand  = (nActuators-1)*Samp;
% Ratio between the PSF fov and the AO-corrected field
nT   = max([1,round(fov/aoBand)]);
% Number of pixels to define the OTF to get a Nyquist PSF-sampling
nOtf     = 2*floor(dk*nT/2);
% Number of pixels to define the pupil accordingly
nPup = nOtf/2;


