function trs = restoreNIRC2Image(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
%1\ Getting data ID
date = trs.date;
path_calib = trs. path_calib;
path_im = trs.path_im;

if strcmp(date,'20130801')
    badPixMap          = fitsread([path_calib,'supermask2013_modified.fits']);
    badPixMap(169,179) = 1;
    badPixMap(39,707)  = 1;
    flat               = fitsread([path_calib,'flat_kp2017.fits']);
else
    badPixMap = fitsread([path_calib,'supermask2017.fits']);
    flat      = fitsread([path_calib,'flat_kp2017.fits']);
    badPixMap(289,136) = 1;
    badPixMap(290,136) = 1;
    badPixMap(289,137) = 1;
    badPixMap(290,137) = 1;
    badPixMap(694,860) = 1;
    badPixMap(578,924) = 1;
    badPixMap(720,859) = 1;
    badPixMap(411,762) = 1;
    badPixMap(449,747) = 1;
end

% Fits info
list      = trs.fitsHdr(:,1);
val       = trs.fitsHdr(:,2);
% Wavelength and PSF sampling
trs.cam.wavelength      = cell2mat(val(strcmp(list,'CENWAVE')))*1e-6;
trs.cam.samp     = 1.5*(trs.cam.wavelength/1.6455e-6);
trs.tel.airmass  = cell2mat(val(strcmp(list,'AIRMASS')));

% Post-processing
im0     = fitsread(path_im);
N = size(im0,1);
if N == 1024
    bMap = badPixMap;
    fMap = flat;
elseif N == 512 %cropped mode (not binning)
    bMap = tools.crop(badPixMap,N);
    fMap = tools.crop(flat,N);
end
trs.cam.frame  = tools.processImage(im0,0,fMap,bMap,trs.cam.samp,...
    'fovInPixel',trs.cam.resolution,'masking',false,'rebin',0,'tfccd',false,...
    'thresholding',-Inf);

% Get Strehl value
 [trs.cam.SR,trs.cam.dSR]    = tools.getStrehl(trs.cam.frame ,trs.tel.pupil,trs.cam.samp);
 [trs.cam.FWHM,trs.cam.dFWHM]= tools.getFWHM(trs.cam.frame,trs.cam.pixelScale,8);
end

