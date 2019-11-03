function psfr = updateReconstructionPrerequisites(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.parse(psfr,varargin{:});

%1\ Update the sampling
[psfr.otf.dk,psfr.otf.nTimes,psfr.trs.tel.resolution,psfr.otf.nOtf] = defineSampling(psfr.trs.dm.nActuators,psfr.trs.cam.samp,psfr.psf.fov);
            
%2\ Update the instrumental features size
psfr.res.static.stat_map_interp = tools.interpolate(psfr.trs.tel.static_map,psfr.trs.tel.resolution);
psfr.trs.tel.pupil = double(logical(psfr.res.static.stat_map_interp));       

% Get the Pupil OTF
P   = tools.enlargePupil(psfr.trs.tel.pupil ,2);
psfr.otf.otfDL = fftshift(tools.fftCorrel(P,P));
psfr.otf.otfDL = psfr.otf.otfDL /max(psfr.otf.otfDL(:));

%3\ Update the influence function
psfr = updateInfluenceFunction(psfr);