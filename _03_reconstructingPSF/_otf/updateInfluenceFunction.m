function psfr = updateInfluenceFunction(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.parse(psfr,varargin{:});

% Update the influence function
bif          = xineticsInfluenceFunction(psfr.trs.dm.pitch);
dmSq     = deformableMirror(psfr.trs.dm.nActuators,'modes',bif,'resolution',psfr.trs.tel.resolution);
psfr.trs.mat.dmIF_hr = tools.idlToMatlabIndex(dmSq.modes.modes,psfr.trs.dm.nActuators,true(psfr.trs.dm.nActuators),'modes');
psfr.trs.mat.dmIF_inv_hr = pinv(full(psfr.trs.mat.dmIF_hr));