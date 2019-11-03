%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the tip-tilt excluded residual phase

INPUT VARS
 psfr          :: The PRIME top-level class
 nZernMode :: (optinional) list of Zernike modes (ex: 4:6 for focus + astigmatisms terms) to enable modal gains retrieval
OUTPUT VARS
 sf_2D             :: bi-dimensional phase structure function map (Toeplitz) of the aliasing error (psfr.nOtf x psfr.nOtf)
Cho                  :: Covariance matrix of the tip-tilt excluded residual error in the acturators space (nActu^2 x nActu^2)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [sf_2D,Cho] = computeResidualPhaseStructureFunction(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.addParameter('nZernMode',[],@isnumeric);
inputs.parse(psfr,varargin{:});
psfr.trs.mat.nZernMode = inputs.Results.nZernMode;

%1\ Get the covariance matrix in the actuators space
Cho  = psfr.trs.rec.res*psfr.trs.rec.res'/psfr.trs.wfs.nExp;
%2\ Get the phase structure function in the actuators space
Du     = (2*pi/psfr.trs.cam.wavelength)^2*(Cho - psfr.res.noise.Cn_ho);
%3\ Get the point-wise phase structure function assuming stationarity
[~,sf_2D]= modes2Otf(Du,psfr.trs.mat.dmIF_hr,psfr.trs.tel.pupil,psfr.otf.nOtf,'method','Vii');

%4\ Get modal phase structure functions
%if the user provides a list of Zernike modes, this code will calculate the residual phase structure function for each mode +
%    the remaining one

if ~isempty(psfr.trs.mat.nZernMode)
    %4.1 Defining the modal basis, here Zernike
    nZ = numel(psfr.trs.mat.nZernMode);
    z  = zernike(psfr.trs.mat.nZernMode,'resolution',psfr.trs.dm.nActuators);
    %4.2 Instantiating
    sfFit_2D_full = sf_2D;
    n = size(sf_2D,1);
    sf_2D = zeros(n,n,nZ+1);
    sf_2D(:,:,end) =  sfFit_2D_full;
    psfr.trs.mat.Hz = cell(1,nZ);
    %4.3 Get the DM spatial filter
    bif     = xineticsInfluenceFunction(psfr.trs.dm.pitch); % specific to Keck
    dm       = deformableMirror(psfr.trs.dm.nActuators,'modes',bif,'resolution',psfr.trs.dm.nActuators);
    dmFilter= pinv(full(dm.modes.modes),1e-1);
    %4.4 loop on considered modes.
    for k=1:nZ
        %4.4.1 Calculate the filtering matrix
        zz = dmFilter*squeeze(z.modes(:,k));
        psfr.trs.mat.Hz{k} = zz*pinv(zz); % keep only the kth mode
        %4.4.2 Get the phase structure functions in the actuators space
        Du     = (2*pi/psfr.trs.cam.wavelength)^2*psfr.trs.mat.Hz{k}*(Cho - psfr.res.noise.Cn_ho)*psfr.trs.mat.Hz{k}';
        %4.4.3 Get the point-wise phase structure funciton assuming stationarity
        [~,sf_2D(:,:,k)]= modes2Otf(Du,psfr.trs.mat.dmIF_hr,psfr.trs.tel.pupil,psfr.otf.nOtf,'method','Vii');
    end
    %4.5 Get the remaining residual SF subtracted from zernike modes.
    sf_2D(:,:,k+1) = psfr.sf.Dho - squeeze(sum(sf_2D(:,:,1:end-1),3));
end

