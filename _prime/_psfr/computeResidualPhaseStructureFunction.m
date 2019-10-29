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
inputs.addRequired('psfr',@(x) isa(x,'prime'));
inputs.addParameter('nZernMode',[],@isnumeric);
inputs.parse(psfr,varargin{:});
nZernMode = inputs.Results.nZernMode;

%1\ Get the covariance matrix in the actuators space
psfr.trs.waveFront= bsxfun(@minus,psfr.trs.waveFront,mean(psfr.trs.waveFront,2));
Cho  = psfr.trs.waveFront*psfr.trs.waveFront'/psfr.trs.nExp;
%2\ Get the phase structure function in the actuators space
Du     = (2*pi/psfr.wvl)^2*(Cho - psfr.Cn_ho);
%3\ Get the point-wise phase structure function assuming stationarity
[~,sf_2D]= tools.modes2Otf(Du,psfr.trs.dmIF_hr,psfr.tel.pupil,psfr.nOtf,1,'method','Vii');

%4\ Get modal phase structure functions
% if the user provides a list of Zernike modes, this code will calculate the residual phase structure function for each mode + 
% the remaining one 

psfr.nGainsHO = 1; % if 1 the gain is common to each mode
if ~isempty(nZernMode)
    %4.1 Defining the modal basis, here Zernike
    nZ = numel(nZernMode);
    z  = zernike(nZernMode,'resolution',psfr.nActu);    
    psfr.nGainsHO = nZ+1;
    %4.2 Instantiating
    sfFit_2D_full = sf_2D;
    n = size(sf_2D,1);
    sf_2D = zeros(n,n,nZ+1);
    sf_2D(:,:,end) =  sfFit_2D_full;
    psfr.Hz = cell(1,nZ);    
    %4.3 Get the DM spatial filter
    bif     = xineticsInfluenceFunction(psfr.pitch); % specific to Keck
    dm       = deformableMirror(psfr.nActu,'modes',bif,'resolution',psfr.nActu);
    dmFilter= pinv(full(dm.modes.modes),1e-1);       
    %4.4 loop on considered modes.
    for k=1:nZ
        %4.4.1 Calculate the filtering matrix        
        zz = dmFilter*squeeze(z.modes(:,k));
        psfr.Hz{k} = zz*pinv(zz); % keep only the kth mode
        %4.4.2 Get the phase structure functions in the actuators space
        Du     = (2*pi/psfr.wvl)^2*psfr.Hz{k}*(Cho - psfr.Cn_ho)*psfr.Hz{k}';
        %4.4.3 Get the point-wise phase structure funciton assuming stationarity
        [~,sf_2D(:,:,k)]= tools.modes2Otf(Du,psfr.trs.dmIF_hr,psfr.tel.pupil,psfr.nOtf,1,'method','Vii');              
    end     
    %4.5 Get the remaining residual SF subtracted from zernike modes.
    sf_2D(:,:,k+1) = psfr.Dho - squeeze(sum(sf_2D(:,:,1:end-1),3));                 
end

