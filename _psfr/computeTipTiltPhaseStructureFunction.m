%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the aliasing phase

INPUT VARS
 psfr          :: The PRIME top-level class

OUTPUT VARS
 sf_2D             :: bi-dimensional phase structure function map (Toeplitz) of the residual tip-tilt (psfr.nOtf x psfr.nOtf)
Ctt                  :: Covariance matrix of the tip-tilt error (2 x 2)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [sf_2D,Ctt] = computeTipTiltPhaseStructureFunction(psfr)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'prime'));
inputs.parse(psfr);

%1\ Discretizing the pupil domain
x          = linspace(-0.5,0.5,psfr.nOtf);
[X,Y]    = meshgrid(x);

%2\ Get the tip-tilt covariance matrix
psfr.trs.tipTilt = bsxfun(@minus,psfr.trs.tipTilt,mean(psfr.trs.tipTilt,2));
Ctt    = psfr.trs.tipTilt*psfr.trs.tipTilt'/size(psfr.trs.tipTilt,2);

%3\ Get the tip-tilt residual phase
Du_tt       = (2*pi/psfr.wvl_ref)^2*(Ctt - psfr.Cn_tt);
sf_2D   = Du_tt(1,1)*X.^2 + Du_tt(2,2)*Y.^2 + Du_tt(1,2)*X.*Y +Du_tt(2,1).*Y'.*X';
