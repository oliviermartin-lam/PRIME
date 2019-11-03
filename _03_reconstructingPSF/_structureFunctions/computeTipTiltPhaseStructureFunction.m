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
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.parse(psfr);

%1\ Discretizing the pupil domain
x          = linspace(-0.5,0.5,psfr.otf.nOtf);
[X,Y]    = meshgrid(x);

%2\ Get the tip-tilt covariance matrix
Ctt    = psfr.trs.tipTilt.com*psfr.trs.tipTilt.com'/psfr.trs.tipTilt.nExp;

%3\ Get the tip-tilt residual phase
Du_tt       = (2*pi/psfr.trs.cam.wavelength)^2*(Ctt - psfr.res.noise.Cn_tt);
sf_2D   = Du_tt(1,1)*X.^2 + Du_tt(2,2)*Y.^2 + Du_tt(1,2)*X.*Y +Du_tt(2,1).*Y'.*X';
