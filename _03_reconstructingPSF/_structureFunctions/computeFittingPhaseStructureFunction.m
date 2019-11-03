%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the
aliasing phase

INPUT VARS
 psfr          :: The PRIME top-level class
 
OUTPUT VARS
 sf_2D             :: bi-dimensional phase structure function map (Toeplitz) of the fitting error (psfr.nOtf x psfr.nOtf)
psd                  :: Normalized (r0 = 1 m) Power spectrum density of the fitting error (psfr.nOtf x psfr.nOtf)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}
function [sf_2D,psd] = computeFittingPhaseStructureFunction(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.addParameter('aoPattern','circle',@ischar)
inputs.parse(psfr,varargin{:});
aoPattern = inputs.Results.aoPattern;

%1\ Parsing inputs
r0 = psfr.res.seeing.r0;
L0 = psfr.res.seeing.L0;
nActu = psfr.trs.dm.nActuators;
d = psfr.trs.dm.pitch;
nT = psfr.otf.nTimes;

%2\ Define the frequency space
nK   = psfr.otf.nOtf;
[kx,ky] = freqspace(nK,'meshgrid');
kc  = 1/(2*d);
kx = kx*kc*nT;
ky = ky*kc*nT;
k   = hypot(kx,ky);
dk  = kc/nActu; % Pixel scale

%3\ Define the atmospheric phase PSD
cst = r0^(-5/3)*(24*gamma(6/5)/5)^(5/6)*(gamma(11/6)^2/(2*pi^(11/3)));
psd = cst*(k.^2 + 1/L0.^2).^(-11/6);

%4\ Filtering the AO-controlled area
if strcmp(aoPattern,'circle')
    idx = k<=kc;
else
    idx = abs(kx)<=kc & abs(ky)<=kc;
end
psd(idx) = 0;

%5\ Define the atmospheric phase PSD
sf_2D = tools.cov2sf(tools.psd2cov(psd,dk));


end


