%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the aliasing phase

INPUT VARS
 psfr          :: The PRIME top-level class

OUTPUT VARS
otfDL              :: Optical transfer function of the telescope pupil only (psfr.nOtf x psfr.nOtf)
otfStat            :: Optical transfer function of the telescope pupil including the static phase psfr.static_map (psfr.nOtf x psfr.nOtf)
Created by      :: O. Beltramo-Martin - ONERA/LAM
Creation date  :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [otfDL,otfStat] = computeStaticOpticalTransferFunction(psfr)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'prime'));
inputs.parse(psfr);

%1\ Interpolating the telescope pupil and static map regarding psfr.nOtf
psfr.tel.resolution = floor(psfr.nOtf/2);
psfr.stat_map_interp = tools.interpolate(psfr.trs.static_map,floor(psfr.nOtf/2));
psfr.tel.pupil = double(logical(psfr.stat_map_interp));
phi_stat   = tools.enlargePupil(psfr.stat_map_interp*1e-9*2*pi/psfr.wvl_ref ,2);

%2\ Get the Telescope OTF (PSF Nyquist Sampling)
P   = tools.enlargePupil(psfr.tel.pupil ,2);
otfDL = fftshift(tools.fftCorrel(P,P));
otfDL =otfDL /max(otfDL(:));

%3\ Telescope + static aberrations OTF (PSF Nyquist Sampling)
E = P.*exp(1i*phi_stat);
otfStat = fftshift(tools.fftCorrel(E,E));
otfStat =otfStat /max(otfStat(:));
