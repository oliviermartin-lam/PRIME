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

function otfStat = computeStaticOpticalTransferFunction(psfr)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.parse(psfr);

% Get Telescope + static aberrations OTF (PSF Nyquist Sampling)
P             =  tools.enlargePupil(psfr.trs.tel.pupil,2);
phi_stat   = tools.enlargePupil(psfr.res.static.stat_map_interp*1e-9*2*pi/psfr.trs.cam.wavelength ,2);
E = P.*exp(1i*phi_stat);
otfStat = fftshift(tools.fftCorrel(E,E));
otfStat =otfStat /max(otfStat(:));
