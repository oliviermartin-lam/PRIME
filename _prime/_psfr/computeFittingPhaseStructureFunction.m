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
            inputs.addRequired('psfr',@(x) isa(x,'prime'));           
            inputs.parse(psfr,varargin{:});

	    %1\ Parsing inputs
            r0 = psfr.atm.r0;
            L0 = psfr.atm.L0;
            nActu = psfr.nActu;
            d = psfr.pitch;
            nT = psfr.nTimes;
            
            %2\ Define the frequency space
            nK   = (2*nActu-1)*nT;
            [kx,ky] = freqspace(nK,'meshgrid');
            kc  = 1/(2*d);
            kx = kx*kc*nT;
            ky = ky*kc*nT;
            k   = hypot(kx,ky);
            dk  = kc/nActu; % Pixel scale
            
            %3\ Define the atmospheric phase PSD            
            cst = (24*gamma(6/5)/5)^(5/6)*(gamma(11/6)^2/(2*pi^(11/3)));
            psd = cst*(k.^2 + 1/L0.^2).^(-11/6);
            
            %4\ Filtering the AO-controlled area
            idx = k<=kc;
            %idx = abs(kx)<=kc & abs(ky)<=kc;
            psd(idx) = 0;
            
            %5\ Define the atmospheric phase PSD            
            sf_2D = tools.cov2sf(tools.psd2cov(psd,dk));
                        
            %if nargout > 2                
            %    %5\ Compute the atmospheric phase covariance map
            %    covMap  = real(tools.psd2cov(psd,dk));                
            %    %6\ Obtain the Toeplitz point-wise covariance matrix rd^2
            %    covMat  = tools.covMap2Matrix(covMap,nActu*nT-1,nActu*nT-1); % not working                
            %    %7\ Get the phase structure function
            %    sf_4D = 2*(diag(diag(covMat)) - covMat);
            % end
end
            

