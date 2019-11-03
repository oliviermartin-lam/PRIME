function bestFit = updateResults(bestFit,beta,fRes,J)
            inputs = inputParser;
            inputs.addRequired('bestFit',@(x) isa(x,'bestFitting'));
            inputs.addRequired('beta',@isnumeric);
            inputs.addRequired('fRes',@isnumeric);
            inputs.addRequired('J',@isnumeric);
            inputs.parse(bestFit,beta,fRes,J);
            
            
            %1\ Compensate photometry and background for the normalization factor
            beta(1:bestFit.nStars) = beta(1:bestFit.nStars).*bestFit.normFactor;
            beta(end) = beta(end).*bestFit.normFactor;
            
            %2\ Get precision of estimates
            dx = diff(nlparci(real(beta),fRes,'jacobian',J),1,2);
            for i=1:length(dx)
                dx(i) = diff(nlparci(real(beta(i)),fRes,'jacobian',J(:,i)),1,2);
            end
            dx = reshape(dx,1,[]);
            bestFit.x_final = beta;
            bestFit.x_prec  = dx;
            
            %3\ Unpacking AO parameters
            wvl_ratio =  bestFit.psfr_.trs.cam.wavelength/500e-9;
            d = bestFit.nStars*3;
            % Seeing estimation
            idx         = bestFit.idxR0 +d;
            bestFit.r0_fit  = sum(beta(idx))^(-3/5);
            bestFit.r0_prec = 3/5*sum(beta(idx))^(-8/5)*sqrt(sum(dx(idx).^2));
            % Cn2 estimation
            idx         = bestFit.idxCn2 +d;
            bestFit.Cn2_fit = beta(idx);
            bestFit.Cn2_prec= dx(idx);
            if isempty(bestFit.Cn2_fit)
                bestFit.Cn2_fit = bestFit.r0_fit^(-5/3);
                bestFit.Cn2_prec= 5/3*bestFit.r0_prec*bestFit.r0_fit^(-8/3);
            end
            % Retrieved atmosphere
            if numel(bestFit.Cn2_fit) == 1
                hl = 0;
            else
                hl          = [bestFit.atm.layer.altitude];
            end
            
            nL  = length(bestFit.Cn2_fit*wvl_ratio^2);
            r0  = sum(bestFit.Cn2_fit*wvl_ratio^2)^(-3/5);
            fl  = bestFit.Cn2_fit/sum(bestFit.Cn2_fit);
           
            %bestFit.atm_fit = atmosphere(photometry.V0,r0,bestFit.atm.L0,'fractionnalR0',fl,...
              %  'altitude',hl,'layeredL0',bestFit.atm.L0,'windSpeed',wS,'windDirection',wD);
            
            % AO parameters
            bestFit.xao_fit = beta([bestFit.idxDho,bestFit.idxDtt,bestFit.idxDal]+3*bestFit.nStars);
            bestFit.xao_prec= dx([bestFit.idxDho,bestFit.idxDtt,bestFit.idxDal]+3*bestFit.nStars);
            
            %4\ Unpacking Stellar parameters
            bestFit.bg_fit           = beta(end);
            bestFit.xstars_fit       = beta(1:d);
            bestFit.xstars_prec      = dx(1:d);
            bestFit.catalog_fit.id   = 1:bestFit.nStars;
            bestFit.catalog_fit.RA   = bestFit.xstars_fit(bestFit.nStars+1:bestFit.nStars*2)*bestFit.psfr_.trs.cam.pixelScale;
            bestFit.catalog_fit.DEC  = bestFit.xstars_fit(2*bestFit.nStars+1:bestFit.nStars*3)*bestFit.psfr_.trs.cam.pixelScale;
            bestFit.catalog_fit.FLUX = bestFit.xstars_fit(1:bestFit.nStars);
            bestFit.catalog_fit.dRA  = bestFit.xstars_prec(bestFit.nStars+1:bestFit.nStars*2)*bestFit.psfr_.trs.cam.pixelScale;
            bestFit.catalog_fit.dDEC = bestFit.xstars_prec(2*bestFit.nStars+1:bestFit.nStars*3)*bestFit.psfr_.trs.cam.pixelScale;
            bestFit.catalog_fit.dFLUX= bestFit.xstars_prec(1:bestFit.nStars);
            bestFit.catalog_fit.dMAG = 2.5*bestFit.catalog_fit.dFLUX./bestFit.catalog_fit.FLUX/log(10);
            
            %5\ Grab the best-fitted PSF + 3-sigma PSF
            wMap        = bestFit.weightMap;
            bestFit.weightMap  = 1;
            bestFit.im.rec = imageModel(beta,bestFit.xdata,bestFit);
            if all(bestFit.x_prec <= max(abs(bestFit.lbounds),abs(bestFit.ubounds)))
                im_3sig_m   = bestFit.modelFUN(bestFit.x_final-bestFit.x_prec,bestFit.xdata);
                psf_3sig_m  = bestFit.psf_fit;
                im_3sig_p   = bestFit.modelFUN(bestFit.x_final+bestFit.x_prec,bestFit.xdata);
                psf_3sig_p  = bestFit.psf_fit;
                bestFit.im.rec_3sig = [im_3sig_m,im_3sig_p];
                bestFit.psf.rec_3sig= [psf_3sig_m,psf_3sig_p];
            end
            bestFit.weightMap= wMap;
            
        end
        
        