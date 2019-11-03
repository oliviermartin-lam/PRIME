function bestFit = estimatePsfMetrics(bestFit)
            
            %6\ PSF statistics
            im           = bestFit.im_sky;
            dx           = bestFit.x_prec;
            bestFit.eqm_fit  = tools.getFVU(im,bestFit.im.rec);
            bestFit.eqm_init = tools.getFVU(im,bestFit.im.rec_init);
            if ~isempty(bestFit.im_3sig)
                eqm_p3s      = tools.getFVU(im,bestFit.im.rec_3sig(:,1:end/2,:)- dx(end));
                eqm_m3s      = tools.getFVU(im,bestFit.im.rec_3sig(:,end/2+1:end,:)+ dx(end));
                bestFit.eqm_3sig = sqrt((0.5*(eqm_p3s + eqm_m3s)).^2 - bestFit.eqm_fit.^2);
            end
            
            % Image statistics
            P = bestFit.psfr_.trs.tel.pupil;
            samp = bestFit.psfr_.trs.cam.samp;
            ps = bestFit.psfr_.trs.cam.pixelScale;
            
           
            [bestFit.SR_init,bestFit.dSR_init]  = tools.getStrehl(bestFit.im.rec_init,P,samp);
            [bestFit.SR_fit,bestFit.dSR_fit]    = tools.getStrehl(bestFit.im.rec,P,samp);            
            [bestFit.FWHM_init,bestFit.dFWHM_init]= tools.getFWHM(bestFit.im.rec_init,ps,8);
            [bestFit.FWHM_fit,bestFit.dFWHM_fit]= tools.getFWHM(bestFit.im.rec,ps,8);
        end
        
        