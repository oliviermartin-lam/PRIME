function bestFit = fittingSetup(bestFit,varargin)
            inputs = inputParser;
            inputs.addRequired('bestFit',@(x) isa(x,'bestFitting'));
            inputs.addParameter('aoinit',[],@isnumeric);
            inputs.addParameter('aobounds',[],@isnumeric);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,true],@islogical);
            inputs.addParameter('nZernMode',[],@isnumeric);
            inputs.parse(bestFit,varargin{:});
            
            nZernMode_history = bestFit.psfr_.trs.mat.nZernMode;
            bestFit.nZernMode= inputs.Results.nZernMode;
            fitCn2       = inputs.Results.fitCn2;
            fitR0        = inputs.Results.fitR0;
            fitGains     = inputs.Results.fitGains;
            aoinit_      = inputs.Results.aoinit;
            aobounds_    = inputs.Results.aobounds;
            
            
            %1\ Recalculate the Residual phase structure function if the
            %user opts for a modal gain retrieval            
            if isempty(bestFit.nZernMode)
                bestFit.psfr_.sf.Dho_z = bestFit.psfr_.sf.Dho;
                bestFit.psfr_.cov.Cho_z = bestFit.psfr_.cov.Cho;
                bestFit.nGainsHO = 1;
            elseif  (~isempty(bestFit.nZernMode)) && (numel(bestFit.nZernMode) ~= numel(nZernMode_history) || (any(bestFit.nZernMode ~= nZernMode_history)))
                [bestFit.psfr_.sf.Dho_z,bestFit.psfr_.cov.Cho_z] = computeResidualPhaseStructureFunction(bestFit.psfr_,'nZernMode',bestFit.nZernMode);
                bestFit.nGainsHO = numel(bestFit.nZernMode)+1; % if 1 the gain is common to each mode
            end
            
            %2\ Define initial guesses
            if (~isempty(aoinit_)) && (~isempty(aobounds_))
                ao_init = aoinit_;
                ao_lb   = aobounds_(1,:);
                ao_ub   = aobounds_(2,:);
            else                                
                bestFit.idxCn2 = [];Cn2_init  = []; Cn2_lb    = [];
                Cn2_ub    = [];bestFit.idxR0  = [];                
                gain_init = [];gain_lb   = [];gain_ub   = [];                                
                nInit    = 1;       
                
                %2.1 Cn2 and r0 init                              
                if fitCn2
                    bestFit.idxCn2= nInit:nInit+bestFit.atm.nLayer-1;
                    bestFit.idxR0 = bestFit.idxCn2;
                    nInit     = nInit+bestFit.atm.nLayer;
                    Cn2_init  = [bestFit.atm.layer.fractionnalR0]*bestFit.atm.r0^(-5/3)*(bestFit.atm.wavelength/bestFit.wvl_ref)^2;
                    Cn2_lb    = 5^(-5/3)*ones(1,bestFit.atm.nLayer);
                    Cn2_ub    = 0.05^(-5/3)*ones(1,bestFit.atm.nLayer);
                else
                    if fitR0
                        bestFit.idxR0 = nInit;
                        nInit     = nInit+1;
                        Cn2_init  = bestFit.psfr_.trs.atm.r0^(-5/3)*(bestFit.psfr_.trs.atm.wavelength/bestFit.psfr_.trs.cam.wavelength)^2;
                        Cn2_lb    = 5^(-5/3);
                        Cn2_ub    = 0.01^(-5/3);
                    end
                end
                
                %2.2. HODM gain              
                if fitGains(1)
                    bestFit.idxDho = nInit:nInit+numel(bestFit.nZernMode);
                    nInit          = nInit+length(bestFit.idxDho);
                    gain_init    = ones(1,numel(bestFit.idxDho));
                    gain_lb      = [gain_lb,zeros(1,numel(bestFit.idxDho))];
                    gain_ub     = [gain_ub,1e2*ones(1,numel(bestFit.idxDho))];               
                end                
                
                %2.3. TTDM gain
                if fitGains(2)
                    bestFit.idxDtt   = nInit;
                    nInit        = nInit+1;
                    gain_init    = [gain_init,1];
                    gain_lb      = [gain_lb,0];
                    gain_ub      = [gain_ub,1e2];              
                end
                
                %2.4. Aliasing model gain
                if fitGains(3)
                    bestFit.idxDal = nInit;
                    gain_init    = [gain_init,1];
                    gain_lb      = [gain_lb,0];
                    gain_ub      = [gain_ub,1e2];                
                end
                
                %2.5 Concatenating vectors
                ao_init   = [Cn2_init,gain_init];
                ao_lb     = [Cn2_lb,gain_lb];
                ao_ub     = [Cn2_ub,gain_ub];
            end
            
            %3\  Define initial guess on stellar parameters
            uMax = 2; % PSF position may change by +/- 2 pixels
            photo_init = ones(1,bestFit.nStars)/bestFit.nStars;            
            pos_init   = [bestFit.psfr_.trs.src.x,bestFit.psfr_.trs.src.y]*constants.radian2mas/bestFit.psfr_.trs.cam.pixelScale;
            bg_init    = 0;
            photo_b   = [zeros(1,bestFit.nStars);2*ones(1,bestFit.nStars)];
            pos_b      = uMax*ones(1,2*bestFit.nStars);
            bg_b       = 5*std(reshape(bestFit.im_sky,bestFit.fov_sky^2,[]),1);
            
            %4\ Concatenating PSF and stellar parameters
            bestFit.x_init = [photo_init,pos_init,ao_init,bg_init];
            bestFit.lbounds= [photo_b(1,:),pos_init-pos_b,ao_lb,-bg_b];
            bestFit.ubounds= [photo_b(2,:),pos_init+pos_b,ao_ub,bg_b];
            
            %5\ Define the additional useful data
            bestFit.xdata = {bestFit.idxCn2,bestFit.idxR0,bestFit.idxDho,bestFit.idxDtt,bestFit.idxDal};
        end