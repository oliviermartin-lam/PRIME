classdef prime < handle
    
    properties
        % --------------------- Observation
        trs;                    %Telemetry class
        im_sky;                  % Observations
        fov_sky;                 % Image fov in pixel
        nStars;                  % Number of stars per frame
        nPSFs;                   % Number of PSFs in the 3D cube
        wvl;                        % Imaging wavelength;
        wvl_ref;
        airmass;                 % Airmass value
        ron;                     % Read out noise in ADU
        SR_sky;                  % Image Strehl ratio [0-1]
        dSR_sky;                 % Uncertainty of the image Strehl ratio [0-1]
        FWHM_sky;                % Image full width at half maximum in mas
        dFWHM_sky;               % Uncertainty of the FWHM
        
        % --------------------- Sub-classes
        tel;                     % telescope class
        atm;                     % atmosphere class
        src;                     % Sources
        gs;                      % Guide star
        ttgs;                    % tip-tilt star
        fao;                     % Spatial frequency adaptive optics
        % --------------------- System configuration
        stat_map_interp;
        stat_map_interp_zer;
        var_n;                   % HO WFS noise variance
        var_n_tt;                % TT WFS noise variance
        r0_tel;
        fwhm_tel;
        dfwhm_tel;
        dr0_tel;
        L0_tel;
        dL0_tel;
        % --------------------- OTFs/Model
        % Sampling
        dk;
        nActu;
        pitch;
        psInMas;
        nTimes;
        fov_fit;
        Samp;
        nOtf;
        nGainsHO;
        % Covariance/PSD
        Cho;
        Cho_z;
        Ctt;
        Cn_ho;
        Cn_tt;
        Cal;
        psdFit;
        Dho;
        Dho_z;
        Dtt;
        Dal;
        Dfit;
        Dani_l;
        Dani;
        nZernMode;
        zer;
        Hz=1;
        rmZerFromPhase = 1;
        
        % OTF
        otfDL;
        otfStat;
        otfAO;
        otfDM;
        otfTT;
        otfFit;
        otfAl;
        otfOnAxis;
        Kani;
        otfShannon;
        otfDen;
        mskOtf;
        % --------------------- Fitting process setup
        % PSF model Initial guess
        modelFUN;
        psf_init;
        psf_ext;
        im_init;
        x_init;
        lbounds;
        ubounds;
        fitOption;
        ydata;
        xdata;
        weightMap;
        normFactor;
        idxR0;
        idxCn2;
        idxDho;
        idxDtt;
        idxDal;        
        x_final;
        x_prec;
        beta_;
        J_;
        fRes_;
        % Stars Initial guess
        nParamStars;            % Number of stars parameters
        p_nStars;               % Number of PSFs. Can Set to one if the PSF does not spatially vary in the FOV
        initStars;              % Stars parameters initial guess
        xStars;yStars;          % Astrometry initial guess
        fluxStars;              % Photometry initial guess
        % Time counter
        t_inst;
        t_fit;
        t_mod;
        % --------------------- PSF outputs
        psf_fit;
        im_fit;
        psf_3sig;
        im_3sig;
        eqm_init;
        eqm_fit;
        eqm_3sig;
        bg_fit;
        SR_fit;
        dSR_fit;
        FWHM_fit;
        dFWHM_fit;
        psdAO;
        SR_init;
        dSR_init;
        FWHM_init;
        dFWHM_init;
        % --------------------- Stellar parameters outputs
        xstars_fit;
        xstars_prec;
        catalog_fit;
        % --------------------- Atmosphere parameters outputs
        atm_fit;
        r0_fit;
        r0_prec;
        Cn2_fit;
        Cn2_prec;
        % --------------------- AO parameters outputs
        xao_fit;
        xao_prec;
        map_fit;
        % --------------------- Error Breakdown
        wfe;
        % --------------------- Flags
        flagToeplitz       = true;
        flagAniso          = false;
        flagNoiseMethod    = 'autocorrelation';  %'RTF'
        flagAnisoMethod    = 'oomao';            %'FOURIER','FLICKER'
        flagAliasingMethod = 'PSD';              %'COV'
        flagInitPSFR       = false;
    end
    
    methods
        
        function obj = prime(trs,varargin)
            inputs = inputParser;
            inputs.addRequired('trs',@(x) isa(x,'telemetry'));
            inputs.addParameter('fov_fit',[],@isnumeric);
            inputs.addParameter('monochromaticBand',10e-9,@isnumeric);
            inputs.addParameter('psf_ext',[],@isnumeric);
            inputs.addParameter('flagToeplitz',false,@islogical);
            inputs.addParameter('flagNoiseMethod','autocorrelation',@ischar);
            inputs.addParameter('flagAnisoMethod','oomao',@ischar);
            inputs.addParameter('flagAliasingMethod','PSD',@ischar);
            inputs.parse(trs,varargin{:});
            
            % Grab inputs
            obj.trs = trs;
            obj.ron                = inputs.Results.trs.parm.cam.ron;
            obj.flagToeplitz       = inputs.Results.flagToeplitz;
            obj.fov_fit            = inputs.Results.fov_fit;
            obj.flagNoiseMethod    = inputs.Results.flagNoiseMethod;
            obj.flagAnisoMethod    = inputs.Results.flagAnisoMethod;
            obj.flagAliasingMethod = inputs.Results.flagAliasingMethod;
            obj.psf_ext            = inputs.Results.psf_ext;
            
            if isempty(obj.fov_fit)
                obj.fov_fit = trs.parm.cam.resolution+4;
            end
            
            %1\ Instantiate the system/image model configuration
            tic();
            obj = obj.instantiateTopLevelClasses();
            t0 = toc();
            
            %2\ Instantiate the PSF reconstruction
            tic;
            obj = obj.forwardPSFR(obj.wvl,obj.fov_fit,trs.parm.cam.pixelScale,...
                'flagToeplitz',obj.flagToeplitz,'flagAnisoMethod',obj.flagAnisoMethod);
            t1 = toc();
            
            obj.t_inst = t0+t1;
            %3\ Prompt messages
            fprintf('-----------------------------------------\n');
            fprintf('Time for class instantiation :\t %.3g s\n',obj.t_inst);
            fprintf('-----------------------------------------\n');
        end
        
        %% %%%%%%%%%%%%%%%%%%%%
        %                 INSTANTIATION                    %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = instantiateTopLevelClasses(obj)
            
            fprintf('Instantiate the system...');
            
            parm_ = obj.trs.parm;
            
            %1\ Detector
            obj.wvl = obj.trs.wvl;
            
            %2\ Scientific source
            photoSci    = photometry_kasp(obj.wvl,0,0);
            [zS,aS]     = utilities.arcsec2polar(parm_.sci.x,parm_.sci.y);
            obj.src     = source('zenith',zS,'azimuth',aS,...
                'wavelength',photoSci);
            obj.nStars  = numel(obj.src);
            
            %3\ Atmosphere
            obj.atm   = atmosphere(parm_.atm.photometry,parm_.atm.r0,mean(parm_.atm.L0),...
                'layeredL0',parm_.atm.L0,'fractionnalR0',parm_.atm.fractionalR0,...
                'altitude',parm_.atm.altitude,'windSpeed',parm_.atm.windSpeed,...
                'windDirection',parm_.atm.windDirection);
            
            %4\ Telescope
            obj.dk = 2*parm_.dm.nActuators+1;
            obj.tel = telescope(parm_.tel.D,'resolution',obj.dk,...
                'obstructionRatio',parm_.tel.obstructionRatio);
            obj.airmass = 1/cos(obj.atm.zenithAngle);
            
            %5\ DM
            obj.pitch = parm_.dm.pitch;
            obj.nActu = parm_.dm.nActuators;
            
            %6\ Guide stars
            if isfield(parm_,'lGs')
                [zS,aS]     = utilities.arcsec2polar(parm_.lGs.x,parm_.lGs.y);
                obj.gs      = source('zenith',zS,'azimuth',aS,...
                    'wavelength',parm_.lGs.photometry,'height',parm_.lGs.height);
                
                [zS,aS]     = utilities.arcsec2polar(parm_.nGs.x,parm_.nGs.y);
                obj.ttgs    = source('zenith',zS,'azimuth',aS,...
                    'wavelength',parm_.nGs.photometry);
            else
                [zS,aS]     = utilities.arcsec2polar(parm_.nGs.x,parm_.nGs.y);
                obj.gs      = source('zenith',zS,'azimuth',aS,...
                    'wavelength',parm_.nGs.photometry);
                obj.ttgs    = [];
            end
            
            fprintf('...Done\n');
        end
        
        %% %%%%%%%%%%%%%%%%%%%%
        %            PSF RECONSTRUCTION               %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = forwardPSFR(obj,wvl,fov,pixelScale,varargin)
            inputs = inputParser;
            inputs.addRequired('wvl',@isnumeric);
            inputs.addRequired('fov',@isnumeric);
            inputs.addRequired('pixelScale',@isnumeric );
            inputs.addParameter('Samp',obj.trs.Samp,@isnumeric);
            inputs.addParameter('flagToeplitz',false,@islogical);
            inputs.addParameter('flagAnisoMethod','oomao',@ischar);
            inputs.parse(wvl,fov,pixelScale,varargin{:});
            
            obj.flagToeplitz    = inputs.Results.flagToeplitz;
            obj.flagAnisoMethod = inputs.Results.flagAnisoMethod;
            
            %1\ Manage sampling
            obj.wvl_ref = wvl; obj.fov_fit  = fov; obj.psInMas  = pixelScale;
            obj.dk          = 2*obj.nActu - 1;
            obj.Samp     = inputs.Results.Samp;
            aoBand       = obj.tel.D*obj.Samp/obj.pitch;
            obj.nTimes   = max([1,round(obj.fov_fit/aoBand)]);
            obj.nOtf     = floor(obj.dk*obj.nTimes);
            obj.tel.resolution = floor(obj.nOtf/2);
            % Update the influence function
            bif          = xineticsInfluenceFunction(obj.pitch);
            dmSq     = deformableMirror(obj.nActu,'modes',bif,'resolution',obj.tel.resolution);
            obj.trs.dmIF_hr = tools.idlToMatlabIndex(dmSq.modes.modes,obj.nActu,true(obj.nActu),'modes');
            obj.trs.dmIF_inv_hr = pinv(full(obj.trs.dmIF_hr));
            
            %2\ Noise covariance matrices estimation
            fprintf('Estimating noise covariance matrices\n');
            [obj.Cn_ho,obj.Cn_tt,obj.var_n,obj.var_n_tt] = estimateNoiseCovarianceFromTelemetry(obj.trs,'method',obj.flagNoiseMethod);
            
            %\2 Seeing estimation
            fprintf('Estimating the seeing');
            [obj.r0_tel,obj.L0_tel,obj.fwhm_tel,obj.dr0_tel,obj.dL0_tel,obj.dfwhm_tel] = ...
                estimateSeeingFromTelemetry(obj,'mskPhase',obj.trs.validActu);
            obj.atm.r0 = obj.r0_tel;
            obj.atm.L0 = obj.L0_tel;
            
            %3\ Diffraction-limit OTF - Nyquist sampling
            [obj.otfDL, obj.otfStat] = computeStaticOpticalTransferFunction(obj);
            
            %4\ Normalized Fitting SF
            [obj.Dfit,obj.psdFit] = computeFittingPhaseStructureFunction(obj);
            
            %5\ Normalized Alasing SF
            [obj.Dal,obj.Cal] = computeAliasingPhaseStructureFunction(obj);
            
            %6\ AO residual SF
            [obj.Dho_z,obj.Cho_z] = computeResidualPhaseStructureFunction(obj);
            obj.Dho = obj.Dho_z;
            obj.Cho = obj.Cho_z;
            
            %7\ Tip-tilt SF
            [obj.Dtt,obj.Ctt] = computeTipTiltPhaseStructureFunction(obj);
            
            %8\ Anisoplanatism
            [obj.Dani_l,obj.Dani,obj.Kani] = computeAnisoplanatismPhaseStructureFunction(obj);
            
            %9\ Reconstruct the PSF
            obj.psf_init = obj.psfModel([],{[],[],[],[],[],[],[]});
            
        end
        
        function [psf] = psfModel(obj,x,xdata)
            
            psf = zeros(obj.fov_fit,obj.fov_fit,obj.nStars);
            
            % ----------------- Check inputs
            
            % Initial values
            gHO  = ones(1,obj.nGainsHO); gTT  = 1;gAl  = 1;
            r053 = obj.atm.r0^(-5/3);Cn2=[];
            
            obj.idxCn2 = xdata{1};
            obj.idxR0 = xdata{2};
            obj.idxDho = xdata{3};
            obj.idxDtt   = xdata{4};
            obj.idxDal = xdata{5};
            
            if ~isempty(obj.idxR0)
                r053= sum(x(obj.idxR0));
            end
            if ~isempty(obj.idxCn2)
                Cn2 = x(obj.idxCn2);
                r053 = sum(Cn2);
            end
            if ~isempty(obj.idxDho)
                gHO = x(obj.idxDho);
            end
            if ~isempty(obj.idxDal)
                gAl = x(obj.idxDal);
            end
            if ~isempty(obj.idxDtt)
                gTT = x(obj.idxDtt);
            end
            
            %2\ Get the on-axis AO residual phase SF
            Don = r053*(obj.Dfit+ obj.Dal*gAl) + sum(bsxfun(@times,obj.Dho_z,reshape(gHO,1,1,[])), 3) + obj.Dtt*gTT;%*meter2rad^2;
            
            for iSrc = 1:obj.nStars
                %3\ Anisoplanatism
                if numel(Cn2) > 1 && obj.flagAniso
                    %Cn2_wvl = Cn2*wvl_ratio^2;
                    obj.Dani   = sum(bsxfun(@times, squeeze(obj.Dani_l(:,:,:,iSrc)), reshape(Cn2,1,1,[])), 3);
                else
                    obj.Dani = 0;
                end
                
                %4\ OTF Multiplication
                obj.otfShannon = obj.otfStat.*exp(-0.5*(Don + obj.Dani));
                %5\ PSF calculation
                psf_ij = tools.otfShannon2psf(obj.otfShannon,obj.Samp(iSrc),obj.fov_fit);
                S = sum(psf_ij(:));
                psf(:,:,iSrc) = psf_ij/S;
            end                       
        end
        
        function [im, J_im] = multiplePSFModel(obj,x,xdata)
            
            % Sort parameters
            xS = x(1:3*obj.nStars);
            xPSF = x(1+3*obj.nStars:end-1);
            xBg = x(end);
            
            % Get the PSF for demanded field positions and wavelengths
            obj.psf_fit = obj.psfModel(xPSF,xdata);
            
            % Manage the inputs
            nSrc = length(xS)/3;
            fSrc = xS(1:nSrc);
            xSrc = xS(1+nSrc:2*nSrc);
            ySrc = xS(1+2*nSrc:3*nSrc);
            dX   = [xSrc;ySrc];
            
            im = zeros(obj.fov_fit,obj.fov_fit);
            for iSrc = 1:obj.nStars
                % Stack PSFs - Residual chromatism effects are not included yet
                psfi = obj.psf_fit(:,:,iSrc);
                psfi = psfi/sum(psfi(:));
                % Translate PSFs at sources position
                im  = im + fSrc(iSrc)*tools.translateImage(psfi,dX(:,iSrc));
            end
            im = im + xBg;
            if obj.fov_sky ~=obj.fov_fit
                im = tools.crop(im,obj.fov_sky);
            end
            im(obj.im_sky == 0) = 0;
            im = im.*obj.weightMap;
            
            %%%%%%%%%%%%%%%%%%%%%
            %JACOBIAN CALCULATION
            
            if nargout >1
                
                %1\ Jacobian with respect to photometry and astrometry
                %Those ones are ok
                
                nX = numel(xS) + numel(xPSF) + numel(xBg);
                J_im = zeros(obj.fov_sky^2,nX);
                for iSrc = 1:nSrc
                    % Photometry
                    J_im(:,1) = im(:)/fSrc;
                    % Astrometry
                    [u,v] = freqspace(size(im),'meshgrid');
                    phasor = exp(-1i*pi*(v*dX(1) + u*dX(2)));
                    otf = tools.psf2otf(im);
                    otf = otf/max(otf(:));
                    J_im(:,2)  = reshape(real(fftshift(ifft2(fftshift(-1i*pi*v.*otf)))).*obj.weightMap,obj.fov_sky^2,1);
                    J_im(:,3) = reshape(real(fftshift(ifft2(fftshift(-1i*pi*u.*otf)))).*obj.weightMap,obj.fov_sky^2,1);
                end
                
                %2\ Jacobian with respect to PSF parameters
                %These ones are not fully ok
                if ~isempty(obj.idxDal)
                    gAl = xPSF(obj.idxDal);
                else
                    gAl = 1;
                end
                if ~isempty(obj.idxR0)
                    r053= sum(xPSF(obj.idxR0));
                else
                    r053 = obj.atm.r0^(-5/3);
                end
                
                % I'm not sure that the following is a proper way to get
                % the derivative
                n=1;
                if ~isempty(obj.idxR0)
                    otf_1 = tools.interpolateOtf(padarray(-0.5*(obj.Dfit+ obj.Dal*gAl).*obj.otfShannon,round((obj.Samp-1)*[obj.fov_fit,obj.fov_fit]/2),'both'),obj.fov_sky);
                    J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,obj.fov_sky^2,[]);
                    n = n+1;
                end
                if ~isempty(obj.idxDho)
                    for k=1:obj.nGainsHO
                        otf_1 = tools.interpolateOtf(padarray(-0.5*obj.Dho_z(:,:,k).*obj.otfShannon,round((obj.Samp-1)*[obj.fov_fit,obj.fov_fit]/2),'both'),obj.fov_sky);
                        J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,obj.fov_sky^2,[]);
                        n = n+1;
                    end
                end
                if ~isempty(obj.idxDtt)
                    otf_1 = tools.interpolateOtf(padarray(-0.5*obj.Dtt.*obj.otfShannon,round((obj.Samp-1)*[obj.fov_fit,obj.fov_fit]/2),'both'),obj.fov_sky);
                    J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,obj.fov_sky^2,[]);
                    n = n+1;
                end
                if ~isempty(obj.idxDal)
                    otf_1 = tools.interpolateOtf(padarray(-0.5*r053*obj.Dal.*obj.otfShannon,round((obj.Samp-1)*[obj.fov_fit,obj.fov_fit]/2),'both'),obj.fov_sky);
                    J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*obj.weightMap,obj.fov_sky^2,[]);
                    n = n+1;
                end
                % Jacobian with respect to background
                %This one is ok
                J_im(:,3*nSrc+n:end) = obj.weightMap(:);
            end
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%
        %               MODEL-FITTING                      %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = fittingSetup(obj,varargin)
            inputs = inputParser;
            inputs.addParameter('aoinit',[],@isnumeric);
            inputs.addParameter('aobounds',[],@isnumeric);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,true],@islogical);
            inputs.addParameter('nZernMode',[],@isnumeric);
            inputs.parse(varargin{:});
            
            nZernMode_history = obj.nZernMode;
            obj.nZernMode= inputs.Results.nZernMode;
            fitCn2       = inputs.Results.fitCn2;
            fitR0        = inputs.Results.fitR0;
            fitGains     = inputs.Results.fitGains;
            aoinit_      = inputs.Results.aoinit;
            aobounds_    = inputs.Results.aobounds;
            
            
            %1\ Recalculate the Residual phase structure function if the
            %user opts for a modal gain retrieval
            if isempty(obj.nZernMode)
                obj.Dho_z = obj.Dho;
                obj.Cho_z = obj.Cho;
                obj.nGainsHO = 1;
            elseif  (~isempty(obj.nZernMode)) && (numel(obj.nZernMode) ~= numel(nZernMode_history) || (any(obj.nZernMode ~= nZernMode_history)))
                [obj.Dho_z,obj.Cho_z] = computeResidualPhaseStructureFunction(obj,'nZernMode',obj.nZernMode);
            end
            
            %2\ Define initial guesses
            if (~isempty(aoinit_)) && (~isempty(aobounds_))
                ao_init = aoinit_;
                ao_lb   = aobounds_(1,:);
                ao_ub   = aobounds_(2,:);
            else                                
                obj.idxCn2 = [];Cn2_init  = []; Cn2_lb    = [];
                Cn2_ub    = [];obj.idxR0  = [];                
                gain_init = [];gain_lb   = [];gain_ub   = [];                                
                nInit    = 1;       
                
                %2.1 Cn2 and r0 init                              
                if fitCn2
                    obj.idxCn2= nInit:nInit+obj.atm.nLayer-1;
                    obj.idxR0 = obj.idxCn2;
                    nInit     = nInit+obj.atm.nLayer;
                    Cn2_init  = [obj.atm.layer.fractionnalR0]*obj.atm.r0^(-5/3)*(obj.atm.wavelength/obj.wvl_ref)^2;
                    Cn2_lb    = 5^(-5/3)*ones(1,obj.atm.nLayer);
                    Cn2_ub    = 0.05^(-5/3)*ones(1,obj.atm.nLayer);
                else
                    if fitR0
                        obj.idxR0 = nInit;
                        nInit     = nInit+1;
                        Cn2_init  = obj.atm.r0^(-5/3)*(obj.atm.wavelength/obj.wvl_ref)^2;
                        Cn2_lb    = 5^(-5/3);
                        Cn2_ub    = 0.01^(-5/3);
                    end
                end
                
                %2.2. HODM gain              
                if fitGains(1)
                    obj.idxDho = nInit:nInit+numel(obj.nZernMode);
                    nInit          = nInit+length(obj.idxDho);
                    gain_init    = ones(1,numel(obj.idxDho));
                    gain_lb      = [gain_lb,zeros(1,numel(obj.idxDho))];
                    gain_ub     = [gain_ub,1e2*ones(1,numel(obj.idxDho))];               
                end                
                
                %2.3. TTDM gain
                if fitGains(2)
                    obj.idxDtt   = nInit;
                    nInit        = nInit+1;
                    gain_init    = [gain_init,1];
                    gain_lb      = [gain_lb,0];
                    gain_ub      = [gain_ub,1e2];              
                end
                
                %2.4. Aliasing model gain
                if fitGains(3)
                    obj.idxDal = nInit;
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
            photo_init = ones(1,obj.nStars)/obj.nStars;
            [xS,yS]    = pol2cart([obj.src.azimuth],[obj.src.zenith]);
            pos_init   = [yS,xS]*constants.radian2mas/obj.psInMas;
            bg_init    = 0;
            photo_b   = [zeros(1,obj.nStars);2*ones(1,obj.nStars)];
            pos_b      = uMax*ones(1,2*obj.nStars);
            bg_b       = 5*std(reshape(obj.im_sky,obj.fov_sky^2,[]),1);
            
            %4\ Concatenating PSF and stellar parameters
            obj.x_init = [photo_init,pos_init,ao_init,bg_init];
            obj.lbounds= [photo_b(1,:),pos_init-pos_b,ao_lb,-bg_b];
            obj.ubounds= [photo_b(2,:),pos_init+pos_b,ao_ub,bg_b];
            
            %5\ Define the additional useful data
            obj.xdata = {obj.idxCn2,obj.idxR0,obj.idxDho,obj.idxDtt,obj.idxDal};
        end
        
        function obj = bestFitting(obj,im_sky,varargin)
            inputs = inputParser;
            inputs.addRequired('im_sky',@isnumeric);
            inputs.addParameter('aoinit',[],@isnumeric);
            inputs.addParameter('aobounds',[],@isnumeric);
            inputs.addParameter('MaxIter',100,@isnumeric);
            inputs.addParameter('TolX',1e-10,@isnumeric);
            inputs.addParameter('TolFun',1e-10,@isnumeric);
            inputs.addParameter('MaxFunEvals',1e3,@isnumeric);
            inputs.addParameter('InitDamping',1,@isnumeric);
            inputs.addParameter('display','iter',@ischar);
            inputs.addParameter('weighting',false,@islogical);
            inputs.addParameter('flagModel','PSF-R',@ischar);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,true],@islogical);
            inputs.addParameter('nZernMode',[],@isnumeric);
            inputs.addParameter('flagNoStellar',false,@islogical);
            inputs.addParameter('flagJacobian',false,@islogical);
            inputs.parse(im_sky,varargin{:});
            
            % Parse inputs
            obj.im_sky   = im_sky;
            MaxIter      = inputs.Results.MaxIter;
            TolX         = inputs.Results.TolX;
            TolFun       = inputs.Results.TolFun;
            MaxFunEvals  = inputs.Results.MaxFunEvals;
            InitDamping  = inputs.Results.InitDamping;
            display      = inputs.Results.display;
            weighting    = inputs.Results.weighting;
            flagJacobian = inputs.Results.flagJacobian;
            
            %1\ Check the frame dimensions
            obj.nStars  = numel(obj.src);
            obj.fov_sky  = size(obj.im_sky,1);
            
            %2\ Normalize the observation and define the weight matrix
            % Normalization
            obj.normFactor = sum(obj.im_sky(:));
            obj.ydata = obj.im_sky/obj.normFactor;
            [~,~,obj.ron] = tools.getFlux(obj.im_sky);
            % Weighting matrix
            obj.weightMap         = 1./sqrt((max(obj.im_sky,0)+obj.ron^2));
            
            if weighting
                obj.weightMap = obj.weightMap.*(obj.ydata>0);
            else
                obj.weightMap = ones(size(obj.ydata)).*(obj.ydata>0);
            end
            obj.ydata = obj.ydata.*obj.weightMap;
            
            %3\ Initial guess and bounds
            aoinit_   = inputs.Results.aoinit;
            aobounds_ = inputs.Results.aobounds;
            
            obj = obj.fittingSetup('aoinit',aoinit_,'aobounds',aobounds_,'fitR0',inputs.Results.fitR0,...
                'fitCn2',inputs.Results.fitCn2,'fitGains',inputs.Results.fitGains,...
                'nZernMode',inputs.Results.nZernMode);
            
            %4\ Define options
            obj.fitOption = optimoptions(@lsqcurvefit,'MaxIter',MaxIter,...
                'TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,...
                'InitDamping',InitDamping,'Display',display,'SpecifyObjectiveGradient',flagJacobian);
            
            %5\ Model definition
            obj.modelFUN = @(x,xdata) obj.multiplePSFModel(x,xdata);
            
            %6\ Calculating the model image using initial guess
            tic;
            wMap         = obj.weightMap;
            obj.weightMap= 1;
            obj.im_init  = obj.modelFUN(obj.x_init,obj.xdata);
            obj.weightMap= wMap;
            obj.t_mod    = toc();
            
            %7\ Non-linear least-squares minimization
            tic;
            [beta,~,fRes,~,~,~,J] = lsqcurvefit(@obj.multiplePSFModel,obj.x_init,obj.xdata,...
                obj.ydata,obj.lbounds,obj.ubounds,obj.fitOption);
            obj.t_fit = toc();
            
            obj.beta_ = beta;
            obj.J_    = J;
            obj.fRes_ = fRes;
            
            %8\ Unpacking results + uncertainties
            obj = obj.updateResults(beta,fRes,J);
            
            %9\ Get PSF statistics
            obj = obj.getPSFstatistics();
            
            fprintf('-----------------------------------------\n');
            fprintf('Time for PSFR instantiation :\t %.3g s\n',obj.t_inst);
            fprintf('Time for computing one PSF :\t %.3g s\n',obj.t_mod);
            fprintf('Time for PSF best-fitting :\t %.3g s\n',obj.t_fit);
            fprintf('Time for the whole process :\t %.3g s\n',obj.t_inst+obj.t_fit+obj.t_mod);
            fprintf('-----------------------------------------\n');
            
        end
        
        %% %%%%%%%%%%%%%%%%%%%%
        %                     UPDATING                        %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = updateResults(obj,beta,fRes,J)
            inputs = inputParser;
            inputs.addRequired('beta',@isnumeric);
            inputs.addRequired('fRes',@isnumeric);
            inputs.addRequired('J',@isnumeric);
            inputs.parse(beta,fRes,J);
            
            %1\ Compensate photometry and background for the normalization factor
            beta(1:obj.nStars) = beta(1:obj.nStars).*obj.normFactor;
            beta(end) = beta(end).*obj.normFactor;
            
            %2\ Get precision of estimates
            dx = diff(nlparci(real(beta),fRes,'jacobian',J),1,2);
            for i=1:length(dx)
                dx(i) = diff(nlparci(real(beta(i)),fRes,'jacobian',J(:,i)),1,2);
            end
            dx = reshape(dx,1,[]);
            obj.x_final = beta;
            obj.x_prec  = dx;
            
            %3\ Unpacking AO parameters
            wvl_ratio =  obj.wvl/500e-9;
            d = obj.nStars*3;
            % Seeing estimation
            idx         = obj.idxR0 +d;
            obj.r0_fit  = sum(beta(idx))^(-3/5);
            obj.r0_prec = 3/5*sum(beta(idx))^(-8/5)*sqrt(sum(dx(idx).^2));
            % Cn2 estimation
            idx         = obj.idxCn2 +d;
            obj.Cn2_fit = beta(idx);
            obj.Cn2_prec= dx(idx);
            if isempty(obj.Cn2_fit)
                obj.Cn2_fit = obj.r0_fit^(-5/3);
                obj.Cn2_prec= 5/3*obj.r0_prec*obj.r0_fit^(-8/3);
            end
            % Retrieved atmosphere
            if numel(obj.Cn2_fit) == 1
                hl = 0;
            else
                hl          = [obj.atm.layer.altitude];
            end
            
            nL  = length(obj.Cn2_fit*wvl_ratio^2);
            r0  = sum(obj.Cn2_fit*wvl_ratio^2)^(-3/5);
            fl  = obj.Cn2_fit/sum(obj.Cn2_fit);
            wS  = 10*ones(1,nL);
            wD  = zeros(1,nL);
            obj.atm_fit = atmosphere(photometry.V0,r0,obj.atm.L0,'fractionnalR0',fl,...
                'altitude',hl,'layeredL0',obj.atm.L0,'windSpeed',wS,'windDirection',wD);
            
            % AO parameters
            obj.xao_fit = beta([obj.idxDho,obj.idxDtt,obj.idxDal]+3*obj.nStars);
            obj.xao_prec= dx([obj.idxDho,obj.idxDtt,obj.idxDal]+3*obj.nStars);
            
            %4\ Unpacking Stellar parameters
            obj.bg_fit           = beta(end);
            obj.xstars_fit       = beta(1:d);
            obj.xstars_prec      = dx(1:d);
            obj.catalog_fit.id   = 1:obj.nStars;
            obj.catalog_fit.RA   = obj.xstars_fit(obj.nStars+1:obj.nStars*2)*obj.psInMas;
            obj.catalog_fit.DEC  = obj.xstars_fit(2*obj.nStars+1:obj.nStars*3)*obj.psInMas;
            obj.catalog_fit.FLUX = obj.xstars_fit(1:obj.nStars);
            obj.catalog_fit.dRA  = obj.xstars_prec(obj.nStars+1:obj.nStars*2)*obj.psInMas;
            obj.catalog_fit.dDEC = obj.xstars_prec(2*obj.nStars+1:obj.nStars*3)*obj.psInMas;
            obj.catalog_fit.dFLUX= obj.xstars_prec(1:obj.nStars);
            obj.catalog_fit.dMAG = 2.5*obj.catalog_fit.dFLUX./obj.catalog_fit.FLUX/log(10);
            
            %5\ Grab the best-fitted PSF + 3-sigma PSF
            wMap        = obj.weightMap;
            obj.weightMap  = 1;
            if all(obj.x_prec <= max(abs(obj.lbounds),abs(obj.ubounds)))
                im_3sig_m   = obj.modelFUN(obj.x_final-obj.x_prec,obj.xdata);
                psf_3sig_m  = obj.psf_fit;
                im_3sig_p   = obj.modelFUN(obj.x_final+obj.x_prec,obj.xdata);
                psf_3sig_p  = obj.psf_fit;
                obj.im_3sig = [im_3sig_m,im_3sig_p];
                obj.psf_3sig= [psf_3sig_m,psf_3sig_p];
            end
            obj.im_fit  = obj.modelFUN(obj.x_final,obj.xdata); %psf_fit is updated inside FUN
            obj.weightMap= wMap;
            
        end
        
        function obj = getPSFstatistics(obj)
            
            %6\ PSF statistics
            im           = obj.im_sky;
            dx           = obj.x_prec;
            obj.eqm_fit  = tools.getFVU(im,obj.im_fit);
            obj.eqm_init = tools.getFVU(im,obj.im_init);
            if ~isempty(obj.im_3sig)
                eqm_p3s      = tools.getFVU(im,obj.im_3sig(:,1:end/2,:)- dx(end));
                eqm_m3s      = tools.getFVU(im,obj.im_3sig(:,end/2+1:end,:)+ dx(end));
                obj.eqm_3sig = sqrt((0.5*(eqm_p3s + eqm_m3s)).^2 - obj.eqm_fit.^2);
            end
            
            % Image statistics
            [obj.SR_sky,obj.dSR_sky]    = tools.getStrehl(im,obj.tel.pupil,obj.Samp);
            [obj.SR_init,obj.dSR_init]  = tools.getStrehl(obj.im_init,obj.tel.pupil,obj.Samp);
            [obj.SR_fit,obj.dSR_fit]    = tools.getStrehl(obj.im_fit,obj.tel.pupil,obj.Samp);
            [obj.FWHM_sky,obj.dFWHM_sky]= tools.getFWHM(im,obj.psInMas,8);
            [obj.FWHM_init,obj.dFWHM_init]= tools.getFWHM(obj.im_init,obj.psInMas,8);
            [obj.FWHM_fit,obj.dFWHM_fit]= tools.getFWHM(obj.im_fit,obj.psInMas,8);
        end
        
        function displayResults(obj,im,fov)
            
            if nargin <3
                fov = size(im,1);
            end
            %1\ 1D plot of PSF
            fov_im     = size(im,1);
            if fov < fov_im
                im = tools.crop(im,fov);
                imfit = tools.crop(obj.im_fit,fov);
            else
                fov = fov_im;
                imfit = obj.im_fit;
            end
            im(im<0)=0;
            n       = floor(fov/2);
            x       = linspace(0,obj.psInMas*fov/2,n);
            xfull   = obj.psInMas*fov/2*linspace(-1,1,fov);
            
            h = figure;
            subplot(1,3,1)
            semilogy(x,radial(im,n,n),'b-');hold on;
            semilogy(x,radial(imfit,n,n),'r--');
            semilogy(x,abs(radial(imfit-im,n,n)),'k:');
            legend({'Sky image','Best-fitted model','Residual'},'interpreter','latex','FontSize',16,'Location','southwest');
            set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
            tools.makeAxisSquared(h);
            subplot(1,3,2)
            semilogy(xfull,im(:,n),'b-');hold on;
            semilogy(xfull,imfit(:,n),'r--');
            semilogy(xfull,abs(imfit(:,n)-im(:,n)),'k:');
            legend({'Sky image','Best-fitted model','Residual'},'interpreter','latex','FontSize',16,'Location','southwest');
            set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
            tools.makeAxisSquared(h);
            subplot(1,3,3)
            semilogy(xfull,im(n,:),'b-');hold on;
            semilogy(xfull,imfit(n,:),'r--');
            semilogy(xfull,abs(imfit(n,:)-im(n,:)),'k:');
            legend({'Sky image','Best-fitted model','Residual'},'interpreter','latex','FontSize',16,'Location','southwest');
            set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
            tools.makeAxisSquared(h);
            
            %2\ 1D plot of OTF
            otf_sky = abs(fftshift(fft2(im)));
            otf_sky = otf_sky/max(otf_sky(:));
            otf_fit = abs(fftshift(fft2(imfit)));
            otf_fit = otf_fit/max(otf_fit(:));
            
            u       = linspace(0,max(obj.Samp,1),fov/2);
            h = figure;
            semilogy(u,(radial(otf_sky)),'b--');hold on;
            semilogy(u,(radial(otf_fit)),'r-');
            semilogy(u,abs(radial(otf_fit-otf_sky)),'k:');
            xlim([0,1]);
            legend({'Sky OTF','Best-fitted model','Residual'},'interpreter','latex','FontSize',18);
            set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
            tools.makeAxisSquared(h);
            
            %3\ 2D PSF;
            P = [im,imfit,imfit-im];
            Plog = log10(abs(P));
            if sign(max(Plog(:))) ==1
                mn = 0.2*max(Plog(:));
                mx = 0.9*max(Plog(:));
            else
                mn = 5*max(Plog(:));
                mx = 1.1*min(Plog(:));
            end
            
            h = figure;
            imagesc(Plog,[mn,mx]);
            set(gca,'XTick',[],'YTick',[]);
        end
    end
end




