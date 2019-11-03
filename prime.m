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
        flagNoiseMethod    = 'autocorrelation';  %'rtf'
        flagAnisoMethod    = 'oomao';            %'FOURIER','FLICKER'        
    end
    
    methods
        
        function obj = prime(trs,varargin)
            inputs = inputParser;
            inputs.addRequired('trs',@(x) isa(x,'telemetry'));
            inputs.addParameter('fov_fit',[],@isnumeric);
            inputs.addParameter('monochromaticBand',10e-9,@isnumeric);
            inputs.addParameter('psf_ext',[],@isnumeric);
            inputs.addParameter('flagToeplitz',true,@islogical);
            inputs.addParameter('flagNoiseMethod','autocorrelation',@ischar);
            inputs.addParameter('flagAnisoMethod','flicker',@ischar);
            inputs.parse(trs,varargin{:});
            
            % Grab inputs
            obj.trs = trs;            
            obj.flagToeplitz       = inputs.Results.flagToeplitz;
            obj.fov_fit            = inputs.Results.fov_fit;
            obj.flagNoiseMethod    = inputs.Results.flagNoiseMethod;
            obj.flagAnisoMethod    = inputs.Results.flagAnisoMethod;
            
            if isempty(obj.fov_fit)
                obj.fov_fit = trs.cam.resolution+4;
            end
            
            
            %1\ Instantiate the PSF reconstruction
            tic;
            obj = obj.forwardPSFR(obj.wvl,obj.fov_fit,trs.cam.pixelScale,...
                'flagToeplitz',obj.flagToeplitz,'flagAnisoMethod',obj.flagAnisoMethod);
             obj.t_inst = toc();
                       
            %2\ Prompt messages
            fprintf('-----------------------------------------\n');
            fprintf('Time for class instantiation :\t %.3g s\n',obj.t_inst);
            fprintf('-----------------------------------------\n');
        end
        
        %% %%%%%%%%%%%%%%%%%%%%
        %                 INSTANTIATION                    %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% %%%%%%%%%%%%%%%%%%%%
        %            PSF RECONSTRUCTION               %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% %%%%%%%%%%%%%%%%%%%%
        %               MODEL-FITTING                      %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
       
        %% %%%%%%%%%%%%%%%%%%%%
        %                     UPDATING                        %
        %%%%%%%%%%%%%%%%%%%%%%%%
        
         end
end




