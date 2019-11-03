classdef bestFitting < handle
    
    properties
        % --------------------- Fitting process setup
        % PSF model Initial guess
        psf;
        psfr_;
        im_sky;
        nStars;
        fov_sky;
        ron;
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
        nZernMode;        
        nGainsHO;
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
    end
    
    methods
        function obj = bestFitting(psfr_,varargin)
            inputs = inputParser;
            inputs.addRequired('psfr_',@(x) isa(x,'psfr'));
            inputs.addParameter('aoinit',[],@isnumeric);
            inputs.addParameter('aobounds',[],@isnumeric);
            inputs.addParameter('MaxIter',100,@isnumeric);
            inputs.addParameter('TolX',1e-10,@isnumeric);
            inputs.addParameter('TolFun',1e-10,@isnumeric);
            inputs.addParameter('MaxFunEvals',1e3,@isnumeric);
            inputs.addParameter('InitDamping',1,@isnumeric);
            inputs.addParameter('display','iter',@ischar);
            inputs.addParameter('weighting',false,@islogical);
            inputs.addParameter('fitR0',true,@islogical);
            inputs.addParameter('fitCn2',false,@islogical);
            inputs.addParameter('fitGains',[true,true,true],@islogical);
            inputs.addParameter('nZernMode',[],@isnumeric);
            inputs.addParameter('flagJacobian',false,@islogical);
            inputs.parse(psfr_,varargin{:});
            
            % Parse inputs            
            obj.psfr_ = psfr_;
            obj.im_sky   = psfr_.trs.cam.frame;
            MaxIter      = inputs.Results.MaxIter;
            TolX         = inputs.Results.TolX;
            TolFun       = inputs.Results.TolFun;
            MaxFunEvals  = inputs.Results.MaxFunEvals;
            InitDamping  = inputs.Results.InitDamping;
            display      = inputs.Results.display;
            weighting    = inputs.Results.weighting;
            flagJacobian = inputs.Results.flagJacobian;
            
            %1\ Check the frame dimensions
            obj.nStars  = psfr_.trs.src.nSrc;
            obj.fov_sky  = size(obj.im_sky,1);
            
            %2\ Normalize the observation and define the weight matrix
            % Normalization
            obj.normFactor = sum(obj.im_sky(:));
            obj.ydata = obj.im_sky/obj.normFactor;
            [~,~,obj.ron] = tools.getFlux(obj.im_sky);
            % Weighting matrix
            obj.weightMap         = 1./sqrt((max(obj.ydata,0)+obj.ron^2));
            
            if weighting
                obj.weightMap = obj.weightMap.*(obj.ydata>0);
            else
                obj.weightMap = ones(size(obj.ydata)).*(obj.ydata>0);
            end
            obj.ydata = obj.ydata.*obj.weightMap;
            
            %3\ Initial guess and bounds
            aoinit_   = inputs.Results.aoinit;
            aobounds_ = inputs.Results.aobounds;
            
            obj = fittingSetup(obj,'aoinit',aoinit_,'aobounds',aobounds_,'fitR0',inputs.Results.fitR0,...
                'fitCn2',inputs.Results.fitCn2,'fitGains',inputs.Results.fitGains,...
                'nZernMode',inputs.Results.nZernMode);
            
            %4\ Define options
            obj.fitOption = optimoptions(@lsqcurvefit,'MaxIter',MaxIter,...
                'TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,...
                'InitDamping',InitDamping,'Display',display,'SpecifyObjectiveGradient',flagJacobian);
            
            %5\ Model definition
            obj.modelFUN = @(x,xdata) imageModel(x,xdata,obj);
            
            %6\ Calculating the model image using initial guess
            tic;
            wMap         = obj.weightMap;
            obj.weightMap= 1;
            obj.im_init  = obj.modelFUN(obj.x_init,obj.xdata);
            obj.weightMap= wMap;
            obj.t_mod    = toc();
            
            %7\ Non-linear least-squares minimization
            tic;
            [beta,~,fRes,~,~,~,J] = lsqcurvefit(obj.modelFUN,obj.x_init,obj.xdata,...
                obj.ydata,obj.lbounds,obj.ubounds,obj.fitOption);
            obj.t_fit = toc();
            
            obj.beta_ = beta;
            obj.J_    = J;
            obj.fRes_ = fRes;
            
            %8\ Unpacking results + uncertainties
            obj = updateResults(beta,fRes,J);
            
            %9\ Get PSF statistics
            %obj = obj.getPSFstatistics();
            
            %             fprintf('-----------------------------------------\n');
            %             fprintf('Time for PSFR instantiation :\t %.3g s\n',obj.t_inst);
            %             fprintf('Time for computing one PSF :\t %.3g s\n',obj.t_mod);
            %             fprintf('Time for PSF best-fitting :\t %.3g s\n',obj.t_fit);
            %             fprintf('Time for the whole process :\t %.3g s\n',obj.t_inst+obj.t_fit+obj.t_mod);
            %             fprintf('-----------------------------------------\n');
            
        end
    end
    
end
