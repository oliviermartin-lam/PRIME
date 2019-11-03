classdef psfr < handle
    
    properties
        % --------------------- Observation
        trs;    %Telemetry class
        % --------------------- Results
        res = [];
        otf = [];
        cov = [];
        sf = [];
        psf = [];    % Reconstructed 2D PSF
        flags = [];
        % --------------------- OTFs/Model   
        % Covariance/PSD       
%         nZernMode;
%         zer;
%         Hz=1;
%         rmZerFromPhase = 1;        
    end
    
    methods
        
        function obj = psfr(trs,varargin)
            inputs = inputParser;
            inputs.addRequired('trs',@(x) isa(x,'telemetry'));
            inputs.addParameter('fov',trs.cam.resolution+4,@isnumeric);
            inputs.addParameter('flagToeplitz',true,@islogical);
            inputs.addParameter('flagNoiseMethod','autocorrelation',@ischar);
            inputs.addParameter('flagAnisoMethod','oomao',@ischar);
            inputs.addParameter('flagAoPattern','circle',@ischar);
            inputs.parse(trs,varargin{:});
            
            %1\ Parsing inputs
            obj.trs = inputs.Results.trs;            
            obj.psf.fov            = inputs.Results.fov;
            obj.flags.toeplitz       = inputs.Results.flagToeplitz;
            obj.flags.noiseMethod    = inputs.Results.flagNoiseMethod;
            obj.flags.anisoMethod    = inputs.Results.flagAnisoMethod;
            obj.flags.aoPattern    = inputs.Results.flagAoPattern;
            
            tic();            
            %1\ Update the reconstruction prerequisites
            obj = updateReconstructionPrerequisites(obj);
                        
            %2\ Reconstruct the PSF
            obj = forwardPSFR(obj,'flagToeplitz',obj.flags.toeplitz,'flagNoiseMethod',obj.flags.noiseMethod,'flagAnisoMethod',obj.flags.anisoMethod);
            t_inst = toc();
                       
            %3\ Prompt messages
            fprintf('-----------------------------------------\n');
            fprintf('Time for class instantiation :\t %.3g s\n',t_inst);
            fprintf('-----------------------------------------\n');
        end
        
       
        
         end
end




