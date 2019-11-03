classdef telemetry < handle
    
    properties (SetObservable=true)
        
        % ---------------------------- PATHS ---------------------------- %
        date;
        path_data;
        path_im;
        path_calib;
        path_profile;
        fitsHdr;
        % ------------------- SYSTEM/OBSERVING CONFIGURATION -------------------- %
        aoMode;
        sysTime;
        % Structures (not OOMAO class)
        src;                                                            % structure   containing science source info
        ngs;                                                            % structure  containing NGS info
        lgs;                                                              % structure  containing LGS info
        tel;                                                              % structure containing telescope info
        atm;                                                            % structure containing the atmosphere parameters
        wfs;                                                            %structure containing HO WFS info
        tipTilt;                                                        %structure containing TT info
        dm;                                                             %structure containing HO DM info
        cam;                                                           % structure containing Science camera info
        rec;                                                            % structure containing reconstructed wavefront
        mat;                                                           % structure containing system matrices
        holoop;                                                       % structure containing HO loop configuration
        ttloop;                                                        % structure containing TT loop configuration              
    end
    
    methods
        
        function obj = telemetry(date,path_data,path_im,path_calib,varargin)
            inputs = inputParser;           
            inputs.addRequired('date', @ischar);
            inputs.addRequired('path_data', @ischar);
            inputs.addRequired('path_im', @ischar);
            inputs.addRequired('path_calib', @ischar);
            inputs.addParameter('path_profile', [],@ischar);
            inputs.parse(date,path_data,path_im,path_calib,varargin{:});
            
            %1\ Parsing inputs
            obj.date         = inputs.Results.date;
            obj.path_data    = inputs.Results.path_data;
            obj.path_im      = inputs.Results.path_im;
            obj.path_calib   = inputs.Results.path_calib;
            obj.path_profile = inputs.Results.path_profile;
            
            % 2\ Retrieve the fits header
            tmp = fitsinfo(obj.path_im);
            obj.fitsHdr =tmp.PrimaryData.Keywords;
            
            %2\ Initialize system structures
            obj = initializeStructures(obj);
            
            %3\ Restoring AO telemetry            
            obj = restoreKeckTelemetry(obj);
                       
            %4\ Model temporal transfer functions
            obj = modelTransferFunctions(obj);
            
            %5\ Restoring and processing NIRC2 images
            obj = restoreNIRC2Image(obj);
            
            %6\ Getting MASS/DIMM data         
            obj =restoreMassDimmMaunaKea(obj);           
            
        end                        
    end
end

