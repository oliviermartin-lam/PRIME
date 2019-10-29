classdef photometry_kasp < handle
    properties (SetObservable=true)
        wavelength % [micron]
        bandwidth  % [micron]
        zeroPoint  % [ph/s]       
    end
    
    methods
        
        %% Constructor
        function obj = photometry_kasp(w,bw,zp)
            obj.wavelength = w;
            obj.bandwidth  = bw;
            obj.zeroPoint  = zp/368;
        end                      
        %% Get nPhoton
        function out = nPhoton(obj,magnitude)
            out = obj.zeroPoint*10^(-0.4*magnitude);
        end
        %% Get magnitude
        function out = magnitude(obj,nPhoton)
            out = -2.5*log10(nPhoton/obj.zeroPoint);
        end
        
        function out = rdivide(obj,val)
            %% WAVELENGTH DIVISION
            %
            % out = obj./val divides the wavelength by val
            
            out = obj.wavelength./val;
        end
        
        function out = mrdivide(obj,val)
            %% WAVELENGTH DIVISION
            %
            % out = obj/val divides the wavelength by val
            
            out = rdivide(obj,val);
        end
    end
end