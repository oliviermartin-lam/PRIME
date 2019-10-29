classdef massProfiler < handle    
    properties (SetObservable=true)
        file;
        year;
        month;
        day;
        hour; %% UT, 24 hour based (int)
        minute;
        second;
        timeInHours; %% UT, 24 hour based, float
        r0;
        seeing;
        wavelength;
        profs;
        altitude;
    end    
    methods
        function obj = massProfiler(massfile,varargin)
            inputs = inputParser;
            inputs.addRequired('dimmfile', @ischar);
            inputs.addParameter('wavelength',500e-9, @isnumeric);
            inputs.parse(massfile,varargin{:});
            
            obj.file        = massfile;
            obj.wavelength = inputs.Results.wavelength;
                        
            %['year', 'month', 'day', 'hour', 'minute', 'second', 
            %'cn2dh_05', 'cn2dh_1', 'cn2dh_2', 'cn2dh_4', 'cn2dh_8', 
            %'cn2dh_16', 'seeing'])
            
            T = csvread(massfile);            
            % Date and time are in HST
            obj.year        = T(:,1);
            obj.month       = T(:,2);
            obj.day         = T(:,3);
            obj.hour        = T(:,4);
            obj.minute      = T(:,5);
            obj.second      = T(:,6);            
            obj.profs       = T(:,7:end-1);
            obj.seeing      = T(:,end);
            obj.altitude    = [0.5,1,2,4,8,16]*1e3;
            % Convert from Hawaii-Aleutian Standard Time (HST) to UT
            obj.hour        = obj.hour + 10;
            idx             = find(obj.hour >= 24);
            obj.day(idx(1)) = obj.day(idx(1)) + 1;
            obj.hour(idx(1))= obj.hour(idx(1)) - 24;
            obj.r0          = 1e-2*0.98*obj.wavelength*constants.radian2arcsec./obj.seeing; % in m
            obj.timeInHours = obj.hour+(obj.minute/60.0)+(obj.second/3600.0);            
        end           
    end
end