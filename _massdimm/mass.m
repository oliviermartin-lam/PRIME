classdef mass < handle    
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
        free_seeing;
        wavelength;
    end    
    methods
        function obj = mass(massfile,varargin)
            inputs = inputParser;
            inputs.addRequired('dimmfile', @ischar);
            inputs.addParameter('wavelength',500e-9, @isnumeric);
            inputs.parse(massfile,varargin{:});
            
            obj.file        = massfile;
            obj.wavelength = inputs.Results.wavelength;
            %('year','month','day','hour','minute','second','seeing')%
            T = csvread(massfile);
            % Date and time are in HST
            obj.year        = T(:,1);
            obj.month       = T(:,2);
            obj.day         = T(:,3);
            obj.hour        = T(:,4);
            obj.minute      = T(:,5);
            obj.second      = T(:,6);
            obj.free_seeing = T(:,7);
            % Convert from HST to UT
            obj.hour        = obj.hour + 10;
            idx             = find(obj.hour >= 24);
            obj.day(idx(1)) = obj.day(idx(1)) + 1;
            obj.hour(idx(1))= obj.hour(idx(1)) - 24;
            obj.r0          = 1e-2*0.98*obj.wavelength*constants.radian2arcsec./obj.free_seeing; % in m
            obj.timeInHours = obj.hour+(obj.minute/60.0)+(obj.second/3600.0);
        end       
    end
end