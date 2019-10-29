classdef profiler < handle
   properties (SetObservable=true) 
       % External Profiler
       utDate;
       saveTo;
       DIMM;
       MASS;
       MASSProfiler;
       % Merged MASS/DIMM information
       wishedTime;
       closestTimeDIMM;
       closestTimeMASS;
       closestTimeMASSPROF;  
       airmass;
       r0;
       seeing;
       alt;
       seeingAlt;
       seeing0;
       cn2h;
       fl;                     
       % Internal profiler
       intProf;
       D;
       wavelength;
       nMin;
       nMax;
       noiseMethod;
       opticalGain;
       fitL0;
       best;
       r0int;
       L0int;
       dr0int;
       dL0int;   
       nBest;
       mskPhase;
   end
    
    methods
        function obj = profiler(sys,varargin)
            
            inputs = inputParser;
            inputs.addRequired('sys',@(x) isa(x,'aoSystem'));              
            inputs.addParameter('utDate','',@ischar);
            inputs.addParameter('saveTo','',@ischar);
            inputs.addParameter('wishedTime',0,@isnumeric);
            inputs.addParameter('intProf',true,@islogical);
            inputs.addParameter('D',sys.tel.D,@isnumeric);
            inputs.addParameter('wavelength',sys.atm.wavelength,@isnumeric);
            inputs.addParameter('nMin',4,@isnumeric);
            inputs.addParameter('nMax',100,@isnumeric);
            inputs.addParameter('noiseMethod',0,@isnumeric);
            inputs.addParameter('fitL0',true,@islogical);
            inputs.addParameter('best',false,@islogical);
            inputs.addParameter('mskPhase',1:size(sys.coefsHO,1),@isnumeric);
            inputs.addParameter('airmass',1,@isnumeric);
            inputs.addParameter('opticalGain',1,@isnumeric);
            inputs.parse(sys,varargin{:});
            
            obj.utDate      = inputs.Results.utDate;
            obj.saveTo      = inputs.Results.saveTo;
            obj.intProf     = inputs.Results.intProf;
            obj.wishedTime  = inputs.Results.wishedTime;
            obj.D           = inputs.Results.D;
            obj.wavelength  = inputs.Results.wavelength;
            obj.nMin        = inputs.Results.nMin;
            obj.nMax        = inputs.Results.nMax;
            obj.noiseMethod = inputs.Results.noiseMethod;
            obj.fitL0       = inputs.Results.fitL0;
            obj.best        = inputs.Results.best;
            obj.mskPhase    = inputs.Results.mskPhase;
            obj.airmass     = inputs.Results.airmass;
            obj.opticalGain = inputs.Results.opticalGain;
            % ------------ Instanciating external profiler
            if ~strcmp(obj.saveTo,'') && ~strcmp(obj.utDate,'')                
                [dimmfile,massfile,proffile] = obj.fetchData(obj.utDate,obj.saveTo);
                
                if ~strcmp(dimmfile,'')
                    obj.DIMM                 = dimm(dimmfile);
                    [iB,obj.closestTimeDIMM] = obj.indexTime(obj.DIMM.timeInHours,obj.wishedTime);                
                    obj.seeing               = obj.airmass^(3/5)*obj.DIMM.seeing(iB);   
                    obj.r0                   = 0.976*3600*180/pi*obj.DIMM.wavelength/obj.seeing;
                end
                if ~strcmp(massfile,'')
                    obj.MASS                 = mass(massfile);
                    [iB,obj.closestTimeMASS] = obj.indexTime(obj.MASS.timeInHours,obj.wishedTime);
                    obj.seeingAlt            = obj.airmass^(3/5)*obj.MASS.free_seeing(iB); 
                end
                if ~strcmp(proffile,'')
                    obj.MASSProfiler             = massProfiler(proffile);
                    [iB,obj.closestTimeMASSPROF] = obj.indexTime(obj.MASSProfiler.timeInHours,obj.wishedTime);
                    % Definig the fractionnal r0
                    obj.alt     = [0 obj.MASSProfiler.altitude];                     
                    obj.seeing0 = (obj.seeing^(5/3) - obj.seeingAlt^(5/3))^(3/5)/obj.airmass^(3/5);
                    obj.seeing  = (obj.seeing0^(5/3) + obj.seeingAlt^(5/3))^(3/5);
                    mu0         = (0.976*constants.radian2arcsec)^(5/3)*0.423*4*pi*pi/obj.MASSProfiler.wavelength^(1/3);
                    cn20        = obj.seeing0^(5/3)/mu0;
                    obj.cn2h    = [cn20 obj.MASSProfiler.profs(iB,:)];
                    obj.fl      = obj.cn2h/sum(obj.cn2h(:));        
                end                        
            end
            % ------------ Estimating the r0/L0 from DM commands
            if obj.intProf
                u       = sys.coefsHO*obj.opticalGain;
                dmModes = sys.dm.modes.modes;
                [obj.r0int,obj.L0int,obj.dr0int,obj.dL0int,obj.nBest,sys.Cnn] = ...
                    psfrTools.getr0L0fromDMcommands(dmModes,u,sys.tel,obj.wavelength,...
                    obj.nMin,obj.nMax,obj.fitL0,obj.best,obj.mskPhase);
            end
        end
    end
    
    methods (Static)
    
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---------------- MASS/DIMM ------------------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function [dimmLoc,massLoc,profLoc] = fetchData(utDate, saveTo)
      %{
        Saves massdimm files to directory specified. The output files
        will have the format:

        <utDate>.mass.dat
        <utDate>.dimm.dat
        <utDate>.massprof.dat

        Parameters
        ----------
        utDate : str
        The string UT date yyyymmdd in the format such as 20170113 for 2017-01-13        
        saveTo : str
        The directory where the retrieved MASS DIMM profiles will be stored.
     %}
           
           % Create the saving folder
           fprintf('Saving MASS/DIMM data to directory %s\n',saveTo);
           saveDimm = strjoin({saveTo,'DIMM/'},'');
           saveMass = strjoin({saveTo,'MASS/'},'');
           saveProf = strjoin({saveTo,'MASSPROF/'},'');
           if ~isdir(saveTo)
               mkdir(saveTo);
           end
           if ~isdir(saveDimm)
               mkdir(saveDimm);
           end
           if ~isdir(saveMass)
               mkdir(saveMass);
           end
           if ~isdir(saveProf)
               mkdir(saveProf);
           end
           
           urlRoot = 'http://mkwc.ifa.hawaii.edu/current/seeing/';
           
           % ------- Save the DIMM file
           dimmFile     = strjoin({utDate,'.dimm.dat'},'');
           dimmLoc      = strjoin({saveDimm,dimmFile},'');
           if exist(dimmLoc,'file')
               fprintf('DIMM file already existing in your directory\n');
           else
               url          = strjoin({urlRoot,'dimm/',dimmFile},'');
               [tmp,status] = urlread(url);               
               if status
                   dlmwrite(dimmLoc,tmp,'');
                   fprintf('Write the DIMM file in %s\n',dimmLoc);
               else
                   fprintf('No DIMM data avalaible at that date !\n');
                   dimmLoc = '';
               end
           end
           
           % ----- Save the MASS file
           massFile     = strjoin({utDate,'.mass.dat'},'');
           massLoc      = strjoin({saveMass,massFile},'');
           if exist(massLoc,'file')
               fprintf('MASS file already existing in your directory\n');
           else
               url          = strjoin({urlRoot,'mass/',massFile},'');
               [tmp,status] = urlread(url);               
               if status
                   dlmwrite(massLoc,tmp,'');
                   fprintf('Write the MASS file in %s\n',massLoc);
               else
                   fprintf('No MASS data avalaible at that date !\n');
                   massLoc = '';
               end
           end
           
           % ------- Save the MASS profile
           massproFile = strjoin({utDate,'.masspro.dat'},'');
           profLoc     = strjoin({saveProf,massproFile},'');
           if exist(profLoc,'file')
               fprintf('MASSPROFILE file already existing in your directory\n');
           else
               url         = strjoin({urlRoot,'masspro/',massproFile},'');
               [tmp,status]= urlread(url);
              
               if status 
                   dlmwrite(profLoc,tmp,'');
                   fprintf('Write the MASSPROFILE file in %s\n',profLoc);
               else
                   fprintf('No MASSPROFILE data avalaible at that date !\n');
                   profLoc = '';
               end
           end
       end
       function [closestIndex,closestDate] =  indexTime(timeList, t0)
           %{
                Fetch the closest row of data for a specified time (UT) in
                hour.
           %}
           delta = abs(timeList - t0);
           [~, closestIndex] = min(delta);
           if delta(closestIndex) > 1.0
               fprintf('Could not find data close to %d\n', t0);
               closestDate = 0;
           else
               closestDate = timeList(closestIndex);
           end           
       end            
       %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%---------------- INTERNAL PROFILER ------------------------
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       
    end
    
    
    
end