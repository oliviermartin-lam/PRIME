function trs = restoreMassDimmMaunaKea(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));

%1\ Getting data ID
date = trs.date;
path_im = trs.path_im;
path_profile = trs.path_profile;

%2\ Retrieve the system time
list  = trs.fitsHdr(:,1);
val  = trs.fitsHdr(:,2);
[hi,mi,si] = string2Date(cell2mat(val(strcmp(list,'EXPSTART'))));
[hf,mf,sf] = string2Date(cell2mat(val(strcmp(list,'EXPSTOP'))));
trs.sysTime   = 0.5*(hi+hf + (mi+mf)/60 + (si+sf)/3600);

%3\ Retrieve the closest files
[dimmfile,massfile,proffile] = profiler.fetchData(date,path_profile);

%4\ Read files
% DIMM data
if ~isempty(dimmfile) || ~isempty(massfile) || ~isempty(proffile)
    DIMM       = dimm(dimmfile);
    [iB,~]     = profiler.indexTime(DIMM.timeInHours,trs.sysTime);
    seeingDIMM = DIMM.seeing(iB);
    % Altitude seeing
    MASS       = mass(massfile);
    [iB,~]     = profiler.indexTime(MASS.timeInHours,trs.sysTime);
    seeingAlt  = MASS.free_seeing(iB);
    
    if seeingAlt > seeingDIMM
        fprintf('Altitude seeing greater than total seeing!!!\n')
        fl  = [0.517 0.119 0.063 0.061 0.105 0.081 0.054];
        alt = [0 0.5 1 2 4 8 16]*1e3;
    else
        % MASS profile
        MASSProfiler = massProfiler(proffile);
        %Ground Cn2
        seeing0 = (seeingDIMM^(5/3) - seeingAlt^(5/3))^(3/5);
        mu0     = (0.976*constants.radian2arcsec)^(5/3)*0.423*4*pi*pi/MASSProfiler.wavelength^(1/3);
        cn20    = seeing0^(5/3)/mu0;
        % total Cn2
        [iB,~] = profiler.indexTime(MASSProfiler.timeInHours,sysTime);
        alt     = [0 MASSProfiler.altitude];
        cn2h    = [cn20 MASSProfiler.profs(iB,:)];
        fl      = cn2h/sum(cn2h(:));
    end
    r0 = 0.976*3600*180/pi*DIMM.wavelength/seeingDIMM;
else
    r0  = 0.16; %median Value at MaunaKea
    fl    = [0.517 0.119 0.063 0.061 0.105 0.081 0.054];
    alt  =  [0 .5 1 2 4 8 16]*1e3;
end

%5\ Update the atmosphere structure
trs.atm.wavelength = 500e-9;
trs.atm.r0 = r0;
trs.atm.L0 = 25;
trs.atm.weights = fl;
trs.atm.heights = alt;
trs.atm.Cn2 = r0^(-5/3)*fl;
trs.atm.windSpeed = [6.8 6.9 7.1 7.5 10.0 26.9 18.5];
trs.atm.windDirection = [0 pi/2 pi/4 3*pi/2 -pi/4 pi/6 -pi/4];
trs.atm.nLayer = length(fl);

function [h,m,s] = string2Date(expStr)
h        = str2double(expStr(1:2));
m        = str2double(expStr(4:5));
s        = str2double(expStr(7:end));

