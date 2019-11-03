function trs = restoreKeckTelemetry(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));

%% 1\ Getting data ID
date = trs.date;
path_data = trs.path_data;
path_calib = trs. path_calib;
path_im = trs.path_im;

%% 2\ Restore telemetry data
trsData = restore_idl('filename',path_data);

%% 3\ Restore fits header 
list  = trs.fitsHdr(:,1);
val   = trs.fitsHdr(:,2);
trs.aoMode = 'NGS';
if isfield(trsData,'B')
    trs.aoMode = 'LGS';
end

%% 4\ Get AO control loop data                 
%4.1. Get slopes in pixels unit
trs.wfs.slopes   = double(trsData.A.OFFSETCENTROID);
trs.wfs.nSl     = size(trs.wfs.slopes,1);
trs.wfs.nExp    = size(trs.wfs.slopes,2);

%4.2. Get DMs commands in OPD units
trs.dm.com = double(trsData.A.DMCOMMAND)*trs.dm.volt2meter;
trs.dm.nCom       = size(trs.dm.com,1);

%4.3. Get tip-tilt measurements and conversion into OPD
if ~isfield(trsData,'B')
    trs.tipTilt.slopes  = flipud(double(trsData.A.RESIDUALWAVEFRONT(trs.dm.nCom+1:trs.dm.nCom+2,:))); %angle in arcsec
    trs.tipTilt.com = double(trsData.A.TTCOMMANDS);
else
    trs.tipTilt.slopes = flipud(double(trsData.B.DTTCENTROIDS));
    trs.tipTilt.com = double(trsData.B.DTTCOMMANDS);
end
trs.tipTilt.slopes= trs.tipTilt.tilt2meter*trs.tipTilt.slopes;
trs.tipTilt.slopes= bsxfun(@minus,trs.tipTilt.slopes,mean(trs.tipTilt.slopes,2));
trs.tipTilt.com= trs.tipTilt.tilt2meter*trs.tipTilt.com;
trs.tipTilt.com = bsxfun(@minus,trs.tipTilt.com,mean(trs.tipTilt.com,2));
trs.tipTilt.nExp = size(trs.tipTilt.slopes,2);

%% 5\ Get system matrices and reconstructed wavefront

%5.1\ Get DM commands reconstructors from slopes
MC              = reshape(double(trsData.RX),trs.wfs.nSl ,trs.dm.nCom+3,trsData.NREC); %command matrix
trs.mat.R    = trs.dm.volt2meter*MC(:,1:trs.dm.nCom,:)';
trs.mat.Rtt = trs.dm.volt2meter*MC(:,trs.dm.nCom+1:trs.dm.nCom+2,:)';

%5.2 Influence DM functions
bif      = xineticsInfluenceFunction(trs.dm.pitch); % to be replaced
dmSq     = deformableMirror(trs.dm.nActuators,'modes',bif,'resolution',2*trs.dm.nActuators-1); % to be replaced
trs.mat.dmIF       = dmSq.modes.modes;
trs.mat.dmIF_inv = pinv(full(trs.mat.dmIF));
trs.mat.Hdm        = trs.mat.dmIF*trs.mat.dmIF_inv;

%5.3\ Get the reconstructed wavefront in OPD and in the actuators space
trs.rec.res    = trs.dm.volt2meter*double(trsData.A.RESIDUALWAVEFRONT(1:trs.dm.nCom,:));
trs.rec.res    = bsxfun(@minus,trs.rec.res,mean(trs.rec.res,2));
trs.rec.focus = double(trsData.A.RESIDUALWAVEFRONT(end,:));
trs.rec.focus = bsxfun(@minus,trs.rec.focus,mean(trs.rec.focus,2));

% fill vector to get 21x21 actuators
trs.dm.validActuators     = logical(load([path_calib,'keckValidActuatorMap.txt']));
u = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:) = trs.rec.res;
trs.rec.res = u;
u = zeros(trs.dm.nActuators^2,trs.wfs.nExp);
u(trs.dm.validActuators,:) = trs.dm.com;
trs.dm.com = u;

%% 5\ Get instrumental features

%5.1. Static aberrations from phase diversity calibration
if strcmp(date,'20130801')
    trs.tel.static_map = fitsread([path_calib,'phasemaps_2013.fits'])*1.6455e3/2/pi;
    trs.tel.static_map = squeeze(trs.tel.static_map(:,:,2));
else
    pc       = [path_calib,'phase_diversity_results_',trs.date,'_average_PD1.sav'];
    cal      = restore_idl('filename',pc);
    ncpaMap1 = double(cal.PHASE);
    ncpaWvl  = double(cal.LAMBDA)*1e-6;
    pc(end-4)= '2';
    cal      = restore_idl('filename',pc);
    ncpaMap  = double(cal.PHASE)*0.5 + 0.5*ncpaMap1;
    trs.tel.static_map = tools.rotateIm(ncpaMap,90)*ncpaWvl*1e9/2/pi;
end

%5.2. Pupil
trs.tel.pupil = double(logical(trs.tel.static_map));
trs.tel.resolution = size(trs.tel.pupil,1); 
trs.tel.pixelscale = trs.tel.Ddm/trs.tel.resolution;  
%High resolution pupil
%trs.tel.pupil = fitsread([trs.path_calib,'keckPupil.fits']);

%% 6\ Get the loop status and model transfer function
%6.1. Delays
wssmprg           = str2num(cell2mat(val(strcmp(list,'WSSMPRG'))));
[trs.holoop.lat,trs.ttloop.lat] = estimateLoopDelay(wssmprg,trs.aoMode);
%6.2. Frequency
trs.holoop.freq = 1/(100e-9*mean(diff(trsData.A.TIMESTAMP)));
trs.ttloop.freq =trs.holoop.freq;
if isfield(trsData,'B')
    trs.ttloop.freq =  1/(100e-9*mean(diff(trsData.B.TIMESTAMP)));   
end
%6.3. RTC controller HO loop
trs.holoop.gain = double(trsData.DM_SERVO(1));
trs.holoop.tf.num=  (double(trsData.DM_SERVO(1:4))');
trs.holoop.tf.den= (double(trsData.DM_SERVO(5:end))');
trs.ttloop.gain = double(trsData.DT_SERVO(1));
trs.ttloop.tf.num=  (double(trsData.DT_SERVO(1:4))');
trs.ttloop.tf.den= (double(trsData.DT_SERVO(5:end))');

end
