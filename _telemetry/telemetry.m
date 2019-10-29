classdef telemetry < handle
    
    properties (SetObservable=true)
        
        % ---------------------------- PATHS ---------------------------- %
        aoSys;
        parm;
        date;
        path_data;
        path_im;
        path_calib;
        path_profile;
        % -------------------------- TELEMETRY -------------------------- %
        slopes;
        nSl;                     % Number of controlled slopes
        nExp;                    % Number of frames
        tipTilt;
        hodm_pos;
        nCom;                    % Number of controlled commands
        ittm_pos;
        focus;
        MI;
        R;
        Rtt;
        Rtomo;
        waveFront;
        validActu;
        % ------------------------ SCIENCE IMAGE ------------------------ %
        im_sky;
        photometry;
        resolution;
        airmass=1;
        wvl;
        psInMas;
        Samp;
        % ------------------------- LOOP STATUS ------------------------- %
        holoop_gain
        holoop_freq ;
        holoop_lat;
        holoop_rtf;
        holoop_ntf;
        holoop_pn;
        ttloop_gain;
        ttloop_freq;
        ttloop_lat;
        ttloop_rtf;
        ttloop_ntf;
        ttloop_pn
        % ------------------------- CALIBRATION ------------------------- %
        static_map;
        dmIF_hr;
        dmIF;
        dmIF_inv;
        dmIF_inv_hr;
        Hdm;
        psf_ncpa;
        otf_ncpa;
        % --------------------------- EXTERNAL -------------------------- %
        Cn2_ext; %fractional weight + altitude at zenith
        r0_ext;  % r0 @500nm and at zenith
    end
    
    methods
        
        function obj = telemetry(aoSys,date,path_data,path_im,path_calib,varargin)
            inputs = inputParser;
            inputs.addRequired('aoSys', @(x) isa(x,'struct') | isa(x,'aoSystem'));
            inputs.addRequired('date', @ischar);
            inputs.addRequired('path_data', @ischar);
            inputs.addRequired('path_im', @ischar);
            inputs.addRequired('path_calib', @ischar);
            inputs.addParameter('path_profile', [],@ischar);
            inputs.parse(aoSys,date,path_data,path_im,path_calib,varargin{:});
            
            %Inputs
            obj.aoSys       = aoSys;
            obj.date         = inputs.Results.date;
            obj.path_data    = inputs.Results.path_data;
            obj.path_im      = inputs.Results.path_im;
            obj.path_calib   = inputs.Results.path_calib;
            obj.path_profile = inputs.Results.path_profile;
            
            
            if isstruct(aoSys)
                obj.parm = aoSys;                
                obj = obj.restoreKeckTelemetry(obj.date,obj.path_data,obj.path_calib);
                obj = obj.restoreNIRC2Image(obj.date,obj.path_im,obj.path_calib);
                
                if ~isempty(obj.path_profile)
                    obj = obj.restoreMassDimmMaunaKea(obj.date,obj.path_im,obj.path_profile);
                end
            else
                obj.parm = aoSys.parm;
                obj = restoreSimulationTelemetry();
            end
        end
        
        function obj = restoreKeckTelemetry(obj,date,path_data,path_calib)
            
            %1\ Restore telemetry
            trs      = restore_idl('filename',path_data);
            
            %2\ Get AO control loop data
            D              = 9.98;
            volt2meter   = 0.5e-6;%4095e-6;
            obj.slopes   = double(trs.A.OFFSETCENTROID);
                        
            if ismatrix(obj.slopes)
                obj.nSl     = size(obj.slopes,1);
                obj.nExp    = size(obj.slopes,2);
            else
                obj.nSl     = size(obj.slopes,1);
                obj.nWfs    = size(obj.slopes,2);
                obj.nExp    = size(obj.slopes,3);
            end
            
            obj.hodm_pos = double(trs.A.DMCOMMAND);
            obj.hodm_pos = obj.hodm_pos*volt2meter;
            obj.nCom       = size(obj.hodm_pos,1);
            
            %3\ Get matrices
            MC           = reshape(double(trs.RX),obj.nSl ,obj.nCom+3,trs.NREC); %command matrix
            MC           = mean(MC,3);
            obj.R        = volt2meter*MC(:,1:obj.nCom)';
            obj.Rtt      = volt2meter*MC(:,obj.nCom+1:obj.nCom+2)';
            obj.focus    = 2.6*flipud(double(trs.A.RESIDUALWAVEFRONT(end,:)));
            
            if ~isfield(trs,'B')
                obj.tipTilt  = 2.6*flipud(double(trs.A.RESIDUALWAVEFRONT(obj.nCom+1:obj.nCom+2,:)));
                obj.ittm_pos = double(trs.A.TTCOMMANDS);
            else
                obj.tipTilt  = 0.26*flipud(double(trs.B.DTTCENTROIDS));
                obj.ittm_pos = 0.26*double(trs.B.DTTCOMMANDS);
            end
            
            % Reconstruction
            obj.waveFront= volt2meter*double(trs.A.RESIDUALWAVEFRONT(1:obj.nCom,:));
            obj.waveFront = bsxfun(@minus,obj.waveFront,mean(obj.waveFront,2));
            
            obj.ittm_pos= D*tan(constants.arcsec2radian*obj.ittm_pos);
            obj.ittm_pos= bsxfun(@minus,obj.ittm_pos,mean(obj.ittm_pos,2));            
            obj.tipTilt = D*tan(constants.arcsec2radian*obj.tipTilt);
            obj.tipTilt = bsxfun(@minus,obj.tipTilt,mean(obj.tipTilt,2));
            
            %4\ Get calibration data
            %4.1. Static map
            if strcmp(date,'20130801')
                obj.static_map = fitsread([path_calib,'phasemaps_2013.fits'])*1.6455e3/2/pi;
                obj.static_map = squeeze(obj.static_map(:,:,2));
            else
                pc       = [path_calib,'phase_diversity_results_',obj.date,'_average_PD1.sav'];
                cal      = restore_idl('filename',pc);
                ncpaMap1 = double(cal.PHASE);
                ncpaWvl  = double(cal.LAMBDA)*1e-6;
                pc(end-4)= '2';
                cal      = restore_idl('filename',pc);
                ncpaMap  = double(cal.PHASE)*0.5 + 0.5*ncpaMap1;
                obj.static_map = tools.rotateIm(ncpaMap,90)*ncpaWvl*1e9/2/pi;
            end
            
            %4.2 Influence DM functions
            bif      = xineticsInfluenceFunction(0.5625);
            nActu1D  = 21;
            nRes     = 2*nActu1D-1;
            dmSq     = deformableMirror(nActu1D,'modes',bif,'resolution',nRes);
            valDM    = logical(load([path_calib,'keckValidActuatorMap.txt']));
            obj.dmIF = dmSq.modes.modes;%keckTools.idlToMatlabIndex(dmSq.modes.modes,nActu1D,true(nActu1D),'modes');
            obj.dmIF_inv = pinv(full(obj.dmIF));
            obj.Hdm      = obj.dmIF*obj.dmIF_inv;
            
            obj.validActu = valDM;
            u = zeros(nActu1D^2,obj.nExp);
            u(valDM,:) = obj.waveFront;
            obj.waveFront = u;
            u = zeros(nActu1D^2,obj.nExp);
            u(valDM,:) = obj.hodm_pos;
            obj.hodm_pos = u;
            
            %5\ Get the loop status
            %5.1 Loop parameters
            obj.holoop_freq = 1/(100e-9*mean(diff(trs.A.TIMESTAMP)));
            obj.holoop_gain = double(trs.DM_SERVO(1));
            obj.holoop_lat  = 2.13e-3; %Vand Dam SPIE 2004
            obj.ttloop_freq = obj.holoop_freq*size(obj.tipTilt,2)/size(obj.slopes,2);
            obj.ttloop_gain = double(trs.DT_SERVO(1));
            obj.ttloop_lat  = 1.65e-3; %Vand Dam SPIE 2004
            
            %5.2 Rejection transfer function
            % Integrator tf
            ho_num  = (double(trs.DM_SERVO(1:4))');
            ho_den  = [1 (double(trs.DM_SERVO(5:end))')];
            h_servo = tf(ho_num,ho_den,1/obj.holoop_freq);
            
            % Lag tf
            delay  = obj.holoop_lat*obj.holoop_freq;
            delta  = (delay - floor(delay));
            if delay > 1
                h_lag = tf('z',1/obj.holoop_freq)^(-1) + delta*tf('z',1/obj.holoop_freq)^(-2);
            else
                h_lag = 1 + delta*tf('z',1/obj.holoop_freq)^(-1);
            end
            obj.nExp = size(obj.slopes,2);
            floc           = logspace(-3,log10(0.5*obj.holoop_freq),obj.nExp/2);
            obj.holoop_rtf = squeeze(bode(1/(1+h_servo*h_lag),2*pi*floc))';
            
            %Noise transfer function
            obj.holoop_ntf = squeeze(bode(h_servo*h_lag/(1+h_servo*h_lag),2*pi*floc));
            obj.holoop_pn  = (trapz(floc,abs(obj.holoop_ntf).^2)*2/obj.holoop_freq);
                        
            % TT Integrator tf
            tt_num  = (double(trs.DT_SERVO(1:4))');
            tt_den  = [1 (double(trs.DT_SERVO(5:end))')];
            ht_servo = tf(tt_num,tt_den,1/obj.ttloop_freq);
            
            % TT Lag tf
            delay  = obj.ttloop_lat*obj.ttloop_freq;
            delta  = (delay - floor(delay));
            if delay > 1
                ht_lag = tf('z',1/obj.ttloop_freq)^(-1) + delta*tf('z',1/obj.ttloop_freq)^(-2);
            else
                ht_lag = 1 + delta*tf('z',1/obj.ttloop_freq)^(-1);
            end
            nF = size(obj.tipTilt,2);
            floc           = logspace(-3,log10(0.5*obj.ttloop_freq),nF/2);
            obj.ttloop_rtf = squeeze(bode(1/(1+ht_servo*ht_lag),2*pi*floc))';
            %Noise transfer function
            obj.ttloop_ntf = squeeze(bode(ht_servo*ht_lag/(1+ht_servo*ht_lag),2*pi*floc));
            obj.ttloop_pn  = (trapz(floc,abs(obj.ttloop_ntf).^2)*2/obj.ttloop_freq);
            
        end
        
        function obj = restoreNIRC2Image(obj,date,path_im,path_calib)
            
            if strcmp(date,'20130801')
                badPixMap          = fitsread([path_calib,'supermask2013_modified.fits']);
                badPixMap(169,179) = 1;
                badPixMap(39,707)  = 1;
                flat               = fitsread([path_calib,'flat_kp2017.fits']);
            else
                badPixMap = fitsread([path_calib,'supermask2017.fits']);
                flat      = fitsread([path_calib,'flat_kp2017.fits']);
                badPixMap(289,136) = 1;
                badPixMap(290,136) = 1;
                badPixMap(289,137) = 1;
                badPixMap(290,137) = 1;
                badPixMap(694,860) = 1;
                badPixMap(578,924) = 1;
                badPixMap(720,859) = 1;
                badPixMap(411,762) = 1;
                badPixMap(449,747) = 1;
            end
            
            % Fits info
            tmp       = fitsinfo(path_im);
            psfHeader = tmp.PrimaryData.Keywords;
            list      = psfHeader(:,1);
            val       = psfHeader(:,2);
            % Wavelength and PSF sampling
            obj.wvl      = cell2mat(val(strcmp(list,'CENWAVE')))*1e-6;
            wvlMax       = cell2mat(val(strcmp(list,'MAXWAVE')))*1e-6;
            obj.photometry = photometry_kasp(obj.wvl,2*(wvlMax-obj.wvl),1e12);
            obj.Samp     = 1.5*(obj.wvl/1.6455e-6);
            obj.psInMas = 9.94;
            %obj.Samp     = constants.radian2mas*obj.wvl/8.992/2/9.94;
            obj.airmass  = cell2mat(val(strcmp(list,'AIRMASS')));
            % Post-processing
            im0     = fitsread(path_im);
            N = size(im0,1);
            if N == 1024
                bMap = badPixMap;
                fMap = flat;
            elseif N == 512 %cropped mode (not binning)
                bMap = tools.crop(badPixMap,N);
                fMap = tools.crop(flat,N);
            end
            obj.im_sky  = tools.processImage(im0,0,fMap,bMap,obj.Samp,...
                'fovInPixel',obj.parm.cam.resolution,'masking',false,'rebin',0,'tfccd',false,...
                'thresholding',-Inf);
        end
        
        function obj = restoreMassDimmMaunaKea(obj,date,path_im,path_profile)
            
            tmp       = fitsinfo(path_im);
            psfHeader = tmp.PrimaryData.Keywords;
            list      = psfHeader(:,1);
            val       = psfHeader(:,2);
            expStart  = cell2mat(val(strcmp(list,'EXPSTART')));
            hi        = str2double(expStart(1:2));
            mi        = str2double(expStart(4:5));
            si        = str2double(expStart(7:end));
            expStop   = cell2mat(val(strcmp(list,'EXPSTOP')));
            hf        = str2double(expStop(1:2));
            mf        = str2double(expStop(4:5));
            sf        = str2double(expStop(7:end));
            sysTime   = 0.5*(hi+hf + (mi+mf)/60 + (si+sf)/3600);
            
            [dimmfile,massfile,proffile] = profiler.fetchData(date,path_profile);
            % DIMM data
            if ~isempty(dimmfile) || ~isempty(massfile) || ~isempty(proffile)
                DIMM       = dimm(dimmfile);
                [iB,~]     = profiler.indexTime(DIMM.timeInHours,sysTime);
                seeingDIMM = DIMM.seeing(iB);
                % Altitude seeing
                MASS       = mass(massfile);
                [iB,~]     = profiler.indexTime(MASS.timeInHours,sysTime);
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
                obj.r0_ext  = 0.976*3600*180/pi*DIMM.wavelength/seeingDIMM;
                obj.Cn2_ext = [fl;alt];
            else
                obj.r0_ext  = 0;
                obj.Cn2_ext = [zeros(1,7);zeros(1,7)];
            end
        end
        
        function obj = restoreSimulationTelemetry(obj)
            %1\ WFS slopes
            obj.slopes  = squeeze(obj.aoSys.loopData.slopes);
            obj.nSl     = size(obj.slopes,1);
            obj.nExp    = size(obj.slopes,2);
            obj.tipTilt = squeeze(obj.aoSys.loopData.tiptilt);
            %2\DM commands
            obj.hodm_pos= squeeze(obj.aoSys.loopData.dmcom);
            obj.ittm_pos= squeeze(obj.aoSys.loopData.tiltCom);
            obj.nCom    = size(obj.hodm_pos,1);
            %3\ Reconstructors
            obj.R       = obj.aoSys.wfs2dm.M;
            obj.Rtt     = obj.aoSys.matrices.slopes2Tilt;
            %4\ Reconstructed Wavefront
            obj.waveFront = obj.R*squeeze(obj.slopes);
            %5\ Calibration
            %5.1. Static map in nm
            obj.static_map = trs.ncpa.map;
            
            %5.2 Influence DM functions in opd/volt
            obj.dmIF = obj.aoSys.dm.modes.modes;
            
%             if size(obj.dmIF,1)~=obj.dk^2
%                 fprintf('\nInterpolate the influence functions...');
%                 tmp  = zeros(obj.dk,obj.dk,obj.nCom);
%                 nPh  = sqrt(size(obj.dmIF,1));
%                 if2d = reshape(full(obj.dmIF),nPh,nPh,[]);
%                 for i=1:obj.nCom
%                     tmp(:,:,i) = tools.interpolate(if2d(:,:,i),obj.dk);
%                 end
%                 obj.dmIF = tmp/max(tmp(:));
%                 obj.dmIF = reshape(obj.dmIF,obj.dk*obj.dk,obj.nCom);
%             end
%             
%             obj.dmIF_inv = pinv(full(obj.dmIF));
%             obj.Hdm      = obj.dmIF*obj.dmIF_inv;
            
            %6\ Loop status
            %6.1. HO loop
            obj.holoop_gain = obj.aoSys.loopStatus.ho.gain;
            obj.holoop_freq = 1/obj.aoSys.tel.samplingTime/obj.aoSys.loopStatus.ho.frameRate;
            obj.holoop_lat  = obj.aoSys.loopStatus.ho.latency;
            
            %6.2 Rejection transfer function
            % Integrator tf
            ho_num  = [-obj.holoop_gain 0];
            ho_den  = [1 -1];
            h_servo = tf(ho_num,ho_den,1/obj.holoop_freq);
            
            % Lag tf
            delay  = obj.holoop_lat*obj.holoop_freq;
            delta  = (delay - floor(delay));
            h_lag  = tf('z',1/obj.holoop_freq)^(-1) + delta*tf('z',1/obj.holoop_freq)^(-2);
            nF     = size(obj.waveFront,2);
            floc           = logspace(-3,log10(0.5*obj.holoop_freq),nF/2);
            obj.holoop_rtf = squeeze(bode(1/(1+h_servo*h_lag),2*pi*floc))';
            
            %6.3 Noise transfer function
            obj.holoop_ntf = squeeze(bode(h_servo*h_lag/(1+h_servo*h_lag),2*pi*floc))';
            obj.holoop_pn  = (trapz(floc,abs(obj.holoop_ntf).^2)*2/obj.holoop_freq);
            
            
            %6.4 Tip-tilt transfer function
            if isfield(trs.loopStatus,'tt')
                obj.ttloop_gain = obj.aoSys.loopStatus.tt.gain;
                obj.ttloop_freq = 1/obj.aoSys.tel.samplingTime/obj.aoSys.loopStatus.tt.frameRate;
                obj.ttloop_lat  = trs.loopStatus.tt.latency;
                
                % TT Integrator tf
                tt_num  = [-obj.ttloop_gain 0];
                tt_den  = [1 -1];
                ht_servo = tf(tt_num,tt_den,1/obj.ttloop_freq);
                
                % TT Lag tf
                delay  = obj.ttloop_lat*obj.ttloop_freq;
                delta  = (delay - floor(delay));
                if delay > 1
                    ht_lag = tf('z',1/obj.ttloop_freq)^(-1) + delta*tf('z',1/obj.ttloop_freq)^(-2);
                else
                    ht_lag = 1 + delta*tf('z',1/obj.ttloop_freq)^(-1);
                end
                nF = size(obj.tipTilt,2);
                floc           = logspace(-3,log10(0.5*obj.ttloop_freq),nF/2);
                obj.ttloop_rtf = squeeze(bode(1/(1+ht_servo*ht_lag),2*pi*floc))';
                %Noise transfer function
                obj.ttloop_ntf = squeeze(bode(ht_servo*ht_lag/(1+ht_servo*ht_lag),2*pi*floc))';
                obj.ttloop_pn  = (trapz(floc,abs(obj.ttloop_ntf).^2)*2/obj.ttloop_freq);
            else
                obj.ttloop_gain = obj.holoop_gain;
                obj.ttloop_lat  = obj.holoop_lat;
                obj.ttloop_freq = obj.holoop_freq;
                obj.ttloop_rtf  = obj.holoop_rtf;
                obj.ttloop_ntf  = obj.holoop_ntf;
                obj.ttloop_pn   = obj.holoop_pn;
            end
            obj.airmass = 1;
        end
    end
end

