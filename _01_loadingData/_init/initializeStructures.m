function trs = initializeStructures(trs)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));

%1\ Sources
trs.src = [];
trs.src.x = 0;
trs.src.y = 0;
trs.src.height = Inf;
trs.src.nSrc = numel(trs.src.x);

trs.lgs = [];
trs.lgs.height = 90e3;
trs.lgs.x = 0;                                                     % LGS x position in the field in arcsec
trs.lgs.y = 0;                                                     % LGS y position in the field in arcsec

trs.ngs = [];
trs.ngs.x = 0;                                                     % NGS x position in the field in arcsec
trs.ngs.y = 0;                                                     % NGS y position in the field in arcsec
trs.ngs.height = Inf;
%1\ Atmosphere structure
trs.atm = [];
trs.atm.wavelength = [];                                    % Wavelength on which chromatic parameters are defined
trs.atm.r0 = [];                                                   % Fried's parameter in meter
trs.atm.L0 = [];                                                  % Outer scale in meter
trs.atm.nLayer = [];                                            % number of layers
trs.atm.weights = [];                                          % Fractional weights (sum = 1)
trs.atm.heights = [];                                           % Layers' height in meters
trs.atm.Cn2 = [];                                                % Cn2 in m^(-5/3)
trs.atm.windSpeed = [];                                     % wind speed profile in m/s
trs.atm.windDirection = [];                                 % wind direction in rad

%1\ Telescope structure
trs.tel = [];
trs.tel.Dhex   = 11.88;                                       %maximal distance inside the haxogonal pupil
trs.tel.Dcircle = 9.98;                                         % equivalent diameter with a similar area
trs.tel.Ddm    = 11.25;                                       % DM actuators pitch multiplied by number of actuators
trs.tel.cobs    = 0.2356;                                     % central obstruction ratio
trs.tel.pupil   = [];                                               %pupil shape
trs.tel.resolution = [];                                        % \# pixel within the pupil
trs.tel.pixelscale = [];                                         % pixel scale in m^-1
trs.tel.static_map = [];                                        % static map in nm - updated in restoreKeckTelemetry
trs.tel.airmass = [];                                             % Telescope airmass

%2\ HO WFS structure
trs.wfs = [];
trs.wfs.slopes = [];                                              % WFS slopes in pixel units;
trs.wfs.nSl = 608;                                               % Number of slopes measurements within the pupil (x and y)
trs.wfs.nSl_c = [];                                               % Number of controlled slopes
trs.wfs.nExp = [];
trs.wfs.theta = 0;                                               % DM pupil rotation


%3\ HO DM structure
trs.dm = [];
trs.dm.com = [];                                                 %DM commands in volts
trs.dm.volt2meter = 0.4095e-6;                         %conversion factor from volts to meter OPD
trs.dm.nActuators = 21;                                    % 1D Number of actuators                                            
trs.dm.nCom = 349;                                          % Number of total actuators within the pupil
trs.dm.pitch = 0.5625;                                       % DM actuator pitch in meter
trs.dm.validActuators = [];

%4\ TT structure
trs.tipTilt = [];
trs.tipTilt.slopes = [];                                         % WFS slopes in pixel units;
trs.tipTilt.com = [];                                             % TT DM commands
trs.tipTilt.tilt2meter   = 12.68e-6;                      %arcsec of tilt to OPD over the Keckpupil 
trs.tipTilt.x = 0;                                                 % TT NGS  x position in the field in arcsec
trs.tipTilt.y = 0;                                                 % TT NGS  y position in the field in arcsec
trs.tipTilt.nExp = [];

%5\ System matrices
trs.mat = [];
trs.mat.R = [];                                                    % Command matrix
trs.mat.Rtt = [];                                                 % Tip-tilt reconstructor
trs.mat.dmIF_hr = [];
trs.mat.dmIF = [];
trs.mat.dmIF_inv = [];
trs.mat.dmIF_inv_hr = [];
trs.mat.Hdm = [];
trs.mat.Hz = [];
trs.mat.nZernMode = [];
%6\ Wavefront reconstruction
trs.rec = [];
trs.rec.focus = [];                                                % Reconstructed focus in the actuator space
trs.rec.res = [];                                                   % TT-excluded residual wavefront

%7\ HO Loop configuration
trs.holoop = [];
trs.holoop.gain = [];                                           % HO loop gain
trs.holoop.freq = [];                                            % HO loop frequency in Hz
trs.holoop.lat= [];                                               % HO loop latency in s
trs.holoop.tf.freq= [];                                         % Temporal frequencies vector (max = trs.holoop.freq /2)
trs.holoop.tf.wfs= [];                                          % WFS transfer function
trs.holoop.tf.lag= [];                                           % Latency transfer function
trs.holoop.tf.num= [];                                         % Numerator of the servo transfer function
trs.holoop.tf.den= [];                                          % Denumerator of the servo transfer function
trs.holoop.tf.servo= [];                                        % Controller transfer function
trs.holoop.tf.ol= [];                                             % Loop-loop transfer function
trs.holoop.tf.ctf= [];                                            % closed-loop transfer function
trs.holoop.tf.rtf= [];                                            % rejection transfer function
trs.holoop.tf.ntf= [];                                           % noise transfer function
trs.holoop.tf.pn= [];                                           % Noise propagation factor

%8\ HO Loop configuration
trs.ttloop = [];
trs.ttloop.gain = [];
trs.ttloop.freq = [];
trs.ttloop.lat= [];
trs.holoop.tf.freq= [];
trs.ttloop.tf.wfs= [];
trs.ttloop.tf.lag= [];
trs.ttloop.tf.num= [];
trs.ttloop.tf.den= [];
trs.ttloop.tf.servo= [];
trs.ttloop.tf.ol= [];
trs.ttloop.tf.ctf= [];
trs.ttloop.tf.rtf= [];
trs.ttloop.tf.ntf= [];
trs.ttloop.tf.pn= [];       
        
%9\ Camera structure
trs.cam = [];
trs.cam.frame = [];                                             % on-sky image
trs.cam.expTime = [];                                         % exposure time
trs.cam.resolution = 168;                                   % field of view in pixels
trs.cam.pixelScale = 9.94;                                  %pixel scale in mas
trs.cam.fov = trs.cam.resolution*9.94/1e3;        %Field of view in arcsec
trs.cam.wavelength = [];                                    % Central wavelength
trs.cam.samp = [];                                              % PSF sampling. Samp=1 -> Shannon sampling
trs.cam.SR = [];        
trs.cam.dSR = [];        
trs.cam.FWHM = [];        
trs.cam.dFWHM = [];  
