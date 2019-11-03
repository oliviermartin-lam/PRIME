%{
------------HEADER-----------------
Objective          ::  Select the valid actuators regarding the outer and
inner pupils size
INPUT VARS
dmModes         :: the acutators commands
dmCom            :: Outer and inner pupil diameter in meter
tel                :: DM actuators pitch
Cn                    :: Number of 1D actuators over the full pupil x/y-axis
fitL0
best
initR0
D1
D2
badModesList
aoMode
nMin
nMax

OUTPUT VARS
r0            :: The estimated r0
L0            :: The estimated L0
dr0            :: The 3-sigma precision on the r0
dL0            :: The 3-sigma precision on the L0
validZernikeModes   :: The j-index vector of Zernikes modes used for fitting
z_coefs_mod            :: The corresponding fitted coefficients
z_coefs_var            :: The empirical Zernike modes vaariance
(denoising/dealiasing and removed from bad modes)

Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function varargout = getr0L0FromDMcommands(dmModes,dmCom,tel,Cn,varargin)
inputs = inputParser;
inputs.addRequired('dmModes',@isnumeric);
inputs.addRequired('dmCom',@isnumeric );
inputs.addRequired('tel',@isstruct);
inputs.addRequired('Cn',@isnumeric);
inputs.addParameter('fitL0',false,@islogical);
inputs.addParameter('best',false,@islogical);
inputs.addParameter('initR0',[],@isnumeric);
inputs.addParameter('D1',0,@isnumeric);
inputs.addParameter('D2',0,@isnumeric);
inputs.addParameter('badModesList',[],@isnumeric);
inputs.addParameter('aoMode','NGS',@ischar);
inputs.addParameter('nMin',4,@isnumeric);
inputs.addParameter('nMax',120,@isnumeric);
inputs.parse(dmModes,dmCom,tel,Cn,varargin{:});

%1\ Parsing inputs
fitL0    = inputs.Results.fitL0;
flagBest = inputs.Results.best;
initR0  = inputs.Results.initR0;
D1       = inputs.Results.D1;
D2       = inputs.Results.D2;
nMin     = inputs.Results.nMin;
nMax     = inputs.Results.nMax;

%2\ Estimating the Zernike modes rms values
%2.1 Defining the Zernike modes over the real telescope pupil
nRes   = sqrt(size(dmModes,1));
P         = tel.pupil;
validZernikeModes = nMin:nMax;
zern_tel  = zernike(validZernikeModes,'resolution',nRes,'D',D1);
z_modes = zern_tel.modes;

%2.2 Defining the Zernike modes over the truncated pupil
if D1 ~=0
    % truncated pupil definition
    x1D     = (-1+mod(nRes,2)/nRes:2/nRes:1-1/nRes)*D1/2;
    [x2D,y2D]  = meshgrid(x1D);    
    [~,r2D] = cart2pol(x2D,y2D);    
    mask     = (r2D<=D1/2) & (r2D>D2/2);
    P_trunc = P.*mask;
    % define Zernike modes over the truncated pupil
    zern_trunc   = zernike(validZernikeModes,'resolution',nRes,'pupil',P_trunc,'D',D1);    
    % Calculate the projection matrix
    zTrunc_modes= zern_trunc.modes;
    proj   = (zTrunc_modes' * zTrunc_modes)\zTrunc_modes'*z_modes;
else
    proj = 1;
    D1 = tel.Dcircle;
end

%3\ Process the telemetry to estimate the Zernike modes
%3.1 Actuators commands to Zernike
u2z      = pinv(full(z_modes))*dmModes;
z_coefs = u2z*dmCom;
%3.2. Denoising
z_coefs_var   = std(z_coefs,[],2).^2 - diag(u2z*Cn*u2z');
%3.3. Dealiasing
%polynomial model of the variance excess due to the aliasing
pp = [0.99893514401310313 7.7066823200766521e-4,...
    -1.8233364176012401e-4 1.9037601227012146e-5 -2.7620673576089771e-7];
[nz,mz] = nmOrder(validZernikeModes);
aliasExcess = 0;
for i=1:numel(pp)
    aliasExcess = aliasExcess+pp(i).*nz.^(i-1);
end
z_coefs_var = z_coefs_var./aliasExcess';

%3.4. Bad modes list
nValid = numel(validZernikeModes);
badModesList = ones(1,nValid);
if strcmp(inputs.Results.aoMode,'NGS')
    wbm = [[4:6],[87,88 101,102 74,75, 62,63, 51 41 33 25 116 117 115]];
else
    wbm = [[4:6],[7,8,11,16,17,18,22,29,30,37,38,...
        44,46,47,56,60,62,64,66,67,68,79,81,82,84,86,...
        88,90,92,93,94,95,106,110,112,114,116]];
end
badModesList(intersect(validZernikeModes,wbm)) = 0;
badModesList(mz <=1) = 0;
badModesList = logical(badModesList);
% Removing bad modes
validZernikeModes = validZernikeModes(badModesList);
z_coefs_var = z_coefs_var(badModesList);


%4\ Model-fitting the Zernike modes
if ~flagBest % Single fitting over the validZernikeModes modes    
    [r0,L0,dr0,dL0] = fitZernikeVariance(z_coefs_var,validZernikeModes,D1,'fitL0',fitL0,'initR0',initR0);
    iBest = 1;
else % Multiple fitting and keep the most precise one   
    % instantiation
    nM   = nMin:nMax/4; iN   = length(nM);
    r0   = zeros(1,iN);dr0  = zeros(1,iN);
    L0   = zeros(1,iN);dL0  = zeros(1,iN);
    % loop over
    for i=1:iN
        [r0(i),L0(i),dr0(i),dL0(i)] = fitZernikeVariance(z_coefs_var(i:end),validZernikeModes(i:end),D1,'fitL0',fitL0,'initR0',initR0);
    end
    
    % Measure the joint precision over adjusted parameters
    eps   = hypot(dr0./r0,dL0./L0);
    iBest = find(eps == min(eps));
    % Keep the most precise result
    r0    = r0(iBest);dr0   = dr0(iBest);
    L0    = L0(iBest); dL0   = dL0(iBest);       
end

%5\ Gathering outputs
varargout{1} = r0;
varargout{2} = L0;
varargout{3} = dr0;
varargout{4} = dL0;
varargout{5} = validZernikeModes(iBest:end);
varargout{6} = zernikeVarianceModel([D1/r0,D1/L0],validZernikeModes(iBest:end));
varargout{7} =z_coefs_var(iBest:end);


