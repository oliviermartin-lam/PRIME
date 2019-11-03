%{
------------HEADER-----------------
Objective          ::  Select the valid actuators regarding the outer and
inner pupils size
INPUT VARS
uin              :: the actuators commands
D1,D2          :: Outer and inner pupil diameter in meter
nX,nY          :: (optional) Number of 1D actuators over the full pupil x/y-axis
pitch          :: (optional)  DM actuators pitch

OUTPUT VARS
uout            :: The truncated actuators command vector
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function uout = selectActuators(uin,D1,D2,varargin)
inputs = inputParser;
inputs.addRequired('wvf',@isnumeric);
inputs.addRequired('D1',@isnumeric);
inputs.addRequired('D2',@isnumeric);
inputs.addParameter('nX',21,@isnumeric);
inputs.addParameter('nY',21,@isnumeric);
inputs.addParameter('pitch',0.5625,@isnumeric);
inputs.parse(uin,D1,D2,varargin{:});

%1\ Define the actuators grid
nX = inputs.Results.nX;
nY = inputs.Results.nY;
actuPitch = inputs.Results.pitch;
pScale = (nX-1)*actuPitch/2;
%1D DM actuators position for which x1D(nX/2+1) = 0 and y1D(nY/2+1) = 0;
x1D     = (-1+mod(nX,2)/nX:2/nX:1-1/nX)*pScale;
y1D     = (-1+mod(nY,2)/nY:2/nY:1-1/nY)*pScale; 
% 2D actuators position
[x2D,y2D]  = meshgrid(x1D,y1D);
% Radius and angle
[~,r2D] = cart2pol(x2D,y2D);

%2\ Define DM actuators masks
maskDM349 = r2D <= 11.25/2;
maskDMtrun = r2D < D1/2 & r2D > D2/2;
nActu_trunc = nnz(maskDMtrun);
%3\ Define the index of the 320 actuators to be selected
uin = uin(std(uin,[],2)~=0,:);
k349 = 0;
ktrunc = 0;
idxTrunc = zeros(1,nActu_trunc);
for iAct = 1:nX
    for jAct = 1:nY
        if maskDM349(iAct,jAct)
            k349 = k349+1;
            if maskDMtrun(iAct,jAct)
                ktrunc = ktrunc + 1;
                idxTrunc(ktrunc) = k349;
            end
        end
    end
end
%4\ Select the actuators 
uout = zeros(nX^2,size(uin,2));
id = find(maskDMtrun);
for iAct = 1:nActu_trunc
    uout(id(iAct),:) = uin(idxTrunc(iAct),:);
end
   