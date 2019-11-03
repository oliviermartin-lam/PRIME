%{
------------HEADER-----------------
Objective          :: Fit a model of Zernike coefficients variance over empirical measurements
 
INPUT VARS
z_coefs_var         :: Zernike coefficients variance
validZernikeModes            :: j-index of Zernike modes considered during
the fitting
D                :: Pupil outer diameter in meter
Cn                    :: Number of 1D actuators over the full pupil x/y-axis
fitL0
initR0
initL0

OUTPUT VARS
r0            :: The estimated r0
L0            :: The estimated L0
dr0            :: The 3-sigma precision on the r0
dL0            :: The 3-sigma precision on the L0


Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function varargout = fitZernikeVariance(z_coefs_var,validZernikeModes,D,varargin)
inputs = inputParser;
inputs.addRequired('z_coefs_var',@isnumeric);
inputs.addRequired('validZernikeModes',@isnumeric);
inputs.addRequired('D',@isnumeric);
inputs.addParameter('fitL0',false,@islogical);
inputs.addParameter('initR0',0.16,@isnumeric);
inputs.addParameter('initL0',25,@isnumeric);
inputs.parse(z_coefs_var,validZernikeModes,D,varargin{:});

r0init = inputs.Results.initR0;
L0init = inputs.Results.initL0;
fitL0   = inputs.Results.fitL0;

%Fitting procedure
FUN = @(x,xdata) (zernikeVarianceModel(x,xdata));
X0  = D./[r0init,L0init];
ub  = D./[0.01,D];
lb  = D./[1,100];

if ~fitL0
    X0 = X0(1);
    ub = ub(1);
    lb = lb(1);
end

opt = optimoptions(@lsqcurvefit,'MaxIter',1e2,'TolFun',1e-12,...
    'TolX',1e-12,'MaxFunEvals',3e2);

[beta,~,R,~,~,~,J] = lsqcurvefit(FUN,X0,validZernikeModes,z_coefs_var,lb,ub,opt);

%[beta,R,J] = nlinfit(nMin:nMax,c0,FUN,X0);
varargout{1}  = abs(D/beta(1));
if ~fitL0
    varargout{2}  = Inf;
else
    varargout{2}  = abs(D/beta(2));
end

if isreal(beta)
    tmp        = diff(nlparci(beta,R,'jacobian',J),1,2);
    % Outputs
    varargout{3} = D*tmp(1)/beta(1)^2;
    if ~fitL0
        varargout{4} = 0;
    else
        varargout{4} = D*tmp(2)/beta(2)^2;
    end
else
    varargout{3} = Inf;
    varargout{4} = Inf;
end

