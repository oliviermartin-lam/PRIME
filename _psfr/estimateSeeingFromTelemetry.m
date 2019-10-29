%{
------------HEADER-----------------
Objective          ::  Estimate the noise covariance matrix

INPUT VARS
psfr              :: The PRIME top-level class
nMin, nMax  :: (optional) 'Min and Max Zernike modes order the wavefront is decomposed onto (defaut: 4, 100)
fitL0             :: (optional) fit the L0 as well if set to true (default: true)
best              :: (optional) Repeat the fitting process for different values of nMin from nMin to nMax/4 and pick up the most precise one
wvl_ref         :: (optional) Wavelength the r0 is given at.
mskPhase     :: (optional) Phase mask to reject non-valid actuators.

OUTPUT VARS
r0               :: The estimated line-of-sight r0 in meters at wvl
L0               :: The estimated L0 in meters
fwhm          :: The estimated PSF FWHM in arcsec
dr0             :: The 3-sigma precision on the r0 estimation
dL0             :: The 3-sigma precision on the L0 estimation
dfwhm        :: The 3-sigma precision on the fwhm estimation

Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}


function [r0,L0,fwhm,dr0,dL0,dfwhm] = estimateSeeingFromTelemetry(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'prime'));
inputs.addParameter('nMin',4,@isnumeric);
inputs.addParameter('nMax',100,@isnumeric);
inputs.addParameter('fitL0',true,@islogical);
inputs.addParameter('best',true,@islogical);
inputs.addParameter('wvl_ref',0.5e-6,@isnumeric);
inputs.addParameter('mskPhase',psfr.trs.validActu,@islogical);
inputs.parse(psfr,varargin{:});

%1\ Parsing inputs
nMax = inputs.Results.nMax;
nMin = inputs.Results.nMin;
fitL0 = inputs.Results.fitL0;
wvl_ref = inputs.Results.wvl_ref;
mskPhase = inputs.Results.mskPhase;

%2\ Getting the reconstructed wavefront and noise covariance at 500 nm
u = psfr.trs.hodm_pos;
w = psfr.trs.waveFront;
dt = psfr.trs.holoop_lat*psfr.trs.holoop_freq;
dt = mod(dt,2);
wvf = (w + (dt.*circshift(u,-3,2) + (1-dt).*circshift(u,-2,2)) )*(2*pi/wvl_ref); % To be reviewed
%wvf = u*(2*pi/wvl_ref);
Cn    = psfr.Cn_ho*(2*pi/wvl_ref)^2;

%3\ Estimating the r0/outer scale
[r0,L0,dr0,dL0] =  getr0FromTelemetry(psfr.trs.dmIF_hr,wvf,psfr.tel,Cn,...
    'best',true,'fitL0',fitL0, 'nMin',nMin,'nMax',nMax,'mskPhase',mskPhase);

fprintf('r_0 estimated from the telemetry: %.3g cm\n',1e2*r0);

%4\ Estimating the atmosphere-limited PSF FWHM
k1 = 1.03*wvl_ref*constants.radian2arcsec;
fwhm = k1/r0;
dfwhm= k1*dr0/r0^2;

if L0 ~=0
    k2     = 2.813;
    k3     = 0.356;
    sqfact = sqrt(1-k2*(r0/L0)^k3);
    fwhm = fwhm*sqfact;
    
    dfdr01 = dr0/r0*fwhm;
    dfdr02 = 0.5*k1*dr0/r0/sqfact*k2*k3*L0^(-k3)*r0^(k3-1);
    dfdr0 = hypot(dfdr01,dfdr02);        
    dfdL0 = 0.5*k1*dL0/r0/sqfact*k2*k3*r0^(k3)*L0^(-k3-1);    
    dfwhm=hypot(dfdr0,dfdL0);
end


function varargout = getr0FromTelemetry(dmModes,waveFront,tel,Cn,varargin)
inputs = inputParser;
inputs.addRequired('dmModes',@isnumeric);
inputs.addRequired('waveFront',@isnumeric );
inputs.addRequired('tel',@(x) isa(x,'telescope'));
inputs.addRequired('Cn',@isnumeric);
inputs.addParameter('nMin',4,@isnumeric);
inputs.addParameter('nMax',100,@isnumeric);
inputs.addParameter('mskPhase',true(1,size(dmModes,2)),@islogical);
inputs.addParameter('fitL0',false,@islogical);
inputs.addParameter('best',false,@islogical);
inputs.parse(dmModes,waveFront,tel,Cn,varargin{:});

nMin     = inputs.Results.nMin;
nMax     = inputs.Results.nMax;
mskPhase = inputs.Results.mskPhase;
fitL0    = inputs.Results.fitL0;
best     = inputs.Results.best;

% Derive the zernike reconstructor
nRes   = sqrt(size(dmModes,1));%tel.resolution;
zern   = zernike(nMin:nMax,'resolution',nRes);
zernP  = zernike(nMin:nMax,'resolution',nRes,'pupil',tel.pupil);
% Project the modes onto the telescope pupil
zModes = zern.modes;
zpModes= zernP.modes;
proj   = (zpModes' * zpModes)\zpModes'*zModes;
% DM commands to zernike
u2z    = proj*pinv(full(zModes))*dmModes(:,mskPhase);
z      = u2z*waveFront(mskPhase,:);

%Zernike variance distribution
c0   = std(z,[],2).^2 - diag(u2z*Cn(mskPhase,mskPhase)*u2z');

if ~best
    [r0,L0,dr0,dL0] = fitZernikeVariance(tel.D,nMin,nMax,c0,fitL0);
    iBest           = 1;
    varargout{5}    = nMin;
else
    nM   = nMin:nMax/4;
    iN   = length(nM);
    r0   = zeros(1,iN);
    L0   = zeros(1,iN);
    dr0  = zeros(1,iN);
    dL0  = zeros(1,iN);
    for i=1:iN
        [r0(i),L0(i),dr0(i),dL0(i)] = fitZernikeVariance(tel.D,...
            nM(i),nMax,c0(i:end),fitL0);
    end
    eps   = hypot(dr0./r0,dL0./L0);
    iBest = find(eps == min(eps));
    r0    = r0(iBest);
    L0    = L0(iBest);
    dr0   = dr0(iBest);
    dL0   = dL0(iBest);
    varargout{5} = nM(iBest);
end

varargout{1} = r0;
varargout{2} = L0;
varargout{3} = dr0;
varargout{4} = dL0;

if nargout > 5
    varargout{6} = covNoise;
end

function varargout = fitZernikeVariance(D,nMin,nMax,c0,fitL0)
%Fitting procedure
FUN = @(x,xdata) (zernikeVarianceModel(x,xdata));
X0  = D./[0.16,50];
ub  = D./[0.01,D];
lb  = D./[1,100];

if ~fitL0
    X0 = X0(1);
    ub = ub(1);
    lb = lb(1);
end

opt = optimoptions(@lsqcurvefit,'MaxIter',1e2,'TolFun',1e-12,...
    'TolX',1e-12,'MaxFunEvals',3e2);

[beta,~,R,~,~,~,J] = lsqcurvefit(FUN,X0,nMin:nMax,c0,lb,ub,opt);

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


function out = zernikeVarianceModel(x,xdata)
%% ZERNIKEVARIANCE Zernike coefficients variance
%
% out = variance(modes,atmosphere) computes the
% variance of Zernike coefficients from the modes and the
% atmosphere object
%
% out = variance(zernike,atmosphere) computes the
% variance of Zernike coefficients from the Zernike polynomials
% object and the atmosphere object
%
% See also zernike, atmosphere

dr0     = abs(x(1));
if length(x) == 1
    dL0  = 0;
else
    dL0 = abs(x(2));
end
jv      = xdata;
[nv,mv] = nmOrder(jv);
nv0     = nv;
index   = diff(nv)~=0;
jv      = [jv(index) jv(end)];
mv      = [mv(index) mv(end)];
nv      = [nv(index) nv(end)];
nf      = length(nv);
out     = zeros(length(jv),1);

for cpt = 1:nf
    j = jv(cpt);
    n = nv(cpt);
    m = mv(cpt);
    out(nv0==n,1) = zernCovCoef(dr0,dL0,j,j,n,m,n,m);
end


function out = zernCovCoef(dr0,dL0,i,j,ni,mi,nj,mj)
if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
    if dL0==0
        if i==1 && j==1
            out = Inf;
        else
            out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                (dr0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                newGamma(-5./6+(ni+nj)./2,...
                [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
        end
    else
        out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
            (dr0./dL0).^(5./3)./dL0.^2.*...
            sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
            UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*dL0);
    end
else
    out = 0;
end

function out = newGamma(a,b)
% NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
% Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
% out = newGamma(a,b)

out = prod(gamma(a))./prod(gamma(b));

function out = UnParamEx4q2(mu,alpha,beta,p,a)
% UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
% of R. Conan (Modelisation des effets de l'echelle externe de coherence
% spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
% en Astronomie, University of Nice-Sophia Antipolis, October 2000)
% http://www-astro.unice.fr/GSM/Bibliography.html#thesis

a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
b1 = [1+alpha+beta 1+alpha 1+beta];
a2 = [(1-mu)./2+p 1+p p];
b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];

out = (1./(2.*sqrt(pi).*gamma(p))).*(...
    newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
    pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
    newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
    pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
function out = pochammerSeries(p,q,a,b,z,tol,nmax)
% POCHAMMERSERIES Computes power series in Pochammer notation
% pochammerSeries(p,q,a,b,z)
% pochammerSeries(p,q,a,b,z,tol)
% pochammerSeries(p,q,a,b,z,[],nmax)
% pochammerSeries(p,q,a,b,z,tol,nmax)

if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
    
    if p==length(a) && q==length(b)
        
        switch nargin
            case 6
                nmax = 1e3;
            case 7
                if isempty(tol)
                    tol = 1e-6;
                end
            otherwise
                tol = 1e-6;
                nmax = 1e3;
        end
        
        out = zeros(size(z));
        
        indz = find(z==0);
        if ~isempty(indz)
            out(indz) = 1;
        end
        
        indnz = find(z~=0);
        if ~isempty(indnz)
            z = z(indnz);
            ck = 1;
            step = Inf;
            k = 0;
            som = ck;
            while (k<=nmax) && (step>tol)
                ckp1 = prod(a+k).*z.*ck./prod(b+k);
                step = abs(abs(ck)-abs(ckp1));
                som = som + ckp1;
                k = k+1;
                ck = ckp1;
            end
            if step>tol
                warning('pochammerSeries','Maximum iteration reached before convergence')
            end
            out(indnz) = som;
        end
        
    else
        error('p and q must be the same length than vectors a and b, respectively')
        
    end
    
else
    error('This generalized hypergeometric function doesn''t converge')
end


function [n,m] = nmOrder(i)
% [n,m] = nmOrder(i,)
%Give the radial and azimutal order of the ith zernike
n = floor((-1.+sqrt(8*(i-1)+1))/2);
p = (i-(n.*(n+1))/2);
k = mod(n,2);
m = floor((p+k)/2)*2 - k;
