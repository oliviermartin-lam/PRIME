%{
------------HEADER-----------------
Objective          :: Provide a model of the Zernike coefficients variance
 
INPUT VARS
x              :: the input parameters [D/r0,D/L0]
xdata       :: the Zernike coefficients j-index

OUTPUT VARS
zern_var           :: The Zernike coefficients variance
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 11/01/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function zern_var = zernikeVarianceModel(x,xdata)
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
zern_var     = zeros(length(jv),1);

for cpt = 1:nf
    j = jv(cpt);
    n = nv(cpt);
    m = mv(cpt);
    zern_var(nv0==n,1) = zernCovCoef(dr0,dL0,j,j,n,m,n,m);
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


