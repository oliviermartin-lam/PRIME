%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the aliasing phase

INPUT VARS
 psfr          :: The PRIME top-level class
 
OUTPUT VARS
 sf_2D             :: bi-dimensional phase structure function map (Toeplitz) of the aliasing error (psfr.nOtf x psfr.nOtf)
Cal                  :: Covariance matrix of the aliasing error in the acturators space (nActu^2 x nActu^2)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [sf_2D,Cal] = computeAliasingPhaseStructureFunction(psfr,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'psfr'));
inputs.addParameter('aoPattern','circle',@ischar)
inputs.parse(psfr,varargin{:});
aoPattern = inputs.Results.aoPattern;

%1\ Parsing inputs
atm = psfr.trs.atm;
D = psfr.trs.tel.Ddm;
nActu = psfr.trs.dm.nActuators;
nT = psfr.otf.nTimes;
ts = 1/psfr.trs.holoop.freq;
td = psfr.trs.holoop.lat;
varn = psfr.res.noise.varn_ho*(2*pi/psfr.trs.cam.wavelength)^2;
g = psfr.trs.holoop.gain;


%2\ Define the frequency space
nK   = 2*nActu-1;
[kx,ky] = freqspace(nK,'meshgrid');
kc  = 0.5*(nActu-1)/D;
kx = kx*kc + 1e-7;
ky = ky*kc + 1e-7;
k   = hypot(kx,ky);
dk  = 2*kc/nK; % Pixel scale


%3\ Define the atmospheric phase PSD
fr0 = atm.weights;
L0 = psfr.res.seeing.L0;
r0 = psfr.res.seeing.r0;
[vlx,vly] = pol2cart(atm.windDirection,atm.windSpeed);
cst = (24*gamma(6/5)/5)^(5/6)*(gamma(11/6)^2/(2*pi^(11/3)));
Wphi = r0^(-5./3)*cst*(k.^2 + 1/L0.^2).^(-11./6);
                     
%4\ Calculating the WFS reconstructor
Wn = varn/(2*kc)^2;
d = D/(nActu-1);
Sx = 1i*2*pi*kx*d ;
Sy = 1i*2*pi*ky*d ;
AvRec     = sinc(d*kx).*sinc(d*ky).*exp(1i*pi*d*(kx+ky));
SxAvRec = Sx.*AvRec;
SyAvRec = Sy.*AvRec;
gPSD = abs(SxAvRec).^2 + abs(SyAvRec).^2 + Wn./Wphi;
Rx = conj(SxAvRec)./gPSD;
Ry = conj(SyAvRec)./gPSD;
Rx(ceil((nK+1)/2),ceil((nK+1)/2)) = 0;
Ry(ceil((nK+1)/2),ceil((nK+1)/2)) = 0;

%5\ Calculating the AO controller trandfer function
%5.1 Spatializing the transfer function
h1          = zeros(nK);
for kLayer = 1:atm.nLayer
    % Define the temporal frequency
    fi = -vlx(kLayer)*kx + -vly(kLayer)*ky;
    idx = abs(fi) <1e-7;
    fi(idx) = 1e-8.*sign(fi(idx));
    % Get the AO closed-loop transfer function
    tfwfs    = wfsTransferFunction(fi,psfr.trs.holoop.freq);
    tflag = delayTransferFunction(fi,psfr.trs.holoop.freq,psfr.trs.holoop.lat*psfr.trs.holoop.freq);
    tfservo = servoTransferFunction(fi,psfr.trs.holoop.freq,psfr.trs.holoop.tf.num,psfr.trs.holoop.tf.den);
    tfol = tfwfs.*tfservo.*tflag;
    h1 = h1 + fr0(kLayer)*tfol./(1+tfol);
end

%5.2 Get the PSD
psd = 0;
w = 2*1i*pi*d;
for mi = -nT:nT %loop on x-axis spatial frequencies
    for ni = -nT:nT %loop on y-axis spatial frequencies
        if mi~=0 || ni~=0           
            % Get the shifted atmospheric spectrum
	    km = kx - mi/d;
            kn = ky - ni/d;
            W_mn = (hypot(km,kn).^2 + 1/L0^2).^(-11/6).*pistonRemoval(D,d,kx,mi,ni);           
            % Spatial filtering
            Q = (Rx.*w.*km + Ry.*w.*kn) .* (sinc(d*km).*sinc(d*kn));
            % Temporal filter
            avr = 0;
            for kLayer = 1:atm.nLayer
                avr = avr + fr0(kLayer)*...
                    (sinc(km*vlx(kLayer)*ts).* sinc(kn*vly(kLayer)*ts).*...
                    exp(1i*2*pi*km*vlx(kLayer)*td).*exp(1i*2*pi*kn*vly(kLayer)*td).*h1);
            end
            q =  Q.* avr;
            psd = psd + W_mn .* (q.*conj(q));
        end
    end
end
psd(isnan(psd)) = 0;
%5.4 Filter the uncorrected part
if strcmp(aoPattern,'circle')
    idx = k<=kc;
else
    idx = abs(kx)<=kc & abs(ky)<=kc;
end
psd = real(cst*psd.*pistonRemoval(D,d,kx,0,0)).*idx;
psd_pad = tools.enlargePupil(psd,nT);
covMap = real(tools.psd2cov(psd_pad,dk));
sf_2D =  tools.cov2sf(covMap);

if size(sf_2D,1)~=psfr.otf.nOtf
    sf_2D = tools.interpolateOtf(sf_2D,psfr.nOtf);
end

if nargout > 1
    Cal = tools.covMap2Matrix(real(tools.psd2cov(psd,dk)),round(nK/2),round(nK/2));
end

function PR = pistonRemoval(D,d,f,mi,ni)
if ~isvector(f)
    f = f(1,:);
end
besselargx     = pi*D/d*(f -mi/d) * d; % here I should not divide by dx because the frequency vector is already in proper units
besselargy     = pi*D/d*(f -ni/d) * d;
[BX, BY]       = meshgrid(besselargx,besselargy);
fltr           = 2*besselj(1,sqrt(BX.^2+BY.^2))./(sqrt(BX.^2+BY.^2));
fltr(isnan(fltr)) = 1;
PR = 1-abs(fltr).^2;
