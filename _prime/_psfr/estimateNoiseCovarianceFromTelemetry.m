%{
------------HEADER-----------------
Objective          ::  Estimate the noise covariance matrix

INPUT VARS
 trs          :: The telemetry top-level class
method    :: (optional) 'method': 'autocorrelation', 'interpolation' or 'rtf'

OUTPUT VARS
Cn                   :: The HO WFS noise covariance matrix (psfr.nActu^2 x psfr.nActu^2)
Cn_tt               :: The TT WFS noise covariance matrix (2 x 2)
var_n               :: The HO WFS noise variance in meter^2
var_n_tt           :: The TT WFS noise variance in meter^2
Created by       :: O. Beltramo-Martin - ONERA/LAM
Creation date   :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [Cn,Cn_tt,var_n,var_n_tt] = estimateNoiseCovarianceFromTelemetry(trs,varargin)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addParameter('method','autocorrelation',@ischar);
inputs.parse(trs,varargin{:});
method = inputs.Results.method;

%1\ HO WFS noise covariance matrix
if strcmp(method,'rtf')
	wvf = trs.waveFront;
else
	u = trs.hodm_pos;
	w = trs.waveFront;
	dt = trs.holoop_lat*trs.holoop_freq;
	dt = mod(dt,2);
	wvf = (w + (dt.*circshift(u,-3,2) + (1-dt).*circshift(u,-2,2)));% To be reviewed
    wvf = u;
end

Cn = getNoiseCovariance(wvf,'method',method,'rtf',trs.holoop_rtf);
var_n = trace(Cn)/size(Cn,1);

%2\ TT WFS noise covariance matrix
if strcmp(method,'rtf')
	wvftt = trs.tipTilt;
else
	utt = trs.ittm_pos;
	wtt = trs.tipTilt;
	dt = trs.ttloop_lat*trs.ttloop_freq;
	dt = mod(dt,2);
	wvftt = (wtt + (dt.*circshift(utt,-3,2) + (1-dt).*circshift(utt,-2,2)) ); % To be reviewed
    wvftt = utt;
end
Cn_tt = getNoiseCovariance(wvftt,'method',method,'rtf',trs.ttloop_rtf);
var_n_tt = trace(Cn_tt);


function Cnn = getNoiseCovariance(s,varargin)
inputs = inputParser;
inputs.addRequired('s',@isnumeric);
inputs.addParameter('method','autocorrelation',@ischar);
inputs.addParameter('nfit',1,@isnumeric);
inputs.addParameter('nshift',1,@isnumeric);
inputs.addParameter('rtf',1,@isnumeric);
inputs.parse(s,varargin{:});

method = inputs.Results.method;
nfit = inputs.Results.nfit;
nshift = inputs.Results.nshift;
rtf = inputs.Results.rtf;

%1\ Mean removal
s       = squeeze(s);
s       = bsxfun(@minus,s,mean(s,2));
[nS,nF] = size(s);
Cnn     = zeros(nS);

%2\ Covariance derivation
if strcmp(method,'interpolation') % Fitting the first samples of the auto-correlation function with a polynomial model
    delay   = linspace(0,1,nfit+2);
    for i=1:nS
        g      = (ifft(fft(s(i,:)).*conj(fft(s(i,:))))/nF);
        mx     = max(g(:));
        fit    = polyfit(delay(2:end),g(2:nfit+2),nfit);
        yfit   = 0*delay;
        for k =1:nfit+1
            yfit   = yfit + fit(k)*delay.^(nfit-k+1);
        end
        Cnn(i,i) = mx - yfit(1);
    end
    Cnn = diag(Cnn);
    
elseif strcmp(method,'autocorrelation') % Deriving the temporal cross-correlation difference over 2xnshift
    % If the turbulence is frozen over the  WFS frame rate, the difference of the autocorrelation and the 1-frame shifted correlation
    % is the noise variance.
    ds_n  = s - circshift(s,[0,-nshift]);
    ds_p  = s - circshift(s,[0,nshift]);
    Cnn = 0.5*(s*ds_n' + s*ds_p')/nF;             
elseif strcmp(method,'rtf') % Adjusting the noise level through the noise rejection transfer function model
    fftS = fft(s,[],2)/nF;
    fftS = fftS(:,1:floor(end/2))./RTF;
    cfftS= conj(fftS);
    nn = round(nF/10);
    
    for i=1:nS
        % Get the cross PSD
        crossPSD  = abs(bsxfun(@times,fftS(i,:),cfftS(i:end,:)));
        % Estimate the noise plateau
        Cnn(i,i:end)  = mean(crossPSD(:,end-nn:end),2);
    end
    Cnn = transpose(Cnn) + Cnn - diag(diag(Cnn));
end
