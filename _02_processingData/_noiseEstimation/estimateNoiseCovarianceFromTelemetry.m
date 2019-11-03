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

function res = estimateNoiseCovarianceFromTelemetry(trs,varargin)
inputs = inputParser;
inputs.addRequired('trs',@(x) isa(x,'telemetry'));
inputs.addParameter('method','autocorrelation',@ischar);
inputs.parse(trs,varargin{:});
method = inputs.Results.method;

res = [];
%1\ HO WFS noise covariance matrix
res.Cn_ho = getNoiseCovariance(trs.dm.com,'method',method,'rtf',trs.holoop.tf.ctf./trs.holoop.tf.wfs,'FSAMP',trs.holoop.freq);
validInput = std(trs.dm.com,[],2)~=0;
res.varn_ho = trace(res.Cn_ho(validInput,validInput))/nnz(validInput);

%2\ TT WFS noise covariance matrix
res.Cn_tt = getNoiseCovariance(trs.tipTilt.com,'method',method,'rtf',trs.ttloop.tf.ctf./trs.ttloop.tf.wfs,'FSAMP',trs.ttloop.freq);
res.varn_tt = trace(res.Cn_tt);


function Cnn = getNoiseCovariance(s,varargin)
inputs = inputParser;
inputs.addRequired('s',@isnumeric);
inputs.addParameter('FSAMP',1e3,@isnumeric);
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
validInput = std(s,[],2)~=0;
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
    fftS = fft(s,[],2); % we verified that std(s(k,:))^2 = sum(abs(fftS(k,:)).^2): Parseval
    fftS = fftS(:,1:floor(end/2))./rtf;
    % truncate to keep high-spatial frequencies only
    sub_dim = floor(9*nF/2/10);
    fftS = fftS(:,sub_dim:end);
    cfftS= conj(fftS);
       
    for i=1:nS 
        if validInput(i)
            % Get the cross PSD
            crossPSD  = bsxfun(@times,fftS(i,:),cfftS(i:end,:));
            % Estimate the noise plateau
            Cnn(i,i:end)  = median(real(crossPSD),2);
        end
    end
    Cnn = transpose(Cnn) + Cnn - diag(diag(Cnn));    
    % Normalization
    Cnn = Cnn/nF;
end
