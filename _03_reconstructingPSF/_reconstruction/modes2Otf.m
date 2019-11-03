function [otf,dphi] = modes2Otf(Cmm,modes,pupil,npsf,varargin)
inputs = inputParser;
inputs.addRequired('Cmm',@isnumeric);
inputs.addRequired('modes',@isnumeric);
inputs.addRequired('pupil',@isnumeric);
inputs.addRequired('npsf',@isnumeric);
inputs.addParameter('D',11.25,@isnumeric);
inputs.addParameter('dk',43,@isnumeric);
inputs.addParameter('method','Vii',@ischar);
inputs.parse(Cmm,modes,pupil,npsf,varargin{:});

method       = inputs.Results.method;
D = inputs.Results.D;
dk = inputs.Results.dk;

if any(Cmm(:))
            
    switch method
        
        case 'Vii'
            % Autocorrelation of the pupil expressed in pupil
            nPx        = sqrt(size(modes,1));
            pupExtended= fftshift(tools.enlargePupil(double(pupil),2));
            fftPup     = fft2(pupExtended);
            conjPupFft = conj(fftPup);
            G          = fftshift(real(fft2(fftPup.*conjPupFft)));
            % Defining the inverse
            den        = zeros(size(G));
            msk        = G./max(G(:)) > 1e-7;
            den(msk)   = 1./G(msk);
            % Diagonalizing the Cvv matrix
            [U,S]   = svd(Cmm);
            s       = diag(S);
            nModes  = length(s);
            M       = modes * U;
            %loop on actuators            
            buf = zeros(size(pupExtended));            
            for k=1:nModes
                Mk   = reshape(M(:,k),nPx,nPx);
                Mk   = fftshift(tools.enlargePupil(Mk,2));
                % Vii computation
                Vk   = real(fft2(Mk.^2.*pupExtended).*conjPupFft) - abs(fft2(Mk.*pupExtended)).^2;
                % Summing modes into dm basis
                buf  = buf + s(k) .* Vk;
            end            
            dphi     = den.*fftshift(real(fft2(2*buf)));
            otf      = G.*exp(-0.5*dphi);
            
        case 'Uij'
            % Autocorrelation of the pupil expressed in pupil
            nPx        = sqrt(size(modes,1));
            pupExtended= fftshift(tools.enlargePupil(double(pupil),2));
            fftPup     = fft2(pupExtended);
            conjPupFft = conj(fftPup);
            G          = fftshift(real(fft2(fftPup.*conjPupFft)));
            % Defining the inverse
            den        = zeros(size(G));
            msk        = G./max(G(:)) > 1e-7;
            den(msk)   = 1./G(msk);
            nm   = size(modes,2);
            dphi = 0*pupExtended;
            
            %Double loops on modes
            for i=1:nm
                Mi = reshape(modes(:,i),nPx,nPx);
                Mi = tools.enlargePupil(Mi,2);
                for j=1:i
                    %Getting modes + interpolation to twice resolution
                    Mj    = reshape(modes(:,j),nPx,nPx);
                    Mj   = tools.enlargePupil(Mj,2);
                    term1 = real(fft2(Mi.*Mj.*pupExtended).*conjPupFft);
                    term2 = real(fft2(Mi.*pupExtended).*conj(fft2(Mj.*pupExtended)));
                    % Uij computation
                    Uij   = real(ifft2(term1-term2));
                    %Summing terms
                    fact = double(i~=j) + 1;
                    dphi = dphi + fact*Cmm(i,j).*Uij;
                end
            end
            dphi = fftshift(dphi).*den.*msk;
            otf  = G.*exp(-0.5*dphi);
            
        case 'zonal'            
            % Grabbing the valid actuators positions in meters           
            loc  = pointWiseLocation(D,D/(dk-1),true(dk));
            %OTF sampling
            nPh  = round(D/dp+1);
            %nU1d = 2*nPh;
            u1D  = (-nPh:1:nPh-1)*dp;
            nU1d = length(u1D);
            % Determining couples of point with the same separation
            [shiftX,shiftY] = mkotf_indpts(nU1d,nPh,u1D,loc,dp);
            % WF Amplitude
            amp0 = ones(nnz(idxValid),1);
            
            %% Long-exposure OTF
            otf = mkotf(shiftX,shiftY,nU1d,amp0,dp,-0.5*Cmm);
            otf_dl = mkotf(shiftX,shiftY,nU1d,amp0,dp,0*Cmm);
            %Interpolation
            if size(otf,1) ~=npsf
                otf      = tools.interpolateOtf(otf,npsf);
                otf_dl = tools.interpolateOtf(otf_dl,npsf);
            end
            otf = otf./otf_dl.*(otf_dl/max(otf_dl(:))>1e-5);
            otf      = otf./max(otf(:));
            dphi = -2*log(otf);
    end
else
    % Diffraction-limit case
    G    = G./max(G(:));
    otf  = G;
    dphi = 0*G;
end
% Interpolation of the OTF => determination of the PSF fov
otf = otf.*(G>1e-5);
otf = otf./max(otf(:));
if size(otf,1) ~= npsf
    otf = tools.interpolateOtf(otf,npsf);
    otf = otf./max(otf(:));
    dphi= tools.interpolateOtf(dphi,npsf);
end

%Zonal calculation
function out = pointWiseLocation(D,dp,idxValid)
% Defining the point-wise locations
xloc                 = -D/2:dp:D/2;
[actuLocX, actuLocY] = meshgrid(xloc);
actuLocX             = actuLocX(idxValid);
actuLocY             = actuLocY(idxValid);
out                  = [actuLocX(:), actuLocY(:)];


function otf = mkotf(indptsc,indptsc2,nU1d,ampl,dp,C_phi)

%Instantiation
otf         = zeros(nU1d);
C_phip_diag = exp(diag(C_phi));
C_phipi     = exp(-2*C_phi);
C_phi_size2 = size(C_phi,2);

for iu=1:nU1d
    for ju=1:nU1d
        indpts  = indptsc{iu, ju};
        indpts2 = indptsc2{iu, ju};
        
        if isempty(indpts)
            otf(iu,ju) = 0;
        else
            myarg      = C_phip_diag(indpts2).*C_phip_diag(indpts)...
                .*C_phipi(C_phi_size2*(indpts-1) + indpts2);
            kernel     = (conj(ampl(indpts2)) .* ampl(indpts))' * myarg;
            otf(iu,ju) = kernel*dp^2;
        end % if isempty(mypts)
    end %ju
end %iu

dc  = sum(abs(ampl).^2) * dp^2;
otf = otf/dc;


function [indptsc,indptsc2] = mkotf_indpts(nU1d,nPh,u1D,loc,dp)

% index pts in a 3x bigger array
locInPitch = loc/dp;
minLoc     = min(locInPitch);
ninloc     = size(loc,1);
nc         = 3*nPh;
n          = nc-nPh;
n1         = (n-mod(n,2))/2 + 1 + mod(nPh,2);
minLoc2    = minLoc -(n1-1);
loc2       = round(locInPitch - ones(ninloc,1)*minLoc2+1);
%embed the loc2 inside ncxnc array.
indx_emb       = loc2(:,1)+(loc2(:,2)-1)*nc;
mask           = zeros(nc);
mask(indx_emb) = 1;
indptsc        = cell(nU1d);
indptsc2       = cell(nU1d);

for iu=1:nU1d
    for ju=1:nU1d
        u2D        = [u1D(iu)  u1D(ju)];
        uInPitch   = u2D/dp;
        iniloc_sh2 = loc2 + ones(ninloc,1)*round(uInPitch);
        indxsh     = iniloc_sh2(:,1) + (iniloc_sh2(:,2)-1)*nc;
        
        %index of points in iniloc_sh2 that are intersect with iniloc_sh
        %indpts is the overlapping points in the shifted array
        indptsc{ju,iu}  = find(mask(indxsh));
        mask2           = zeros(nc);
        mask2(indxsh)   = 1;
        indptsc2{ju,iu} = find(mask2(indx_emb));
    end
end
