function [im, J_im] = imageModel(x,xdata,bestFit)  

%% 1\ Split parameters
xS = x(1:3*bestFit.nStars);
xPSF = x(1+3*bestFit.nStars:end-1);
xBg = x(end);
                          
%% 2\ Grab inputs
%2.1. PSF parameters init
gHO  = ones(1,bestFit.nGainsHO); gTT  = 1;gAl  = 1;
r053 = bestFit.psfr_.trs.atm.r0^(-5/3);Cn2=[];
fov_fit = bestFit.psfr_.psf.fov;
fov_sky = bestFit.psfr_.trs.cam.resolution;
samp = bestFit.psfr_.trs.cam.samp;
sf = bestFit.psfr_.sf;
otf = bestFit.psfr_.otf;

%2.2. Fitted parameters indexes
bestFit.idxCn2 = xdata{1};
bestFit.idxR0 = xdata{2};
bestFit.idxDho = xdata{3};
bestFit.idxDtt   = xdata{4};
bestFit.idxDal = xdata{5};
%2.3. Update parameters
if ~isempty(bestFit.idxR0)
    r053= sum(xPSF(bestFit.idxR0));
end
if ~isempty(bestFit.idxCn2)
    Cn2 = xPSF(bestFit.idxCn2);
    r053 = sum(Cn2);
end
if ~isempty(bestFit.idxDho)
    gHO = xPSF(bestFit.idxDho);
end
if ~isempty(bestFit.idxDal)
    gAl = xPSF(bestFit.idxDal);
end
if ~isempty(bestFit.idxDtt)
    gTT = xPSF(bestFit.idxDtt);
end
%2.4. Stellar parameters
nSrc = length(xS)/3;
fSrc = xS(1:nSrc);
xSrc = xS(1+nSrc:2*nSrc);
ySrc = xS(1+2*nSrc:3*nSrc);
dX   = [xSrc;ySrc];
                        
%% 3\ Get the image model accordingly the input parameters
bestFit.psf.rec = zeros(fov_fit,fov_fit,bestFit.nStars);
im  = zeros(fov_fit,fov_fit);
%3.1 Phase structure function on-axis
Don = r053/bestFit.psfr_.res.seeing.r0^(-5/3)*(sf.Dfit+ sf.Dal*gAl) + sum(bsxfun(@times,sf.Dho_z,reshape(gHO,1,1,[])), 3) + sf.Dtt*gTT;%*meter2rad^2;
sf.Dani = 0;
for iSrc = 1:bestFit.nStars
    %3.2 Anisoplanatism
    if numel(Cn2) > 1 && bestFit.flagAniso        
        sf.Dani   = sum(bsxfun(@times, squeeze(sf.Dani_l(:,:,:,iSrc)), reshape(Cn2,1,1,[])), 3);  
    end
    
    %3.3 OTF Multiplication
    otf.otfShannon = otf.otfStat.*exp(-0.5*(Don + sf.Dani));
    
    %3.4 PSF calculation
    psf_ij = tools.otfShannon2psf(otf.otfShannon,samp(iSrc),fov_fit);
    bestFit.psf.image(:,:,iSrc) = psf_ij/sum(psf_ij(:));
    
    %3.5 Include photometry and astrometry
    im  = im + fSrc(iSrc)*tools.translateImage(bestFit.psf.image(:,:,iSrc) ,dX(:,iSrc));
end
               
%3.6. Include background                       
im = im + xBg;

%3.7. Cropping
if fov_sky ~= fov_fit
    im = tools.crop(im,fov_sky);
end

% 3.8. Weighting
im(bestFit.im_sky == 0) = 0;
im = im.*bestFit.weightMap;
            
%% 4.\ JACOBIAN CALCULATION

if nargout >1
    
    %1\ Jacobian with respect to photometry and astrometry
    %Those ones are ok
    
    nX = numel(xS) + numel(xPSF) + numel(xBg);
    J_im = zeros(fov_sky^2,nX);
    for iSrc = 1:nSrc
        % Photometry
        J_im(:,1) = im(:)/fSrc;
        % Astrometry
        [u,v] = freqspace(size(im),'meshgrid');
        phasor = exp(-1i*pi*(v*dX(1) + u*dX(2)));
        otf_sky = tools.psf2otf(im);
        otf_sky = otf_sky/max(otf_sky(:));
        J_im(:,2)  = reshape(real(fftshift(ifft2(fftshift(-1i*pi*v.*otf_sky)))).*bestFit.weightMap,fov_sky^2,1);
        J_im(:,3) = reshape(real(fftshift(ifft2(fftshift(-1i*pi*u.*otf_sky)))).*bestFit.weightMap,fov_sky^2,1);
    end
    
    %2\ Jacobian with respect to PSF parameters
    %These ones are not fully ok
    if ~isempty(bestFit.idxDal)
        gAl = xPSF(bestFit.idxDal);
    else
        gAl = 1;
    end
    if ~isempty(bestFit.idxR0)
        r053= sum(xPSF(bestFit.idxR0));
    else
        r053 = bestFit.psfr_.trs.atm.r0^(-5/3);
    end
    
    % I'm not sure that the following is a proper way to get
    % the derivative
    n=1;
    if ~isempty(bestFit.idxR0)
        otf_1 = tools.interpolateOtf(padarray(-0.5*(sf.Dfit+ sf.Dal*gAl).*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
        J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*bestFit.weightMap,fov_sky^2,[]);
        n = n+1;
    end
    
    if ~isempty(bestFit.idxDho)
        for k=1:bestFit.nGainsHO
            otf_1 = tools.interpolateOtf(padarray(-0.5*sf.Dho_z(:,:,k).*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
            J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*bestFit.weightMap,fov_sky^2,[]);
            n = n+1;
        end
    end
    if ~isempty(bestFit.idxDtt)
        otf_1 = tools.interpolateOtf(padarray(-0.5*sf.Dtt.*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
        J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*bestFit.weightMap,fov_sky^2,[]);
        n = n+1;
    end
    if ~isempty(bestFit.idxDal)
        otf_1 = tools.interpolateOtf(padarray(-0.5*r053*sf.Dal.*otf.otfShannon,round((samp-1)*[fov_fit,fov_fit]/2),'both'),fov_sky);
        J_im(:,3*nSrc+n)    = reshape(real(fftshift(ifft2(fftshift(otf_1.*phasor)))).*bestFit.weightMap,fov_sky^2,[]);
        n = n+1;
    end
    % Jacobian with respect to background
    %This one is ok
    J_im(:,3*nSrc+n:end) = bestFit.weightMap(:);
end                  
        