function displayResults(obj,fov)

im = obj.psfr_.trs.cam.frame;
if nargin <3
    fov = size(im,1);
end
%1\ 1D plot of PSF
fov_im     = size(im,1);
if fov < fov_im
    im = tools.crop(im,fov);
    imfit = tools.crop(obj.im.rec,fov);
else
    fov = fov_im;
    imfit = obj.im.rec;
end
im(im<0)=0;
n       = floor(fov/2);
x       = linspace(0,obj.psfr_.trs.cam.pixelScale*fov/2,n);
xfull   = obj.psfr_.trs.cam.pixelScale*fov/2*linspace(-1,1,fov);

h = figure;
subplot(1,3,1)
semilogy(x,radial(im,n,n),'b-');hold on;
semilogy(x,radial(imfit,n,n),'r--');
semilogy(x,abs(radial(imfit-im,n,n)),'k:');
legend({'Sky image','Best-fitted model','Residual'},'interpreter','latex','FontSize',16,'Location','southwest');
set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
tools.makeAxisSquared(h);
subplot(1,3,2)
semilogy(xfull,im(:,n),'b-');hold on;
semilogy(xfull,imfit(:,n),'r--');
semilogy(xfull,abs(imfit(:,n)-im(:,n)),'k:');
legend({'Sky image','Best-fitted model','Residual'},'interpreter','latex','FontSize',16,'Location','southwest');
set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
tools.makeAxisSquared(h);
subplot(1,3,3)
semilogy(xfull,im(n,:),'b-');hold on;
semilogy(xfull,imfit(n,:),'r--');
semilogy(xfull,abs(imfit(n,:)-im(n,:)),'k:');
legend({'Sky image','Best-fitted model','Residual'},'interpreter','latex','FontSize',16,'Location','southwest');
set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
tools.makeAxisSquared(h);

%2\ 1D plot of OTF
otf_sky = abs(fftshift(fft2(im)));
otf_sky = otf_sky/max(otf_sky(:));
otf_fit = abs(fftshift(fft2(imfit)));
otf_fit = otf_fit/max(otf_fit(:));

u       = linspace(0,max(obj.psfr_.trs.cam.samp,1),fov/2);
h = figure;
semilogy(u,(radial(otf_sky)),'b--');hold on;
semilogy(u,(radial(otf_fit)),'r-');
semilogy(u,abs(radial(otf_fit-otf_sky)),'k:');
xlim([0,1]);
legend({'Sky OTF','Best-fitted model','Residual'},'interpreter','latex','FontSize',18);
set(gca,'FontSize',18,'FontName','cmr12','TickLabelInterpreter','latex' );
tools.makeAxisSquared(h);

%3\ 2D PSF;
P = [im,imfit,imfit-im];
Plog = log10(abs(P));
if sign(max(Plog(:))) ==1
    mn = 0.2*max(Plog(:));
    mx = 0.9*max(Plog(:));
else
    mn = 5*max(Plog(:));
    mx = 1.1*min(Plog(:));
end

h = figure;
imagesc(Plog,[mn,mx]);
set(gca,'XTick',[],'YTick',[]);
end
