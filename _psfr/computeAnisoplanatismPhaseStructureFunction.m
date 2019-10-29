%{
------------HEADER-----------------
Objective          ::  Compute bi-dimensional phase structure function of the anisoplanatism error

INPUT VARS
 psfr          :: The PRIME top-level class
 
OUTPUT VARS
 Dani_l             :: bi-dimensional layer-by-layer phase structure function map  (psfr.nOtf x psfr.nOtf x psfr.atm.nLayer x psfr.nStars)
Dani                  :: bi-dimensional integrated phase structure function map  (psfr.nOtf x psfr.nOtf x psfr.nStars)
Kani                  :: Anisoplanatism spatial filter based on psfr.atm (psfr.nOtf x psfr.nOtf x psfr.nStars)
Created by         :: O. Beltramo-Martin - ONERA/LAM
Creation date      :: 10/04/2019
                      
Change Record:     ::
------------HEADER END----------------
%}

function [Dani_l,Dani,Kani] = computeAnisoplanatismPhaseStructureFunction(psfr)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'prime'));
inputs.parse(psfr);



%1\ Verify if there is anisoplanatism
Dani_l = 0;
Dani = 0;
Kani = 1;
if numel(psfr.gs) > 1 || numel(psfr.ttgs) >1
    psfr.flagAniso = false;
else
    thereisAngularAniso = any(any([psfr.src.directionVector]~=[psfr.gs.directionVector]));
    thereisFocalAniso   = any(any([psfr.src.height]~=[psfr.gs.height]));
    thereisAnisokinetism= (~isempty(psfr.ttgs)) && (any(any([psfr.src.directionVector]~=psfr.ttgs.directionVector)));
    psfr.flagAniso = thereisAngularAniso + thereisFocalAniso + thereisAnisokinetism;
end


if psfr.flagAniso
    psfr.atm.wavelength = psfr.src(1).photometry;
    %2\ Get the anisoplanatism phase structure function due to angular and focal anisoplanatism
    Dani_l= instantiateAnisoplanatism(psfr,psfr.gs);   
    Dani  = squeeze(sum(bsxfun(@times, psfr.Dani_l, reshape(psfr.atm.r0^(-5/3)*[psfr.atm.layer.fractionnalR0],1,1,[])), 3));
    
    if thereisAnisokinetism
    %3\ Get the tip-tilt anisoplanatism phase structure function
        DaniTT_l  = instantiateAnisoplanatism(psfr,psfr.ttgs,'isTT',true);
        Dani  = squeeze(Dani+ sum(bsxfun(@times, DaniTT_l, reshape(psfr.atm.r0^(-5/3)*[psfr.atm.layer.fractionnalR0],1,1,[])), 3));
        Dani_l= Dani_l  + DaniTT_l;
    end
    
    %4\ Get the anisoplanatism spatial filter
    if psfr.flagToeplitz
        Kani = exp(-0.5*psfr.Dani);
    else
        idx       = true(psfr.dk);
        dp        = psfr.tel.D/(psfr.dk-1);
        psfr.otfDen= tools.zonalCovarianceToOtf(Dani*0,psfr.nOtf,psfr.tel.D,dp,idx);
        psfr.mskOtf= psfr.otfDen > 1e-6;
        Kani  = tools.zonalCovarianceToOtf(Dani,psfr.nOtf,psfr.tel.D,dp,idx);
        Kani(psfr.mskOtf) = Kani(psfr.mskOtf)./psfr.otfDen(psfr.mskOtf);
    end
end

function Dani_l = instantiateAnisoplanatism(psfr,gs,varargin)
inputs = inputParser;
inputs.addRequired('psfr',@(x) isa(x,'prime'));
inputs.addRequired('gs',@(x) isa(x,'source'));
inputs.addParameter('isTT',false,@islogical);
inputs.parse(psfr,gs,varargin{:});


fprintf('Instantiate the anisoplanatism model...');

%1\ Defining the spatial filters
D       = psfr.tel.D;
npt  = psfr.dk;
isTT = inputs.Results.isTT;
if gs.height ~= Inf
    zern   = zernike(2:3,'resolution',npt);
    TT     = zern.modes;
    Hfilter= eye(npt^2) - TT*pinv(TT);
    Hfilter= Hfilter*psfr.trs.Hdm;
elseif gs.height == Inf && isTT
    x      = (-1+1/npt:2/npt:1-1/npt);
    [X,Y]  = meshgrid(x,x);
    TT     = [X(:),Y(:)];
    Hfilter= TT*pinv(TT);
else
    Hfilter      = psfr.Hdm;
end

%2\ SF Calculation
if psfr.flagToeplitz
    Dani_l = zeros(psfr.nOtf,psfr.nOtf,psfr.atm.nLayer,psfr.nStars);
else
    Dani_l  = zeros(psfr.dk^2,psfr.dk^2,psfr.atm.nLayer,psfr.nStars);
end

if strcmpi(psfr.flagAnisoMethod,'FLICKER') % Use the Flicker's 2008 report to derive the SF
    % Get inputs
    f0      = 2*pi/psfr.atm.L0;
    % Phase sample locations in the pupil
    x       = -D/2:D/(psfr.dk-1):D/2;
    [x1,y1] = meshgrid(x);
    X1      = (x1(:)*ones(1,psfr.dk^2))';
    Y1      = repmat(y1,[psfr.dk,psfr.dk]);
    % Samples separation in the pupil
    rhoX    = bsxfun(@minus,X1,x1(:));
    rhoY    = bsxfun(@minus,Y1,y1(:)');
    % Instantiation
    Ialpha  = @(x,y) tools.mcDonald(f0*hypot(x,y));
    I0      = Ialpha(0,0);
    I1      = Ialpha(rhoX,rhoY);
    cte     = 0.12184*0.06*(2*pi)^2;
    
    % Anisoplanatism Structure Function
    for iSrc = 1:numel(psfr.src)
        sep = psfr.src(iSrc).directionVector - gs.directionVector;
        thx = sep(1);
        thy = sep(2);
        for l = 1:psfr.atm.nLayer
            zl   = psfr.atm.layer(l).altitude;
            if gs.height == Inf  || isempty(gs.height)
                I2    = Ialpha(rhoX+zl*thx,rhoY+zl*thy);
                I3    = Ialpha(zl*thx,zl*thy);
                I4    = Ialpha(rhoX-zl*thx,rhoY-zl*thy);
                tmp   = I0 - 2*I1 + I2 - 2*I3  + I4;
            else
                g     = zl/gs.height;
                I2    = Ialpha(rhoX*(1-g),rhoY*(1-g));
                I3    = Ialpha(rhoX-g*X1+zl*thx,rhoY-g*Y1+zl*thy);
                I4    = Ialpha(g*X1-zl*thx,g*Y1-zl*thy);
                I5    = Ialpha(g*(rhoX-X1)-zl*thx,g*(rhoY-Y1)-zl*thy);
                I6    = Ialpha((1-g)*rhoX+g*X1-zl*thx,(1-g)*rhoY+g*Y1-zl*thy);
                tmp   = 2.*I0 - I1 - I2 + I3 - I4 - I5 + I6;
            end
            
            if psfr.flagToeplitz
                tmp = tools.covMatrix2Map(tmp,psfr.dk,psfr.dk);
                tmp = tools.interpolateOtf(tmp,psfr.nOtf);
            end
            Dani_l(:,:,l,iSrc)  = Dani_l(:,:,l,iSrc) + cte*psfr.atm.L0^(5/3)*Hfilter*tmp*Hfilter';
        end
    end
    
elseif strcmpi(psfr.flagAnisoMethod,'FOURIER') % Model the anisoplanatism into spatial fequencies domain
    
    fao = spatialFrequencyAdaptiveOptics(psfr.tel,psfr.atm,psfr.nActu,...
        psfr.var_n*(2*pi/psfr.wvl)^2,psfr.holoop_gain,1/psfr.holoop_freq,...
        psfr.holoop_lat,psfr.dk,'rigaut',1,'nTimes',psfr.nTimes);
    close();
    % Redefine the sources
    [xS,yS]     = pol2cart([psfr.src.azimuth],[psfr.src.zenith]);
    [xG,yG]     = pol2cart([psfr.gs.azimuth],[psfr.gs.zenith]);
    [aS,zS]     = cart2pol(xS-xG,yS-yG);
    srcOff      = source('wavelength',psfr.src(1).photometry,'azimuth',aS,'zenith',zS);
    fao.src = srcOff;
    atmRef      = tools.duplicateAtmosphere(psfr.atm);
    % Loop on stars
    for iSrc = 1:numel(psfr.src)
        fao.atm.wavelength = srcOff(iSrc).wavelength;
        for l=1:psfr.atm.nLayer
            subAtm    = psfr.atm.slab(l);
            if subAtm.layer.altitude
                % Get the aniso PSD
                fao.atm = subAtm;
                psdAniso = fao.anisoplanatismPSD(fao.fx,fao.fy,iSrc);
                % Get the SF
                sfAniso  = fftshift(tools.cov2sf(tools.psd2cov(psdAniso,2*fao.fc/psfr.dk)));
                % Normalize
                f0 = (psfr.atm.r0^(-5/3)*psfr.atm.layer(l).fractionnalR0);
                Dani_l(:,:,l,iSrc) = tools.interpolateOtf(sfAniso/f0,psfr.nOtf);
            end
        end
    end
else % Use the native OOMAO routines in phaseStats to calculate the anisoplanatism SF.
    %Layers-by-layers covariance matrices
    cov = @(dk,D,atm,gs) phaseStats.spatioAngularCovarianceMatrix(dk,D,psfr.atm,gs);
    cross = @(dk,D,atm,src,gs) phaseStats.spatioAngularCovarianceMatrix(dk,D,psfr.atm,src,'srcCC',gs);        
    thereIsFocAniso = gs.height~=Inf;
    
    for l=1:psfr.atm.nLayer
        subAtm    = psfr.atm.slab(l);
        if subAtm.layer.altitude % Check that the layer is not at the ground
            f0 = (psfr.atm.r0^(-5/3)*psfr.atm.layer(l).fractionnalR0);
            % Atmospheric covariance matrix in the guide star direction
            Cgg = cov(psfr.dk,D,subAtm,gs);
            
            for iSrc=1:psfr.nStars
                thereIsAngAniso = any(psfr.src(iSrc).directionVector - gs.directionVector);
                % Phase covariance matrix
                if psfr.src(iSrc).height ~= gs(iSrc).height
                    %LGS case: Css ~= Csg
                    Css = cov(psfr.dk,D,subAtm,psfr.src(iSrc));
                else
                    %NGS case: Css = Csg
                    Css = Cgg;
                end
                % Cross covariance matrix
                if thereIsAngAniso || thereIsFocAniso
                    Cgs = cross(psfr.dk,D,subAtm,psfr.src(iSrc),gs);
                    Cgs = Cgs{1};
                else
                    Cgs = 0;
                end
                % Anisoplanatic covariance matrix
                Cani = Hfilter*(Css + Cgg - Cgs - Cgs')*Hfilter'/f0;
                if psfr.flagToeplitz
                    Cani = tools.covMatrix2Map(Cani,psfr.dk,psfr.dk);
                    Cani = tools.interpolateOtf(Cani,psfr.nOtf);
                    Dani_l(:,:,l,iSrc) = 2*(max(Cani(:)) - Cani);
                else                
                    Dani_l(:,:,l,iSrc) = 2*(diag(Cani) - Cani);
                end
            end
        end
    end
end


