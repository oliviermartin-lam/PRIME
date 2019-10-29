clear all;close all;


%% -------------- 2013 NGS DATA
date = '20130801';
idxO = {'n0004','n0005','n0006','n0007','n0008','n0009','n0010','n0011','n0012',...
    'n0013','n0050','n0051','n0052','n0053','n0054','n0065','n0066','n0067',...
    'n0068','n0069','n0085','n0086','n0070','n0071','n0072','n0073','n0074',...
    'n0087','n0088','n0089','n0090','n0092','n0093','n0094','n0095','n0096'...
    ,'n0097','n0098','n0099','n0100','n0101','n0102','n0103','n0104','n0105'};
trs_suffix   = 'NGS';
nObj = numel(idxO);

% Instantiating output tables
SR     = zeros(2,3,nObj);
FWHM   = zeros(2,3,nObj);
EQM    = zeros(2,nObj);
wfe    = zeros(11,nObj);
r0_fit = zeros(2,nObj);

% loop on files
for iObj = 1:nObj
    % Run PSFR
    obj_name     = cell2mat(idxO(iObj));
    scriptKeck_PSFR
    % Buffering results
    %1\ PSF statistics
    SR(1,:,iObj)   = [psfr.SR_sky,psfr.SR_init,psfr.SR_fit];
    SR(2,:,iObj)   = [psfr.dSR_sky,psfr.dSR_init,psfr.dSR_fit];
    FWHM(1,:,iObj) = [psfr.FWHM_sky,psfr.FWHM_init,psfr.FWHM_fit];
    FWHM(2,:,iObj) = [psfr.dFWHM_sky,psfr.dFWHM_init,psfr.dFWHM_fit];
    EQM(:,iObj)    = [psfr.eqm_init,psfr.eqm_fit]*1e2;
    %2\ Error breakdown
    wfe(:,iObj)    = struct2array(psfr.wfe);
    %3\ Parameters analysis
    r0_fit(:,iObj) = [psfr.r0_fit,psfr.r0_prec];
    
end

res.SR   = SR;
res.FWHM = FWHM;
res.wfe  = wfe;
res.r0   = r0_fit;
save(['/home/omartin/Projects/PSFR/KECK_PSFR/Results/Sky/',date,'/keck_psfr_',date,'_',trs_suffix,'.mat'],'res');

%% DISPLAY
h = figure;
errorbar(squeeze(SR(1,1,:)),squeeze(SR(1,2,:)),-squeeze(SR(2,2,:)),squeeze(SR(2,2,:)),...
    -squeeze(SR(2,1,:)),squeeze(SR(2,1,:)),'bs','MarkerFaceColor','b','MarkerSize',5);
hold on
errorbar(squeeze(SR(1,1,:)),squeeze(SR(1,3,:)),-squeeze(SR(2,3,:)),squeeze(SR(2,3,:)),...
    -squeeze(SR(2,1,:)),squeeze(SR(2,1,:)),'ro','MarkerFaceColor','r','MarkerSize',5);
plot([0,1],[0,1],'--');


h = figure;
errorbar(squeeze(FWHM(1,1,:)),squeeze(FWHM(1,2,:)),-squeeze(FWHM(2,2,:)),squeeze(FWHM(2,2,:)),...
    -squeeze(FWHM(2,1,:)),squeeze(FWHM(2,1,:)),'bs','MarkerFaceColor','b','MarkerSize',5);
hold on
errorbar(squeeze(FWHM(1,1,:)),squeeze(FWHM(1,3,:)),-squeeze(FWHM(2,3,:)),squeeze(FWHM(2,3,:)),...
    -squeeze(FWHM(2,1,:)),squeeze(FWHM(2,1,:)),'ro','MarkerFaceColor','r','MarkerSize',5);
plot([0,120],[0,120],'--');

