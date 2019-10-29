%Must be user-defined
path_root = '/run/media/omartin/HDD_OBM/KECK_DATA/';
%%
if false
    date         = '20130801';
    obj_name     = 'n0084';
    trs_suffix   = 'NGS';
end

path_data    = [path_root,date,'/TRS/',obj_name,'_',trs_suffix,'_trs.sav'];
path_im      = [path_root,date,'/NIRC2/',obj_name,'.fits'];
path_calib   = [path_root,'CALIBRATION/'];
path_profile = [path_root,'MASSDIMM/'];

%% GRAB TELEMETRY
if strcmp(trs_suffix,'NGS')
    parFileKeck_NGS();
else
    parFileKeck_LGS();
end
telem = telemetry('KECK',date,path_data,path_im,path_calib,'path_profile',path_profile,...
    'resolution',parm.cam.resolution);

%% DEFINE ATM
if ~isempty(telem.r0_ext)
    fl  = telem.Cn2_ext(1,:);
    alt = telem.Cn2_ext(2,:);
    wS  = [6.8 6.9 7.1 7.5 10.0 26.9 18.5];
    wD  = [0 pi/2 pi/4 3*pi/2 -pi/4 pi/6 -pi/4];
    %Removw very wek layers
    idl = fl>0.01;
    fl  = fl(idl)/sum(fl(idl));
    alt = alt(idl);
    wS  = wS(idl);
    wD  = wD(idl);
    % Update the parameters list
    parm.atm.photometry    = photometry.V0;
    parm.atm.r0            = telem.r0_ext;
    parm.atm.fractionalR0  = fl;
    parm.atm.altitude      = alt;
    parm.atm.windSpeed     = wS;
    parm.atm.windDirection = wD;
    if strcmp(trs_suffix,'LGS')
        parm.lGs.height    = parm.lGs.height*telem.airmass;
    end
end
    %% PSFR
psfr = prime(parm,telem,'dk',43,'flagToeplitz',true);
psfr = psfr.bestFitting(telem.im_sky,'weighting',false,'fitGains',[true,true,false],'fitCn2',false,'nZernMode',[4:6]);
psfr.errorBreakDown();