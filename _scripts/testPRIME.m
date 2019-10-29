clear all;
close all;
%% MANAGE WORKING PATHS
path_workspace = '/home/omartin/Projects/KVS/CODES/';
addpath(genpath(path_workspace));
%% DEFINE DATA PATH
path_root = '/run/media/omartin/HDD_OBM/KECK_DATA/';
date          = '20130801';
obj_name  = 'n0004';
trs_suffix   = 'NGS';
path_data    = [path_root,date,'/TRS/',obj_name,'_',trs_suffix,'_trs.sav'];
path_im      = [path_root,date,'/NIRC2/',obj_name,'.fits'];
path_calib   = [path_root,'CALIBRATION/'];
path_profile = [path_root,'MASSDIMM/'];

%% INITIALIZE THE TRS CLASS
parFileKeck_NGS();
trs = telemetry(parm,date,path_data,path_im,path_calib,'path_profile',path_profile);

%% INITIALIZE THE PRIME CLASS
psfr = prime(trs);

%% MODEL ADJUSTMENT
psfr = psfr.bestFitting(trs.im_sky,'weighting',true,'fitGains',[true,true,true],'fitCn2',false,'flagJacobian',true,'nZernMode',[4:6]);
psfr.displayResults(psfr.im_sky);

%% IMAGE-ASSISTED AO ERROR BREAKDOWN
errorBreakDown(psfr);
