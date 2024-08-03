%% Combine rsfMRI connectivity matrices for different projects (i.e., 3-min zach scans and 10-min Soyun whole brain)

%% Step 1 - set directories
HomeDir = '/Volumes/yasssamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/GraphTheory/';
CONNDir_zachScans = '/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/BEACoN_SALSA_n34_preprocessing/results/firstlevel/RRC_02_BNA151_05mmMotionThreshold'
CONNDir_soyunScans = '/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/BEACon_SALSA_n4_SoyunWholeBrain_preprocessing/results/firstlevel/RRC_02_BNA141_05mmMotionThreshold' %typo in folder name, should be 151
CONNDir_yassaScans = '/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/BEACoN_SALSA_n2_YassaWholeBrain_preprocessing/results/firstlevel/RRC_02_BNA151_05mMotionThreshold' %type here too whoops
%% Step 2 - Load in files and save subjects matrices with identifier for scan type

load([CONNDir_zachScans,'/resultsROI_Condition001.mat']); % FC matrix here
zachScanSubjects = Z; % subjects with the Zach scan resting state
zachDOF = DOF;
zachSE = SE;
clear DOF names names2 regressors SE xyz Z; % keep workspace tidy

load([CONNDir_soyunScans,'/resultsROI_Condition001.mat']); 
soyunDOF = DOF;
soyunSE = SE;
soyunScanSubjects = Z; % subjects with the soyun scan resting state
clear SE DOF; 

load([CONNDir_yassaScans,'/resultsROI_Condition001.mat']); 
yassaDOF = DOF;
yassaSE = SE;
yassaScanSubjects = Z; % subjects with the soyun scan resting state
clear SE DOF; 

%% Step 3 - combine the three 3D matrices containing the connectivity matrices for each scan type
CombinedScans = cat(3, zachScanSubjects, soyunScanSubjects, yassaScanSubjects);

%% Step 4 - create subject list
cd(HomeDir);
sub_info = readtable('BEACoN_SALSA_n40_sublist.csv');
sublist = sub_info.beacon_id;

%% Step5 - Save combined file
cd(HomeDir);
save('resultsROI_Condition001_n40_CombinedZachSoyunYassaScans_05mmMotionThreshold.mat','CombinedScans', 'sublist', 'zachDOF', 'zachSE','soyunDOF','soyunSE','yassaDOF', 'yassaSE','names', 'names2', 'regressors', 'xyz');

