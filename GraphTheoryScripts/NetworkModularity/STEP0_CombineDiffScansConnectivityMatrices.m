%Bianca Leonard
%btleonar@uci.edu
%adapted from Miranda Chappel-Farley, PhD GraphTheoryScripts
%% Combine rsfMRI connectivity matrices for Conte Center 1.0 Data
clear all; 
%% Step 1 - set directories
HomeDir = '/Users/bianca/Mirror/GitHub/PVTNetworkStats/GraphTheoryScripts';
CONNDir = '/Users/bianca/Mirror/Yassa_Lab/Projects/Conte Center/CC_1p0_Data/FirstLevelDataCopy/ConteCenter/Longitudinal_Study_8_1_23/first_level_results/RRC_01_semipartial'

%% Step 2 - Load in files and save subjects matrices with identifier for scan type

load([CONNDir,'/resultsROI_Condition001.mat']); % FC matrix here

%% Step 3 - combine the three 3D matrices containing the connectivity matrices for each scan type
CombinedScans = Z;

%% Step 4 - create subject list
cd(HomeDir);
sub_info = readtable('CC1p0_ROI_CONN_n178_sublist.csv');
sublist = sub_info.nsubid;
connlist = sub_info.conn_id; 

%% Step5 - Save combined file
cd(HomeDir);
save('resultsROI_Condition001.mat','CombinedScans', 'sublist', 'connlist', 'SE','DOF','names', 'names2', 'regressors', 'xyz');

