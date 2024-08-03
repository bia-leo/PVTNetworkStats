%% Network Visualization

% @author: Miranda Chappel-Farley; mgchappe@uci.edu | mgchappelfarley@gmail.com
% 2/7/2022

% This script is adapted from the Univeristy of Utah R25 Advanced
% Statistical Methods in Neuroimaging and Genetics Course. 
% Please cite accordingly. 
%% Load in directories
clear all; clc;

HomeDir = '/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/GraphTheory';
FCDir = [HomeDir,'/FC_matrix_bivariate'];
BinDir = [FCDir, '/WeightedMatrix/'];
CommDir = [BinDir, 'BCT_output_Step3_NetworkModularity_Weighted/'];
FigureDir = [HomeDir, '/Figures/CommunityAssignmentMatrices/'];

cd(HomeDir);
sub_info = readtable('BEACoN_SALSA_n40_sublist.csv');
sub_info.to_exclude(42) = 1;
subsToExclude = sub_info.to_exclude;

%sublist = sub_info.beacon_id;

cd(FCDir);
sublist = load('FCmatrix_151ROIs_N40.mat', 'sublist');
sublist = sublist.sublist;
numSubs = length(sublist);

% Select diretory with input data (community affiliation vectors)
cd(CommDir);
% Set generic file path
comm_file_path = fullfile(CommDir, '*FinalLouvainPartitions.mat');
% Create object with all the matrix files
files = dir(comm_file_path);
clear comm_file_path

%Load Weighted Matrices
cd(BinDir);
mat_file_path = fullfile(BinDir, '*BCTweighted.mat');
matFiles = dir(mat_file_path);

% Let's load in the first file so we can automatically set the number of nodes so we don't have to do it manually.
cd(CommDir);
load(files(1).name, 'numNodes', 'gamma');
levelsOfGamma = length(gamma);

%% Create individual subject community affiliation vector matrix figures
for i = 1:length(files)
    cd(CommDir);
    fname = matFiles(i).name;
    sub = fname(1:4);
    fprintf(1, 'Generating Affiliation Matrix for %s\n', sub);
    load(files(i).name, 'Dfinal'); %load in the final partition matrix for each level of gamma
    communities = Dfinal; 
    cd(BinDir);
    load(matFiles(i).name, 'A_sym_norm');% load in the weighted connectvity matrix
    CIJctx = A_sym_norm;
    cd(FigureDir);
    for i_gamma = 1:levelsOfGamma % for each level of gamma
        Ci = communities(:,i_gamma); % community affiliation vector
        [gx,gy,idxsort] = grid_communities(Ci);
        f = figure('units','inches','position',[2,2,4,4]);
        imagesc(CIJctx(idxsort,idxsort));
        hold on; 
        plot(gx,gy,'r', 'linewidth', 2);
        colormap(parula(length(unique(CIJctx))));
        colorbar;
        xlabel('Nodes (ROIs)'); ylabel('Nodes (ROIs)');
        axis square;
        filename = [sub, '_communityMatrix_gamma', num2str(gamma(i_gamma)), '.png'];
        saveas(f, filename);
        clear Ci f filename
        close all;
    end
    
end

fprintf('\n Community Affiliation Matrix Figures generated and saved \n');
%% Determine whether hipp and amygdala are always grouped in the same community
cd(BinDir);
load('NodeRoles_n40.mat', 'ROIlabels');
hipp = find(contains(ROIlabels.regionLabel, 'Hipp'));
amg = find(contains(ROIlabels.regionLabel, 'Amyg'));
amy_hipp_idx = vertcat(amg, hipp);
cd(CommDir);

% pre-allocate matrix to hold categorizations as to whether hippocampus and amygdala are
% within or not within the same community
hippAssignment = zeros(numSubs, levelsOfGamma);
amyAssignment = zeros(numSubs, levelsOfGamma);
modAssignment = zeros(numSubs, levelsOfGamma);

for i = 1:length(files)
    fname = files(i).name;
    sub = fname(1:4);
    fprintf(1, 'Checking hippocampal and amygdala community assignment for %s\n', sub);
    load(files(i).name, 'Dfinal'); %load in the final partition matrix for each level of gamma
    comms = Dfinal;
    for i_gamma = 1:levelsOfGamma
        fprintf(1, 'Checking gamma level %s\n', num2str(i_gamma));
        amg_hipp_comms = comms(amy_hipp_idx, i_gamma); % lists out the modules at that level of gamma
        if range(amg_hipp_comms(1:4)) ~= 0 % if the communitity affiliation vectors are not all the same, range will be >1
            hippAssignment(i,i_gamma) = 1; % if hipp ROIs arent all in same mod, put a 1
        else
            hippAssignment(i,i_gamma) = 0; % otherwise put a 0
        if range(amg_hipp_comms(5:8)) ~= 0 % if amygdala ROIS arent all in same mod
            amyAssignment(i,i_gamma) = 1; % put a 1
        else
            amyAssignment(i,i_gamma) = 0; % otherwise put 0
        end
        end
        % now to compare the two vectors 
        % take the mode of each vector so we know whether, more often or
        % not,thee nodes are assigned to the same modules
        if length(unique(amg_hipp_comms(1:4))) ~= length(unique(amg_hipp_comms(5:8)))
            modAssignment(i, i_gamma) = 1;
        else
            modAssignment(i,i_gamma) = 0;
            
        end   
    end
    clear fname sub comms 
end 

% Calculate the mode of the community assignment for each to use for
% analyses

%pre-allocate
mode_modAssignment = zeros(numSubs,1);
mode_hippAssignment = zeros(numSubs,1);
mode_amyAssignment = zeros(numSubs,1);

for i = 1:numSubs
    mode_modAssignment(i) = mode(modAssignment(i,:));
    mode_hippAssignment(i) = mode(hippAssignment(i,:));
    mode_amyAssignment(i) = mode(amyAssignment(i,:));
end

% Add the mode to the dataframes
modAssignment(:,levelsOfGamma + 1) = mode_modAssignment;
hippAssignment(:,levelsOfGamma + 1) = mode_hippAssignment;
amyAssignment(:,levelsOfGamma + 1) = mode_amyAssignment;

modAssignment(:,levelsOfGamma + 2) = sublist; % This contains the final values for analyses on whether hipp and amg are in the same community 
hippAssignment(:,levelsOfGamma + 2) = sublist;
amyAssignment(:,levelsOfGamma + 2) = sublist;
%%
cd(BinDir);
save([BinDir,'HippAmyg_CommunityAssignment_n40.mat'],'sublist', 'numSubs','ROIlabels', 'levelsOfGamma', 'numNodes', 'modAssignment', 'amyAssignment', 'hippAssignment', 'mode_modAssignment', 'mode_hippAssignment', 'mode_amyAssignment');
fprintf('Hipp Amygdala Assigment Matrix saved in .mat file: HippAmyg_CommunityAssignment_n40.mat ');
%% Write to csv files for analyses - MANUALLY CHANGE COLUMN NAMES
cd(HomeDir);
table = array2table(modAssignment);
colnames = ["agreement_gamma0", "agreement_gamma0_05", "agreement_gamma0_1", "agreement_gamma0_15", "agreement_gamma0_20", "agreement_gamma0_25", "agreement_gamma0_30", "mode_HippAmyCommAgreement", "beacon_id"];
table.Properties.VariableNames = colnames;
writetable(table, 'HippocampusAmygdalaCommunityAssignmentAgreement_n40_7gamma.csv');



