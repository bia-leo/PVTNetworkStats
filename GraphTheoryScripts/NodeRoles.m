% Miranda Chappel-Farley
% The purpose of this script is to calculate several node-specific mesasures includeing:
% Node Strength

% Eigenvector Centrality - characterizes the global spread of an event
% happening at a particular node

% Participation coefficient - strenght of a nodes connections within and
% between communiities

% Node stregnth- measure of the influence of a node 

% Betweenness Centrality

%% Load in directories
clear all; clc;

HomeDir = '/Users/bianca/Mirror/GitHub/PVTNetworkStats/GraphTheoryScripts';
FCDir = [HomeDir,'/FC_matrix_bivariate'];
BinDir = [FCDir, '/WeightedMatrix/'];
CommDir = [BinDir, 'BCT_output_Step3_NetworkModularity_Weighted/'];

cd(HomeDir);
sub_info = readtable('CC1p0_ROI_CONN_n178_sublist.csv');
% Read the CSV file into a table
include_data = readtable('/Users/bianca/Mirror/GitHub/PVTNetworkStats/REST_Incude_Var.csv');

% Extract the first column
subsToInclude = include_data{:, 1};

% Reverse the values (0s become 1s and 1s become 0s)
subsToExclude = 1 - subsToInclude;

cd(FCDir);
sublist = load('FCmatrix_11ROIs_N178.mat', 'sublist');
sublist = sublist.sublist;
numSubs = length(sublist);

% Select diretory with input data (symmetric matrices)
%sym_matrix_dir = BinDir;
cd(FCDir);
sym_matrix_dir = FCDir;
% Set generic file path
sym_file_path = fullfile(sym_matrix_dir, '*BCTweighted.mat');
% Create object with all the matrix files
files = dir(sym_file_path);
clear sym_file_path

% Let's load in the first file so we can automatically set the number of nodes so we don't have to do it manually.
load(files(1).name, 'A_sym_norm');
adj_matrix = A_sym_norm;
clear A_sym_norm;
numNodes = size(adj_matrix, 1);
clear adj_matrix ;
%% Calculate Node Strengths
% -----
% Node strength (weighted network) - the sum of connectivity weights of the
% edges attached to each node - can compute positive and negative strength
% Node strength is the sum of weights of links connected to the node.
% -----
cd(BinDir);
TotalNodeWeight = zeros(numNodes,2,numSubs); 
PosStrength = zeros(numNodes, numSubs);
NegStrength = zeros(numNodes, numSubs);
NormStrength = zeros(numNodes, numSubs);
% col1 = Spos - total positive weight
% col2 = Sneg - total negative weight


for i = 1:length(files)
    if i == 128 || 154%skip file 128 for now 
        continue;
    end
    fprintf(1, 'Compiling all subject matrices: now reading in and calculating node strengths for %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load the weighted matrix
    A = A_sym_norm; 
    [Spos, Sneg, vpos, vneg] = strengths_und_sign(A);
    PosStrength(:,i) = Spos;
    NegStrength(:,i) = Sneg;
    Snorm = Spos - ((Sneg/Spos + Sneg).*Sneg);
    % normalized node strength, Rubinov & Sporns 2011
    % s* = s+ - (s-/s+ + s-)s-
    NodeStrengths(:,1, i) = vpos;
    NodeStrengths(:,2, i) = vneg;
    clear A 
end

fprintf('Node strengths calculated and compiled for all subjects: saved in NodeStrengths \n');
%% Plot the distributions of node strengths
figure
subplot(3,1,1)
histogram(PosStrength(:,:))
xlabel('Nodal Strength for Positive Weights')
ylabel('Count')
subplot(3,1,2)
histogram(NegStrength(:,:))
xlabel('Nodal Strength for Negative Weights')
ylabel('Count')
subplot(3,1,3)
histogram(NormStrength(:,2))
xlabel('Normalized Nodal Strength')
ylabel('Count')

figure
subplot(2,1,1)
histogram(TotalNodeWeight(:,1,:))
xlabel('Total Positive Weight')
ylabel('Count')
subplot(2,1,2)
histogram(TotalNodeWeight(:,2,:))
xlabel('Total Negative Weight')
ylabel('Count')


%% Eigenvector Centrality
%Eigenvector centrality accounts for the quantity and quality of a nodes
%degree -- can compute with weighted networks, either add 1 to negative
%weights, but then strong negative correlation will be given lower weight
%than a weak pos corr. If you take abs val. a strong neg correlation will
%have equal weighting as a strong pos. 

% choosing to remap by adding 1, since the meaning behind neg correlation
% is unclear, but don't want to disregard them entirely or give them teh
% same weight. (Lohmann 2010: Eigenvector centrality mapping for analyzing
% connectivity patterns in fmri data of the human brain)

evCentrality = zeros(numSubs, numNodes);
B = ones(numNodes, numNodes); % create a matrix of ones with the same size as our adj matrix
for i = 1:length(files)
    if i == 128 || 154%skip file 128 for now 
        continue;
    end
    fprintf(1, 'Compiling all subject matrices: now reading in and calculating eigenvector centrality for %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load the weighted matrix
    A = A_sym_norm; 
    Apos = A + B; % add 1 to everything so no negatives, but not entirely ignoring them
    v = eigenvector_centrality_und(Apos);
    evCentrality(i,:) = v;
end

histogram(evCentrality(:,:)) % chceking distribution
xlabel('Eigenvector Centrality')
ylabel('Count')
fprintf('Eigenvector Centrality calculated and compiled for all subjects: saved in evCentrality \n');
%% Participation coefficient - requires module assignment
mod_dir = CommDir;
% Set file path to the partitions
cd(CommDir);
mod_file_path = fullfile(mod_dir, '*FinalLouvainPartitions.mat');
% Create object with all the matrix files
mod_files = dir(mod_file_path);
load(mod_files(1).name, 'gamma');
% pre-alllocate matrix to hold participation coefficients for each level of
% gamma
levelsOfgamma = length(gamma);

Pos_pCoef = zeros(numNodes,levelsOfgamma, numSubs);
Neg_pCoef = zeros(numNodes,levelsOfgamma, numSubs);

for i = 1:length(mod_files)
    fprintf(1, 'Calculating Participation Coefficient for %s\n', mod_files(i).name);
    load(mod_files(i).name); %load in the final partition matrix for each level of gamma
    cd(BinDir);
    load(files(i).name) % load in the weighted connectvity matrix
    cd(CommDir);
    communities = Dfinal; % the community affiliation vectors at each level of gamma
    W = A_sym_norm; % undirected connection matrix with pos and neg weights
    for i_gamma = 1:levelsOfgamma
        Ci = communities(:,i_gamma); % community affiliation vector
        [Ppos, Pneg] = participation_coef_sign(W,Ci);
        Pos_pCoef(:,i_gamma,i) = Ppos;
        Neg_pCoef(:,i_gamma,i) = Pneg;
    end
    
end

figure
subplot(2,1,1)
histogram(Pos_pCoef(:,1));
xlabel('Participation Coeff for Positive Weights')
ylabel('Count')
subplot(2,1,2)
histogram(Neg_pCoef(:,1)); % chceking distribution
xlabel('Participation Coeff for Negative Weights')
ylabel('Count')



fprintf('Positive and Negative Participation Coeffs calculated and compiled for all subjects: saved in Pos_pCoef and Neg_pCoef \n');
%% Calculate correlatio coefficient for the participation coefficient across levels of gamma for analyses - Pos_pCoef
%collapse across third dimension (num of subjects)

%Skipping this for now, going to go based on individual ROIs

pCoef = reshape(Pos_pCoef, numNodes*numSubs,levelsOfgamma);
pCoef_scaled = pCoef * 10000 ; 
pCoef_scaled_round = round(pCoef_scaled); %NMI only takes integers so we are multiplying by 100 then rounding

nmiVals = zeros(2,6);

nmiVals(1,1) = 12;
nmiVals(2,1) = nmi(pCoef_scaled_round(:,1), pCoef_scaled_round(:,2));

nmiVals(1,2) = 13;
nmiVals(2,2) = nmi(pCoef_scaled_round(:,1), pCoef_scaled_round(:,3));

nmiVals(1,3) = 14;
nmiVals(2,3) = nmi(pCoef_scaled_round(:,1), pCoef_scaled_round(:,4));

nmiVals(1,4) = 23;
nmiVals(2,4) = nmi(pCoef_scaled_round(:,2), pCoef_scaled_round(:,3));

nmiVals(1,5) = 24;
nmiVals(2,5) = nmi(pCoef_scaled_round(:,2), pCoef_scaled_round(:,4));

nmiVals(1,6) = 34;
nmiVals(2,6) = nmi(pCoef_scaled_round(:,3), pCoef_scaled_round(:,4));
 
% Gamma levels 1 and 4 are most dissimilar, so I will include those in my
% analyses
%% Load in ROI labels for analyses
%maskDir = '/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/AnalysisMask/n40_AnalysisMask'
%cd(maskDir);
%ROIlabels = readtable('n40_finalMask_ROINumsAndLabels.csv');
%% SAVE to mat file 
save([BinDir,'NodeRoles_n178.mat'],'numNodes','sublist', 'numSubs', 'PosStrength', 'NegStrength', 'NormStrength', 'TotalNodeWeight', 'evCentrality', 'Pos_pCoef', 'Neg_pCoef' , 'names');

%save recalculated pos participation coefficient - EC did not need to be
%recalc with new range of gamma
%save([BinDir,'ParticipationCoefficients_n40.mat'],'numNodes','sublist', 'numSubs', 'gamma', 'levelsOfgamma','Pos_pCoef', 'Neg_pCoef' , 'ROIlabels');
fprintf('Participation Coefficient .mat file generated');
%% Pull values for specific ROIs
% Subcortical ROIs in mask
%hipp = find(contains(ROIlabels.regionLabel, 'Hipp'));
%amg = find(contains(ROIlabels.regionLabel, 'Amyg'));

% Frontal ROIs in Mask
%ifs = find(contains(ROIlabels.regionLabel, 'IFS'));
%a44 = find(contains(ROIlabels.regionLabel, 'A44'));
%a45 = find(contains(ROIlabels.regionLabel, 'A45'));
%a14m = find(contains(ROIlabels.regionLabel, 'A14m'));

%vmPFC ROIs in mask
%A13, A11, A10
%a13 = find(contains(ROIlabels.regionLabel, 'A13'));
%a11 = find(contains(ROIlabels.regionLabel, 'A11'));
%a10 = find(contains(ROIlabels.regionLabel, 'A10'));

%Control ROIs in mask: Occipital lobe
%vmPos = find(contains(ROIlabels.regionLabel, 'vmPOS'));

%Entorhinal Cortex and temporal lobe control (parahippocampal)
%a35 = find(contains(ROIlabels.regionLabel, 'A35')); %parahippocampal
%a28 = find(contains(ROIlabels.regionLabel, 'A28')); %entorhinal cortex

%Compile all indices of ROIS%
%ROIS = vertcat(hipp, amg, ifs, a44, a45, a14m, vmPos, a13, a11, a10, a35, a28);
%ROIS = sort(ROIS); 

centrality_ROIS = evCentrality(:,ROIS);
Spos_ROIS = PosStrength(ROIS,:)';
Pcoef_ROIS = Pos_pCoef(ROIS,:,:);

%grab the least similar calc participation coefficients according to NMI
Pcoef_ROIS_gamma1 = squeeze(Pcoef_ROIS(:,1,:))'; 
Pcoef_ROIS_gamma4 = squeeze(Pcoef_ROIS(:,4,:))'; 

%pick median gamma value - median is 0.0996
Pcoef_ROIS_medianGamma = squeeze(Pcoef_ROIS(:,3,:))';
selectedROIS = ROIlabels(ROIS,:); %table of selected ROIS with labels
%% SAVE to mat file 

labels = selectedROIS.regionLabel'; % pull labels

centrality = array2table(centrality_ROIS);
centrality.Properties.VariableNames = labels; 
centrality.beacon_id = sublist;

Spos = array2table(Spos_ROIS);
Spos.Properties.VariableNames = labels;
Spos.beacon_id = sublist; 

Pcoef_mg = array2table(Pcoef_ROIS_medianGamma);
Pcoef_mg.Properties.VariableNames = labels;
Pcoef_mg.beacon_id = sublist;

% Pcoef_g1 = array2table(Pcoef_ROIS_gamma1);
% Pcoef_g1.Properties.VariableNames = labels;
% Pcoef_g1.beacon_id = sublist;
% 
% Pcoef_g4 = array2table(Pcoef_ROIS_gamma4);
% Pcoef_g4.Properties.VariableNames = labels;
% Pcoef_g4.beacon_id = sublist;

%To save just participation coeff
%save([BinDir,'ParticipationCoeff_37ROIs_n40.mat'],'sublist', 'Pcoef_mg', 'selectedROIS');

save([BinDir,'NodeRoles_37ROIs_n40.mat'],'sublist', 'numSubs','centrality', 'Spos', 'Pcoef_mg', 'selectedROIS');
fprintf('Node Roles for selected ROIS saved in .mat file: NodeRoles_37ROIs_n40.mat ');
%% Write to csv files for analyses
cd(HomeDir);
% convert total group Eff values to csv
writetable(centrality, 'EigenVectorCentrality_n40_37ROIs.csv');
writetable(Spos, 'PositiveNodeStrength_n40_37ROIs.csv');
writetable(Pcoef_mg, 'ParticipationCoefficient_medianGamma_n40_37ROIs.csv');

%% Betweenness Centrality -- a measure of information flow. Determines how central a vertex is for information propogation
% Data must be converted to connection-length matrix. Nodes with a high
% betweenness centrality participate in a large number of shortest paths.
% 
% The input matrix must be a connection-length matrix, typically
% obtained via a mapping from weight to length. For instance, in a
% weighted correlation network higher correlations are more naturally
% interpreted as shorter distances and the input matrix should
% consequently be some inverse of the connectivity matrix. 
% 
% Betweenness centrality may be normalised to the range [0,1] as
% BC/[(N-1)(N-2)], where N is the number of nodes in the network.

betweenness_centrality = zeros(numSubs, numNodes); % pre-allocate a matrix for our metrics
cd(BinDir);
for i = 1:length(files)
    fprintf(1, 'Compiling all subject matrices: now reading in and calculating betweenness centrality for %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load the weighted matrix
    A_sym_pos = max(A_sym_norm,0); % keeps everything zero and above (i.e., drops negative corrs)
    A_sym_pos_len = weight_conversion(A_sym_pos, "lengths"); %convert corrs to lengths. high corr = short length, low corr = long length
    bc = betweenness_wei(A_sym_pos_len); % calculate betweenness centrality on the connection-length metrix
    bc_norm = bc/((numNodes-1)*(numNodes-2)); %normalize to [0,1]
    betweenness_centrality(i,:) = bc_norm;
end
clear i;


betweenness_ROIs = betweenness_centrality(:,ROIS); %pull the ROIS of interest
bc_table = array2table(betweenness_centrality); %make it a table
bc_table.Properties.VariableNames = names;
bc_table.nsubid = sublist; 
%bc_table.ROI = names; 

%add labels and sub list
histogram(betweenness_centrality) % chceking distribution
xlabel('Betweenness Centrality')
ylabel('Count')

cd(HomeDir);
writetable(bc_table, 'BetweennessCentrality_n178_11ROIs.csv');



