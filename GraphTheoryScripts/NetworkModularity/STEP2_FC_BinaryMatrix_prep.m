%% Format Weighted Matrices for BCT Graph Theory Analyses

% Originally written by Jenna Adams to calculate weighted matrices but
% adapted by Miranda Chappel-Farley to calculate Binary adjacency matrices

% 2/9/2022

% The only parameters to change are the directories and subject list in Step 1

%% Notes to Myself:
% going to read up on community_louvain to determine whether or not to keep
% negative edge weights (losing a lot of data to drop them).

%% Step 1. Set directories and load subject list mat file.

HomeDir = '/Volumes/yassamri/BEACoN/ProcessedData/rsfMRI/For_CONN/Beacon_Miranda/BEACoN_SALSA_N24/GraphTheory';
FCDir = [HomeDir, '/FC_matrix_bivariate/'];
BinDir = [FCDir, 'BinarizedMatrix/'];

FCmatrix = load('FCmatrix_151ROIs_N40.mat');
sublist = FCmatrix.sublist;

fprintf(' \n Step 1 Complete- directories set and subject list loaded');

%% Step 2. Initialize zeros matrix to hold information about negative edge weights, symmetrize matrix, save.
negatives_summary = zeros(length(sublist),2);

fprintf(' \n Step 2 Complete- Negatives Summary Matrix Initialized');
%% Step 3. loop through all the subject files
for i = 1:length(sublist)
    subject = num2str(sublist{i});
    disp(['Running ', subject]);
    load([FCDir,subject,'_FCmatrix_160ROIs.mat'])
    % make A into our functional connectivity matrix
    A = FCmatrix;
    clear FCmatrix;
    num_rois = length(A);
    % Make matrix symmetric by summing it with its transpose and dividing all elements in half. 
    % This is essentially taking the average of the corresponding elements across the diagonal, 
    % and setting each element to the average. The result is a symmetric matrix, where each element
    % in the average of the two beta-interaction terms for eachs et of  ROI-pairs
    A_sym = (A + A') / 2;
    % Ensure symmetry and set NaNs to Zero     
    A_sym = weight_conversion(A_sym, 'autofix'); 

    % Count how many negatives and save
    negatives_orig = (length(nonzeros(A_sym(A_sym<0)))/2); % divide by 2 b/c of lower half of triangle
    negatives_orig_perc = negatives_orig/((length(A)*length(A))/2); % calculate the percentage of all that were negatives
    negatives_summary(i,1) = negatives_orig; %saves number of negatives as a group summary stat
    negatives_summary(i,2) = negatives_orig_perc; % saves percentage as group summary stat

    % Remove negative values
    A_sym_pos = max(A_sym,0); % keeps everything zero and above (i.e., drops negative corrs)


    % apply costs - calculate metrics later across costs to look for stability
    costs = 0.05:0.05:0.5; % 10 costs, ranging from .5 to 50 in increments of 0.05

    % create 3D adjacency matrix for subject across costs
    adj_matrix = zeros(num_rois, num_rois, length(costs));

    % keep track of total number of connections per cost -each row is 0.05 increase in cost
    connections_per_cost = zeros(10,1); 
    
    % Keep track of network density per cost
    density_per_cost = zeros(10,1);

    % time to loop through costs

    for cost = 1:length(costs) % for each cost in costs 
        fprintf(1, 'Creating %f cost adjacency matrix\n', costs(cost));
        adj_sym = threshold_proportional(A_sym_pos, costs(cost)); % threshold on cost (e.g., @ 0.05 only the top 5% of edges are left)
        adj_sym_bin = weight_conversion(adj_sym, 'binarize'); %binarize based on cost
        adj_matrix(:,:,cost) = adj_sym_bin; % add the cost-adjusted matrix to 3D matrix
        connections_per_cost(cost,1) = length(nonzeros(adj_sym_bin));% calculate how remaining connections at cost
        density_per_cost(cost,1) = density_und(adj_sym_bin);
    
    end 
     %Save Binarized 3D Matrix and other structures for each Subject
     save([BinDir,subject,'_FCmatrix_BCTbinarized.mat'],'A','A_sym','A_sym_pos','names','negatives_summary','density_per_cost','adj_sym', 'adj_sym_bin', 'adj_matrix', 'connections_per_cost');
     
end 
    
save([BinDir, 'NegativeEdges_Summary.mat'], 'negatives_summary');
fprintf('\n Step 3 Complete- Symmetrical Binary Adjacency Matrices Calculated at Costs for Each Subject Saved');