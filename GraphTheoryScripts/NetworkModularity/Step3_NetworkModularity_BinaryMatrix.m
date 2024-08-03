%% Calculate network modularity across costs

% This script was originally written by Jesse Yaros, then adapted by Jenna Adams, and finally
% adapted by Miranda Chappel-Farley for Binarized Adjacency Matrices to run over a series of 
% costs and resolution parameters. Lastly to compare to a uniform null model at the suggestion
% of Rick Betzel, as default null model in community_louvain is not suitable for correlational data. 

% Script is based on methods from Lacichinetti & Fortuanto 2012 and Cohen & D'Espocito 2016. 

%% Step 1- Set directories
HomeDir = '/Volumes/yassamri/BEACoN/ProcessedData/rsfMRI/For_CONN/Beacon_Miranda/BEACoN_SALSA_N24/GraphTheory/';
FCDir = [HomeDir,'/FC_matrix_bivariate/'];
BinDir = [FCDir, 'BinarizedMatrix/'];

sublist = [1001 1004 1006 1007 1008 1009 1010 1012 1013 1021 1023 1032 1033 1034 1036 1041 1042 1046 1047 1048 1054 1059 1065 1069];
sublist = num2cell(sublist);
numSubs = length(sublist);

% Create folder to store node partitions for each subject*condition  matrix
partitionDir = [HomeDir,'/BCT_output_Step3_NetworkModularity/'];
% mkdir(partitionDir);
fprintf('Step 1 complete- necessary directories set and subject list created. \n NOTE: Adjust free parameters in Step 2 before running. \n');

%% Step 2- Load in files, select total number of partitions

%  -------FREE PARAMETERS- ADJUST ACCORDINGLY -------
% Set number of Louvain partitions to create - currently this is based on Cohen & D'Espocito
numPartitions = 150; % select  number of desired partitions
gamma = 0.1:0.1:1.5; % 15 possible values of gamma, ranging from .1 to 1.5 in increments of 0.1
%  -------------------------------------------------

% Select diretory with input data (symmetric matrices)
sym_matrix_dir = BinDir;
% Set generic file path
sym_file_path = fullfile(sym_matrix_dir, '*BCTbinarized.mat');
% Create object with all the matrix files
files = dir(sym_file_path);
clear sym_file_path

% Let's load in the first file so we can automatically set the number of nodes so we don't have to do it manually.
load(files(1).name, 'adj_matrix');
numNodes = size(adj_matrix, 1);

% For a given matrix, generate an empty matrix P of possible Louvain node partitions across costs.
% P will have dimensions = (node# x partition# x cost#)
numCosts = size(adj_matrix, 3);
P = zeros(numNodes, numPartitions, numCosts);

% For a given matrix, generate an empty matrix Qstat of possible Q values across partitions and costs.
% Q stat will have dimensions = (partition# x cost#)
Qstat = zeros(numPartitions, numCosts);
clear adj_matrix % we can clear adj_matrix since we load each file later in the loop in step 3

% Pre-allocate empty cells with proper dimensions to hold each subjects
% partitions and Q values across costs while varying the level of gamma
levelsOfGamma = length(gamma);

% Pre-allocate empty matrix with #ROIx#ROI dimensions for each cost
% Pre-allocate cells with proper dimensions to hold the agreement matrices for each subject for each cost/level of gamma
D = zeros(numNodes, numNodes, numCosts);
AgreementMatrices = num2cell(zeros(numSubs,levelsOfGamma)); 


% pre-allocate
% for each subject (row), create an empty cell that will contain a matrix of adjancency
% matrices (150 partitions at each cost or Q at each cost) at a given level of gamma (column)
for subj=1:numSubs
    for lGamma = 1:levelsOfGamma
        VaryingAdjMatrices{subj,lGamma}  = num2cell(zeros(numNodes, numPartitions, numCosts)); %holds partitions
        QstatMatrix{subj,lGamma} = num2cell(zeros(numPartitions, numCosts())); % holds Q 
    end
end 
clear subj lGamma;

fprintf('Step 2 Complete- Files loaded, parameters set, cells and matrices pre-allocated for step 2 \n');
 
%% Step 3- Loop through the files and calculate Modularity and generate module structure for each partition for each cost

% Cycle through each file
for index = 1:length(files)
  % Assign file name to variable  
  fileName = files(index).name;
  % Attach path to file name
  fullFileName = fullfile(sym_matrix_dir, fileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  % Load file into matrix A
  load(fullFileName,'adj_matrix');
  A = adj_matrix;
  
  % We need to iterate through different resolution parameters. When
  % gamma (resolution parameter) < 1, larger communities are resolved,
  % whereas gamma > 1 yeilds more communities containing fewer nodes. One
  % method for selecting the resolution parameter is to  report community
  % structure at the value of gamma at which partitions are most similar to
  % each other (Bassett et al. 2013 & Sporns & Betzel 2016). This can be
  % computed using normalized mutual information (use the z-score). As a
  % benchmarch, Cohen & D'Esposito 2016 used gamma = 1.25. According to
  % Fortunato & Barthelemy 2007 " a value of Q larger than 0.3–0.4 is a clear 
  % indication that the subgraphs of the corresponding partition are
  % modules."
  
  % ----------------------- 
  % need to add rick betzels uniform null model into community_louvain
  % once this is working: ci = community_louvain(ones(N),[],[],(A - gamma).*~eye(length(A)));
  % Larger values of gamma yeild smaller communities
  ----------------------- 
  
  for i = 1:levelsOfGamma % for each level of gamma
     for iter = 1:numPartitions % for each individual partition
      for cost = 1:numCosts % at each cost
        [M,Q1] = community_louvain(A(:,:,cost),gamma(i),[],'modularity'); % calculate partition and modularity at a given level of gamma 
        P(:,iter,cost) = M; % save partition into P for partision # and cost at a given level of gamma
        Qstat(iter,cost) = Q1; % save Q into Qstat for each partition and cost at a given level of gamma
        VaryingAdjMatrices{index,i} = P; % add p into VaryingAdjMatrices so we can see how community detection changes across gamma for each subject
        QstatMatrix{index,i} = Qstat; % add Qstat into QstatMatrix so we can see how Q changes across gamma for each subject
      end 
     end 
  end
  
  clear i iter cost;
  
  % --------- Left off here --------- 
  
  % Now need to generate an agreement matrix for each cost at a given value
  % of gamma (each cell in VaryingAdjMatrices, which each has 10 layers)
  % create cell of similar size
  
  % Need to create a consensus matrix across each layer of 150 partitions.
  % Then, determine the optimal cost and level of gamma based on the
  % conensus matrix at each cost. Cohen & D'Espocito thresholded each cell
  % in the consensus matrix to 0.5, so that the agreement of all pairs that
  % were not assigned to the same matrix at least 50% of the time were set
  % to 0. Then, they created a consensus partition by running the Louvain
  % algorithm on the agree ment matrix 100 times until a single
  % representative partition was obtained.
  
  for i = 1:levelsOfGamma
    for j = 1:cost
      matrix = VaryingAdjMatrices{index,i}; % pulling one column and not entire structure within cell
      D(:,:,j) = agreement(matrix(:,:,j)); % put agreement matrix for given cost in D
      AgreementMatrices{index,i} = D; % ERROR IS HERE
    end
  end 
  
  % size(VaryingAdjMatrices{7,14})
  
  % ----------------------- 
  
  % Create agreement matrix D whose elements indicate the number of times 
  % any two vertices were assigned to the same class. Should have same
  % dimensions as original connectivity matrix (ie numNodes x numNodes).
  % This function is from the BCT toolbox
  D = agreement(P); % this is for one subject
  
  % Convert agreement matrix D to proportion of agreement matrix Dp .
  % This line converts counts in the agreement matrix to probabilities. 
  % Each element in Dp will have a value equal to or less than the total 
  % number of partitions. We do this by dividing matrix D by the scalar
  % total number of partitions. This converts each element into a proportion.
  Dp = D/numPartitions;  % this is for one subject
  
  
  %Use BCT consensus function to run the agreement matrix through the
  %louvain algorithm again. Set threshold to .3 Lancichinetti & Fortuando 
  %2012 recommend a low threshold --.3 or below -- see figure S2. 
  %One thing to consider -- consensus_und implements louvain with gamma =
  %1. (the default)For now I am not altering function, but I don't know if
  %it is weird that the intial partitions were create with gamme of 1.25,
  % while the concensus partitions are created with gamma set to 1....
  
  %Also note, I followed Cohen and D'Espocito for gamma selection in the
  %original louvain partitions. However I am deviating from them here. They
  %used tau/threhsold = .5 which is a bit high based on Lancichinett and
  %Fortuanto 2012 recommendations
  
  final_clusters = consensus_und(Dp,.3,100);  % this is for one subject
  
  subject = char(sublist(index));  % this is for one subject
  
  %Save partitions and final clusters
  save([partitionDir,subject,'_LouvainPartitions.mat'],'P','D','Dp','final_clusters');
  
  % P is each partition (150), D is agreement matrix, Dp is percentage of
  % of modules assigned to same module across particition, 'final_clusters'
  % is the optimal clustering. 

end

save([partitionDir, 'ModularityStatistic_LouvainPartitions.mat'], 'Qstat');
%Qstat is the modulairty for each partition

%% to plot modularity distributions
figure
histogram(Qstat(1,:))
hold on
histogram(Qstat(2,:))
histogram(Qstat(3,:))
histogram(Qstat(4,:))
histogram(Qstat(5,:))
histogram(Qstat(6,:))
histogram(Qstat(7,:))
histogram(Qstat(8,:))
histogram(Qstat(9,:))
histogram(Qstat(10,:))
histogram(Qstat(11,:))
histogram(Qstat(12,:))
histogram(Qstat(13,:))
histogram(Qstat(14,:))
histogram(Qstat(15,:))
histogram(Qstat(16,:))
histogram(Qstat(17,:))
histogram(Qstat(18,:))
legend('SALSA002', 'SALSA008', 'SALSA009', 'SALSA010', 'SALSA011', 'SALSA013', 'SALSA014', 'SALSA015', 'SALSA021', 'SALSA027', 'SALSA028', 'SALSA030', 'SALSA032', 'SALSA033', 'SALSA038', 'SALSA046', 'SALSA050', 'SALSA052')
%PC = participation_coef(A, final_clusters);

% need to calculate Q for varying measures of gamma (Resolution parameter)
% to find optimal modularity value for each subject
%% to find most frequent modularity value over the 150 partitions and save
FinalQ = array2table(zeros(length(sublist), 2));
FinalQ.FinalQ1 = sublist.';

for i = [1:length(sublist)];
    Qmode = mode(Qstat(i,:));
    FinalQ.FinalQ2(i) = Qmode;
end

save([partitionDir, 'FinalModularityStatistic_ModeLouvainPartition']);