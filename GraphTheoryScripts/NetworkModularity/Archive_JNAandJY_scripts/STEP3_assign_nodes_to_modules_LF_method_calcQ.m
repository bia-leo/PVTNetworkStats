%Use this script to assign nodes to modules. 
%This script requires installation and use of functions from the Brain
%Connectivity Toolbox
%Script is based on methods from Lacichinetti & Fortuanto 2012 and Cohen &
%D'Espocito 2016. 


HomeDir = '/Volumes/yassapublic/Users/Miranda/SALSA/data/imaging/BEACoN_SALSA18/R21PrelimData';
FCDir = [HomeDir,'/FC_matrix_bivariate/'];

conn_values = readtable('/Volumes/yassapublic/Users/Miranda/SALSA/data/imaging/BEACoN_SALSA18/R21PrelimData/r21_restingstate_BF_BNA_MTL_connectivityValues.xlsx');
sublist = conn_values.salsaID.';

%Create folder to store node partitions for each subject*condition  matrix
partitionDir = [HomeDir,'/BCT_output_bivariate/Step3_assign_nodes_to_modules/'];
%mkdir(partitionDir);

%select diretory with input data (symmetric matrices)
sym_matrix_dir = FCDir;
%set generic file path
sym_file_path = fullfile(sym_matrix_dir, '*BCTweighted.mat');
%create object with all the matrix files
files = dir(sym_file_path);

%For matrix A, generate a matrix L of possible Louvain node partitions.
%L will have dimensions = (node# x partition#)

%set number of Louvain partitions to create - this is based on Cohen &
%D'Espocito
numPartitions = 150; % select  number
%set number of rows using node number
numNodes = 115;
%initialize empty matrix P with proper dimensions
P = zeros(numNodes, numPartitions);
% initialize empty matrix Qstat with proper dimensions - needs to be max
% and pick partition with that Q
Qstat = zeros(length(sublist), numPartitions);

%cycle through each file
for index = 1:length(files)
  %assign file name to variable  
  fileName = files(index).name;
  %attach path to file name
  fullFileName = fullfile(sym_matrix_dir, fileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  %load file into matrix A
  load(fullFileName,'A_sym_norm_pos');
  A = A_sym_norm_pos;
  
 
  %Iterate as many times as numPartitions. 
  %For each iteration call community_lovain function and set the
  %corresponding columnNumber in matrix P to equal the new partition. 
  %When complete, each column of P will be a different partition. 
  
  %Chose gamma = 1.25 which Coehn & D'Espocito also used. Forms about 4
  %clusters. Not sure how to decide how high resolution we want...
  for iter = 1:numPartitions
      [M,Q1] = community_louvain(A,0.8,[],'modularity');
      P(:,iter) = M;
      Qstat(index,iter) = Q1; 
  end
  
  %Create agreement matrix D whose elements indicate the number of times 
  %any two vertices were assigned to the same class. Should have same
  %dimensions as original connectivity matrix (ie numNodes x numNodes).
  %This function is from the BCT toolbox
  D = agreement(P); % this is for one subject
  
  %Convert agreement matrix D to proportion of agreement matrix Dp .
  %This line converts counts in the agreement matrix to probabilities. 
  %Each element in Dp will have a value equal to or less than the total 
  %number of partitions. We do this by dividing matrix D by the scalar
  %total number of partitions. This converts each element into a proportion.
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

save([partitionDir, 'FinalModularityStatistic_ModeLouvainPartitions.mat'], 'FinalQ');