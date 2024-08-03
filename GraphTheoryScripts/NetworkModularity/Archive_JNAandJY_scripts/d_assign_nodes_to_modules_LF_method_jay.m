%Use this script to assign nodes to modules. This script requires installation and use of functions from the Brain
%Connectivity Toolbox. The script is based on methods from Lacichinetti & Fortuanto 2012 and Cohen & D'Espocito 2016. 

%Add Brain Connectivity Toolbox
addpath(genpath('/Volumes/yassapublic/Users/Miranda/SALSA/code/BCT'));
ResultsDir = '/Volumes/yassapublic/Users/Miranda/SALSA/data/imaging/BEACoN_SALSA18/R21PrelimData/FC_matrix_bivariate'

%Create a folder to store node partitions for each subject*condition matrix
if not(isfolder('STEP2_d_assign_nodes_to_modules'))
    mkdir('STEP2_d_assign_nodes_to_modules')
end

%select directory to save output to
partitionDir = uigetdir

%select diretory with input data (symmetric matrices)
sym_matrix_dir = uigetdir
%set generic file path
sym_file_path = fullfile(sym_matrix_dir, '*symmetric.csv')
%create object with all the matrix files
files = dir(sym_file_path)

%cycle through each file
for index = 1:length(files)
  %assign file name to variable  
  fileName = files(index).name;
  %attach path to file name
  fullFileName = fullfile(sym_matrix_dir, fileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  %load file into matrix A
  A = readmatrix(fullFileName);
  
  %For matrix A, generate a matrix L of possible Louvain node partitions.
  %L will have dimensions = (node# x partition#)
  
  %set number of Louvain partitions to create
  numPartitions = 150 ;
  %set number of rows using node number
  numNodes = length(A);
  %initialize empty matrix P with proper dimensions
  P = zeros(numNodes, numPartitions);
 
  %Iterate as many times as numPartitions. 
  %For each iteration call community_lovain function and set the
  %corresponding columnNumber in matrix P to equal the new partition. 
  %When complete, each column of P will be a different partition. 
  
  %Chose gamma = 1.25 which Coehn & D'Espocito also used. Forms about 4
  %clusters. Not sure how to decide how high resolution we want...
  for iter = 1:numPartitions
      P(:,iter) = community_louvain(A,1.25,[],'negative_sym');
  end
  
  %Create agreement matrix D whose elements indicate the number of times 
  %any two vertices were assigned to the same class. Should have same
  %dimensions as original connectivity matrix (ie numNodes x numNodes).
  %This function is from the BCT toolbox
  D = agreement(P);
  
  %Convert agreement matrix D to proportion of agreement matrix Dp .
  %This line converts counts in the agreement matrix to probabilities. 
  %Each element in Dp will have a value equal to or less than the total 
  %number of partitions. We do this by dividing matrix D by the scalar
  %total number of partitions. This converts each element into a proportion.
  Dp = D/numPartitions;
  
  
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
  
  final_clusters_LF_3 = consensus_und(Dp,.3,100);  
  final_clusters_LF_5 = consensus_und(Dp,.5,100);  
  
  
  
  
  %Create file name for matrix of partitions
  %%partitionMatrix_FileName = fullfile(partitionDir, erase(fileName,"_symmetric.csv") + "_LouvainPartitions.csv");
  %Export matrix to csv
  %%writematrix(P, partitionMatrix_FileName);
end

%SURPRISED THIS GENERATES ONLY THREE MODULES... ISNT THAT WAY TOO SMALL?
%SHOULD I READ BRAINNETOME TO SEE IF IF IT HAD DIFFERENT 'NETWRORKSS...."
%figure out how many functional networks brainnetome found, if any, and try
%to adjust gamma to be consistent with that? 3 modules seems incredible
%small...
