%% Calculate Local and Global Efficiency 
% Author: Miranda Chappel-Farley, 10/6/2022
% The purpose of this script is to calculate Global Efficiency from our
% weighted matrices. To do this, we first drop negative correlation
% coefficients, then rescale corrs to lengths. Then we derive a length
% matrix to calculate global efficiency. 

%% Step 1- Set directories
% Make sure you are in the directory: /Volumes/yassamri/BEACoN/ProcessedData/rsfMRI/For_CONN/Beacon_Miranda/BEACoN_SALSA_N24/GraphTheory/FC_matrix_bivariate/WeightedMatrix
clear all; clc;

HomeDir = '/Volumes/yassamri3/SALSA_SleepStudy/BEACoN_SALSA_N40/GraphTheory';
FCDir = [HomeDir,'/FC_matrix_bivariate'];
BinDir = [FCDir, '/WeightedMatrix/'];
% Create folder to store node partitions for each subject*condition  matrix
partitionDir = OutputDir;

cd(HomeDir);
sub_info = readtable('BEACoN_SALSA_n40_sublist.csv');
subsToExclude = sub_info.to_exclude;
%sublist = sub_info.beacon_id;

cd(FCDir);
sublist = load('FCmatrix_151ROIs_N40.mat', 'sublist');
sublist = sublist.sublist;
numSubs = length(sublist);
cd(BinDir);

%% Step 2- Set paths and number of nodes 

% Select diretory with input data (symmetric matrices)
sym_matrix_dir = BinDir;
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

%% Step 3 - Compile weighted length matrices and calculate GE
% initialize empty matrix to contain all of the processed matrices for each
% subject
mat = zeros(size(adj_matrix,1),size(adj_matrix,1),length(files));
GlobalEff_CharPath = zeros(numSubs, 3); % pre-allocate to hold GE and lambda values


for i = 1:length(files)
    fprintf(1, 'Compiling all subject matrices: now reading in and calculating GE for %s\n', files(i).name);
    load(files(i).name, 'A_sym_norm'); %load the weighted matrix
    A_sym_pos = max(A_sym_norm,0); % keeps everything zero and above (i.e., drops negative corrs)
    A_sym_pos_len = weight_conversion(A_sym_pos, "lengths"); %convert corrs to lengths. high corr = short length
    [D,B] = distance_wei(A_sym_pos_len); % D - distance (shortest weighted path matrix), B- number of edges in shortest weighted path matrix
    % Characteristic path length is defined here as the mean shortest path length between all pairs of nodes,
    mat(:,:,i) = D; % put the matrix into the respective subject slice
    [lambda, GE] = charpath(D); % lambda - network characteristic path length, GE - network global eff, 
    GlobalEff(i,2) = GE;
    GlobalEff(i,3) = lambda;
end
clear i;

GlobalEff(:,1) = sublist;
%% Plot Global Efficiency
figure
histogram(GlobalEff(:,2)) % GE
figure
histogram(GlobalEff(:,3)) % lambda


cd(HomeDir);
% convert total group Eff values to csv
csvwrite('GroupGloablEffValues_n40.csv', GlobalEff);


save([BinDir,'PositiveWeightedLengthMatrix_N40_151ROIs.mat'],'data','varnames')
%% Small worldness

% small-worldness

% assortativity coefficeint --- split participants onthis, assortativity_wei



