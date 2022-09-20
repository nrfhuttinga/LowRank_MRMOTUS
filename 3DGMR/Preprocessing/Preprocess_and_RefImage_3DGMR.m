%% Script to preprocess raw k-space data acquired and reconstruct a reference image
%% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

clear all;

%% Set paths

lowrank_mrmotuspath = '/nfs/arch11/researchData/USER/nhutting/code/LowRank_MRMOTUS/';
cd(lowrank_mrmotuspath)
addpath(genpath(pwd));

%% Load data

data_folder  = [lowrank_mrmotuspath,'Data/3DGMR/'];
load([data_folder,'DataStruct_raw.mat']);
ref_im_parameters_3DGMR;

%% Preprocess and reconstruct reference image
ref_im_parameters.parallel_reconstruction = 0;
[DataStruct,HighresReferenceImage]=Preprocess_and_RefImage(DataStruct,ref_im_parameters);
% Save results in the new DataStruct_processed struct that will be used in the MR-MOTUS recons
disp('Saving...')
save([data_folder,'DataStruct_processed.mat'],'DataStruct');
save([data_folder,'HighresReferenceImage.mat'],'HighresReferenceImage');
disp('Done!')
