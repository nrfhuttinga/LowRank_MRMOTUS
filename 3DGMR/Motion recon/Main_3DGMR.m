%% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.


restoredefaultpath;addpath(genpath('../LowRank_MRMOTUS'));
close all;
clear all;
rng(1)

%% 1) load parameters and data

disp('=== Loading parameters and data ===');
% load parameters and data
Parameters_3Dt_RespMotion
load(param_struct.DataStruct_path)

%% 2) Reconstruct motion-fields

[Phi,Psi,MRMOTUS_recon,export_suffix]=LowRank_MRMOTUS_Recon(DataStruct,param_struct);
param_struct.export_suffix = export_suffix;

%% 3) Visualization of results

disp('=== Postprocessing... ===')
% Constructing motion-fields from low-rank components:
% [Note: may take a large amount of RAM in 3D+t]
mf = reshape(Phi*Psi.',[],MRMOTUS_recon.NumberOfSpatialDims,MRMOTUS_recon.NumberOfDynamics);
MRMOTUS_PostProcessing(mf, param_struct);

%%

disp('=== DONE! ===')

