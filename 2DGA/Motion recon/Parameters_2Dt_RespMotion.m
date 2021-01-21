% Parameters for 2D+t respiratory motion reconstructions. Copy and edit this file 
% for your own reconstructions
%
% Niek Huttinga - UMC Utrecht - 2020


% Specify where data can be found, and recons should be exported
base_path                                   = ['/nfs/arch11/researchData/USER/nhutting/code/LowRank_MRMOTUS/'];
    export_folder                           = [base_path,'2DGA/Exports/'];
    highres_referenceimage_path             = [base_path,'Data/2DGA/HighresReferenceImage.mat'];
    DataStruct_path                         = [base_path,'Data/2DGA/DataStruct_processed.mat'];

% Specify which data to extract for the recon
param_struct.NumberOfDynamics               = 100;          % Total dynamics use for reconstructions
param_struct.BeginReadoutIdx                = 1;            % From this readout and onwards the data will be used for the reconstruction
param_struct.ReadoutsPerDynamic             = 20;           % number of readouts per dynamic
param_struct.RespResolvedReconstruction     = 0;            % 0 - time-resolved recon; 1: resp-resolved recon


% Specify reconstruction parameters
param_struct.lambda_det                     = 0*1;          % Jacobian determinant regularization parameter
param_struct.lambda_TV                      = 1e-4;         % TV regularization parameter
param_struct.eps_TV                         = 1e-2;         % epsilon for smooth TV approximation
param_struct.ParallelComputationFlag        = 0;            % parallel computations or not [0 / 1]
    param_struct.NumberOfThreads            = 4;            % Number of threads used in parallel pool
param_struct.NumberOfReconIterations        = 100;           % Number of iterations in reconstruction
param_struct.lbfgs_termination_threshold    = 10e11;         % Threshold for reconstructions (higher threshold = lower accuracy = earlier stopping)
param_struct.PreconditionParam              = 1/2;          % Power to which to raise the precondition matrix (don't touch)
param_struct.NumberOfComponents             = 2;            % Number of components in the low-rank motion model

    
% Parameters for the spline bases
param_struct.NumberOfSpatialSplines         = [30].';                                                                       % Number of splines per spatial dimension
param_struct.NumberOfTemporalSplines        = round(param_struct.ReadoutsPerDynamic*param_struct.NumberOfDynamics/800*5);   % Number of splines over the temporal dimension

% Visualization flags
param_struct.VisualizationFlag              = 1;            % Flag to specify visualizations during LBFGS reconstructions
param_struct.HighresVisualizationFlag       = 1;            % Flag to specify whether to upscale motion-fields after the reconstruction for visualization