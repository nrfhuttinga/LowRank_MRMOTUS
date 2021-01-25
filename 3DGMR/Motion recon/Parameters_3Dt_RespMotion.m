% Parameters for 3D+t respiratory motion reconstructions. Copy and edit this file 
% for your own reconstructions
%
% Niek Huttinga - UMC Utrecht - 2020



% Specify where data can be found, and recons should be exported
base_path                                   = ['../LowRank_MRMOTUS/'];
    export_folder                           = [base_path,'3DGMR/Exports/'];
    highres_referenceimage_path             = [base_path,'Data/3DGMR/HighresReferenceImage.mat'];
    DataStruct_path                         = [base_path,'Data/3DGMR/DataStruct_processed.mat'];

% Specify which data to extract for the recon
param_struct.RespResolvedReconstruction     = 1;            % 0 - time-resolved recon; 1: resp-resolved recon
param_struct.NumberOfDynamics               = 10;           % Total dynamics use for reconstructions (= number of respiratory phases in resp-resolved recon)
param_struct.BeginReadoutIdx                = 1;            % From this readout and onwards the data will be used for the reconstruction
param_struct.ReadoutsPerDynamic             = 20;           % Number of readouts per dynamic (will only be used for time-resolved recon)


% Specify reconstruction parameters
param_struct.lambda_det                     = 0;            % Jacobian determinant regularization parameter
param_struct.lambda_TV                      = 5e-7;         % TV regularization parameter
param_struct.eps_TV                         = 1e1;          % epsilon for smooth TV approximation
param_struct.ParallelComputationFlag        = 1;            % Parallel computations or not [0 / 1]
    param_struct.NumberOfThreads            = 4;            % Number of threads used in parallel pool
param_struct.NumberOfReconIterations        = 60;           % Number of iterations in reconstruction
param_struct.lbfgs_termination_threshold    = 5e11;         % Threshold for reconstructions (higher threshold = lower accuracy = earlier stopping)
param_struct.PreconditionParam              = 1/2;          % Power to which to raise the precondition matrix (don't touch)
param_struct.NumberOfComponents             = 1;            % Number of components in the low-rank motion model

    
% Parameters for the spline bases
param_struct.NumberOfSpatialSplines         = [20].';                                                                       % Number of splines per spatial dimension
param_struct.NumberOfTemporalSplines        = round(param_struct.ReadoutsPerDynamic*param_struct.NumberOfDynamics/800*5);   % Number of splines over the temporal dimension

% Visualization flags
param_struct.VisualizationFlag              = 1;            % Flag to specify visualizations during LBFGS reconstructions
param_struct.HighresVisualizationFlag       = 1;            % Flag to specify whether to upscale motion-fields after the reconstruction for visualization