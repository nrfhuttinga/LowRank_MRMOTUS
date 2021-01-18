


% Specify where data can be found, and recons should be exported
base_path                                   = ['/nfs/arch11/researchData/USER/nhutting/code/ForLowRankPaper/2DGA/'];
    export_folder                           = [base_path,'Exports/'];
    highres_referenceimage_path             = [base_path,'Data/HighresReferenceImage.mat'];
    DataStruct_path                         = [base_path,'Data/DataStruct_processed.mat'];

% Specify which data to extract for the recon
param_struct.NumberOfDynamics               = 100;
param_struct.BeginReadoutIdx                = 1;  
param_struct.ReadoutsPerDynamic             = 20;   % used in time-resolved recon only
param_struct.RespResolvedReconstruction     = 0;    % 0 - time-resolved recon; 1: resp-resolved recon

% Specify reconstruction parameters
param_struct.lambda_det                     = 0*1;
param_struct.lambda_TV                      = 8e-5;
param_struct.eps_TV                         = 1e-2;
param_struct.ParallelComputationFlag        = 0;                          % parallel computations or not [0 / 1]
    param_struct.NumberOfThreads            = 8;                          % number of threads used in parallel pool
    param_struct.ParallelRegularizationFlag = 1;
param_struct.NumberOfReconIterations        = 40;                         % number of iterations in reconstruction
param_struct.lbfgs_termination_threshold    = 5e11;
param_struct.ReconstructionAlgorithm        = 'LBFGS';
param_struct.PreconditionParam              = 1/2;
param_struct.LowRankReconstruction          = 1;
    param_struct.NumberOfComponents         = 3;                                           % number of components in the low-rank motion model

    
% Parameters for the spline bases
param_struct.NumberOfSpatialSplines         = [40].';                                           % number of splines per spatial dimension
param_struct.NumberOfTemporalSplines        = round(param_struct.ReadoutsPerDynamic*param_struct.NumberOfDynamics/800*5); % number of splines over the temporal dimension


param_struct.VisualizationFlag              = 1;
param_struct.HighresVisualizationFlag       = 1;