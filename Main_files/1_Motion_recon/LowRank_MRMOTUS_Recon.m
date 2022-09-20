function [Phi,Psi,MRMOTUS_recon, param_struct]=LowRank_MRMOTUS_Recon(DataStruct_processed,param_struct)
%% Function to reconstruct motion-field components 'Phi' and 'Psi', given preprocessed 'DataStruct' and parameters in 'param_struct'
%
% DataStruct_processed should contain the following fields:
%
%       .RawKspaceData                          -> [ #readoutsamples    #readouts                  ]                    % Coil-combined single-channel k-space data
%       .Coordinates                            -> [ 2/3                #readoutsamples  #readouts ]                    % K-space trajectory coordinates
%       [.SelfNavigator.SurrogateSignal]        -> [ #readouts ]                                                        % Respiratory surrogate signal, you can provide your own, required for resp-resolved mode
%
% See "/LowRank_MRMOTUS/3DGMR/Motion recon/Parameters_3Dt_RespMotion.m" and "/LowRank_MRMOTUS/2DGA/Motion recon/Parameters_2Dt_RespMotion.m" for all required parameters and explanation.
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.


%% Calibrate and sort data according to parameters


[DataStruct_processed.ReferenceImage,DataStruct_processed.RawKspaceData] = CalibrateReferenceAndKdata(DataStruct_processed.ReferenceImage,DataStruct_processed.RawKspaceData,DataStruct_processed.Coordinates,param_struct);

[DataStruct_processed, param_struct] = SortData(DataStruct_processed,param_struct);

% export_suffix = param_struct.export_suffix;

    
DataStruct_processed.Coordinates = demax(DataStruct_processed.Coordinates)/2;

%% Init MR MOTUS operator

MRMOTUS_recon           = MRMOTUS_Operator((DataStruct_processed.ReferenceImage(:)),permute(DataStruct_processed.Coordinates,[2 1 3]),param_struct);
% MRMOTUS_recon.RegularizationFlag = 0;

%% Actual reconstructions

% Handle to evaluate forward model and gradients at iterate 'x', this is
% required for lbfgs
f_handle = @(x) MRMOTUS_recon.forward_and_gradient_lowrank(x,DataStruct_processed.RawKspaceData);

      

        
disp('=== Running LBFGS-B Reconstructions ===')
% set LBFGS options
clearvars options
options.x0      = MRMOTUS_recon.SolutionVariables_init(:);
options.lb      = options.x0*Inf;
options.ub      = options.x0*Inf;
options.factr   = param_struct.lbfgs_termination_threshold;
options.maxIts  = param_struct.NumberOfReconIterations;
options.m       = 6;
options.plotting= param_struct.VisualizationFlag;
options.errFcn  = {@(x) x(:),@(x) toc};

% Run recons
tic
[dvf,~,info]=lbfgsb(f_handle,options.lb,options.ub,options);



% Export results
disp('+Saving some reconstruction results...');
save([param_struct.export_folder,'dvf',param_struct.export_suffix,'.mat'],'dvf','-v7.3')
save([param_struct.export_folder,'recon_info',param_struct.export_suffix,'.mat'],'info');

disp('+Constructing low-rank motion-field components from coefficients')
[Phi,~,Psi,DVF_ops] = MRMOTUS_recon.ExpandMotionfieldCoefficients(dvf);



end
