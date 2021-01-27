function [Phi,Psi,MRMOTUS_recon,export_suffix]=LowRank_MRMOTUS_Recon(DataStruct_processed,param_struct)
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


%% calibrate and cut data according to parameters



[DataStruct_processed.ReferenceImage,DataStruct_processed.RawKspaceData] = CalibrateReferenceAndKdata(DataStruct_processed.ReferenceImage,DataStruct_processed.RawKspaceData,round(size(DataStruct_processed.RawKspaceData,1)/2+1));


DataStruct_processed.RawKspaceData                    = DataStruct_processed.RawKspaceData(:,param_struct.BeginReadoutIdx:end,:,:);
DataStruct_processed.Coordinates                      = DataStruct_processed.Coordinates(:,:,param_struct.BeginReadoutIdx:end);
DataStruct_processed.SelfNavigator.SurrogateSignal    = DataStruct_processed.SelfNavigator.SurrogateSignal(param_struct.BeginReadoutIdx:end);

 % Check trajectory
disp('=== Plotting trajectory ===');
figure;PlotTrajectory(DataStruct_processed.Coordinates(:,:,1:100))


%% sort the data in case of respiratory-resolved reconstruction

disp('=== Sorting data ===');

if param_struct.RespResolvedReconstruction
    % perform phase binning
    binning_pars.surrogate_signal           = DataStruct_processed.SelfNavigator.SurrogateSignal;
    binning_pars.binning_strategy           = 'phase';
    binning_pars.thresh                     = 0.005;
    binning_pars.resp_phases                = param_struct.NumberOfDynamics;
    binning_pars.return_extreme_phases      = 0;

    % sort the data based on the respiratory phase
    [sorting_indices,phase]                 = RespiratoryBinning(binning_pars);
    DataStruct_processed.RawKspaceData                = DataStruct_processed.RawKspaceData(:,sorting_indices,:,:);
    DataStruct_processed.Coordinates                  = DataStruct_processed.Coordinates(:,:,sorting_indices,:,:);

    % force some other parameters for the resp resolved recon
    param_struct.ReadoutsPerDynamic         = floor( size(DataStruct_processed.RawKspaceData,2) / param_struct.NumberOfDynamics );
    param_struct.NumberOfTemporalSplines    = round(4000/800); 
end


%% automatic parameters [don't touch]

NumberOfSpatialDims             = size(DataStruct_processed.Coordinates,1);                                                
spatial_ordering                = [2 1 3];
param_struct.IndicesOnReadout   = 1:size(DataStruct_processed.RawKspaceData,1);

RefImDims = size(DataStruct_processed.ReferenceImage);
svrs = structvars(param_struct);for i=1:size(svrs,1);eval(svrs(i,:));end
export_suffix = generate_export_suffix(NumberOfDynamics,ReadoutsPerDynamic,BeginReadoutIdx,RespResolvedReconstruction,lambda_det,lambda_TV,eps_TV,NumberOfComponents,NumberOfSpatialSplines,NumberOfTemporalSplines,RefImDims);

%% Reshape the snapshot data and kspace coordinates according to specified parameters



ReadoutIndices            = 1+[0:param_struct.ReadoutsPerDynamic*param_struct.NumberOfDynamics-1];
DataStruct_processed.Coordinates    = reshape(DataStruct_processed.Coordinates(spatial_ordering(1:NumberOfSpatialDims),param_struct.IndicesOnReadout,ReadoutIndices),NumberOfSpatialDims,numel(param_struct.IndicesOnReadout)*param_struct.ReadoutsPerDynamic,param_struct.NumberOfDynamics);
DataStruct_processed.RawKspaceData  = double(reshape(DataStruct_processed.RawKspaceData(param_struct.IndicesOnReadout,ReadoutIndices),numel(param_struct.IndicesOnReadout)*param_struct.ReadoutsPerDynamic,param_struct.NumberOfDynamics));


DataStruct_processed.Coordinates = demax(DataStruct_processed.Coordinates)/2;

%% Init MR MOTUS operator
MRMOTUS_recon           = MRMOTUS_Operator((DataStruct_processed.ReferenceImage(:)),permute(DataStruct_processed.Coordinates,[2 1 3]),param_struct);

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
info
b=toc


% Export results
disp('+Saving some reconstruction results...');
save([param_struct.export_folder,'dvf',export_suffix,'.mat'],'dvf','-v7.3')
save([param_struct.export_folder,'recon_info',export_suffix,'.mat'],'info');

disp('+Constructing low-rank motion-field components from coefficients')
[Phi,~,Psi,~] = MRMOTUS_recon.ExpandMotionfieldCoefficients(dvf);



end
