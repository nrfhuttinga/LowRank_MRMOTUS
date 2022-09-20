function [Phi,Psi,MRMOTUS_recon, param_struct]=LowRank_MRMOTUS_Alternating_Recon(DataStruct_processed,param_struct)
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
    
max_traj = max(abs(DataStruct_processed.Coordinates(:)));

DataStruct_processed.Coordinates = DataStruct_processed.Coordinates/max_traj/2;%*(size(DataStruct_processed.ReferenceImage,1)/max_traj);


%% set LBFGS options
clearvars options


options.factr   = param_struct.lbfgs_termination_threshold;
options.maxIts  = param_struct.NumberOfReconIterations;
options.m       = 2;
options.plotting= param_struct.VisualizationFlag;
options.errFcn  = {@(x) x(:),@(x) toc};

ref_image = DataStruct_processed.ReferenceImage(:);

%%
for recon_iter=1:2
    slicer5d(reshape_to_square(ref_image))

    
    %% Init MR MOTUS operator

    MRMOTUS_recon           = MRMOTUS_Operator((ref_image(:)),permute(DataStruct_processed.Coordinates,[2 1 3]),param_struct);
    MRMOTUS_recon.RegularizationFlag = 0;

    %% Actual reconstructions

    % Handle to evaluate forward model and gradients at iterate 'x', this is
    % required for lbfgs
    f_handle = @(x) MRMOTUS_recon.forward_and_gradient_lowrank(x,DataStruct_processed.RawKspaceData);




    disp('=== Running LBFGS-B Reconstructions ===')


   
%     if recon_iter==1
%         param_struct.NumberOfSpatialSplines=72;
%     else
%         param_struct.NumberOfSpatialSplines=16;
%     end
    
    if recon_iter >= 1
        % Run recons
        tic
        options.x0      = MRMOTUS_recon.SolutionVariables_init(:);
        options.lb      = options.x0*Inf;
        options.ub      = options.x0*Inf;
        
        [dvf,~,info]=lbfgsb(f_handle,options.lb,options.ub,options);
    else
        options.x0      = MRMOTUS_recon.SolutionVariables_init(:)*0;
        dvf = options.x0;

        options.lb      = options.x0*Inf;
        options.ub      = options.x0*Inf;
    end
    
     [~,~,~,DVF_ops]=f_handle(dvf);
    
    disp('+ Updating reference image...')
    ref_image=0;
%     max_op = @(x) max(abs(x(:)));
%     max_op = @(x) quantile(abs(x(:)),.9);
    max_op = @(x) norm(x(6, :),'fro');
%     op_scaling=0;

    for i=1:numel(DVF_ops)
        disp(['+     Dynamic ',num2str(i),'/',num2str(numel(DVF_ops))])
        ref_image_i = DVF_ops{i}'*(MRMOTUS_recon.DCF(:,i).*DataStruct_processed.RawKspaceData(:,i));
        ref_image = ref_image + ref_image_i;
    end
    for i=1:numel(DVF_ops)
        disp(['+     Dynamic ',num2str(i),'/',num2str(numel(DVF_ops))])
%         op_scaling = op_scaling + max_op(DVF_ops{i}*ref_image(:)).^2;
        data(:,i)=DVF_ops{i}*ref_image(:);
    end
%     op_scaling = 
    
    
%     ref_image = ref_image * max_op(DataStruct_processed.RawKspaceData)/sqrt(op_scaling);
    ref_image = ref_image * max_op(DataStruct_processed.RawKspaceData)/max_op(data);
    clearvars data;
%     ref_image = ref_im_test;

end

    
% Export results
disp('+Saving some reconstruction results...');
save([param_struct.export_folder,'dvf',param_struct.export_suffix,'.mat'],'dvf','-v7.3')
save([param_struct.export_folder,'recon_info',param_struct.export_suffix,'.mat'],'info');

disp('+Constructing low-rank motion-field components from coefficients')
[Phi,~,Psi,~] = MRMOTUS_recon.ExpandMotionfieldCoefficients(dvf);



end
