restoredefaultpath;addpath(genpath('../LowRank_MRMOTUS')); % TB edit: this is easier for people that download the repo
close all;
clear all;
rng(1)


%% load parameters and data

disp('=== Loading parameters and data ===');

% load parameters
Parameters_2Dt_RespMotion

% Load reference image, self-navigation signal, k-space data and k-space trajectory
load(DataStruct_path)

[DataStruct.ReferenceImage,DataStruct.RawKspaceData] = CalibrateReferenceAndKdata(DataStruct.ReferenceImage,DataStruct.RawKspaceData,round(size(DataStruct.RawKspaceData,1)/2+1));


DataStruct.RawKspaceData                    = DataStruct.RawKspaceData(:,param_struct.BeginReadoutIdx:end,:,:);
DataStruct.Coordinates                      = DataStruct.Coordinates(:,:,param_struct.BeginReadoutIdx:end);
DataStruct.SelfNavigator.SurrogateSignal    = DataStruct.SelfNavigator.SurrogateSignal(param_struct.BeginReadoutIdx:end);


%% check trajectory

disp('=== Plotting trajectory ===');

figure;PlotTrajectory(DataStruct.Coordinates(:,:,1:20))
 


%% sort the data in case of respiratory-resolved reconstruction

disp('=== Sorting data ===');

if param_struct.RespResolvedReconstruction
    % perform phase binning
    binning_pars.surrogate_signal           = DataStruct.SelfNavigator.SurrogateSignal;
    binning_pars.binning_strategy           = 'phase';
    binning_pars.thresh                     = 0.005;
    binning_pars.resp_phases                = param_struct.NumberOfDynamics;
    binning_pars.return_extreme_phases      = 0;

    % sort the data based on the respiratory phase
    [sorting_indices,phase]                 = RespiratoryBinning(binning_pars);
    DataStruct.RawKspaceData                = DataStruct.RawKspaceData(:,sorting_indices,:,:);
    DataStruct.Coordinates                  = DataStruct.Coordinates(:,:,sorting_indices,:,:);

    % force some other parameters for the resp resolved recon
    param_struct.ReadoutsPerDynamic         = floor( size(DataStruct.RawKspaceData,2) / param_struct.NumberOfDynamics );
    param_struct.NumberOfTemporalSplines    = round(4000/800); 
end


%% automatic parameters [don't touch]

NumberOfSpatialDims             = size(DataStruct.Coordinates,1);                                                
spatial_ordering                = [2 1 3];
param_struct.IndicesOnReadout   = 1:size(DataStruct.RawKspaceData,1);

svrs = structvars(param_struct);for i=1:size(svrs,1);eval(svrs(i,:));end
RefImDims = size(DataStruct.ReferenceImage);
export_suffix = generate_export_suffix(NumberOfDynamics,ReadoutsPerDynamic,BeginReadoutIdx,RespResolvedReconstruction,lambda_det,lambda_TV,eps_TV,NumberOfComponents,NumberOfSpatialSplines,NumberOfTemporalSplines,RefImDims);

%% Reshape the snapshot data and kspace coordinates according to specified parameters



ReadoutIndices            = 1+[0:param_struct.ReadoutsPerDynamic*param_struct.NumberOfDynamics-1];
DataStruct.Coordinates    = reshape(DataStruct.Coordinates(spatial_ordering(1:NumberOfSpatialDims),param_struct.IndicesOnReadout,ReadoutIndices),NumberOfSpatialDims,numel(param_struct.IndicesOnReadout)*param_struct.ReadoutsPerDynamic,param_struct.NumberOfDynamics);
DataStruct.RawKspaceData  = double(reshape(DataStruct.RawKspaceData(param_struct.IndicesOnReadout,ReadoutIndices),numel(param_struct.IndicesOnReadout)*param_struct.ReadoutsPerDynamic,param_struct.NumberOfDynamics));
ImDim                     = round((numel(DataStruct.ReferenceImage(:)))^(1/NumberOfSpatialDims));


DataStruct.Coordinates = demax(DataStruct.Coordinates)/2;

%% Init MR MOTUS operator
MRMOTUS_recon           = MRMOTUS_Operator((DataStruct.ReferenceImage(:)),permute(DataStruct.Coordinates,[2 1 3]),param_struct);

%% Actual reconstructions

% Handle to evaluate forward model and gradients at iterate 'x', this is
% required for lbfgs
f_handle = @(x) MRMOTUS_recon.forward_and_gradient_lowrank(x,DataStruct.RawKspaceData);

      

        
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


disp('+Constructing low-rank motion-field components from coefficients')
[Phi,Phi_rshp,Psi,PsiT] = MRMOTUS_recon.ExpandMotionfieldCoefficients(dvf);
disp('+Constructing motion-field from low-rank motion-field components')
% Constructing motion-fields from low-rank components:
% [Note: this may take a large amount of RAM in 3D!]
mf = reshape(Phi*Psi.',[],MRMOTUS_recon.NumberOfSpatialDims,MRMOTUS_recon.NumberOfDynamics);


% Export results
disp('+Saving some reconstruction results...');
save([export_folder,'dvf',export_suffix,'.mat'],'dvf','-v7.3')
save([export_folder,'recon_info',export_suffix,'.mat'],'info');

disp('=== Done! ===')     


%%  ==== Warping reference image with reconstructed motion-fields ======

disp('=== Warping reference image with reconstructed motion-fields ===')


% Load the high resolution reference image for visualizations
if param_struct.HighresVisualizationFlag
    load(highres_referenceimage_path);
else
    HighresReferenceImage = DataStruct.ReferenceImage;
end
N_vis = size(HighresReferenceImage,1);


% Visualize the reference image that will be used for visualizations from
% now on onwards
slicer5d(abs(reshape(single(abs(HighresReferenceImage)),N_vis,N_vis)));


% Actual warping of the high-resolution reference image
result = WarpReferenceImage(HighresReferenceImage,mf);
disp('+Saving warped reference image');
save([export_folder,'result',export_suffix,'.mat'],'dvf','-v7.3')

% Visualize warped reference image results
slicer5d(permute(abs(result),[1 2 3 5 4]));



%% plot determinants


disp('=== Performing validations on Jacobian determinants of motion-fields ===')

% Computing determinants in batches at resolution of high-res reference image
% 1) Compute determinants on low-res motion-fields
% 2) Upscale resulting image to high-res ref image resolution
clearvars det_rc
parfor i=1:size(Psi,1)
    det_rc(:,:,:,i)=imresize(single(DeterminantMotionFields(squeeze(mf(:,:,i)),0)),[N_vis,N_vis]);
end

% Select indices to visualize determinant maps for
[~,minimum_motion_index] = min(sum(abs(Psi),2),[],1);
[~,maximum_motion_index] = max(sum(abs(Psi),2),[],1);
max_ip_det = squeeze(abs(det_rc(:,:,:,maximum_motion_index)));
min_ip_det = squeeze(abs(det_rc(:,:,:,minimum_motion_index)));


% Some visualization parameters [don't touch]
determinant_scale = [0 2];
alpha_ = 0.5;
export = 1;
image = HighresReferenceImage;

% Visualization handles to specify e.g. crops and rotations if necessary
visualization_handle_noabs = @(x) squeeze(x(:,25:end-25,:,:,:));
visualization_handle_abs = @(x) abs(visualization_handle_noabs(x));

ref_mask_path = [get_data_dir(DataStruct_path),'RefMask.mat'];
if exist(ref_mask_path)>0
    load(ref_mask_path);
else
    RefMask = Poly2Binary(image);
    save(ref_mask_path,'RefMask','-v7.3');
end

% Actual visualizations
fig1=figure('Renderer', 'painters');
set_background_black;
set_figure_size(fig1,[0 0 1920 1100]);
ha = tight_subplot(1,2,[.07 -0.40],[.1 .1],[.01 .01]);
PlotOverlayedImage(visualization_handle_abs(image).*visualization_handle_abs(RefMask),visualization_handle_abs(max_ip_det).*visualization_handle_abs(RefMask),alpha_,ha(1),determinant_scale,0.16,1,.02);
PlotOverlayedImage(visualization_handle_abs(image).*visualization_handle_abs(RefMask),visualization_handle_abs(min_ip_det).*visualization_handle_abs(RefMask),alpha_,ha(2),determinant_scale,0.16,0,.02);

% Save visualizations
save_as = [export_folder,'ImageDetOverlayed',export_suffix];
export_fig(save_as,'-png')


%% Motion image overlay 2Dt
disp('=== Overlaying motion-fields on warped reference image... ===')

% Upscale the motion-fields and apply the same visualization handles as for
% the determinant visualization
mf_highres = UpscaleMotionFields(mf,ImDim,N_vis);
mf_highres = visualization_handle_noabs(reshape(mf_highres,N_vis,N_vis,NumberOfSpatialDims,param_struct.NumberOfDynamics));
clearvars mf_new;
max_for_gif = min(param_struct.NumberOfDynamics,800);
for i=1:max_for_gif
    for j=1:size(mf_highres,3)
        mf_new{i}(:,:,j)=mf_highres(:,:,j,i);
    end
end

% Overlay motion-fields as vector-field on warped reference images
[a,b]=MotionImageOverlay_2Dt(visualization_handle_abs((abs(result(:,:,:,1:max_for_gif))).^(4/5)),mf_new,0,3,'g',10,1,0,rot90(visualization_handle_abs(RefMask),3));

% Optionally add some text at the bottom per GIF to specify dynamics
for i=1:param_struct.NumberOfDynamics;text_num{i}=[num2str(round(DataStruct.Sequence.Repetition_time*param_struct.ReadoutsPerDynamic*1000*(i-1))),' ms'];end
text_num = strjust(pad(text_num),'right');
for i=1:param_struct.NumberOfDynamics;text{i} = ['  Time: ',text_num{i}];end

% Export resulting images as GIF
CellsToGif(a,b,param_struct.ReadoutsPerDynamic*DataStruct.Sequence.Repetition_time,[export_folder,'/MotionImageOverlay',export_suffix,'.gif'],text)

disp('=== Done! ===')
