restoredefaultpath;addpath(genpath('/nfs/arch11/researchData/USER/nhutting/code/LowRank_MRMOTUS'));
close all;
clear all;
rng(1)


%% load parameters and data

disp('=== Loading parameters and data ===');

% load parameters
Parameters_3Dt_RespMotion

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

disp('+Done saving');


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
slicer5d(abs(reshape(single(abs(HighresReferenceImage)),N_vis,N_vis,N_vis)));


% Actual warping of the high-resolution reference image
result = WarpReferenceImage(HighresReferenceImage,mf);
disp('+Saving warped reference image');
save([export_folder,'result',export_suffix,'.mat'],'dvf','-v7.3')
disp('+Done saving');

% Visualize warped reference image results
slicer5d(permute(abs(result),[1 2 3 5 4]));



%%  set visualization handles

crop_coronal        = @(x) x(:,10:end);
crop_sagittal       = @(x) x(:,30:end-10);
crop_transverse     = @(x) x(:,10:end);

cor_slice = 67;
sag_slice = 64;
trans_slice = 64;

handle_coronal      = @(x) rot90((squeeze(abs(x(cor_slice,:,:,:)))),1);
handle_sagittal     = @(x) rot90((squeeze(abs(x(:,sag_slice,:,:)))),1);
handle_transverse   = @(x) rot90((squeeze(abs(x(:,:,trans_slice,:)))),0);


%% plot determinants

% Computing determinants in batches at resolution of high-res reference image
% 1) Compute determinants on low-res motion-fields
% 2) Upscale resulting image to high-res ref image resolution
clearvars det_rc
for i=1:size(Psi,1)
    det_rc(:,:,:,i)=imresize3(single(DeterminantMotionFields(mf(:,:,i))),[N_vis,N_vis,N_vis]);
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

ref_mask_path = [get_data_dir(DataStruct_path)];
try
    load([ref_mask_path,'/RefMask.mat']);
catch
    mask_coronal=Poly2Binary(handle_coronal(image));
    mask_sagittal=Poly2Binary(handle_sagittal(image));
    mask_transverse=Poly2Binary(handle_transverse(image));
    save([ref_mask_path,'/RefMask.mat'],'mask_sagittal','mask_transverse','mask_coronal');
end



fig1=figure('Renderer', 'painters');
set_background_black;
set_figure_size(fig1,[0 0 1920 1100]);


ha = tight_subplot(2,3,[.07 -.01],[.1 .1],[.22 .29]);


% coronal
PlotOverlayedImage( crop_coronal(handle_coronal(image).*mask_coronal),crop_coronal(handle_coronal(max_ip_det).*mask_coronal),alpha_,ha(1),determinant_scale,-0.02,1)
PlotOverlayedImage( crop_sagittal(handle_sagittal(image).*mask_sagittal),crop_sagittal(handle_sagittal(max_ip_det).*mask_sagittal),alpha_,ha(2),determinant_scale)
PlotOverlayedImage( crop_transverse(handle_transverse(image).*mask_transverse),crop_transverse(handle_transverse(max_ip_det).*mask_transverse),alpha_,ha(3),determinant_scale)

PlotOverlayedImage( crop_coronal(handle_coronal(image).*mask_coronal),crop_coronal(handle_coronal(min_ip_det).*mask_coronal),alpha_,ha(4),determinant_scale,-0.02,1)
PlotOverlayedImage( crop_sagittal(handle_sagittal(image).*mask_sagittal),crop_sagittal(handle_sagittal(min_ip_det).*mask_sagittal),alpha_,ha(5),determinant_scale)
PlotOverlayedImage( crop_transverse(handle_transverse(image).*mask_transverse),crop_transverse(handle_transverse(min_ip_det).*mask_transverse),alpha_,ha(6),determinant_scale)


% Save visualizations
save_as = [export_folder,'ImageDetOverlayed',export_suffix];
export_fig(save_as,'-png')


%% Motion image overlay 3Dt
disp('=== Overlaying motion-fields on warped reference image... ===')

% Upscale the motion-fields and apply the same visualization handles as for
% the determinant visualization
mf_highres = UpscaleMotionFields(mf,ImDim,N_vis);
clearvars mf_new;
max_for_gif = min(param_struct.NumberOfDynamics,800);




if HighresVisualizationFlag ==1
    display_factor = 3;     % downsampling factor of motion-field for visualization
else
    display_factor = 1;     % downsampling factor of motion-field for visualization
end
threshold       = -10;      % threshold for visualization
color           = 'g';      % color of motion-field
padding         = 10;       % boundary to remove 
scaling         = 1.3;      % scaling for visualizion

% coronal
dimension       = 1;
slice           = cor_slice;
rotations       = 1;
[images_coronal,cm_coronal]=MotionImageOverlay_3Dt(result(:,:,:,:,1:max_for_gif),mf_highres,dimension,slice,threshold,display_factor,color,padding,scaling,rotations,mask_coronal);

% sagittal
dimension       = 2;
slice           = sag_slice;
rotations       = 1;
[images_sagittal,cm_sagittal]=MotionImageOverlay_3Dt(result(:,:,:,:,1:max_for_gif),mf_highres,dimension,slice,threshold,display_factor,color,padding,scaling,rotations,mask_sagittal);

% axial
dimension       = 3;
slice           = trans_slice;
rotations       = 0;
[images_axial,cm_axial]=MotionImageOverlay_3Dt(result(:,:,:,:,1:max_for_gif),mf_highres,dimension,slice,threshold,display_factor,color,padding,scaling,rotations,mask_transverse);

close all;

if databinning
    delay_time = 4/no_dynamics;
else
    delay_time = cones_per_dynamic*U_processed.Sequence.Repetition_time;
end


for i=1:no_dynamics;text_num{i}=[num2str(round(4/20*1000*(i-1))),' ms'];end
text_num = strjust(pad(text_num),'right');
for i=1:no_dynamics;text{i} = ['  Time: ',text_num{i}];end


