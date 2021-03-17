function result=MRMOTUS_PostProcessing( MF, param_struct )

if nargin<2 || (~isfield(param_struct,'DataStruct_path') || ~isfield(param_struct,'highres_referenceimage_path'))
    error('Specify atleast param_struct.DataStruct_path and param_struct.highres_referenceimage_path');
end

if ~isfield(param_struct,'postprocessing')
    param_struct.postprocessing=[];
end


param_struct.postprocessing = set_default(param_struct.postprocessing,'WarpRefImageFlag',           1);
param_struct.postprocessing = set_default(param_struct.postprocessing,'MotionImageOverlayFlag',     1);
param_struct.postprocessing = set_default(param_struct.postprocessing,'JacDeterminantsFlag',        1);
param_struct.postprocessing = set_default(param_struct.postprocessing,'HighresVisualizationFlag',   1);
param_struct.postprocessing = set_default(param_struct.postprocessing,'visualization_handle_noabs', @(x) x);
param_struct.postprocessing = set_default(param_struct.postprocessing,'HighresVisualizationFlag',   1);
param_struct.postprocessing = set_default(param_struct.postprocessing,'flips',                      zeros(1,size(MF,2)));

param_struct = set_default(param_struct,'export_folder',                '');
param_struct = set_default(param_struct,'export_suffix',                '');
param_struct = set_default(param_struct,'RespResolvedReconstruction',   0);
param_struct = set_default(param_struct,'ReadoutsPerDynamic',           31);


NumberOfSpatialDims = size( MF , 2);
no_dyns             = size( MF , 3);


% Load the high resolution reference image for visualizations
if param_struct.postprocessing.HighresVisualizationFlag
    HighresReferenceImage = load(param_struct.highres_referenceimage_path);
    fnames = fieldnames(HighresReferenceImage);
    HighresReferenceImage = HighresReferenceImage.(fnames{:});
else
    load(param_struct.DataStruct_path);
    HighresReferenceImage = DataStruct.ReferenceImage;
end


% flip reference image if required
for i=1:numel(param_struct.postprocessing.flips)
    if param_struct.postprocessing.flips(i)
        HighresReferenceImage = flip(HighresReferenceImage,i);
    end
end
HighresReferenceImage = HighresReferenceImage.^.7;
% flip motion-field if required
MF = FlipMotionField(MF, param_struct.postprocessing.flips);

ImDim     = round((size(MF,1)).^(1/NumberOfSpatialDims));
ImDim_vis = size(HighresReferenceImage,1);


param_struct.postprocessing = set_default(param_struct.postprocessing,'cor_slice',          round(ImDim_vis/2+1));
param_struct.postprocessing = set_default(param_struct.postprocessing,'sag_slice',          round(ImDim_vis/2+1));
param_struct.postprocessing = set_default(param_struct.postprocessing,'trans_slice',        round(ImDim_vis/2+1));
param_struct.postprocessing = set_default(param_struct.postprocessing,'crop_coronal',       @(x) x);
param_struct.postprocessing = set_default(param_struct.postprocessing,'crop_sagittal',      @(x) x);
param_struct.postprocessing = set_default(param_struct.postprocessing,'crop_transverse',    @(x) x);


%%  ==== Warping reference image with reconstructed motion-fields ======


if param_struct.postprocessing.WarpRefImageFlag || param_struct.postprocessing.MotionImageOverlayFlag
    disp('=== Warping reference image with reconstructed motion-fields ===')
    % Visualize the reference image that will be used for visualizations from
    % now on onwards
    slicer5d(abs(reshape_to_square(single(abs(HighresReferenceImage)),NumberOfSpatialDims)));

    % Actual warping of the high-resolution reference image
    result = WarpReferenceImage(HighresReferenceImage,MF);
    disp('+Saving warped reference image');
    try
        save([param_struct.export_folder,'result',param_struct.export_suffix,'.mat'],'result','-v7.3')
    catch
        warning('error in saving, attempting to save in working directory...');
        save(['result',param_struct.export_suffix,'.mat'],'result','-v7.3')
    end
    
    disp('+Done saving');

    % Visualize warped reference image results
    slicer5d(abs(result))
else
    result = [];
end



%%  Set visualization handles

if NumberOfSpatialDims == 3
    param_struct.postprocessing = set_default(param_struct.postprocessing,'handle_coronal',@(x) rot90((squeeze(abs(x(param_struct.postprocessing.cor_slice,:,:,:)))),1));
    param_struct.postprocessing = set_default(param_struct.postprocessing,'handle_sagittal',@(x) rot90((squeeze(abs(x(:,param_struct.postprocessing.sag_slice,:,:)))),1));
    param_struct.postprocessing = set_default(param_struct.postprocessing,'handle_transverse',@(x) rot90((squeeze(abs(x(:,:,param_struct.postprocessing.trans_slice,:)))),0));
else 
    visualization_handle_abs = @(x) abs(param_struct.postprocessing.visualization_handle_noabs(x));
end



%% Load reference image mask for visualization
image_for_vis = HighresReferenceImage;

ref_mask_path = [get_data_dir(param_struct.DataStruct_path)];
if param_struct.postprocessing.JacDeterminantsFlag || param_struct.postprocessing.MotionImageOverlayFlag
try
    load([ref_mask_path,'/RefMask.mat']);
catch
    if NumberOfSpatialDims == 3
        mask_coronal=Poly2Binary(param_struct.postprocessing.handle_coronal(image_for_vis));
        mask_sagittal=Poly2Binary(param_struct.postprocessing.handle_sagittal(image_for_vis));
        mask_transverse=Poly2Binary(param_struct.postprocessing.handle_transverse(image_for_vis));
        save([ref_mask_path,'/RefMask.mat'],'mask_sagittal','mask_transverse','mask_coronal');
    else
        RefMask = Poly2Binary(image_for_vis);
        save(ref_mask_path,'RefMask','-v7.3');
    end
end
end

%% Jacobian determinants

if param_struct.postprocessing.JacDeterminantsFlag
    % Computing determinants in batches at resolution of high-res reference image
    % 1) Compute determinants on low-res motion-fields
    % 2) Upscale resulting image to high-res ref image resolution
    clearvars det_rc
    for i=1:size(MF,3)
        if NumberOfSpatialDims==3
            det_rc(:,:,:,i)=imresize3(single(DeterminantMotionFields(MF(:,:,i))),[ImDim_vis,ImDim_vis,ImDim_vis]);
        else
            det_rc(:,:,:,i)=imresize(single(DeterminantMotionFields(MF(:,:,i))),[ImDim_vis,ImDim_vis]);
        end
    end


    % Select indices to visualize determinant maps for; dynamics with min and max motion in FH direction
    [~,minimum_motion_index] = min(squeeze(sum(MF(:,end,:),1)),[],1);
    [~,maximum_motion_index] = max(squeeze(sum(MF(:,end,:),1)),[],1);
    max_ip_det = squeeze(abs(det_rc(:,:,:,maximum_motion_index)));
    min_ip_det = squeeze(abs(det_rc(:,:,:,minimum_motion_index)));


    % Some visualization parameters [don't touch]
    determinant_scale = [0 2];
    alpha_ = .5;

    if NumberOfSpatialDims == 3
        fig1=figure('Renderer', 'painters');    
        set_background_black;
        set_figure_size(fig1,[0 0 1920 1100]);
        ha = tight_subplot(2,3,[.07 -.01],[.1 .1],[.22 .29]);

        % #1
        PlotOverlayedImage( param_struct.postprocessing.crop_coronal(param_struct.postprocessing.handle_coronal(image_for_vis).*mask_coronal),param_struct.postprocessing.crop_coronal(param_struct.postprocessing.handle_coronal(max_ip_det).*mask_coronal),alpha_,ha(1),determinant_scale,-0.02,1)
        PlotOverlayedImage( param_struct.postprocessing.crop_sagittal(param_struct.postprocessing.handle_sagittal(image_for_vis).*mask_sagittal),param_struct.postprocessing.crop_sagittal(param_struct.postprocessing.handle_sagittal(max_ip_det).*mask_sagittal),alpha_,ha(2),determinant_scale)
        PlotOverlayedImage( param_struct.postprocessing.crop_transverse(param_struct.postprocessing.handle_transverse(image_for_vis).*mask_transverse),param_struct.postprocessing.crop_transverse(param_struct.postprocessing.handle_transverse(max_ip_det).*mask_transverse),alpha_,ha(3),determinant_scale)

        % #2
        PlotOverlayedImage( param_struct.postprocessing.crop_coronal(param_struct.postprocessing.handle_coronal(image_for_vis).*mask_coronal),param_struct.postprocessing.crop_coronal(param_struct.postprocessing.handle_coronal(min_ip_det).*mask_coronal),alpha_,ha(4),determinant_scale,-0.02,1)
        PlotOverlayedImage( param_struct.postprocessing.crop_sagittal(param_struct.postprocessing.handle_sagittal(image_for_vis).*mask_sagittal),param_struct.postprocessing.crop_sagittal(param_struct.postprocessing.handle_sagittal(min_ip_det).*mask_sagittal),alpha_,ha(5),determinant_scale)
        PlotOverlayedImage( param_struct.postprocessing.crop_transverse(param_struct.postprocessing.handle_transverse(image_for_vis).*mask_transverse),param_struct.postprocessing.crop_transverse(param_struct.postprocessing.handle_transverse(min_ip_det).*mask_transverse),alpha_,ha(6),determinant_scale)
    else
        fig1=figure('Renderer', 'painters');
        set_background_black;
        set_figure_size(fig1,[0 0 1920 1100]);
        ha = tight_subplot(1,2,[.07 -0.40],[.1 .1],[.01 .01]);
        
        % #1
        PlotOverlayedImage(visualization_handle_abs(image_for_vis).*visualization_handle_abs(RefMask),visualization_handle_abs(max_ip_det).*visualization_handle_abs(RefMask),alpha_,ha(1),determinant_scale,0.16,1,.02);
        
        % #2
        PlotOverlayedImage(visualization_handle_abs(image_for_vis).*visualization_handle_abs(RefMask),visualization_handle_abs(min_ip_det).*visualization_handle_abs(RefMask),alpha_,ha(2),determinant_scale,0.16,0,.02);

    end
    
    % Save visualizations
    save_as = [param_struct.export_folder,'ImageDetOverlayed',param_struct.export_suffix];
    export_fig(save_as,'-png')
end

%% Motion image overlay

if param_struct.RespResolvedReconstruction
    delay_time = 4/no_dyns;
else
    delay_time = param_struct.ReadoutsPerDynamic*4.4e-3;
end

for i=1:no_dyns;text_num{i}=[num2str(delay_time*1000*(i-1)),' ms'];end
text_num = strjust(pad(text_num),'right');
for i=1:no_dyns;text{i} = ['  Time: ',text_num{i}];end


if param_struct.postprocessing.MotionImageOverlayFlag
    disp('=== Overlaying motion-fields on warped reference image... ===')

    % Upscale the motion-fields and apply the same visualization handles as for
    % the determinant visualization
    MF_highres = UpscaleMotionFields(MF,ImDim,ImDim_vis);
    clearvars MF_new;
    max_for_gif = min(no_dyns,800);

    if param_struct.postprocessing.HighresVisualizationFlag ==1
        display_factor = 3;     % downsampling factor of motion-field for visualization
    else
        display_factor = 1;     % downsampling factor of motion-field for visualization
    end
    
    threshold       = -10;      % threshold for visualization
    color           = 'g';      % color of motion-field
    padding         = 10;       % boundary to remove 
    scaling         = 1.4;      % scaling for visualizion

    if NumberOfSpatialDims == 3
        % coronal
        dimension       = 1;
        slice           = param_struct.postprocessing.cor_slice;
        rotations       = 1;
        [images_coronal,cm_coronal]=MotionImageOverlay_3Dt(result(:,:,:,:,1:max_for_gif).^(4/5),MF_highres,dimension,slice,threshold,display_factor,color,padding,scaling,rotations,rot90(mask_coronal,rotations+3));

        file_name = [param_struct.export_folder,'images_coronal',param_struct.export_suffix,'.gif'];
        CellsToGif(images_coronal,cm_coronal,delay_time,file_name)

        
        % sagittal
        dimension       = 2;
        slice           = param_struct.postprocessing.sag_slice;
        rotations       = 1;
        [images_sagittal,cm_sagittal]=MotionImageOverlay_3Dt(result(:,:,:,:,1:max_for_gif).^(4/5),MF_highres,dimension,slice,threshold,display_factor,color,padding,scaling,rotations,rot90(mask_sagittal,3+rotations));
        
        file_name = [param_struct.export_folder,'images_sagittal',param_struct.export_suffix,'.gif'];
        CellsToGif(images_sagittal,cm_sagittal,delay_time,file_name)

        
        % axial
        dimension       = 3;
        slice           = param_struct.postprocessing.trans_slice;
        rotations       = 0;
        [images_axial,cm_axial]=MotionImageOverlay_3Dt(result(:,:,:,:,1:max_for_gif).^(4/5),-MF_highres,dimension,slice,threshold,display_factor,color,padding,scaling,rotations,rot90(mask_transverse,rotations));

        file_name = [param_struct.export_folder,'images_axial',param_struct.export_suffix,'.gif'];
        CellsToGif(images_axial,cm_axial,delay_time,file_name)
        
        

        

        
        

    else
        
        MF_highres = param_struct.postprocessing.visualization_handle_noabs(reshape(MF_highres,ImDim_vis,ImDim_vis,NumberOfSpatialDims,no_dyns));
        clearvars mf_new;
        max_for_gif = min(no_dyns,800);
        for i=1:max_for_gif
            for j=1:size(MF_highres,3)
                mf_new{i}(:,:,j)=MF_highres(:,:,j,i);
            end
            
            
        end

        % Overlay motion-fields as vector-field on warped reference images
        [a,b]=MotionImageOverlay_2Dt(visualization_handle_abs((abs(result(:,:,:,:,1:max_for_gif))).^(4/5)),mf_new,-10,3,'g',10,1,0,visualization_handle_abs(RefMask));

        
        % Export resulting images as GIF
        CellsToGif(a,b,delay_time,[param_struct.export_folder,'/MotionImageOverlay',param_struct.export_suffix,'.gif'],text)

        close all;
    end
 
    

    


end
