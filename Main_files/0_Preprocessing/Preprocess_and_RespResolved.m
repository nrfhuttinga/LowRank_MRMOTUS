function [DataStruct,HighresReferenceImage,binned_data]=Preprocess_and_RespResolved(DataStruct,ref_im_parameters)
%% Function to preprocess raw k-space data acquired and reconstruct a reference image
%
% DataStruct should contain the following fields:
%
%       .RawKspaceData                          -> [ #readoutsamples    #readouts        1          #coils  ]           % Raw k-space data
%       .Coordinates                            -> [ 2/3                #readoutsamples  #readouts ]                    % K-space trajectory coordinates
%       .Coils
%           .Sensitivities                      -> [ ImDims             ImDims           ImDims     #coils  ]           % Coil sensitivities estimates with espirit
%           .Noise_covariance                   -> [ #coils             #coils ]                                        % Noise covariance matrix between the channels [set to identity if not available]
%       [.SelfNavigator.SurrogateSignal]        -> [ #readouts ]                                                        % Respiratory surrogate signal, you can provide your own
%
% See "/LowRank_MRMOTUS/3DGMR/Preprocessing/ref_im_parameters_2DGA.m" and "/LowRank_MRMOTUS/3DGMR/Preprocessing/ref_im_parameters_3DGMR.m" for all required parameters and explanation.
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

% DataStruct = DataStruct;

if nargin<2
    ref_im_parameters=[];
end

ref_im_parameters = set_default(ref_im_parameters,'centerout_flag',0);
ref_im_parameters = set_default(ref_im_parameters,'readout_downsampling_highres',1);
% ref_im_parameters = set_default(ref_im_parameters,'parallel_reconstruction',0);
ref_im_parameters = set_default(ref_im_parameters,'coil_compression',1);


ref_im_parameters.binning_pars = set_default(ref_im_parameters.binning_pars,'split_in_exhale',1);

if isfield(ref_im_parameters,'surrogate_pars')
    surrogate_pars = ref_im_parameters.surrogate_pars;
end

ncoils = size(DataStruct.RawKspaceData,4);
coil_flag = ncoils>1 


%% Extract image dimensions of reference image from coordinates before we rescale them
indices_on_readouts         = Crop1D(size(DataStruct.RawKspaceData,1),ref_im_parameters.readout_downsampling,ref_im_parameters.centerout_flag);
indices_on_readouts_highres = Crop1D(size(DataStruct.RawKspaceData,1),ref_im_parameters.readout_downsampling_highres,ref_im_parameters.centerout_flag);

ref_im_parameters.ImDims    = make_even(round(max(reshape(sqrt(sum(DataStruct.Coordinates(:,indices_on_readouts,:).^2,1)),[],1),[],1)));
ImDims_highres              = make_even(round(max(reshape(sqrt(sum(DataStruct.Coordinates(:,indices_on_readouts_highres,:).^2,1)),[],1),[],1)));

NumberOfSpatialDims         = size(DataStruct.Coordinates,1);

ref_im_parameters = set_default(ref_im_parameters,'lambda_cc',0);

%% Estimate coil compression coefficients for linear homogeneous coil compression
% As explained in Supporting information 2: "Extension of MR-MOTUS to multi-coil acquisitions"
disp('+Computing coil combination coefficients');

if ~isfield(DataStruct,'Coils')
    DataStruct.Coils.Sensitivities=ones;
end

if coil_flag && ~isfield(DataStruct.Coils,'Noise_covariance') 
    DataStruct.Coils.Noise_covariance = eye(ncoils,ncoils);
end




if ref_im_parameters.coil_compression && coil_flag
    disp('+Computing coil combination coefficients');
    
    disp('  +Reconstructing low-res coil img');
    pars.readout_downsampling=2.25;
    pars.pics_flag=0;
    pars.readout_indices_ref = 1:4000;
    pars.bart.regularization_lambda = 1e-7;
    rtrspla01_config
    RefImage_mc=ReconstructRefImage(DataStruct,pars);
    DataStruct_uncompressed = DataStruct;
    DataStruct.Coils.CompressionCoefficients = HomogeneousCoilCompressionCoefficients(RefImage_mc,eye(size(DataStruct.Coils.Noise_covariance)),0);         
end


%% Extract surrogate from spokes [if required]

if ~isfield(DataStruct,'SelfNavigator')
    if ref_im_parameters.centerout_flag==0
        surrogate_pars.k0_index = size(DataStruct.RawKspaceData,1)/2+1;
    else
        surrogate_pars.k0_index = 1;
    end
    disp('+Extracting self-navigation signal');
%     prs.surrogate_parameters.force_pc = 2;
    
%     prs.lowpass_filtering=0;prs.visualize=1;prs.surrogate_parameters.force_pc = 2;

%     DataStruct    = ExtractSelfnavData(DataStruct,prs);
end



%% Binning

% nt=1000;tsnecoords=tsne(abs(reshape(permute(DataStruct.SelfNavigator.ImageData(Crop1D(232,5,0),1:nt,:,1),[1 4 2 3]),[],nt)).');
% figure;scatter(tsnecoords(:,1),tsnecoords(:,2))

% nclusters=ref_im_parameters.binning_pars.nclusters;
% nt = length(DataStruct.SelfNavigator.ReadoutIndices);
% simba_pars.nclusters = nclusters;
% simba_pars.dim_reduction_method = 'pca';
% clustered_points=simba_test(DataStruct.RawKspaceData(:,DataStruct.SelfNavigator.ReadoutIndices(1:nt),:,1),simba_pars);


% Coil compress if not in PI mode
if ref_im_parameters.coil_compression && coil_flag
    disp('+Performing linear coil compression');
    DataStruct.RawKspaceData = LinearCoilCompression(DataStruct.RawKspaceData,4,DataStruct.Coils.CompressionCoefficients);
%     DataStruct.RawKspaceData = DataStruct.RawKspaceData(:,:,:,1);
end

prs = ref_im_parameters;
prs.lowpass_filtering=0;prs.visualize=1;prs.surrogate_parameters.force_pc = 2;
DataStruct = ExtractSelfnavData(DataStruct,prs);

% DataStruct.SelfNavigator.SurrogateSignal = medfilt1(DataStruct.SelfNavigator.SurrogateSignal,10,[],1);

% == Perform respiratory binning
ref_im_parameters_inex.binning_pars = ref_im_parameters.binning_pars;

ref_im_parameters_inex.binning_pars.surrogate_signal = DataStruct.SelfNavigator.SurrogateSignal;
ref_im_parameters_inex.binning_pars.binning_strategy = 'hybrid';
ref_im_parameters_inex.binning_pars.resp_phases = 10;

ref_im_parameters_inex.binning_pars=set_default(ref_im_parameters_inex.binning_pars,'thresh',0.08);


% 1) separate inhale and exhale
[~,~,inhale_idx,exhale_idx] = RespiratoryBinning(ref_im_parameters_inex.binning_pars);

if ~ref_im_parameters_inex.binning_pars.split_in_exhale
    inhale_idx = [inhale_idx(:);exhale_idx(:)];
    exhale_idx = [];
end


nclusters=ref_im_parameters_inex.binning_pars.nclusters;
[cluster_idx_inhale,b_in,ptcd_in,D_inhale]=kmeans(ref_im_parameters_inex.binning_pars.surrogate_signal([inhale_idx(:)]),nclusters,'MaxIter',300);
% sort the cluster centroids
[~,s_ind_in]=sort(b_in,'descend');
cluster_ordering = [s_ind_in(:)];

if ref_im_parameters_inex.binning_pars.split_in_exhale
    [cluster_idx_exhale,b_ex,ptcd_ex,D_exhale]=kmeans(ref_im_parameters_inex.binning_pars.surrogate_signal([exhale_idx(:)]),nclusters,'MaxIter',300);
    [~,s_ind_ex]=sort(b_ex,'ascend');
    cluster_ordering = [s_ind_in(:);s_ind_ex(:)+nclusters]
    nclusters = nclusters*2;
end




% D_total = cat(2,D_inhale,D_exhale);
% [~,in_ex_cluster]=min(b_in); % exhale cluster in exhale part
% [~,ex_ex_cluster]=min(b_ex); % exhale cluster in inhale part




% order the clusters to ensure smooth transitions from inhale to exhale along cluster dimensions



cmap = colormap(jet(nclusters));
cmap = cmap(randperm(length(cmap)),:);

figure;
% k=1;

%%

pars_csm.overgrid = 2;
pars_csm.nreadouts = size(DataStruct.RawKspaceData,2);
pars_csm.indices_on_readouts = Crop1D(size(DataStruct.RawKspaceData,1),ref_im_parameters.readout_downsampling);

if ~ref_im_parameters.coil_compression && ref_im_parameters.pics_flag
% if ~isfield(DataStruct.Coils,'Sensitivities') || isempty(DataStruct.Coils.Sensitivities)
     DataStruct.Coils.Sensitivities = EstimateCSM(DataStruct.RawKspaceData , DataStruct.Coordinates, pars_csm);
% end
end

k=1;
for i=cluster_ordering.'
    
    
    disp(['Reconstructing cluster ',num2str(k),'/',num2str(nclusters)]);
    
%     respiratory_bin_idx = clustered_points(:,i).';
    
    
    if k>nclusters/2 && ref_im_parameters_inex.binning_pars.split_in_exhale % exhale
        % select N_per_bin spoke sets closest to the cluster centroid
        [~,top_X_points] = sort(D_exhale(:,i-nclusters/2),'ascend');
        N_per_bin = ceil(length(ref_im_parameters_inex.binning_pars.surrogate_signal)/nclusters);
        top_X_points=top_X_points(1:N_per_bin);
        
        respiratory_bin_idx = exhale_idx(top_X_points);

    else %inhale
        % select N_per_bin spoke sets closest to the cluster centroid
        [~,top_X_points] = sort(D_inhale(:,i),'ascend');
        N_per_bin = ceil(length(ref_im_parameters_inex.binning_pars.surrogate_signal)/nclusters);
        top_X_points=top_X_points(1:N_per_bin);
        respiratory_bin_idx = inhale_idx(top_X_points);
    end
    
%     respiratory_bin_idx = inhale_exhale_sorted(find(cluster_idx(:)==i));
    
    plot(ref_im_parameters_inex.binning_pars.surrogate_signal([respiratory_bin_idx]),'Color',cmap(i,:));
    drawnow;
    hold on;




    if DataStruct.Sequence.Self_navigation_interval>1
            readouts_selfnav=1:length(DataStruct.SelfNavigator.SurrogateSignal);
            readouts_ref_binned = readouts_selfnav(respiratory_bin_idx);
            ref_im_parameters.readout_indices_ref = reshape((SelfnavToNormal(DataStruct.SelfNavigator.ReadoutIndices(readouts_ref_binned),DataStruct.Sequence.Self_navigation_interval)),[],1);
            ref_im_parameters.readout_indices_ref(ref_im_parameters.readout_indices_ref>size(DataStruct.RawKspaceData,2))=[];
    else
            ref_im_parameters.readout_indices_ref               = matrix_to_vec(respiratory_bin_idx);
    end
    disp(['Extracted ',num2str(numel(ref_im_parameters.readout_indices_ref)),' spokes with binning']);



    %% Reconstruct two reference images: low-resolution for motion estimations, high-resolution for visualizations
    % DataStruct.Coordinates = DataStruct.Coordinates/2;
    
    
    % Reconstruct low-res ref image
    ReferenceImage(:,:,:,:,k) = ReconstructRefImage( DataStruct, ref_im_parameters);

    DataStruct.SelectedReadoutsReferenceImage{k} = ref_im_parameters.readout_indices_ref;

    k=k+1;
    

end

% figure;imagesc(abs(ReferenceImage));colormap gray; axis image; axis off;

% Reconstruct high-res ref image
%%
highres_ref_im_parameters                       = ref_im_parameters;
highres_ref_im_parameters.readout_downsampling  = ref_im_parameters.readout_downsampling_highres;
highres_ref_im_parameters.ImDims                = ImDims_highres;

if ref_im_parameters.readout_downsampling_highres~=ref_im_parameters.readout_downsampling
    HighresReferenceImage                           = ReconstructRefImage( DataStruct, highres_ref_im_parameters);
else
    HighresReferenceImage = ReferenceImage;
end



% Register both images;
% Otherwise upscaled motion-fields reconstructed with a low-res ref image will not have
% the same geometry as the high res reference image that will be used for visualziations
if ~isfield(ref_im_parameters,'registration') || ref_im_parameters.registration
    
    disp('+Registering high-res reference image to low-res reference image');
    if NumberOfSpatialDims==2
        ReferenceImage_upscaled = imresize(abs(ReferenceImage),size(HighresReferenceImage));
    else
        ReferenceImage_upscaled = imresize3(abs(ReferenceImage),size(HighresReferenceImage));
    end
    [opt,metric]=imregconfig('multimodal');opt.MaximumIterations = 700;
    tform=imregtform(demax(abs(HighresReferenceImage)),demax(ReferenceImage_upscaled),'similarity',opt,metric,'DisplayOptimization',1);
    disp('+Applying registration results');
    if NumberOfSpatialDims==2
        HighresReferenceImage_warped =imwarp(squeeze(real(HighresReferenceImage)),tform,'OutputView',imref2d(size(HighresReferenceImage))) + 1i*imwarp(squeeze(imag(HighresReferenceImage)),tform,'OutputView',imref2d(size(HighresReferenceImage)));
    else
        HighresReferenceImage_warped =imwarp(squeeze(real(HighresReferenceImage)),tform,'OutputView',imref3d(size(HighresReferenceImage))) + 1i*imwarp(squeeze(imag(HighresReferenceImage)),tform,'OutputView',imref3d(size(HighresReferenceImage)));
    end

   

    HighresReferenceImage                                   = HighresReferenceImage_warped;

end
    
 if NumberOfSpatialDims==2
    figure;imagesc(abs(HighresReferenceImage));axis image; axis off; colormap gray;title('High-res reference image');
    figure;imagesc(abs(ReferenceImage));axis image; axis off; colormap gray;title('Low-res reference image');
else
    slicer5d((HighresReferenceImage));
    slicer5d((ReferenceImage));
end

%% Summarize everything in the _processed struct


if ~isfield(DataStruct.SelfNavigator,'LowpassFilterDelay')
    DataStruct.SelfNavigator.LowpassFilterDelay = 0;
end
DataStruct.RawKspaceData                      = DataStruct.RawKspaceData(indices_on_readouts,1: end-DataStruct.SelfNavigator.LowpassFilterDelay,:,:);
DataStruct.Coordinates                        = DataStruct.Coordinates(:,indices_on_readouts,1: end-DataStruct.SelfNavigator.LowpassFilterDelay,:,:);

DataStruct.ReferenceImage                     = ReferenceImage;


end
