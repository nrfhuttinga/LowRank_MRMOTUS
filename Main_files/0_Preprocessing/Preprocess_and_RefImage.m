function [DataStruct,HighresReferenceImage]=Preprocess_and_RefImage(DataStruct,ref_im_parameters)
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

if ~isfield(ref_im_parameters,'centerout_flag')
    ref_im_parameters.centerout_flag = 0;
end

if ~isfield(ref_im_parameters,'readout_downsampling_highres')
    ref_im_parameters.readout_downsampling_highres = 1;
end

if isfield(ref_im_parameters,'surrogate_pars')
    surrogate_pars = ref_im_parameters.surrogate_pars;
end

coil_flag = size(DataStruct.RawKspaceData,4)>1;


if coil_flag && ~isfield(DataStruct.Coils,'Noise_covariance') 
    DataStruct.Coils.Noise_covariance = eye(size(DataStruct.Coils.Sensitivities,4),size(DataStruct.Coils.Sensitivities,4));
end


%% Extract image dimensions of reference image from coordinates before we rescale them
indices_on_readouts         = Crop1D(size(DataStruct.RawKspaceData,1),ref_im_parameters.readout_downsampling,ref_im_parameters.centerout_flag);
indices_on_readouts_highres = Crop1D(size(DataStruct.RawKspaceData,1),ref_im_parameters.readout_downsampling_highres,ref_im_parameters.centerout_flag);

ref_im_parameters.ImDims    = make_even(round(max(reshape(sqrt(sum(DataStruct.Coordinates(:,indices_on_readouts,:).^2,1)),[],1),[],1)));
ImDims_highres              = make_even(round(max(reshape(sqrt(sum(DataStruct.Coordinates(:,indices_on_readouts_highres,:).^2,1)),[],1),[],1)));

NumberOfSpatialDims         = size(DataStruct.Coordinates,1);



%% Estimate coil compression coefficients for linear homogeneous coil compression
% As explained in Supporting information 2: "Extension of MR-MOTUS to multi-coil acquisitions"
disp('+Computing coil combination coefficients');

if coil_flag
    % Set region to optimize the homogeneity on
    sens_target = (sum(abs(DataStruct.Coils.Sensitivities),4)>0)*1;
    sens_target = sens_target(:);

    % Perform compression coefficients
    DataStruct.Coils.CompressionCoefficients = HomogeneousCoilCompressionCoefficients(DataStruct.Coils.Sensitivities,diag(diag(DataStruct.Coils.Noise_covariance)),ref_im_parameters.lambda_cc,sens_target);
end

%% Extract surrogate from spokes [if required]

if ~isfield(DataStruct,'SelfNavigator')
    if ref_im_parameters.centerout_flag==0
        surrogate_pars.k0_index = size(DataStruct.RawKspaceData,1)/2+1;
    else
        surrogate_pars.k0_index = 1;
    end
    disp('+Extracting self-navigation signal');
    
    DataStruct    = ExtractSurrogateSignal(DataStruct,surrogate_pars);
end



%% Binning


% Coil compress if not in PI mode
if ~ref_im_parameters.parallel_reconstruction && coil_flag
    disp('+Performing linear coil compression');
    DataStruct.RawKspaceData = LinearCoilCompression(DataStruct.RawKspaceData,4,DataStruct.Coils.CompressionCoefficients);
end

% Perform respiratory binning
ref_im_parameters.binning_pars.surrogate_signal = DataStruct.SelfNavigator.SurrogateSignal(ref_im_parameters.total_readouts);
[respiratory_bin_idx,~]       = RespiratoryBinning(ref_im_parameters.binning_pars);

if DataStruct.Sequence.Self_navigation_interval>1
    readouts_selfnav=1:length(DataStruct.SelfNavigator.SurrogateSignal);
    readouts_ref_binned = readouts_selfnav(respiratory_bin_idx);
    ref_im_parameters.readout_indices_ref = reshape((SelfnavToNormal(DataStruct.SelfNavigator.ReadoutIndices(readouts_ref_binned),DataStruct.Sequence.Self_navigation_interval)),[],1);
    ref_im_parameters.readout_indices_ref(ref_im_parameters.readout_indices_ref>size(DataStruct.RawKspaceData,2))=[];
else
    ref_im_parameters.readout_indices_ref               = ref_im_parameters.total_readouts(matrix_to_vec(respiratory_bin_idx));
end
disp(['Extracted ',num2str(numel(ref_im_parameters.readout_indices_ref)),' spokes with binning']);



%% Reconstruct two reference images: low-resolution for motion estimations, high-resolution for visualizations
% DataStruct.Coordinates = DataStruct.Coordinates/2;

% Reconstruct low-res ref image
ReferenceImage = ReconstructRefImage( DataStruct, ref_im_parameters);
% figure;imagesc(abs(ReferenceImage));colormap gray; axis image; axis off;

% Reconstruct high-res ref image

highres_ref_im_parameters                       = ref_im_parameters;
highres_ref_im_parameters.readout_downsampling  = ref_im_parameters.readout_downsampling_highres;
highres_ref_im_parameters.ImDims                = ImDims_highres;

if ref_im_parameters.readout_downsampling_highres~=ref_im_parameters.readout_downsampling
    HighresReferenceImage                           = ReconstructRefImage( DataStruct, highres_ref_im_parameters);
else
    HighresReferenceImage = ReferenceImage;
end



%% Register both images;
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
DataStruct.SelectedReadoutsReferenceImage     = ref_im_parameters.readout_indices_ref;
DataStruct.ReferenceImage                     = ReferenceImage;


end
