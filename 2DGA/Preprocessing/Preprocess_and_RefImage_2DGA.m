% Script to preprocess raw k-space data acquired and reconstruct a reference image
%
% Niek Huttinga, UMCU Utrecht, 2020.


clear all;

lowrank_mrmotuspath = '../LowRank_MRMOTUS/';
cd(lowrank_mrmotuspath)
addpath(genpath(pwd));



%% 1) Load data

data_folder  = [lowrank_mrmotuspath,'Data/2DGA/'];
load([data_folder,'DataStruct_raw.mat']);

% DataStruct_raw should contain the following fields:
%
%       .RawKspaceData                          -> [ #readoutsamples    #readouts        1          #coils  ]           % Raw k-space data
%       .Coordinates                            -> [ 2/3                #readoutsamples  #readouts ]                    % K-space trajectory coordinates
%       .Coils
%           .Sensitivities                      -> [ ImDims             ImDims           ImDims     #coils  ]           % Coil sensitivities estimates with espirit
%           .Noise_covariance                   -> [ #coils             #coils ]                                        % Noise covariance matrix between the channels [set to identity if not available]
%       [.SelfNavigator.SurrogateSignal]        -> [ #readouts ]                                                        % Respiratory surrogate signal, you can provide your own


DataStruct_processed = DataStruct;




%% 3) Cut the data according to 'ref_im_parameters.total_readouts'


% Scale coordinates and extract proper readouts
DataStruct_processed.Coordinates      = double(demax(DataStruct_processed.Coordinates(:,:,ref_im_parameters.total_readouts))/2);
DataStruct_processed.RawKspaceData    = DataStruct_processed.RawKspaceData(:,ref_im_parameters.total_readouts,:,:);


%% 3) Estimate coil compression coefficients for linear homogeneous coil compression
% As explained in Supporting information 2: "Extension of MR-MOTUS to multi-coil acquisitions"
disp('+Computing coil combination coefficients');

% Set region to optimize the homogeneity on
sens_target = (sum(abs(DataStruct_processed.Coils.Sensitivities),4)>0)*1;
sens_target = sens_target(:);

% Perform compression coefficients
DataStruct_processed.Coils.CompressionCoefficients = HomogeneousCoilCompressionCoefficients(DataStruct_processed.Coils.Sensitivities,diag(diag(DataStruct_processed.Coils.Noise_covariance)),ref_im_parameters.lambda_cc,sens_target);


%% 4) Extract surrogate from spokes [if required]

if ~isfield(DataStruct_processed,'SelfNavigator')
    surrogate_pars.k0_index = size(DataStruct_processed.RawKspaceData,1)/2+1;
    DataStruct_processed    = ExtractSurrogateSignal(DataStruct_processed,surrogate_pars);
end



%% 5) Binning


% Coil compress if not in PI mode
if ~ref_im_parameters.parallel_reconstruction
    disp('+Performing linear coil compression');
    DataStruct_processed.RawKspaceData = LinearCoilCompression(DataStruct_processed.RawKspaceData,4,DataStruct_processed.Coils.CompressionCoefficients);
end

% Perform respiratory binning
ref_im_parameters.binning_pars.surrogate_signal = DataStruct_processed.SelfNavigator.SurrogateSignal;
[respiratory_bin_idx,~]       = RespiratoryBinning(ref_im_parameters.binning_pars);
readout_indices_ref               = matrix_to_vec(respiratory_bin_idx);


disp(['Extracted ',num2str(numel(readout_indices_ref)),' spokes with binning']);



%% 6) Reconstruct two reference images: low-resolution for motion estimations, high-resolution for visualizations

% Reconstruct low-res ref image
ReferenceImage = ReconstructRefImage( DataStruct_processed, ref_im_parameters);

% Reconstruct high-res ref image
highres_ref_im_parameters                       = ref_im_parameters;
highres_ref_im_parameters.readout_downsampling  = 1;
highres_ref_im_parameters.ImDims                = ImDims_highres;
HighresReferenceImage                           = ReconstructRefImage( DataStruct_processed, highres_ref_im_parameters);


% Register both images;
% Otherwise upscaled motion-fields reconstructed with a low-res ref image will not have
% the same geometry as the high res reference image that will be used for visualziations
disp('+Registering high-res reference image to low-res reference image');
ReferenceImage_upscaled = imresize(abs(ReferenceImage),size(HighresReferenceImage));
[opt,metric]=imregconfig('multimodal');opt.MaximumIterations = 500;
tform=imregtform(demax(abs(HighresReferenceImage)),demax(ReferenceImage_upscaled),'similarity',opt,metric,'DisplayOptimization',1);
disp('+Applying registration results');
HighresReferenceImage_warped =imwarp(squeeze(real(HighresReferenceImage)),tform,'OutputView',imref2d(size(HighresReferenceImage))) + 1i*imwarp(squeeze(imag(HighresReferenceImage)),tform,'OutputView',imref2d(size(HighresReferenceImage)));
slicer5d(HighresReferenceImage_warped);
slicer5d(ReferenceImage);

%% 7) Summarize everything in the _processed struct


HighresReferenceImage                                   = HighresReferenceImage_warped;

DataStruct_processed.RawKspaceData                      = DataStruct_processed.RawKspaceData(indices_on_readouts,1: end-DataStruct_processed.SelfNavigator.LowpassFilterDelay,:,:);
DataStruct_processed.Coordinates                        = DataStruct_processed.Coordinates(:,indices_on_readouts,1: end-DataStruct_processed.SelfNavigator.LowpassFilterDelay,:,:);
DataStruct_processed.SelectedReadoutsReferenceImage     = readout_indices_ref;
DataStruct_processed.ReferenceImage                     = ReferenceImage;

DataStruct = DataStruct_processed;

% Save results in the new DataStruct_processed struct that will be used in the MR-MOTUS recons
disp('Saving...')
save([data_folder,'DataStruct_processed.mat'],'DataStruct');
save([data_folder,'HighresReferenceImage.mat'],'HighresReferenceImage');
disp('Done!')
