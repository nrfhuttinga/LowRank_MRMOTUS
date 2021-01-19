% Script to preprocess raw k-space data acquired and reconstruct a reference image
%
% Copyright Niek Huttinga, UMC Utrecht, 2020, all rights reserved.


clear all;

lowrank_mrmotuspath = '/nfs/arch11/researchData/USER/nhutting/code/LowRank_MRMOTUS/';
cd(lowrank_mrmotuspath)
addpath(genpath(pwd));



%% 1) Load data

data_folder  = [lowrank_mrmotuspath,'Data/3DGMR/'];
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


%% 2) Set parameters

% general parameters
readout_scaling                 = 2.25;                                                 % Factor that determines downsampling the readout direction in order to decrease the resolution of the reference image                
parallel_imaging_reconstruction = 0;                                                    % Can optionally be set to 1 to use a PI-reconstructed high-res reference image for visualizations
recon_overgridding              = 2;                                                    % Overgridding for BART reconstructions, relative to original
vis_overgridding                = 1.1;                                                  % Final overgridding of the reference image, relative to orignal
NumberOfSpatialDims             = size(DataStruct_processed.Coordinates,1);             % Dimensionality of the image [2D/3D]
lambda_cc                       = 0;                                                    % Regularization parameter for homogeneous coil compression 
total_readouts                  = 700:size(DataStruct_processed.RawKspaceData,2);       % Readouts to extract, we typically don't use the first 700 readouts to let the signal settle in a steady state first.

% parameters specify for bart
ref_im_parameters.bart.iterations               = 500;                                  % Number of iterations
ref_im_parameters.bart.regularization_lambda    = .0005;                                % L1-Wavelet regularization parameter
ref_im_parameters.bart.version                  = 6;                                    % BART version [required because there is a slight change in commands over the versions]


% parameters for binning
binning_pars.binning_strategy        = 'phase';                                         % Phase or amplitude binning
binning_pars.thresh                  = 0.005;                                           % Threshold for peak detection algorithm
binning_pars.resp_phases             = 5;                                               % Number of respiratory phases
binning_pars.return_extreme_phase    = 2;                                               % 0 - return all phases, 1 - return inhale, 2 - return exhale




%% 3) Cut the data according to 'total_readouts'

% Extract image dimensions of reference image from coordinates before we rescale them
indices_on_readouts = Crop1D(size(DataStruct_processed.RawKspaceData,1),readout_scaling);
N                   = make_even(round(max(reshape(sqrt(sum(DataStruct_processed.Coordinates(:,indices_on_readouts,:).^2,1)),[],1),[],1)*vis_overgridding));
N_highres           = make_even(round(max(reshape(sqrt(sum(DataStruct_processed.Coordinates(:,:,:).^2,1)),[],1),[],1)*vis_overgridding));

% Scale coordinates and extract proper readouts
DataStruct_processed.Coordinates                        = double(demax(DataStruct_processed.Coordinates(:,:,total_readouts))/2);
DataStruct_processed.RawKspaceData                      = DataStruct_processed.RawKspaceData(:,total_readouts,:,:);


%% 3) Estimate coil compression coefficients for linear homogeneous coil compression

disp('+Computing coil combination coefficients');

% Set region to optimize the homogeneity on
sens_target = (sum(abs(DataStruct_processed.Coils.Sensitivities),4)>0)*1; 
sens_target=sens_target(:);

% Perform compression coefficients
DataStruct_processed.Coils.CompressionCoefficients = HomogeneousCoilCompressionCoefficients(DataStruct_processed.Coils.Sensitivities,diag(diag(DataStruct_processed.Coils.Noise_covariance)),lambda_cc,sens_target);         


%% 4) Extract surrogate from spokes [if required]

if ~isfield(DataStruct_processed,'SelfNavigator')
    surrogate_pars.k0_index = size(DataStruct_processed.RawKspaceData,1)/2+1;
    DataStruct_processed=ExtractSurrogateSignal(DataStruct_processed,surrogate_pars);
end



%% 5) Binning


% Coil compress if not in PI mode
if ~parallel_imaging_reconstruction
    disp('+Performing linear coil compression');
    DataStruct_processed.RawKspaceData = LinearCoilCompression(DataStruct_processed.RawKspaceData,4,DataStruct_processed.Coils.CompressionCoefficients);
end

% Perform respiratory binning
binning_pars.surrogate_signal        = DataStruct_processed.SelfNavigator.SurrogateSignal;
[respiratory_bin_idx,~]=RespiratoryBinning(binning_pars);
ref_im_readouts = matrix_to_vec(respiratory_bin_idx);


disp(['Extracted ',num2str(numel(ref_im_readouts)),' spokes with binning']);



%% 6) Reconstruct two reference images: low-resolution for motion estimations, high-resolution for visualizations

% Reconstruct low-res ref image
ref_im_parameters.parallel_reconstruction   = parallel_imaging_reconstruction;
ref_im_parameters.readout_indices_ref       = ref_im_readouts;
ref_im_parameters.readout_downsampling      = readout_scaling;
ref_im_parameters.ImDims                    = N;
ref_im_parameters.recon_overgridding        = recon_overgridding;
ReferenceImage = ReconstructRefImage( DataStruct_processed, ref_im_parameters);

% Reconstruct high-res ref image
highres_ref_im_parameters = ref_im_parameters;
highres_ref_im_parameters.readout_downsampling = 1;
highres_ref_im_parameters.ImDims = N_highres;
HighresReferenceImage = ReconstructRefImage( DataStruct_processed, highres_ref_im_parameters);


% Register both images;
% Otherwise upscaled motion-fields reconstructed with a low-res ref image will not have the same geometry as the high res reference image that will be used for visualziations
disp('+Registering high-res reference image to low-res reference image');
ReferenceImage_upscaled = imresize3(abs(ReferenceImage),size(HighresReferenceImage));
[opt,metric]=imregconfig('multimodal');opt.MaximumIterations = 500;
tform=imregtform(demax(abs(HighresReferenceImage)),demax(ReferenceImage_upscaled),'similarity',opt,metric,'DisplayOptimization',1);
disp('+Applying registration results');
HighresReferenceImage_warped =imwarp(squeeze(real(HighresReferenceImage)),tform,'OutputView',imref3d(size(HighresReferenceImage))) + 1i*imwarp(squeeze(imag(HighresReferenceImage)),tform,'OutputView',imref3d(size(HighresReferenceImage)));
slicer5d(HighresReferenceImage_warped);
slicer5d(ReferenceImage);
%% 7) Summarize everything in the _processed struct


HighresReferenceImage                                   = HighresReferenceImage_warped;

DataStruct_processed.RawKspaceData                      = DataStruct_processed.RawKspaceData(indices_on_readouts,1: end-DataStruct_processed.SelfNavigator.LowpassFilterDelay,:,:);
DataStruct_processed.Coordinates                        = DataStruct_processed.Coordinates(:,indices_on_readouts,1: end-DataStruct_processed.SelfNavigator.LowpassFilterDelay,:,:);
DataStruct_processed.SelectedReadoutsReferenceImage     = ref_im_readouts;
DataStruct_processed.ReferenceImage                     = ReferenceImage;

DataStruct = DataStruct_processed;

% Save results in the new DataStruct_processed struct that will be used in the MR-MOTUS recons
disp('Saving...')
save([data_folder,'DataStruct_processed.mat'],'DataStruct');
save([data_folder,'HighresReferenceImage.mat'],'HighresReferenceImage');
disp('Done!')
