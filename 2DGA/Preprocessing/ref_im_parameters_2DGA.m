%% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

% general parameters
ref_im_parameters.readout_downsampling          = 2.25;                    	% Factor that determines downsampling the readout direction in order to decrease the resolution of the reference image                
ref_im_parameters.parallel_reconstruction       = 0;                       	% Can optionally be set to 1 to use a PI-reconstructed high-res reference image for visualizations
ref_im_parameters.recon_overgridding            = 2;                       	% Overgridding for BART reconstructions, relative to original
ref_im_parameters.vis_overgridding              = 1.1;                   	% Final overgridding of the reference image, relative to orignal
ref_im_parameters.lambda_cc                     = 1e0;                    	% Regularization parameter for homogeneous coil compression 
ref_im_parameters.total_readouts                = 700:7000;                % Readouts to extract, we typically don't use the first 700 readouts to let the signal settle in a steady state first.

% parameters specify for bart
ref_im_parameters.bart.iterations               = 500;                     	% Number of iterations
ref_im_parameters.bart.regularization_lambda    = .0005;                   	% L1-Wavelet regularization parameter
ref_im_parameters.bart.version                  = 6;                      	% BART version [required because there is a slight change in commands over the versions]


% parameters for binning
ref_im_parameters.binning_pars.binning_strategy        = 'phase';          	% Phase or amplitude binning
ref_im_parameters.binning_pars.thresh                  = 0.005;            	% Threshold for peak detection algorithm
ref_im_parameters.binning_pars.resp_phases             = 5;                	% Number of respiratory phases
ref_im_parameters.binning_pars.return_extreme_phase    = 2;                	
% 0 - return all phases, 1 - return inhale, 2 - return exhale

