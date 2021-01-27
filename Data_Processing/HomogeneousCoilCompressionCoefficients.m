function combinationCoefficients = HomogeneousCoilCompressionCoefficients(coil_sensitivities,noise_covariance,lambda,varargin)
    % Function that returns 'combinationCoefficients' for a linear coil combination that 
    % produces nearly homogeneous sensitivities.
    %
    % Input:
    %   coil_sensitivities      - coil sensitivities [N x N (xN) x N_c]
    %   noise_covariance        - noise covariance matrix [N_c x N_c]
    %   lambda                  - regularization parameter that weighs the influence of the noise covariance matrix
    %   varargin{1}             - sensitivity target [N x N (xN)]; region where to optimze for homogeneity
    %
    % Output:
    %   combinationCoefficients - [1 x 1 x 1 x N_c]
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

    coil_sensitivities_reshaped = double(reshape(coil_sensitivities,[],size(coil_sensitivities,4)));



    if nargin>3
        sensitivity_target = varargin{1};
    else
        sensitivity_target = ones(size(coil_sensitivities_reshaped,1),1)*100;
    end


    R = chol(noise_covariance);
    covar = R'*R;

    combinationCoefficients = pinv(coil_sensitivities_reshaped'*coil_sensitivities_reshaped+lambda*covar,1e-12)*(coil_sensitivities_reshaped'*sensitivity_target) ;
    combinationCoefficients = permute(combinationCoefficients,[2 3 4 1]);


    coil_sens_combined = LinearCoilCompression(coil_sensitivities,4,combinationCoefficients);
    slicer5d(coil_sens_combined)



end