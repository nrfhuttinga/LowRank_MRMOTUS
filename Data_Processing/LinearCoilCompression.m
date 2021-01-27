function compressed_inputs=LinearCoilCompression(inputs_to_compress,compression_dimension,compression_coefficients)
    % Performs linear coil compression with 'compression_coefficients' along
    % 'compression_dimension'
    %
    % Input:
    %   inputs_to_compress          - input that should be compressed
    %   compression_dimension       - dimension to compress 'inputs_to_compress'
    %   compression_coefficients    - coefficients for the linear coil compression (see HomgeneousCoilCompressionCoefficients.m)
    %
    % Output:
    %   compressed_inputs           - 'inputs_to_compress' compressed along 'compression_dimension' 
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
    



    if size(compression_coefficients,compression_dimension)~=size(inputs_to_compress,compression_dimension)
        warning('Input data not consistent in compression dimension. Data already compressed?');
        compressed_inputs = inputs_to_compress;
    else

        compressed_inputs = sum(bsxfun(@times,inputs_to_compress,compression_coefficients),compression_dimension);
    end

end




