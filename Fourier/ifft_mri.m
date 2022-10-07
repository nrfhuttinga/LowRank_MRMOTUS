function output=ifft_mri(data)
    % ifft of 'data' along all dimensions
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

    
    output=ifftshift(ifftn(ifftshift(data)));

end

