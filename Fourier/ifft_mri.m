function output=ifft_mri(data)
    % ifft of 'data' along all dimensions
    %
    % Niek Huttinga, UMC Utrecht, 2020

    output=ifftshift(ifftn(ifftshift(data)));

end

