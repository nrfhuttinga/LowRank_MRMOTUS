function output=fft_mri(data)
    % fft of 'data' along all dimensions
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

    output=fftshift(fftn(fftshift(data)));

end
