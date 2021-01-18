function output=fft_mri(data)
    % fft of 'data' along all dimensions
    %
    % Niek Huttinga, UMC Utrecht, 2020

    output=fftshift(fftn(fftshift(data)));

end
