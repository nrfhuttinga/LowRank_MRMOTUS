function fft_data=ifft_mri_1d(input_data,dimension,npoints)
    % 1d fft of 'input_data' along 'dimension'
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
    
    if nargin==2
        fft_data=fftshift(ifft(fftshift(input_data,dimension),[],dimension),dimension);
    else
        fft_data=fftshift(ifft(fftshift(input_data,dimension),npoints,dimension),dimension);

end
