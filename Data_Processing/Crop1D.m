function cropped_indices=Crop1D(N,scale,centerout_flag)
    % Return indices to remove a factor 'scale' from a readout that originally consists of N points.
    %
    % With centerout_flag=0, function assume k0 point is at round(N/2)+1, for a readout of N samples.
    % With centerou_flag=1, function assume k0 point is at the first sample.
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
    
    if nargin<3
        centerout_flag = 0;
    end
    
    if N==1 && scale==1
        cropped_indices = 1;
    else
        if ~centerout_flag
            cropped_indices = [-round(N/scale/2):round(N/scale/2)-1] + round(N/2)+1;
        else
            cropped_indices = 1:round(N/scale);
        end
    end


end