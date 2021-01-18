function cropped_indices=Crop1D(N,scale)
    % Return indices to remove a factor 'scale' from a readout that originally consists of N points.
    %
    % Note: function assume k0 point is at round(N/2)+1, for a readout of N samples.
    %
    % Niek Huttinga, UMC Utrecht, 2020.
    
    cropped_indices = [-round(N/scale/2):round(N/scale/2)-1] + round(N/2)+1;



end