function slice=SliceData(data,d,si, squeeze_flag)
    % Extract slice with index 'si' from 'data', along dimension 'd'
    %
    % UMC Utrecht - Niek Huttinga - 2020
    %
    % Source: https://nl.mathworks.com/matlabcentral/answers/49904-accessing-certain-dimension-of-multidimensional-array
    
    if nargin< 4
        squeeze_flag = 1;
    end
    
    S=size(data);
    Ar=reshape(data,prod(S(1:d-1)),[], prod(S(d+1:end)));
    slice = reshape( Ar(:,si,:) , [S(1:d-1), 1, S(d+1:end)]);
    
    if squeeze_flag
        slice = squeeze(slice);
    end
    
    
end
