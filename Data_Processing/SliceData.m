function slice=SliceData(data,d,si)
    % Extract slice with index 'si' from 'data', along dimension 'd'
    %
    % UMC Utrecht - Niek Huttinga - 2020
    %
    % Source: https://nl.mathworks.com/matlabcentral/answers/49904-accessing-certain-dimension-of-multidimensional-array
    
    S=size(data);
    Ar=reshape(data,prod(S(1:d-1)),[], prod(S(d+1:end)));
    slice = reshape( Ar(:,si,:) , [S(1:d-1), 1, S(d+1:end)]);
    
    
end
