function size_vector=size_ext(data,dimensions)
    % Extension of the built-in size function.
    % Input: array + dimensions
    % Output: vector with sizes along specified dimensions
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

    if numel(dimensions)>numel(size(data))
        error('Dimension vector larger than number of dimensions, aborting.')
    elseif nargin < 2
        error('Not enough inputs, aborting');
    else
        size_vector = zeros(1,numel(dimensions));
        for i=1:numel(dimensions)
            size_vector(i)=size(data,dimensions(i));
        end
    end

end