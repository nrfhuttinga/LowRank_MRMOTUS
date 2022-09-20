function image=ForceSquareShape(image)
% Checks if 'image' is provided as a square image, rather than a vector. If
% not, then attemps reshaping.
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

    if numel(image)==max(size(image))
        warning('Provided image in vector format, attempting to reshape to square image.');
        try
            image = reshape_to_square(image,3);
            warning('Reshaped image to 3D square');

        catch
            image = reshape_to_square(image,2);
            warning('Reshaped image to 2D square');

        end
    else
        image=image;
    end

end

    