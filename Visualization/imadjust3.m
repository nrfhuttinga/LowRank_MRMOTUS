function [ img ] = imadjust3( img, inLevels, outLevels, gamma, useSingle)
%IMADJUST3 Adjust image intensity values for N-D images (supports gpuArray)
%
%% HTML help: <a href="matlab:web('html/imadjust3.html')">imadjust3</a>.
%
%% Syntax
%   J = IMADJUST3(I)
%   J = IMADJUST3(I, PERCENT)
%   J = IMADJUST3(I, [LOW_IN; HIGH_IN])
%   J = IMADJUST3(I, INLEVEL, [LOW_OUT; HIGH_OUT])
%   J = IMADJUST3(I, INLEVEL, [LOW_OUT; HIGH_OUT])
%   J = IMADJUST3(I, INLEVEL, [LOW_OUT; HIGH_OUT], GAMMA)
%   J = IMADJUST3(I, INLEVEL, [LOW_OUT; HIGH_OUT], GAMMA, USESINGLE)
%   GPUARRAYB = IMADJUST3(GPUARRAYA, ___)
%
%% Description
%   J = IMADJUST3(I) maps the intensity values in a N-D grayscale image I
%   to new values in J such that 1% of data is saturated (Note that
%   imadjust defaults to 2%). This increases the constrast of the output
%   image J.
%   
%   J = IMADJUST3(I, PERCENT) maps the intensity values in I to new values
%   in J such that the given PERCENT percentage of the image is saturated.
%   This increases the constrast of the output image J.
%   
%   J = IMADJUST3(I, INLEVEL, [LOW_OUT; HIGH_OUT]) maps the
%   values in intensity image I to new values in J such that values between
%   boundaries map to values between LOW_OUT and HIGH_OUT. INLEVEL can be a
%   percentage as described above or a vector with [LOW_IN; HIGH_IN]
%   supplied directly. Values below LOW_IN and above HIGH_IN are clipped
%   that is, values below LOW_IN map to LOW_OUT, and those above HIGH_IN
%   map to HIGH_OUT. You can use an empty matrix ([]) for [LOW_IN; HIGH_IN]
%   or for [LOW_OUT; HIGH_OUT] to specify the default of [0 1]. If you omit
%   the argument, [LOW_OUT; HIGH_OUT] defaults to [0 1].
%
%   J = IMADJUST3(I, INLEVEL ,[LOW_OUT; HIGH_OUT], GAMMA) maps the
%   values of I to new values in J as described in the previous syntax.
%   GAMMA specifies the shape of the curve describing the relationship
%   between the values in I and J. If GAMMA is less than 1, the mapping is
%   weighted toward higher (brighter) output values. If GAMMA is greater
%   than 1, the mapping is weighted toward lower (darker) output values. If
%   you omit the argument, GAMMA defaults to 1 (linear mapping).
%
%   J = IMADJUST3(IMADJUST3(I, INLEVEL, [LOW_OUT; HIGH_OUT], GAMMA,
%   USESINGLE) forces single precision in case the input is an integer
%   datatype. This limits memory usage especially when working with
%   GPUARRAY.
%   
%   Note that if HIGH_OUT < LOW_OUT, the output image is reversed, as in a
%   photographic negative.
%   
%   
%% Class Support
%   The input image can be uint8, uint16, int16, double, single or a
%   gpuArray with one of these datatypes underlying. The output image has
%   the same class as the input image.
%   
%% Examples
%   Adjust Contrast of a N-D Grayscale Image
%   Read a low-contrast 4-D grayscale image into the workspace and display
%   a montage of it.
%       vol = load('mri');
%       figure;
%       subplot(1,2,1);
%       montage(vol.D);
%       title('Original image volume');
%   Adjust the contrast of  the image so that 1% of all voxels are
%   saturared and display a montage of it.
%       volAdj = imadjust3(vol.D);
%       subplot(1,2,2);
%       montage(volAdj);
%       title('1% of voxels saturated')
%
%   Adjust Contrast of a N-D Grayscale Image by Saturating a given
%   Percentage of Image Elements.
%   Read a low-contrast 4-D grayscale image into the workspace and display
%   a montage of it.
%       vol = load('mri');
%       figure
%       subplot(1,2,1);
%       montage(vol.D);
%       title('Original image volume');
%   Adjust the contrast of  the image so that 0.1% of all voxels are
%   saturared and display a montage of it.
%       volAdj = imadjust3(vol.D, 0.001);
%       subplot(1,2,2);
%       montage(volAdj);
%       title('0.1% of voxels saturated')
% 
%   Adjust Contrast of a N-D Grayscale Image Specifying Contrast Limits
%   Read a low-contrast 4-D grayscale image into the workspace and display
%   a montage of it.
%       vol = load('mri');
%       figure
%       subplot(1,2,1);
%       montage(vol.D);
%       title('Original image volume');
%   Adjust the contrast of  the image, specifying contrast limits
%       volAdj = imadjust3(vol.D, [0.3 0.7]);
%       subplot(1,2,2);
%       montage(volAdj);
%       title('Specified contrast limits [0.3 0.7]')
%   
%   Adjust Contrast of a N-D Grayscale Image Specifying non-linear Gamma
%   Read a low-contrast 4-D grayscale image into the workspace and display
%   a montage of it.
%       vol = load('mri');
%       figure
%       subplot(1,2,1);
%       montage(vol.D);
%       title('Original image volume');
%   Adjust the contrast of  the image, specifying a gamma value
%       volAdj = imadjust3(vol.D, [], [], 0.5);
%       subplot(1,2,2);
%       montage(volAdj);
%       title('Specified gamma of 0.5')
%
%   Adjust Contrast of a N-D Grayscale Image
%   Read a low-contrast 4-D grayscale image into a gpuArray and display
%   a montage of it.
%       vol = load('mri');
%       D = gpuArray(vol.D);
%       figure
%       subplot(1,2,1);
%       montage(D);
%       title('Original image volume');
%   Adjust the contrast of  the image so that 1% of all voxels are
%   saturared and display a montage of it.
%       DAdj = imadjust3(D);
%       subplot(1,2,2);
%       montage(volAdj);
%       title('1% of voxels saturated, computed on the GPU')
%
%% Input Arguments
%   I -- Image to be adjusted (gpuArray supported)
%       grayscale N-D image
%       Data Types: single | double | int16 | uint8 | uint16 | uint32 |
%       gpuArray (with the previous underlying data types)
%   INLEVEL -- PERCENT or [LOW_IN, HIGH_IN] - Contrast limits
%       0.01 (Default) | scalar between 0 and 1 | [0 1] (Default for empty
%       entry) | two-element numeric vector with values between 0 and 1
%       Contrast limits either as a percentage of pixels to be saturated or
%       as direct lower and upper limits. Values below LOW_IN and above
%       HIGH_IN are clipped; that is, values below LOW_IN map to LOW_OUT,
%       and those above HIGH_IN map to HIGH_OUT. If you specify an empty
%       matrix ([]), imadjust3 uses the default limits [0 1]. If only an
%       image is supplied imadjust3 will saturate 1% of image elements.
%       Data Types: double
%   OUTLEVEL -- Contrast limits for output image
%       [0 1] (Default) | two-elment numeric vector with values between 0
%       and 1
%       Contrast limits for the output image, specified as a two-element
%       numeric vector with values between 0 and 1. Values below low_in and
%       above high_in are clipped; that is, values below low_in map to
%       low_out, and those above high_in map to high_out. If you specify an
%       empty matrix ([]), imadjust3 uses the default limits [0 1].
%       Data Types: double
%   GAMMA -- Shape of the curve describing relationship of input and output
%   values
%       1 (default) | double real scalar
%       Shape of curve describing relationship of input and output values,
%       specified as a numeric value. If gamma is less than 1, imadjust
%       weights the mapping toward higher (brighter) output values. If
%       gamma is greater than 1, imadjust weights the mapping toward lower
%       (darker) output values. If you omit the argument, gamma defaults to
%       1 (linear mapping).
%       Data Types: double
%   USESINGLE -- Scalar bool flag to indicate whether to use single
%   precision for integer based images
%       false (default) | bool scalar
%       In case of integer format images double precision is not always
%       necessary. A true scalar will force imadjust3 to use single
%       precision for these image types thus saving memory.
%       Data Types: logical | numeric
%       
%% Output Arguments
%   J -- Adjusted image (gpuArray supported)
%       grayscale N-D image
%       Adjusted image, returned as a grayscale image. J has the same class
%       as the input image.
%       Data Types: single | double | int16 | uint8 | uint16 | uint32 |
%       gpuArray (with the previous underlying data types)
%
%% See Also
% imadjust
img = demax(abs(single(img)));
% Check if Image Processing Toolbox is installed
if ~license('test', 'image_toolbox')
    error('imadjust3:imageToolboxMissing', ...
        'Imadjust3 requires Matlab Image Processing Toolbox')
end

% Input check and default values for supplied inputs
narginchk(1,5);
if ~exist('useSingle', 'var') || isempty(useSingle)
    useSingle = false;
end
if ~exist('gamma', 'var') || isempty(gamma)
    gamma = 1.;
end
if ~exist('outLevels', 'var') || isempty(outLevels)
    outLevels = [0, 1];
end
if ~exist('inLevels', 'var')
    inLevels = 0.01;
end
if isempty(inLevels)
    inLevels = [0,1];
end

% Check input validitiy
validImClasses = {'single', 'double', 'uint8', 'uint16', 'int16'};
validateattributes(img,{validImClasses{:}, 'gpuArray'}, {'nonempty', 'nonsparse'}, ...
    mfilename, 'Image', 1);
validateattributes(inLevels, {'double'}, {'vector', '>=', 0, '<=', 1}, ...
    mfilename, 'In levels', 2);
validateattributes(outLevels, {'double'}, {'vector', 'numel', 2, '>=', 0, ...
    '<=', 1}, mfilename, 'Out levels', 3);
validateattributes(gamma, {'double'}, {'scalar', 'real', 'nonnegative'}, ...
    mfilename, 'Gamma', 4);
validateattributes(useSingle, {'numeric', 'logical'}, {'scalar', 'binary'}, ...
    mfilename, 'Use single precision', 5);

if numel(inLevels) > 2
    error('imadjust3:levelsMustBe1Or2ElVec', ...
        'Levels must be a 1 or 2-element vector')
end

% Get image class in case of gpuArray
imClass = class(img);
if strcmpi(imClass, 'gpuArray')
    imClass = classUnderlying(img);
    if ~any(strcmp(imClass, validImClasses))
        error('imadjust3:I', ['gpuArray I must have a underlying class of: ', ...
            join(validImClasses, ', ')]);
    end
end

% Check image data type and in case of integer convert to double or single
if strcmpi(imClass, 'double') || strcmpi(imClass, 'single')
    if any(img(:) > 1) || any(img(:) < 0)
        error('imadjust3:I', 'If image type is double or single values must be between 0 and 1')
    end
    classChanged = false;
else
    classChanged = true;
    if useSingle
        img = im2single(img);
    else
        img = im2double(img);
    end
end

lowOut = outLevels(1);
highOut = outLevels(2);

% If instead of discrete levels a certain number of image elements should
% be saturated.
if isscalar(inLevels)
    [lowIn, highIn] = computeLevels(img, inLevels);
else
    lowIn = inLevels(1);
    highIn = inLevels(2);
    if lowIn >= highIn
        error('imadjust3:lowMustBeSmallerThanHigh', ...
            'LOW_IN must be less than HIGH_IN')
    end
end

if abs(lowIn - highIn) < 1e-10
    warning('imadjust3:lowMustBeDifferentThanHigh', ...
        'Resulting levels are too close to eachother, set to [0 1]')
    lowIn = 0;
    highIn = 1;
end

% Compute transform
img =  max(lowIn, min(highIn,img));
img = ((((img - lowIn) ./ (highIn - lowIn)) .^ gamma) .* (highOut - lowOut)) + lowOut;

% clip = img > highOut;
% if any(clip(:))
%     img(clip) = highOut;
% end
% clip = img < lowOut;
% if any(clip(:))
%     img(clip) = lowOut;
% end

% In case the class changed revert change
if classChanged
    img = images.internal.changeClass(imClass,img);
end
end

% Subfunction to compute low and high levels to saturate a given percentage
% of image elements
function [low, high] = computeLevels(img, alpha)
sortedVals = sort(img(isfinite(img)));
n = numel(sortedVals);

% compute values that enclose (1-alpha) percents of all values
low = sortedVals( floor((n-1) * alpha/2) + 1);
high = sortedVals( floor((n-1) * (1-alpha/2)) + 1);


end
