function [dvf,plotImage] = ExtractSlices_3D(Image,motionField,dimension,varargin)
% Inputs: 
%   - Image:        image used to overlay motionField
%   - motionField:  1x3 cell; motion in true x,y,z order (2 1 3 matlab
%                   convention), top is postive, bottom is negative
%   - dimension:    dimension along which to present a slice.
%
% Note:
%   - dim1 = y, dim2 = x, dim3 = z;
%
% Example:
%   cam_man = double(imread('cameraman.tif'));
%   for i=1:3; motion{i} = ones(size(cam_man)); end;
%   dimension = 3;
%   MotionImageOverlay_3D(cam_man,motion,dimension)
%
% Copyright Niek Huttinga, UMC Utrecht, 2020

if dimension > 3
    error('Dimension cannot be larger than 3!')
end



inputSize = size(Image);
if numel(inputSize)==2
    if dimension ~= 3
        error('Cannot visualize third dimension on 2D input image');
    else
        inputSize(3) = 1;
    end
end

if numel(varargin)>0
    display_slice = varargin{1};
else
    display_slice = ceil(inputSize(dimension)/2);
end


if dimension == 3
    I_12 = squeeze(reshape(Image(:, :, display_slice),inputSize(1),inputSize(2)));
    motionField_12_HOR = reshape(motionField{1}(:,:,display_slice),inputSize(1),inputSize(2)); % horizontal in dim1xdim2 = dim2 = x
    motionField_12_VER = reshape(motionField{2}(:,:,display_slice),inputSize(1),inputSize(2)); % vertical in dim1xdim2 = dim1 = y
    plotImage = I_12;
    plotMotion_HOR = motionField_12_HOR;
    plotMotion_VER = motionField_12_VER;
elseif dimension == 2
    I_13 = squeeze(reshape(Image(:, display_slice,:),inputSize(1),inputSize(3))); 
    motionField_13_HOR = reshape(motionField{3}(:, display_slice,:),inputSize(1),inputSize(3)); % horizontal in dim1xdim3 = dim3 = z
    motionField_13_VER = reshape(motionField{2}(:, display_slice,:),inputSize(1),inputSize(3)); % vertical in dim1xdim3 = dim1 = y
    plotImage = I_13;
    plotMotion_HOR = motionField_13_HOR;
    plotMotion_VER = motionField_13_VER;
elseif dimension == 1
    I_23 = squeeze(reshape(Image(display_slice,:,:),inputSize(2),inputSize(3)));
    motionField_23_HOR = reshape(motionField{3}(display_slice,:,:),inputSize(2),inputSize(3)); % horizontal in dim2xdim3 = dim3 = z
    motionField_23_VER = -reshape(motionField{1}(display_slice,:,:),inputSize(2),inputSize(3)); % vertical in dim2xdim3 = dim2 = x
    plotImage = I_23;
    plotMotion_HOR = motionField_23_HOR;
    plotMotion_VER = motionField_23_VER;
end

dvf = cat(3,plotMotion_HOR,plotMotion_VER);



end

