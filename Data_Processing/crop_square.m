function cropped_image = crop_square(im,dim_crop)
% Crop 'dim_crop' from the boundaries of 'im'
%
% Inputs
%   im              - input data N x N (x N)
%   dim_crop        - slices to remove from each dimension [d1 d2 d3]
%
% Outputs
%   cropped_image   - outputs [N-2*d1 x N-2*d2, N-2*d3]
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

im_size = size(im);

if (numel(im_size)~=numel(dim_crop) && dim_crop(end)==0)
    im_size(end+1)=1;
end

if sum(im_size/2<=dim_crop)>0 || (numel(dim_crop)>3 && dim_crop(4)~=0)
   msg = 'Error in input values';
   error(msg)
end


if dim_crop(1)~=0
    im(1:dim_crop(1),:,:)=[];
    im(end-dim_crop(1)+1:end,:,:)=[];
end
if numel(dim_crop) > 1 && dim_crop(2)~=0
    im(:,1:dim_crop(2),:)=[];
    im(:,end-dim_crop(2)+1:end,:)=[];
end
if numel(dim_crop) > 2 && dim_crop(3)~=0
    im(:,:,1:dim_crop(3))=[];
    im(:,:,end-dim_crop(3)+1:end)=[];
end

cropped_image = im;
end
