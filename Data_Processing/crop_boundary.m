function data_out=crop_boundary(data_in,new_size)
% Perform centered crop to resize 'data_in' to 'new_size'; output size is 'new_size'
%
% Input:
%   data_in         - input data N x N x (N)
%   new_size        - new size = [N1_new x N2_new x (N3_new)]
%
% Output:
%   data_out        - cropped data [N1_new x N2_new x (N3_new)]
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

old_size = size(data_in);

for i=1:size(data_in,4)
    

if new_size(1) < old_size(1)
    border_to_crop(1) = (old_size(1)-new_size(1))/2;
else
    border_to_crop(1) = 0;
end

if new_size(2) < old_size(2)
    border_to_crop(2) = (old_size(2)-new_size(2))/2;
else
    border_to_crop(2) = 0;
end


if numel(new_size)<3
    border_to_crop(3) = 0;
else
    if new_size(3) < size(data_in,3) && size(data_in,3)~=1
        border_to_crop(3) = (old_size(3)-new_size(3))/2;
    else
        border_to_crop(3) = 0;
    end
end

data_out(:,:,:,i) = crop_square(squeeze(data_in(:,:,:,i)),border_to_crop(1:3));
end

end

