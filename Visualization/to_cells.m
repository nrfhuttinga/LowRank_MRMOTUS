function [img_out,cm_out]=to_cells(fig_obj,im_capture_flag)

% make snapshot of 'ax_obj' and return cell with image and colormap
%
% Niek Huttinga - UMC Utrecht - 2020

if nargin<2
    im_capture_flag=0;
end

if im_capture_flag
    img_out = getimage(fig_obj);
    cm_out=[];
else
    frame = getframe(fig_obj);
    [imgs] = frame2im(frame);
    [img_out,cm_out] = rgb2ind(imgs,256,'dither');
end


hold off;  
end
