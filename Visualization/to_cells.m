function [img_out,cm_out]=to_cells(fig_obj)
% make snapshot of 'ax_obj' and return cell with image and colormap
%
% Niek Huttinga - UMC Utrecht - 2020

frame = getframe(fig_obj);%,[547,120,891,891]);
[imgs] = frame2im(frame);
hold off;  
[img_out,cm_out] = rgb2ind(imgs,256);