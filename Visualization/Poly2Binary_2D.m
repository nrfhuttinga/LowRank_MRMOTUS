function bin_mask = Poly2Binary_2D(image,shape)

if nargin<2
    shape='poly';
end


figure;imagesc(abs(image)); axis image;set_figure_fullscreen;drawnow;

switch shape
    case 'poly'
        h = impoly();
        mask_coords = getPosition(h);
    case 'rect'
        h = imrect();
        mask_coords = bbox2points(getPosition(h));
end



bin_mask = poly2mask(mask_coords(:,1),mask_coords(:,2),size(image,1),size(image,2));
figure;imagesc(bin_mask);axis image;drawnow;

end