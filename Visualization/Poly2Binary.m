function [mask_total] = Poly2Binary(image)

N = size(image,1);

if numel(size(squeeze(image)))>2
    for i=1:3
        mask_2D = Poly2Binary_2D(squeeze(sum(abs(image),i)));

        switch i
            case 1
                for j=1:N
                    mask_3D{i}(j,:,:)=mask_2D;
                end
            case 2
                for j=1:N
                    mask_3D{i}(:,j,:)=mask_2D;
                end
            case 3 
                for j=1:N
                    mask_3D{i}(:,:,j)=mask_2D;
                end
        end
    end

    mask_total = mask_3D{1}+mask_3D{2}+mask_3D{1}>=3;

else
    
    mask_total = Poly2Binary_2D(squeeze(abs(image)));
end





% mask_coronal_3D=
% mask_sagittal_3D=
% mask_transverse_3D=