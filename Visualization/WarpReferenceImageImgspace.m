function output=WarpReferenceImageImgspace( ReferenceImage , MotionFields )

% Warps reference image 'ReferenceImage' [N_new x N_new (x N_new)] with
% d-dimensional motion-fields 'MotionFields' [N^d x d x T]. 
%
%
% Inputs:
%   ReferenceImage      - [N_new x N_new (xN_new)]  complex image
%   MotionFields        - [N^d x d x T] motion-fields
%
% Output:
%   output: [N_new x N_new x N_new / 1 x 1 x T] image sequence.
%
% Niek Huttinga - 2020 - UMC Utrecht

       
%     if ~isreal(ReferenceImage)
%         warning('Taking absolute value of ReferenceImage...')
%         ReferenceImage = abs(ReferenceImage);
%     end
       
    addpath(genpath('/nfs/rtsan02/userdata/home/nhutting/Documents/git/MRMOTUS/utilities/Undersampling/'));
    
    NumberOfSpatialDims = size(MotionFields,2);
    ReferenceImage = reshape_to_square(ReferenceImage,NumberOfSpatialDims);
    
    

    NumberOfDynamics = size(MotionFields,3);
    N_old = round(size(MotionFields,1).^(1/NumberOfSpatialDims));
    N_new = size(ReferenceImage,1);

    if NumberOfSpatialDims==2
        batch_size = NumberOfDynamics;
    elseif NumberOfSpatialDims == 3
        if N_new <= 80
            batch_size = 8;
        else
            batch_size = 6;
        end
    end

    % ==== warping of reference with highres motion-fields ======
    disp('+     Warping high-res reference with reconstructed motion-fields')
    ReferenceGrid = MRMOTUS_Operator.MakeReferenceGrid(N_new,NumberOfSpatialDims);


    no_batches  = ceil(NumberOfDynamics/batch_size);
    for batch_index = 1:no_batches
        disp(['+Processing batch ',num2str(batch_index),'/',num2str(no_batches),'...'])
        vis_indices = [(batch_index-1)*batch_size+1:min((batch_index-1)*batch_size+batch_size,NumberOfDynamics)];

        if N_old ~= N_new
            mf_highres(:,:,vis_indices) = UpscaleMotionFields(MotionFields(:,:,vis_indices),N_old,N_new);
        else
            mf_highres(:,:,vis_indices) = MotionFields(:,:,vis_indices);
        end
        
        
        
        clearvars resultt
        % warp high res reference image
        if NumberOfSpatialDims==2
            for rr=1:numel(vis_indices)
                new_grid = ReferenceGrid + squeeze(mf_highres(:,:,vis_indices(rr)));
                resultt(:,:,:,:,rr) = crop_boundary_zeroes( single(abs(interp2(reshape_to_square(ReferenceGrid(:,1),NumberOfSpatialDims),reshape_to_square(ReferenceGrid(:,2),NumberOfSpatialDims),ReferenceImage,reshape_to_square(new_grid(:,1),NumberOfSpatialDims),reshape_to_square(new_grid(:,2),NumberOfSpatialDims),'cubic'))),size(ReferenceImage)-2);
            end
        elseif NumberOfSpatialDims==3
            for rr=1:numel(vis_indices)
                new_grid = ReferenceGrid + squeeze(mf_highres(:,:,vis_indices(rr)));
                resultt(:,:,:,:,rr) = single(abs(interp3(reshape_to_square(ReferenceGrid(:,1),NumberOfSpatialDims),reshape_to_square(ReferenceGrid(:,2),NumberOfSpatialDims),reshape_to_square(ReferenceGrid(:,3),NumberOfSpatialDims),ReferenceImage,reshape_to_square(new_grid(:,1),NumberOfSpatialDims),reshape_to_square(new_grid(:,2),NumberOfSpatialDims),reshape_to_square(new_grid(:,3),NumberOfSpatialDims),'cubic')));
                resultt(:,:,:,:,rr) = crop_boundary_zeroes(resultt(:,:,:,:,rr),size(ReferenceImage)-2);
            end
        end
        % store result
        output(:,:,:,:,vis_indices)=resultt;
    end
    
    output(isnan(output))=0;
    det_flag=0;
    if det_flag
        
        output = output .* reshape( DeterminantMotionFields(mf_highres), size(output));
    end
    
end


