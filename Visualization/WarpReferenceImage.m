function output=WarpReferenceImage(ReferenceImage,MotionFields)
% Warps reference image 'ReferenceImage' [N_new x N_new (x N_new)] with
% d-dimensional motion-fields 'MotionFields' [N^d x d x T]. This is done by
% generating [N_new x N_new (x N_new) x T] Cartesian k-space data, followed by
% an inverse FFT.
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


    ReferenceImage = ForceSquareShape(ReferenceImage);

    NumberOfSpatialDims = size(MotionFields,2);
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

        mf_highres(:,:,vis_indices) = UpscaleMotionFields(MotionFields(:,:,vis_indices),N_old,N_new);

        clearvars resultt
        % warp high res reference image
        parfor rr=1:numel(vis_indices)
            resultt(:,:,:,:,rr) = single(abs(ifft_mri(reshape_to_square(MotionFieldOperator(ReferenceGrid/N_new,ReferenceGrid,squeeze(mf_highres(:,:,vis_indices(rr))))*single(abs(ReferenceImage(:))),NumberOfSpatialDims))));
        end

        % store result
        output(:,:,:,:,vis_indices)=resultt;
    end
end


