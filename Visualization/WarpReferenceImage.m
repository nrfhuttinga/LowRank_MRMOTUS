function output=WarpReferenceImage(ReferenceImage,MotionFields,N_new, batch_size,nufft_type)
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
    warning('off','all')

    if nargin<5
        nufft_type='gpunufft';
    end

%     ReferenceImage = ForceSquareShape(ReferenceImage);
    if numel(N_new)==1
        N_new = ones(1,size(MotionFields,2))*N_new;
    end
    
    
    NumberOfSpatialDims = size(MotionFields,2);
    
    if nargin<=2 || isempty(N_new)
        N_new = size(ReferenceImage);
    end
    
%     N_new=N_new;
    
    if sum(N_new~=size(ReferenceImage))>0
        
        if NumberOfSpatialDims==2
            ReferenceImage = imresize(ReferenceImage,N_new);
        else
            ReferenceImage = imresize3(ReferenceImage,N_new);
        end
        
    end
    
    
    NumberOfDynamics = size(MotionFields,3);
    N_old = N_new;%round(size(MotionFields,1).^(1/NumberOfSpatialDims));
    
    if nargin<4
        if NumberOfSpatialDims==2
            batch_size = NumberOfDynamics;
        elseif NumberOfSpatialDims == 3
            if N_new <= 80
                batch_size = 8;
            else
                batch_size = 6;
            end
        end
    end

    % ==== warping of reference with highres motion-fields ======
    disp('+     Warping high-res reference with reconstructed motion-fields')
    ReferenceGrid = MRMOTUS_Operator.MakeReferenceGrid(N_new,NumberOfSpatialDims);


    no_batches  = ceil(NumberOfDynamics/batch_size);
    for batch_index = 1:no_batches
        disp(['+Processing batch ',num2str(batch_index),'/',num2str(no_batches),'...'])
        vis_indices = [(batch_index-1)*batch_size+1:min((batch_index-1)*batch_size+batch_size,NumberOfDynamics)];

        if sum(N_old ~= N_new)==0
            disp('Upscaling DVFs');
            mf_highres = UpscaleMotionFields(MotionFields(:,:,vis_indices),N_old,N_new);
        else
            mf_highres = MotionFields(:,:,vis_indices);
        end
        
        
        
%         clearvars resultt
        % warp high res reference image
       
        
        if strcmp(nufft_type,'gpunufft')
            
            MFO = MotionFieldOperator(ReferenceGrid,ReferenceGrid,mf_highres,nufft_type);
             
            warped_kspace = MFO*single(abs(ReferenceImage(:)));

            for rr=1:numel(vis_indices)
    %             warped_kspace_ = reshape_to_square(MotionFieldOperator(ReferenceGrid/N_new,ReferenceGrid,squeeze(mf_highres(:,:,vis_indices(rr))),nufft_type)*single(abs(ReferenceImage(:))),NumberOfSpatialDims);
                warped_kspace_ = reshape_to_square(warped_kspace(:,rr),NumberOfSpatialDims);
    %             warped_kspace_=crop_boundary(warped_kspace_,size(warped_kspace_)/4*2);

    %             warped_kspace_=crop_boundary(reshape_to_square(warped_kspace(:,rr),NumberOfSpatialDims),size(warped_kspace)/4*2);
    %             window = hanning(size(warped_kspace,1))*hanning(size(warped_kspace,2))'.^4;
    %             warped_kspace = warped_kspace.*window;
                resultt{batch_index}(:,:,:,:,rr) = single(abs(ifft_mri(warped_kspace_)));
            end
        else
            
            for rr=1:numel(vis_indices)
                rr
                MFO = MotionFieldOperator(ReferenceGrid,ReferenceGrid,mf_highres(:,:,rr),nufft_type);
             
                warped_kspace = MFO*single(abs(ReferenceImage(:)));
                warped_kspace_ = reshape(warped_kspace,N_new);
                
                resultt{batch_index}(:,:,:,:,rr) = single(abs(ifft_mri(warped_kspace_)));
            end
        end
            

         % store result
        output(:,:,:,:,vis_indices)=resultt{batch_index};

       
        
    end

%     for batch_index = 1:no_batches
%         disp(['+Processing batch ',num2str(batch_index),'/',num2str(no_batches),'...'])
%         vis_indices = [(batch_index-1)*batch_size+1:min((batch_index-1)*batch_size+batch_size,NumberOfDynamics)]; 
% 
%         
%     end
    
    
end


