function mf_upscaled = UpscaleMotionFields(mf,N_old,N_new)
    % Upscale d-dimensional motion-field mf [N_old^d x d x T] from N_old voxels to N_new
    % voxels per spatial dimension
    %
    % Inputs:
    %   mf          - input motion-fields [N_old^d x d x T]
    %   N_old       - old image dimensions [integer]
    %   N_new       - new image dimensions [integer>N_old]
    %
    % Outputs
    %   mf_upscaled - upscaled motion-fields [N_new^d x d x T]
    %
    % Niek Huttinga, UMC Utrecht, 2020

    mf = mf*(N_new/N_old);
    NumberOfSpatialDims = size(mf,2);
    % interpolate to higher resolution
    for ddd=1:size(mf,3)
        for sss=1:NumberOfSpatialDims
            if NumberOfSpatialDims==2
                mf_upscaled(:,sss,ddd) = single(reshape(imresize(reshape_to_square(mf(:,sss,ddd),NumberOfSpatialDims),ones(1,NumberOfSpatialDims)*N_new,'lanczos3'),[],1));
                sum(reshape(isnan(mf_upscaled(:,sss,ddd)),[],1))
            elseif NumberOfSpatialDims==3
                mf_upscaled(:,sss,ddd) = single(reshape(imresize3(reshape_to_square(mf(:,sss,ddd),NumberOfSpatialDims),ones(1,NumberOfSpatialDims)*N_new,'linear'),[],1));
            end
%             mf_upscaled(:,sss,ddd) = demax(mf_upscaled(:,sss,ddd)-min(mf_upscaled(:,sss,ddd)))*(max(mf(:,sss,ddd))-min(mf(:,sss,ddd)))+min(mf(:,sss,ddd));
        end
    end

    


end
