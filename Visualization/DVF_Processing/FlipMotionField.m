function MF_flipped=FlipMotionField( MF , flips )

    % MF    = [N^d x d x T ]
    % flip  = [ 0/1 ... 0/1] \in \mathbb{R}^d

    
    d = numel(flips);

    ordering = [2 1 3];
    flips_mf = flips(ordering(1:d));
    
    
    if d~=size(MF,2)
        error('Wrong input dimensions: nums in second arg ~= size of second dimension in first arg');
    end

    dims = 1:d;
    N = round((size(MF,1).^(1/d)));
%     MF_flipped = zeros(size(MF));
    
%     for i=dims % loop over motion-field dimension (1st 2nd 3rd)
%         if flips(i) % if we have to flip this dimension, then flip all directions, and all dynamics
%             for k=dims % loop over motion-field direction (x,y,z)
%                 MF_dim = MF(:,k,:);
%                 for j=1:size(MF_dim,3) % loop over all dynamics
%                     MF_dim_t = reshape_to_square(MF_dim(:,:,j),d);
%                     if k~=i % just flip 
%                         MF_flipped(:,k,j) = matrix_to_vec(flip(MF_dim_t,i));
%                     else
%                         % if the motion-field direction is equal to the dimension we flip on, then we need to invert the direction aswell
%                         MF_flipped(:,k,j) = -matrix_to_vec(flip(MF_dim_t,i));
%                     end
%                 end
%             end
%         else
%             MF_flipped(1:N+(i-1)*N,:,:) = MF(1:N+(i-1)*N,:,:);
%         end
%     end

    MF_flipped = reshape(MF,[ones(1,d)*N,d,size(MF,3)]);
    for i=dims
        if flips(i)
            MF_flipped  = flip(MF_flipped,i);
        end
    end
    
    MF_flipped = reshape(MF_flipped,size(MF));
    
    MF_flipped(:,flips_mf>0,:)=-MF_flipped(:,flips_mf>0,:);
    
    
        
                    

              
    
end



