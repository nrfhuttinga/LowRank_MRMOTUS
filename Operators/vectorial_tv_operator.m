function [tv_total,gradient_total] = vectorial_tv_operator(x,eps_TV,D)
    % Compute vectorial total variation and its gradients for all columns in x
    %
    %
    %
    % Inputs:
    %   x           - input over which to compute vectorial TV [voxels spatial_dimensions dynamics]
    %   eps_TV      - epsilon for smooth TV approximation
    %   D           - cell {x,y,z} with finite differences matrices per spatial dimension as element
    %
    % Outputs:
    %   tv_total        - total variation of the input
    %   gradient_total  - corresponding gradients
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
    
    grad_out_size = size(x);


    x = reshape(x,size(x,1),numel(D),[]);
    gradient_vector = cell(numel(x),1);



    
for j=1:size(x,3)
    pointwise_l2_gradient_pooled{j} = 0;
    
for i=1:size(x,2)
    
    
    data = squeeze(x(:,i,j));
    pointwise_l2_gradient{i,j} = 0;
    pointwise_l2_gradient_eps{i,j} = 0;
    
    
for k=1:numel(D)
    gradient_vector{k} = D{k}*data;
    
    pointwise_l2_gradient{i,j} = pointwise_l2_gradient{i,j} + gradient_vector{k}.^2; % pointwise l2 norm squared over k of gradient vector
    

end
    pointwise_l2_gradient_eps{i,j} = sqrt(pointwise_l2_gradient{i,j}+eps_TV); % sqrt( pointwise l2 norm squared + eps ) ~~ abs( pointwise l2 norm )
    
    pointwise_l2_gradient_pooled{j} = pointwise_l2_gradient_pooled{j} + pointwise_l2_gradient{i}+eps_TV; % l2 norm squared over i of pointwise_l2_gradient_eps{i,j} 
    
    

end
    pointwise_l2_gradient_pooled{j} = sqrt(pointwise_l2_gradient_pooled{j}); % sqrt( l2 norm squared over i of pointwise_l2_gradient_eps{i,j} ) = l2 norm over i of pointwise_l2_gradient_eps{i,j}
end


gradient_total=[];
tv_total=0;
for j=1:size(x,3)
    tv_total = tv_total + sum(pointwise_l2_gradient_pooled{j});
    
    for i=1:size(x,2)
        data = squeeze(x(:,i,j));
        gradient_ijk = 0;

        for k=1:numel(D)
            gradient_vector{k} = D{k}*data;

            gradient_ijk = gradient_ijk +  D{k}'* ( gradient_vector{k} ./ pointwise_l2_gradient_eps{i,j} ./ pointwise_l2_gradient_pooled{j} ) ; % THE 2 IS FOR DERIVATIVE OF SQUARED L2 NORM!!
        end
        gradient_total = [gradient_total;gradient_ijk];
    end
end


gradient_total = reshape(gradient_total,grad_out_size);



end