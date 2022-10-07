function B = spline_basis_1D(N_gridpoints,SplineOrder,DerivativeOrder)
% return a basis B with spline_order cubic b-spline functions as columns, 
% x can be a vector


% Example:
%   B = spline_basis_1D(20,6,0);
%   figure;plot(B);title('Spline basis functions');legend();
%   figure;plot(sum(B,2));title('Sum over all basis functions per grid
%   point');
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

if nargin<3 || isempty(DerivativeOrder)
    DerivativeOrder = 0;
end


if DerivativeOrder > 2
    error('Order > 2 not implemented')
end

% leave to 4 for cubic splines, determines the number of control points that influences every point in the representation
control_point_overlaps=4;

% this ensures every spline function has - along a certain dimension - only
% overlap with four other splines along that dimension
number_of_nonoverlapping_splines = (SplineOrder - 1)/control_point_overlaps;

% every spline function is defined on the interval [-2 2], with length 4
spline_fov = number_of_nonoverlapping_splines*4;


% specify the center for each spline function along the fov
spline_centers=linspace(-spline_fov/2,spline_fov/2,SplineOrder);


B = zeros(N_gridpoints,length(spline_centers));
dx = spline_centers(2)-spline_centers(1);
deltax = 1/N_gridpoints;

spline_center_index = 1;
for j=spline_centers
    x1 = linspace(-spline_fov/2,spline_fov/2,N_gridpoints)-j; % 
    B(:,spline_center_index) =  (cubic_spline_1d(x1,DerivativeOrder).').*deltax^DerivativeOrder;%
    
    if spline_center_index==1
        for k=1:2
            x1 = x1 + dx; % 
            B(:,spline_center_index) = B(:,spline_center_index) + (cubic_spline_1d(x1,DerivativeOrder).')*deltax^DerivativeOrder;%.*deltax^derivative_order;
        end
    end
    
    if spline_center_index==length(spline_centers)
        for k=1:2
            x1 = x1-dx; % 
            B(:,spline_center_index) = B(:,spline_center_index) + (cubic_spline_1d(x1,DerivativeOrder).')*deltax^DerivativeOrder;%.*deltax^derivative_order;
        end
    end
    
    spline_center_index = spline_center_index+1;
end




end



    

% B=sparse(B);