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

    % Remove the two outer splines as we do those separately.
    SplineOrder = SplineOrder - 2;
    % leave to 4 for cubic splines, determines the number of control points that influences every point in the representation
    control_point_overlaps=4;

    % this ensures every spline function has - along a certain dimension - only
    % overlap with four other splines along that dimension
    number_of_nonoverlapping_splines = (SplineOrder - 1)/control_point_overlaps;

    % every spline function is defined on the interval [-2 2], with length 4
    spline_fov = number_of_nonoverlapping_splines*4;


    % specify the center for each spline function along the fov
    spline_centers=linspace(-spline_fov/2,spline_fov/2,SplineOrder);

    % center splines + two outer spines
    B = zeros(N_gridpoints,SplineOrder+2);
    dx = spline_centers(2)-spline_centers(1);
    spline_center_index = 2;
    xs = linspace(-spline_fov/2,spline_fov/2,N_gridpoints);
    % do center splines
    for j=spline_centers
        x1 = xs-j; % 
        B(:,spline_center_index) =  cubic_spline_1d(x1,DerivativeOrder).';%.*deltax^derivative_order;
        spline_center_index = spline_center_index + 1;
    end
    
    % do splines at edge
    x_begin = xs-spline_centers(1);
    x_end = xs-spline_centers(end);
    for k=1:2
        x_begin = x_begin + dx;
        x_end   = x_end - dx;
        B(:, 1)   = B(:, 1)   + cubic_spline_1d(x_begin, DerivativeOrder).';
        B(:, end) = B(:, end) + cubic_spline_1d(x_end, DerivativeOrder).';
    end
end



    

% B=sparse(B);
