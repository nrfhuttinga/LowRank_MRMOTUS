function BasisHandle=ConstructBasisHandle(options)
%  Construct a basis handle to analytically compute (partial/mixed) spatial
%  spline derivatives up to 2nd order.
%
%  B=BasisHandle(x_ord,y_ord[,z_ord]) returns a basis B with cubic b-spline
%  basis functions that are differentiated x_ord,y_ord[,z_ord] times in x,y[,z]
%  respectively

%  Inputs: options struct with fields
%      spline_orders    -   Number of splines per dimension, for
%                           multi-resolution spline basis use a vector here  
%      spline_orders_z  -   Optionally, you can specify different spline order(s) in the z direction
%      N                -   Image dimensions
%      spatial_ordering -   Don't change this
%      dimension        -   Number of spatial dimensions [2/3]
%
%
%   Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

if ~isfield(options,'spline_orders_z')
    options.spline_orders_z = options.spline_orders;
end

% options.spatial_ordering = [2 1 3 4];

options=set_default(options,'spatial_ordering',[2 1 3 4]);

if numel(options.N)==1
    options.N = repmat(options.N,1,options.dimension);
end

    for i=1:numel(options.spline_orders)

        if options.dimension == 1 %1D
            
            Bases{i} = @(x_order) (ConstructSplineBasis([options.N(1)],[options.spline_orders(i)],options.dimension,options.spatial_ordering,x_order));
        
        elseif options.dimension == 2 %2D

            Bases{i}    = @(x_order,y_order) (ConstructSplineBasis([options.N(1),options.N(2)],[options.spline_orders(i),options.spline_orders(i)],options.dimension,options.spatial_ordering,[x_order,y_order]));
%             Bases{i}    = @(x_order,y_order) speye(options.N(1)^2,options.N(1)^2);
        
        elseif options.dimension==3 %3D
            
           
            if options.spline_orders(i)==options.N(i)
                Bases{i}    = @(x_order,y_order, z_order) speye(options.N(i)^3,options.N(1)^3);
            else
                Bases{i}    = @(x_order,y_order,z_order) (ConstructSplineBasis([options.N(1),options.N(2),options.N(3)],[options.spline_orders(i),options.spline_orders(i),options.spline_orders_z(i)],options.dimension,options.spatial_ordering,[x_order,y_order,z_order]));
            end
        end
        
        
            
    end
    
    if options.dimension == 1
        BasisHandle = @(x_order) sparse(cell2mat(cellfun(@(f) f(x_order),Bases,'UniformOutput',0)));
    elseif options.dimension == 2
        BasisHandle = @(x_order,y_order) sparse(cell2mat(cellfun(@(f) f(x_order,y_order),Bases,'UniformOutput',0)));
    elseif options.dimension == 3
        BasisHandle = @(x_order,y_order,z_order) sparse(cell2mat(cellfun(@(f) f(x_order,y_order,z_order),Bases,'UniformOutput',0)));
    else
        error('This number of dimensions is not implemented!');
    end
    
