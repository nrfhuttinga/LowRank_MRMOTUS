function [BasisMatrix] = ConstructSplineBasis(N_gridpoints,SplineOrder,varargin)
    % Function to construct a 3D cubic spline basis using kronecker products. 

    % Inputs:
    %   N_gridpoints    - cell with number of grid points per dimension as elements, [x y z] order
    %   splineOrder     - cell with spline order per dimension as elements, [x y z] order
    %   NoDimensions    - number of dimensions to build spline basis on
    %   varargin{1}     - Number of dimensions (2/3d)
    %   varargin{2}     - SpatialOrdering [x y z] order in the input data
    %   varargin{3}     - Order of the spatial derivatives for each dimension (default = 0)
    %
    % Outputs 
    %   BasisMatrix     - Matrix with the spline basis [N^d x prod(splineOrder)]
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

    
    if nargin >= 5
        DerivativeOrder = varargin{3};
        SpatialOrdering = varargin{2};
        NoDimensions = varargin{1};
    elseif nargin >= 4
        SpatialOrdering = varargin{2};
        NoDimensions = varargin{1};
        DerivativeOrder = zeros(1,NoDimensions);
    elseif nargin >= 3
        SpatialOrdering = [2 1 3 4];
        NoDimensions = varargin{1};
        DerivativeOrder = zeros(1,NoDimensions);
    elseif nargin <= 2
        SpatialOrdering = [2 1 3 4];
        NoDimensions = 3;
        DerivativeOrder = zeros(1,NoDimensions);
    end

    spline_function_handle = cell(NoDimensions,1);
    for i=1:NoDimensions
        spline_function_handle{i}=@(N_gridpoints_1D,SplineOrder_1D) spline_basis_1D(N_gridpoints_1D,SplineOrder_1D,DerivativeOrder(SpatialOrdering(i)));

    end



    for k=1:NoDimensions
    % first dimension that is added in the kronecker product is vertical, so we 
    % first add y basis, then x then z. so we swap the orders for this to [2 1 3] 
    % because N_gridpoints and SplineOrder are [x y z]
    %     lattice{SpatialOrdering(k)}=lattice{SpatialOrdering(k)}/SplineOrder{SpatialOrdering(k)};
        if k==1
            BasisMatrix = sparse(spline_function_handle{k}(N_gridpoints(SpatialOrdering(k)),SplineOrder(SpatialOrdering(k))));
        else
            BasisMatrix = sparse(kron(sparse(spline_function_handle{k}(N_gridpoints(SpatialOrdering(k)),SplineOrder(SpatialOrdering(k)))),BasisMatrix));
        end


    end
end

