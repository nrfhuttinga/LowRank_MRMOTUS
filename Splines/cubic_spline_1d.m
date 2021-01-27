function y=cubic_spline_1d(x,derivative_order)
% returns cubic spline
% Inputs:
%   - x: grid point values
%   - derivative_order: derivative of the spline function to return
%
% Output:
%   - y: cubic b-spline differentiated up to derivative_order, evaluated at
%   grid points x.
%
% Sanity checks with gradient function: 
%   1st order derivative: 
%       figure;
%       plot(cubic_spline_1d(linspace(-2,2,20),1));
%       hold on;
%       plot(gradient(cubic_spline_1d(linspace(-2,2,20),0),4/19))
%   2nd order: 
%        figure,plot(cubic_spline_1d(linspace(-2,2,20),2));hold on;plot(gradient(gradient(cubic_spline_1d(linspace(-2,2,20),0),4/19),4/20))
%
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

ind1 = x>=-2 & x<-1;
ind2 = x>=-1 & x< 0;
ind3 = x>=0  & x< 1;
ind4 = x>=1  & x< 2;

dx = (max(x(:))-min(x(:)))/(numel(x)-1);

y=zeros(size(x));

if derivative_order>2
    error('Derivative order > 2 not implemented...')
end

switch derivative_order
    case 0
        y(ind1) = ((x(ind1)+2).^3)/6;
        y(ind2) = (-x(ind2).^3-2*(x(ind2)+1).^3+6*(x(ind2)+1))/6;
        y(ind3) =  (x(ind3).^3+2*(x(ind3)-1).^3-6*(x(ind3)-1))/6;
        y(ind4) = ((2-x(ind4)).^3)/6;
    case 1
%         y(ind1) = ((x(ind1)+2).^2)/2;
%         y(ind2) = -0.5*(x(ind2).*(3.*x(ind2)+4));
%         y(ind3) =  0.5*(x(ind3).*(3.*x(ind3)-4));
%         y(ind4) = -0.5*(x(ind4)-2).^2;
        y       = gradient(cubic_spline_1d(x,0),dx);
    case 2
%         y(ind1) =  x(ind1)+2;
%         y(ind2) = -3*x(ind2)-2;
%         y(ind3) =  3*x(ind3)-2;
%         y(ind4) = 2-x(ind4);
        y       = gradient(cubic_spline_1d(x,1),dx);
end
        
        

end



