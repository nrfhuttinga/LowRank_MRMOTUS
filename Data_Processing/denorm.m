function res = denorm(x)
% divide input data by the norm
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
res=(x/norm(x(:)));

end