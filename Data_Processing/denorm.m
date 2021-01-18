function res = denorm(x)
% divide input data by the norm
%
% Niek Huttinga, UMC Utrecht, 2020
res=(x/norm(x(:)));

end