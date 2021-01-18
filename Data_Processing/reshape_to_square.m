function image=reshape_to_square(image,varargin)
% Reshapes vector 'image' to a square
%
% Inputs
%   image           - input vector N*d x 1
%   varargin{1}     - dimensionality 'd' of the output image (2/3)
%   
% Outputs
%   image           - N x N (xN)
%
% Niek Huttinga, UMC Utrecht, 2020

if nargin>1
    dimension = varargin{1};
else
    dimension = 3;
end

image = reshape(image,round(numel(image).^(1/dimension))*ones(1,dimension));

end
