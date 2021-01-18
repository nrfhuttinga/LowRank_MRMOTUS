function dcf = iterative_dcf_estimation(traj,varargin)
% Iteratively estimate density for arbitrary k-spaces in 3D. This code
% wraps code downloaded from the ISMRM-unbound website and written by Nick 
% Zwart.
% 
% Input 
%   traj        - k-space trajectory [2/3 x #readout_samples x #interleaves]
%   varargin{1} - number of iterations to perform
%   varargin{2} - over sampling factor
% 
% Output
%   dcf         - Density compensation function [#readout_samples x #interleaves]
%
%
% Wrapper by T. Bruijnen, UMC Utrecht.

if size(traj,2)==2 || size(traj,2)==3
    traj = permute(traj,[2 1 3 4 5]);
end


% Check input
kdim=c12d(size(traj));

% Provide number of iter in varargin{1}
if ~isempty(varargin)
    n_iter=varargin{1};
else
    n_iter=6;
end

if nargin > 2
    osf=varargin{2};
else
    osf=2.1;
end

% Image space dimensions
if kdim(4) == 1
    idim(1:3)=ceil(max(abs(matrix_to_vec(traj([1 2],:)))));
else
    %idim(1:2)=ceil(max(abs(matrix_to_vec(traj([1 2],:)))));
    %idim(3)=round(2*max(abs(traj(3,:))));
    idim(1:3)=ceil(max(abs(matrix_to_vec(traj([1 2],:)))));
end

% Scale k-space to [-.5, .5]
traj=.5*traj/max(abs(traj(:)));

% Loop over partitions and compute k-space
for p=1:prod(kdim(5:end)) % Loop over partitions
    dcf(:,:,:,p)=sdc3_MAT(traj(:,:,:,:,p),n_iter,max(idim),0,osf);
end

disp('+     ... Density compensation function iteratively estimated.')

% END
end