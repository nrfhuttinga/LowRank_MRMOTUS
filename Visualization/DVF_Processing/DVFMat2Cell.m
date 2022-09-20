function dvf_cell = DVFMat2Cell( dvf_mat , varargin)

if nargin<2
    dyn_dim = numel(size(dvf_mat));
else
    dyn_dim = varargin{1};
end

for i=1:size(dvf_mat,dyn_dim)
    
    dvf_cell{i}=SliceData(dvf_mat,dyn_dim,i);
    
end

end
