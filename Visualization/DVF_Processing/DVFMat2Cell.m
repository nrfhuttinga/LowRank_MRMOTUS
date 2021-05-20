function dvf_cell = DVFMat2Cell( dvf_mat )

for i=1:size(dvf_mat,5)
    
    dvf_cell{i}=dvf_mat(:,:,:,:,i);
    
end

end
