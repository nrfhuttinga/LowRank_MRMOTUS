function smoothed_dvf=SmoothDVF( dvf , sigma)

if size(dvf,2)==2
    for j=1:size(dvf,3)
    for i=1:size(dvf,2)
        smoothed_dvf(:,i,j) = reshape(imgaussfilt(reshape_to_square(dvf(:,i,j),2),sigma),size(dvf(:,i,j)));
    end
    end
else
    error('not implemented');
end

end
