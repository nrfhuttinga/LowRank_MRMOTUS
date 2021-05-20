function rotatedMF = RotateMotionField2D( motionfield , no_rotations )
% rotate motionfield 'no_no_rotations' counterclockwise

no_rotations = mod(no_rotations,4);

if iscell(motionfield)
    cell_flag = 1;
    motionfield = DVFCell2Mat(motionfield);
%     rotatedMF = DVFMat2Cell(RotateMotionField2D( motionfield, no_rotations));
else
    cell_flag = 0;
end



if size(motionfield,5)==1

    if no_rotations == 1 % counterclockwise
        % x -> y
        % y -> -x


    %         for i=1:dynamics
            mf_new(:,:,1)=motionfield(:,:,2);
            mf_new(:,:,2)=-motionfield(:,:,1);
            rotatedMF=rot90(mf_new,no_rotations);
    %         end



    elseif no_rotations == 2
        % x -> y -> -x
        % y -> -x -> -y

            mf_new(:,:,:,2)=-motionfield(:,:,2);
            mf_new(:,:,:,1)=-motionfield(:,:,1);
            rotatedMF=rot90(mf_new,no_rotations);

    elseif no_rotations == 3
        % x -> y -> -x -> -y
        % y -> -x -> -y -> x

            mf_new(:,:,:,1)=-motionfield(:,:,2);
            mf_new(:,:,:,2)=motionfield(:,:,1);
            rotatedMF=rot90(mf_new,no_rotations);

    elseif no_rotations == 0
        rotatedMF = motionfield;
    end

else

    for i=1:size(motionfield,5)
        rotatedMF(:,:,:,:,i)=RotateMotionField2D(motionfield(:,:,:,:,i),no_rotations);
    end
end
    
if cell_flag
    rotatedMF = DVFMat2Cell(rotatedMF);
end
    

    


end
