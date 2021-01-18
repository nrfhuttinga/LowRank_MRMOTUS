function DET_rec = DeterminantMotionFields(MotionFields,varargin)
    % Function to compute determinants of motion-fields with d spatial dimensions and T dynamics
    %
    % Inputs:
    %   MotionFields        - Motionfields [N^d x d x T ]
    %   varargin{1}         - Visualization flag [0/1]
    %
    % Outputs:
    %   DET_rec             - Jacobian determinants of 'MotionFields' [N x N (x N) x T]
    %
    %
    % Niek Huttinga, UMC Utrecht, 2020.

    dim=size(MotionFields,2);



    if dim==3
        N = round(size(MotionFields,1).^(1/dim));
        [rx,ry,rz]=meshgrid([-N/2:N/2-1],[-N/2:N/2-1],[-N/2:N/2-1]);
        ReferenceGrid = [rx(:),ry(:),rz(:)];

        ds=[1 1 1];
        for i=1:size(MotionFields,3)
            [dx1dx,dx1dy,dx1dz] = gradient( reshape(MotionFields(:,1,i) + ReferenceGrid(:,1) ,N,N,N)*ds(1), ds(1),ds(2),ds(3));
            [dy1dx,dy1dy,dy1dz] = gradient( reshape(MotionFields(:,2,i) + ReferenceGrid(:,2) ,N,N,N)*ds(2), ds(1),ds(2),ds(3));
            [dz1dx,dz1dy,dz1dz] = gradient( reshape(MotionFields(:,3,i) + ReferenceGrid(:,3) ,N,N,N)*ds(3), ds(1),ds(2),ds(3));
            DET_rec(:,:,:,i) =  dx1dx.*(dy1dy.*dz1dz-dy1dz.*dz1dy)- ...
                                dx1dy.*(dy1dx.*dz1dz-dy1dz.*dz1dx)+ ...
                                dx1dz.*(dy1dx.*dz1dy-dy1dy.*dz1dx);
        end

    if nargin>1
        if varargin{1}
            slicer5d(DET_rec);
        end
    end        


    elseif dim==2
        N = round(size(MotionFields,1).^(1/dim));
        [rx,ry]=meshgrid([-N/2:N/2-1],[-N/2:N/2-1]);
        ReferenceGrid = [rx(:),ry(:)];

        ds=[1 1];


        for i=1:size(MotionFields,3)
            [dx1dx,dx1dy] = gradient( reshape(MotionFields(:,1,i),N,N), ds(1),ds(2));
            [dy1dx,dy1dy] = gradient( reshape(MotionFields(:,2,i),N,N), ds(1),ds(2));
            DET_rec(:,:,i) =  (dx1dx+1).*(dy1dy+1)-dx1dy.*dy1dx;
        end
    end


end

