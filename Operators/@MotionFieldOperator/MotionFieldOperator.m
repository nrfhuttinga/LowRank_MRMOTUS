classdef MotionFieldOperator
    % Main warping object
    %
    % Niek Huttinga - UMC Utrecht - 2020
    
    properties 
        KspaceTrajectory
        NewGrid
        TrajSize
        NumberOfSpatialDimensions
        NufftOptions
        AdjointFlag
        NufftFwdHandle
        NufftAdjDvfHandle
        NufftAdjImHandle
    end
    
    methods
        
        % Constructor
        function obj = MotionFieldOperator(KspaceTrajectory,ReferenceGrid,MotionFields,varargin)
            
            % Intialize internal variables
            obj.KspaceTrajectory            = KspaceTrajectory;
            obj.AdjointFlag                 = 0;
            obj.NumberOfSpatialDimensions   = size(ReferenceGrid,2);
            obj.NewGrid = zeros(size(MotionFields));
            obj.TrajSize            = size(obj.KspaceTrajectory,1);
            
            % Apply motion-fields to warp the grid
            for i=1:obj.NumberOfSpatialDimensions
            for j=1:size(MotionFields,5)
                obj.NewGrid(:,i,j) = ReferenceGrid(:,i)+squeeze(MotionFields(:,i,j));
            end
            end
            
            % Initialize FINUFFT options            
            finufft_options.debug                       = 0;
            finufft_options.spread_debug                = 0;
            finufft_options.chkbnds                     = 0;
            finufft_options.nthreads                    = 1;
            finufft_options.fftw                        = 1;
            finufft_options.upsampfac                   = 2;
            finufft_options.kerevalmeth                 = 1;
            finufft_options.eps                         = 1e-6;
            obj.NufftOptions = finufft_options;

            % Create function handles for convenience, you could replace
            % these with your own
            if obj.NumberOfSpatialDimensions == 3
                obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d3(2*pi*(new_x),2*pi*(new_y),2*pi*(new_z),double(squeeze(x)),-1,obj.NufftOptions.eps,kspace_traj(:,1),kspace_traj(:,2),kspace_traj(:,3),obj.NufftOptions);
                obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,new_z,x,dimension)   finufft3d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*kspace_traj(:,3), 2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x  ,1,obj.NufftOptions.eps,new_x,new_y,new_z,obj.NufftOptions);
                obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*kspace_traj(:,3), double(x),1,obj.NufftOptions.eps,new_x,new_y,new_z,obj.NufftOptions);

            elseif obj.NumberOfSpatialDimensions == 2
                obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,x)                  finufft2d3(2*pi*(new_x),2*pi*(new_y),double(squeeze(x)),-1,obj.NufftOptions.eps,kspace_traj(:,1),kspace_traj(:,2),obj.NufftOptions);
                obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,x,dimension)        finufft2d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x,1,obj.NufftOptions.eps,new_x,new_y,obj.NufftOptions);
                obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,x)                  finufft2d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),double(x),1,obj.NufftOptions.eps,new_x,new_y,obj.NufftOptions);
            end
                    
        end % end constructor
        
        function s = mtimes(obj,x) % Apply motionfield operator to the reference image q
            
            if obj.AdjointFlag==0 % Forward mode F(u)*q
                s = zeros(size(obj.KspaceTrajectory,1),size(obj.KspaceTrajectory,3));
                
                if obj.NumberOfSpatialDimensions==3     % 3D
                    for i=1:size(s,3)
                        s(:,i) = obj.NufftFwdHandle(obj.KspaceTrajectory(:,:,i),obj.NewGrid(:,1,i),obj.NewGrid(:,2,i),obj.NewGrid(:,3,i),x(:,1));
                    end
                
                elseif obj.NumberOfSpatialDimensions==2 % 2D
                    for i=1:size(s,2)
                        s(:,i) = obj.NufftFwdHandle(obj.KspaceTrajectory(:,:,i),obj.NewGrid(:,1,i),obj.NewGrid(:,2,i),x(:,1));
                    end
                end
                
                s = s(:);
                
            elseif obj.AdjointFlag==1 % Adjoint mode 1 F_u(u).'*q (w.r.t. motion-fields)

                s = zeros(size(obj.NewGrid));
                x = reshape(x,[],size(s,3));
                if obj.NumberOfSpatialDimensions==3     % 3D
                    for i=1:size(s,2)
                    for jj=1:size(s,3)
                        s(:,i,jj) = obj.NufftAdjDvfHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),obj.NewGrid(:,3,jj),x(:,jj),i);
                    end
                    end
                elseif obj.NumberOfSpatialDimensions==2 % 2D
                    for i=1:size(s,2)
                    for jj=1:size(s,3)
                        s(:,i,jj) = obj.NufftAdjDvfHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),x(:,jj),i);
                    end
                    end
                    
                end
                
            elseif obj.AdjointFlag==2 % Adjoint mode 2 F_q(u)'*q (w.r.t. images)

                s = zeros(size(obj.NewGrid,1),1);
                x = reshape(x,[],size(obj.NewGrid,3));
                if obj.NumberOfSpatialDimensions==3 % 3D
                    for jj=1:size(obj.NewGrid,3)
                        s = s + obj.NufftAdjImHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),obj.NewGrid(:,3,jj),x(:,jj));
                    end
                elseif obj.NumberOfSpatialDimensions==2 % 2D
                    for jj=1:size(obj.NewGrid,3)
                        s = s + obj.NufftAdjImHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),x(:,jj));
                    end
                end
                
                
            end
            
            
               
        end % end mtimes
        
        function obj = transpose(obj) % apply jacobian^T w.r.t. motionfields to a vector
            
            obj.AdjointFlag = 1;

        end % end transpose
        
        function obj = ctranspose(obj) % apply jacobian^T w.r.t. reference image to a vector
            
            obj.AdjointFlag = 2;
            
        end % end ctranspose
        
        
        
    end
    
end

            