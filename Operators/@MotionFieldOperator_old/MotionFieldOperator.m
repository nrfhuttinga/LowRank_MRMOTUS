classdef MotionFieldOperator
    
    properties 
        KspaceTrajectory
        NewGrid
        TrajSize
        NumberOfSpatialDimensions
        NufftOptions
        AdjointFlag
        NufftType
        NufftFwdHandle
        NufftAdjDvfHandle
        NufftAdjImHandle
        G_op
    end
    
    methods
        
        % Constructor
        function obj = MotionFieldOperator(KspaceTrajectory,ReferenceGrid,MotionFields,varargin)
            obj.KspaceTrajectory            = KspaceTrajectory;
            obj.AdjointFlag                 = 0;
            obj.NumberOfSpatialDimensions   = size(ReferenceGrid,2);

            obj.NewGrid = zeros(size(MotionFields));
            for i=1:obj.NumberOfSpatialDimensions
            for j=1:size(MotionFields,3)
                obj.NewGrid(:,i,j) = ReferenceGrid(:,i)+MotionFields(:,i,j);
            end
            end
            
            obj.TrajSize            = size(obj.KspaceTrajectory,1);
            
            if nargin>3
                obj.NufftType = varargin{1};
            else
                obj.NufftType = 'finufft';
            end
            
%             disp(['>>>  NUFFT type: ',obj.NufftType])
            
            
            if strcmp(obj.NufftType,'finufft')
                finufft_options.debug                       = 0;
                finufft_options.spread_debug                = 0;
                finufft_options.chkbnds                     = 0;
                finufft_options.nthreads                    = 1;
                finufft_options.fftw                        = 1;
                finufft_options.upsampfac                   = 1.25;
                finufft_options.kerevalmeth                 = 1;
                finufft_options.eps                         = 1e-2;
                obj.NufftOptions = finufft_options;
                
                if obj.NumberOfSpatialDimensions == 3
                    obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d3(2*pi*(new_x),2*pi*(new_y),2*pi*(new_z),double(squeeze(x)),-1,obj.NufftOptions.eps,kspace_traj(:,1),kspace_traj(:,2),kspace_traj(:,3),obj.NufftOptions);
                    obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,new_z,x,dimension)   finufft3d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*kspace_traj(:,3), 2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x  ,1,obj.NufftOptions.eps,new_x,new_y,new_z,obj.NufftOptions);
                    obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*kspace_traj(:,3), double(x),1,obj.NufftOptions.eps,new_x,new_y,new_z,obj.NufftOptions);

                elseif obj.NumberOfSpatialDimensions == 2
                    obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,x)                  finufft2d3(2*pi*(new_x),2*pi*(new_y),double(squeeze(x)),-1,obj.NufftOptions.eps,kspace_traj(:,1),kspace_traj(:,2),obj.NufftOptions);
                    obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,x,dimension)        finufft2d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x,1,obj.NufftOptions.eps,new_x,new_y,obj.NufftOptions);
                    obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,x)                  finufft2d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),double(x),1,obj.NufftOptions.eps,new_x,new_y,obj.NufftOptions);
                end
                    
            elseif strcmp(obj.NufftType,'gpunufft')

                N=round(size(obj.NewGrid,1)^(1/obj.NumberOfSpatialDimensions));
                
                % type 1 = adjoint

                
                if obj.NumberOfSpatialDimensions == 3
%                     obj.G_op=gpuNUFFT3D(-permute(obj.NewGrid(:,[2 1 3],:,:,:)/N,[2 1 4 5 3]),[],N);
                    obj.G_op=gpuNUFFT2D(-permute(obj.NewGrid(:,[2 1 3],:,:,:)/N,[2 1 4 5 3]),N);
                    
                    obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,new_z,x)                  reshape(obj.G_op'*x,[],size(x,2));
                    obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,new_z,x,dimension)        obj.G_op*(2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x); % pointwise mult with sampling mask first
                    obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,new_z,x)                  obj.G_op*(double(x));

                elseif obj.NumberOfSpatialDimensions == 2
                    obj.G_op=gpuNUFFT2D(-permute(obj.NewGrid(:,[2 1],:,:,:)/N,[2 1 4 5 3]),N);
                    obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,x)                  reshape(obj.G_op'*x,[],size(x,2));
                    obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,x,dimension)        obj.G_op*(2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x); % pointwise mult with sampling mask first
                    obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,x)                  obj.G_op*(double(x));
                end

               
                
            elseif strcmp(obj.NufftType,'finufft_cartesian')

                finufft_options.debug                       = 0;
                finufft_options.spread_debug                = 0;
                finufft_options.chkbnds                     = 0;
                finufft_options.nthreads                    = 1;
                finufft_options.fftw                        = 1;
                finufft_options.upsampfac                   = 1.25;
                finufft_options.kerevalmeth                 = 1;
                finufft_options.eps                         = 1e-2;
                obj.NufftOptions = finufft_options;
                
                N=round(size(obj.NewGrid,1)^(1/obj.NumberOfSpatialDimensions));

                if obj.NumberOfSpatialDimensions == 3
                    
%                     obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d3(2*pi*(new_x),2*pi*(new_y),2*pi*(new_z),double(squeeze(x)),-1,obj.NufftOptions.eps,kspace_traj(:,1),kspace_traj(:,2),kspace_traj(:,3),obj.NufftOptions);
%                     obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,new_z,x,dimension)   finufft3d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*kspace_traj(:,3), 2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x  ,1,obj.NufftOptions.eps,new_x,new_y,new_z,obj.NufftOptions);
%                     obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d3(2*pi*kspace_traj(:,1),2*pi*kspace_traj(:,2),2*pi*kspace_traj(:,3), double(x),1,obj.NufftOptions.eps,new_x,new_y,new_z,obj.NufftOptions);
                    
                    
                    obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d1(2*pi*(new_x)/N,2*pi*(new_y)/N,2*pi*(new_z)/N,double(squeeze(x)),-1,obj.NufftOptions.eps,round(length(kspace_traj(:,1)).^(1/3)),round(length(kspace_traj(:,2)).^(1/3)),round(length(kspace_traj(:,3)).^(1/3)),obj.NufftOptions);
                    obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,new_z,x,dimension)   finufft3d2(2*pi*(new_x)/N,2*pi*(new_y)/N,2*pi*(new_z)/N,1,obj.NufftOptions.eps, permute(reshape_to_square(2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x,3),[2,1,3]),obj.NufftOptions);
                    obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,new_z,x)             finufft3d2(2*pi*(new_x)/N,2*pi*(new_y)/N,2*pi*(new_z)/N,1,obj.NufftOptions.eps,permute(reshape_to_square(double(x),3),[2,1,3]),obj.NufftOptions);
                    
                    

                elseif obj.NumberOfSpatialDimensions == 2
                    obj.NufftFwdHandle      = @(kspace_traj,new_x,new_y,x)                  finufft2d1(2*pi*(new_x)/N,2*pi*(new_y)/N,double(squeeze(x)),-1,obj.NufftOptions.eps,round(length(kspace_traj(:,1)).^(1/2)),round(length(kspace_traj(:,2)).^(1/2)));%,obj.NufftOptions);
                    obj.NufftAdjDvfHandle   = @(kspace_traj,new_x,new_y,x,dimension)        finufft2d2(2*pi*(new_x)/N,2*pi*(new_y)/N,1,obj.NufftOptions.eps, reshape_to_square(2*pi*1i*spdiags(kspace_traj(:,dimension),0,obj.TrajSize,obj.TrajSize)*x,2).',obj.NufftOptions);
                    obj.NufftAdjImHandle    = @(kspace_traj,new_x,new_y,x)                  finufft2d2(2*pi*(new_x)/N,2*pi*(new_y)/N,1,obj.NufftOptions.eps,reshape_to_square(double(x),2).',obj.NufftOptions);
                end
                
            end
                
                        
            
        end % end constructor
        
        function s = mtimes(obj,x) % apply motionfield operator to q_ref
            
            if obj.AdjointFlag==0 % Forward operator F(u)*q
                s = zeros(size(obj.KspaceTrajectory,1),size(obj.KspaceTrajectory,3));
                if obj.NumberOfSpatialDimensions==3
                    
                if strcmp(obj.NufftType,'finufft') || strcmp(obj.NufftType,'finufft_cartesian')
                    for i=1:size(s,3)
                        s(:,i) = reshape(permute(obj.NufftFwdHandle(obj.KspaceTrajectory(:,:,i),obj.NewGrid(:,1,i),obj.NewGrid(:,2,i),obj.NewGrid(:,3,i),x(:,1)),[2,1,3]),[],1);
                    end
                elseif strcmp(obj.NufftType,'gpunufft')
%                     for i=1:size(obj.NewGrid,3)
                        s = reshape(obj.NufftFwdHandle(obj.KspaceTrajectory(:,:,1),obj.NewGrid(:,1,1),obj.NewGrid(:,2,1),obj.NewGrid(:,3,1),x(:,1)),[],size(obj.NewGrid,3));
%                     end
                end
                    
                elseif obj.NumberOfSpatialDimensions==2
                    
                    if strcmp(obj.NufftType,'finufft') || strcmp(obj.NufftType,'finufft_cartesian')
                        for i=1:size(s,2)
                            size(obj.NewGrid(:,:,i))
                            size(obj.KspaceTrajectory(:,:,i))
                            s(:,i) = reshape(obj.NufftFwdHandle(obj.KspaceTrajectory(:,:,i),obj.NewGrid(:,1,i),obj.NewGrid(:,2,i),x(:,1)).',[],1);
                        end
                    elseif strcmp(obj.NufftType,'gpunufft')
%                         for i=1:size(obj.NewGrid,3)
                            s  = reshape(obj.NufftFwdHandle(obj.KspaceTrajectory(:,:,1),obj.NewGrid(:,1,1),obj.NewGrid(:,2,1),x),[],size(obj.NewGrid,3));
%                         end
                    end
                    
                end
                
            elseif obj.AdjointFlag==1 % Adjoint operator F_u(u).'*q (w.r.t. motion)
%                 disp('AdjointFlag 1')

                s = zeros(size(obj.NewGrid));
                if obj.NumberOfSpatialDimensions==3 % 3D

                    if strcmp(obj.NufftType,'finufft') || strcmp(obj.NufftType,'finufft_cartesian')
                        for i=1:size(s,2)
                        for jj=1:size(s,3)
                            s(:,i,jj) = obj.NufftAdjDvfHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),obj.NewGrid(:,3,jj),x,i);
                        end
                        end
                    elseif strcmp(obj.NufftType,'gpunufft')
                        for i=1:size(s,2)
                            s(:,i,:) = obj.NufftAdjDvfHandle(obj.KspaceTrajectory(:,:,1),obj.NewGrid(:,1,1),obj.NewGrid(:,2,1),obj.NewGrid(:,3,1),x,i);
                        end
                    end

                elseif obj.NumberOfSpatialDimensions==2 % 2D
                    
                    if strcmp(obj.NufftType,'finufft') || strcmp(obj.NufftType,'finufft_cartesian')
                        for i=1:size(s,2)
                        for jj=1:size(s,3)
                            s(:,i,jj) = obj.NufftAdjDvfHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),x,i);
                        end
                        end
                    elseif strcmp(obj.NufftType,'gpunufft')
                        for i=1:size(s,2)
                            s(:,i,:) = obj.NufftAdjDvfHandle(obj.KspaceTrajectory(:,:,1),obj.NewGrid(:,1,1),obj.NewGrid(:,2,1),x,i);
                        end
%                         for i=1:size(s,2)
%                         for jj=1:size(s,3)
%                             s(:,i,jj) = obj.NufftAdjDvfHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),x,i);
%                         end
%                         end
                    end
                end
                
                    
                   
                
            elseif obj.AdjointFlag==2 % Adjoint operator F_q(u)'*q (w.r.t. images)
%                 disp('AdjointFlag 2')

                s = zeros(size(obj.NewGrid,1),size(obj.NewGrid,3));
                
                if strcmp(obj.NufftType,'finufft') || strcmp(obj.NufftType,'finufft_cartesian')

                    if obj.NumberOfSpatialDimensions==3 % 3D
                        for jj=1:size(s,2)
                            s = s + obj.NufftAdjImHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),obj.NewGrid(:,3,jj),x(:,jj));
                        end
                    elseif obj.NumberOfSpatialDimensions==2 % 2D
                        for jj=1:size(s,2)
                            s = s + obj.NufftAdjImHandle(obj.KspaceTrajectory(:,:,jj),obj.NewGrid(:,1,jj),obj.NewGrid(:,2,jj),x(:,jj));
                        end
                    end
                elseif strcmp(obj.NufftType,'gpunufft')
                    if obj.NumberOfSpatialDimensions==3 % 3D
                        s = obj.NufftAdjImHandle(obj.KspaceTrajectory(:,:,1),obj.NewGrid(:,1,:),obj.NewGrid(:,2,:),obj.NewGrid(:,3,:),x);
                        s = sum(s,2);
                    elseif obj.NumberOfSpatialDimensions==2 % 2D
                        s = obj.NufftAdjImHandle(obj.KspaceTrajectory(:,:,1),obj.NewGrid(:,1,:),obj.NewGrid(:,2,:),x);
                        s = sum(s,2);
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

            