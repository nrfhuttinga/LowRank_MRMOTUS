classdef MRMOTUS_Operator
    % Main class for MR-MOTUS reconstructions: performs all required initializations and contains
    % forward operators, regularization operators, and all corresponding gradients
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
    
    properties 
        ReferenceImage
        KspaceTrajectory
        ReferenceGrid
        DiagReferenceImage
        NumberOfDynamics
        NumberOfSpatialDims 
        ImDims
        referenceImageMask
        RegularizationFlag
        RegularizationOptions
        DCF
        PreconditionMatrix
        SpatialBasis
        SpatialBasisHandle
        TemporalBasis
        OperatorScaling
        static_transfer_struct
        static_reg_struct
        Ncoefficients_spatial
        Ncoefficients_temporal
        param_struct
        SolutionVariables_init
        CalibrationsPerformed
    end
    
    
    methods
        
        %% Constructor
        function obj = MRMOTUS_Operator(ReferenceImage,KspaceTrajectory,param_struct)
            % Constructor.
            %
            % Inputs:
            %   ReferenceImage          - Complex reference image [N^d x 1]
            %   KspaceTrajectory        - Kspace trajectory for the forward model [#readout_samples x d x interleaves]
            %   param_struct            - All parameters specified in the example parameter files, e.g. 'Parameters_2Dt_RespMotion.m'
            
            disp('=== Initializing MR-MOTUS Operator object ===')

            ReferenceImage = ReferenceImage(:);
            
            % set internal parameters - don't touch
            obj.KspaceTrajectory        = KspaceTrajectory;
            obj                         = obj.set_reference_image(ReferenceImage(:));
            obj.NumberOfDynamics        = size(obj.KspaceTrajectory,3);
            obj.NumberOfSpatialDims     = size(obj.KspaceTrajectory,2);
            obj.ImDims                  = round(numel(ReferenceImage(:,1))^(1/obj.NumberOfSpatialDims));
            obj.ReferenceGrid           = MRMOTUS_Operator.MakeReferenceGrid(obj.ImDims,obj.NumberOfSpatialDims);
            obj.param_struct            = param_struct;
            
            
            disp(['+Detected ',num2str(obj.NumberOfSpatialDims),'D data']);
            
            % initialize all required components
            obj = obj.initialize_bases();
            obj = obj.initialize_solution_variables();
            obj = obj.initialize_regularization();
            obj = obj.initialize_parallel_computation();
            obj = obj.initialize_preconditioning();
            
            
            disp('=== Done! ===');
            
        end % end constructor
        
        
        
        %% Forward and gradient - multi dynamic & low-rank
        function [ObjfuncValue,Gradient,ForwardSignal,DvfOperator] = forward_and_gradient_lowrank(obj,MotionFieldCoefficients,KspaceData)
            % Main function to perform lowrank MR-MOTUS: returns objective funtion value and gradients
            %
            % Inputs
            %   MotionFieldCoefficients{1} = vectorization of Phi_matrix:         [ Component 1, Dim 1 ]
            %                                                                     [      ...      ]  
            %                                                                     [ Component 1, Dim d ]
            %                                                                     [      ...      ]
            %                                                                     [ Component R, Dim d ]
            %         
            %   MotionFieldCoefficients{2} = vectorization of (Psi_matrix)^T:     [ Component 1 ]
            %                                                                     [   ...       ]
            %                                                                     [ Component R ]
            %   or a vertical concatenation of the components above into a
            %   single vector
            %
            %   KspaceData                  - target kspace data to compute residuals [#readoutsamples x #interleaves]
            %
            % Outputs
            %   ObjfuncValue                - Objective function value 0.5*|| F(MotionFieldCoefficients) - KspaceData||_2^2
            %   Gradient                    - Gradient of the objective function w.r.t. MotionFieldCoefficients (see code for size specifications)
            %   ForwardSignal               - F(MotionFieldCoefficients)
            %   DvfOperator                 - Cell with MotionFieldOperator objects for each dynamic
            
            
            
            [Phi,Phi_rshp,Psi,PsiT,SpatialCoefficients_rshp2,MotionFieldCoefficients] = obj.ExpandMotionfieldCoefficients(MotionFieldCoefficients);
            
%             for i=1:size(Phi,2)
%                 Phi_scaling(1,i) = norm(Phi(:,i));
%             end
            Phi_scaling = ones(1,obj.param_struct.NumberOfComponents);
            
            SpatialCoefficients_rshp2 = SpatialCoefficients_rshp2   .* Phi_scaling;
%             Phi_rshp                  = Phi_rshp                    * Phi_scaling;
            Phi                       = Phi                         .* Phi_scaling;
            
            lowrankstart=tic;
            
            if isempty(obj.Ncoefficients_spatial)
                obj.Ncoefficients_spatial=size(obj.SpatialBasis,2);
            end
                
            
            %% Initialize...
            % - Gradient of objective w.r.t. spatial motion-fields, i.e. dE/d(Phi),
            % - Size: N^d x d*R
            GradientPhi                 = zeros([obj.ImDims^obj.NumberOfSpatialDims,obj.NumberOfSpatialDims*obj.param_struct.NumberOfComponents]);  
            % - Gradient w.r.t. temporal profile, i.e. dE/d(Psi^T)
            % - Size: R x M
            GradientPsiT                = zeros([obj.param_struct.NumberOfComponents,obj.NumberOfDynamics]); 
            ForwardSignal               = zeros(size(KspaceData));
            % - Gradient of objective function w.r.t. ALL cubic B-spline coefficients, i.e. [dE/d(alpha);dE/d(beta)]
            % - Size: (d*N_phi + N_psi)*R
            Gradient                    = zeros([obj.Ncoefficients_spatial+size(obj.TemporalBasis,2)*obj.param_struct.NumberOfComponents,1]);
            % - Gradient of objective function w.r.t. spatial cubic B-spline coefficients, i.e. dE/d(alpha)
            % - Size: d*N_phi*R x 1
            GradientPhiCoeff            = zeros(numel(MotionFieldCoefficients{1}),1);
            % - Set to zero to check if regularizations produce correct
            % results
            data_fid_lambda             = 1;
            % - Objective function values: total=data_fid + reg, data fidelity, regularization
            ObjfuncValue                = 0;
            ObjfuncValueDataFid         = zeros(obj.NumberOfDynamics,1);
            ObjfuncValueRegDynamic      = zeros(obj.NumberOfDynamics,1);
                
            
            
            
            % Make sure MATLAB doesn't broadcast the complete object, i.e.
            % don't make a copy for every worker:
            KspaceTrajectory        = obj.KspaceTrajectory;
            PreconditionMatrix      = obj.PreconditionMatrix;
            NumberOfDynamics        = obj.NumberOfDynamics;
            ReferenceGrid           = obj.ReferenceGrid;
            NumberOfComponents      = obj.param_struct.NumberOfComponents;
            NumberOfSpatialDims     = obj.NumberOfSpatialDims;
            RegularizationOptions   = obj.RegularizationOptions;
            referenceImageMask      = obj.referenceImageMask;
            RegularizationFlag      = obj.RegularizationFlag;
            static_transfer_struct  = obj.static_transfer_struct;
            static_reg_struct       = obj.static_reg_struct;
            fwd_and_gradient        = @ MRMOTUS_Operator.forward_and_gradient_singledynamic;
            fwd_grad_reg            = @ MRMOTUS_Operator.forward_and_gradient_regularization_singledynamic;
            ParallelComputationFlag = obj.param_struct.ParallelComputationFlag;
            
            % objective function scalings - don't touch
            scaling_reg             = 1e4/NumberOfDynamics;
            scaling_datafid = 0;
            for i=1:obj.NumberOfDynamics
                scaling_datafid        = scaling_datafid + norm(obj.PreconditionMatrix{i}*KspaceData(:,i))^2;
            end
            scaling_datafid = 1/sqrt(scaling_datafid)*1e4;
            
            
           
            


            
            %% Main computation block
%             rand_dyns = randsample(obj.NumberOfDynamics,floor(obj.NumberOfDynamics/2)).';
            % computations in serie
            if ~obj.param_struct.ParallelComputationFlag
                
                % loop over all dynamics
                for dynamic = 1:obj.NumberOfDynamics 
                    
                    % compute motion-field for this dynamic
                    MotionField_current                         = reshape(Phi*PsiT(:,dynamic),size(ReferenceGrid));
                    
                    % compute obj and gradient value for this dynamic
                    % w.r.t. complete motion-field
                    [ObjfuncValueDynamic,GradientDvfDynamic,ForwardSignal(:,dynamic),DvfOperator{dynamic}] = fwd_and_gradient(MotionField_current,squeeze(KspaceTrajectory(:,:,dynamic)),squeeze(KspaceData(:,dynamic)),PreconditionMatrix{dynamic},static_transfer_struct,ParallelComputationFlag);
                    ObjfuncValueDynamic = ObjfuncValueDynamic     * data_fid_lambda * scaling_datafid;
                    GradientDvfDynamic  = GradientDvfDynamic(:)   * data_fid_lambda * scaling_datafid;
                    
                    % seperate complete gradient to phi and psi components
                    GradientPsiT(:,dynamic)                     = Phi'*GradientDvfDynamic;
                    GradientPhi                                 = GradientPhi + reshape(kron(PsiT(:,dynamic),GradientDvfDynamic),[],NumberOfComponents*NumberOfSpatialDims);%reshape(GradientDVF_t,[],obj.NumberOfSpatialDims)),[],obj.NumberOfSpatialDims*NumberOfComponents),[],1);
                    
                    % update obj. function value
                    ObjfuncValueDataFid(dynamic) = ObjfuncValueDynamic;     
                    
                    
                    if obj.RegularizationFlag
                        
                        % Compute 'effective' spatial coefficients alpha * Psi(:,t)^T for this dynamics
                        TotalCoefficients_rshp = reshape(SpatialCoefficients_rshp2*PsiT(:,dynamic),[],NumberOfSpatialDims);
                        
                        % Compute the regularization function value R(D_t) and gradient dR/dD_t
                        [RegValueDynamic,RegGradientDvfDynamic] = fwd_grad_reg(TotalCoefficients_rshp, RegularizationOptions,static_reg_struct,ParallelComputationFlag); 
                        
                        % Update internal objective function value
                        ObjfuncValueRegDynamic(dynamic) = RegValueDynamic                                                           * scaling_reg;
                        
                        % Compute dR/d(alpha) and dR/d(Psi^T) from dR/D_t
                        GradientPhiCoeff = GradientPhiCoeff + kron(PsiT(:,dynamic),RegGradientDvfDynamic(:))                        * scaling_reg;
                        GradientPsiT(:,dynamic) = GradientPsiT(:,dynamic) + SpatialCoefficients_rshp2.'*RegGradientDvfDynamic(:)    * scaling_reg;
                    end
                    
                end           
                
            % computations in parallel
            else
%                 spe = speye(size(PreconditionMatrix{1},1),size(PreconditionMatrix{1},1));
                a=1;
                % loop over all dynamics
                parallel_reg = 1;
                parfor dynamic = 1:obj.NumberOfDynamics % some notes: matlab transfers all variables in the loop, regardless of whether they are accessed or within if statements
                    
                    
                    % compute motion-field for this dynamic
                    MotionField_current                         = reshape(Phi*PsiT(:,dynamic),size(ReferenceGrid));
                    
                    % compute obj and gradient value for this dynamic
                    % w.r.t. complete motion-field
                    [ObjfuncValueDynamic,GradientDvfDynamic,ForwardSignal(:,dynamic),DvfOperator{dynamic}] = fwd_and_gradient(MotionField_current,squeeze(KspaceTrajectory(:,:,dynamic)),KspaceData(:,dynamic),PreconditionMatrix{dynamic},static_transfer_struct,ParallelComputationFlag);
                    ObjfuncValueDynamic = ObjfuncValueDynamic       * data_fid_lambda * scaling_datafid;
                    GradientDvfDynamic  = GradientDvfDynamic(:)     * data_fid_lambda * scaling_datafid;
                    
                    % seperate complete gradient to phi and psi components
                    GradientPsiT(:,dynamic)                     = Phi'*GradientDvfDynamic;
                    GradientPhi                                 = GradientPhi + reshape(kron(PsiT(:,dynamic),GradientDvfDynamic),[],NumberOfComponents*NumberOfSpatialDims);%reshape(GradientDVF_t,[],obj.NumberOfSpatialDims)),[],obj.NumberOfSpatialDims*NumberOfComponents),[],1);
                    
                    % update obj. function value
                    ObjfuncValueDataFid(dynamic) = ObjfuncValueDynamic;                             
                    
                    if RegularizationFlag && parallel_reg
                        
                        % Compute 'effective' spatial coefficients alpha * Psi(:,t)^T for this dynamics
                        TotalCoefficients_rshp = reshape(SpatialCoefficients_rshp2*PsiT(:,dynamic),[],NumberOfSpatialDims);
                        % Compute the regularization function value R(D_t) and gradient dR/dD_t
                        [RegValueDynamic,RegGradientDvfDynamic] = fwd_grad_reg(TotalCoefficients_rshp, RegularizationOptions,static_reg_struct,ParallelComputationFlag); 
                        
                        % Update internal objective function value
                        ObjfuncValueRegDynamic(dynamic) = RegValueDynamic            * scaling_reg;

                        
                        % Compute dR/d(alpha) and dR/d(Psi^T) from dR/D_t
                        regradient_rshp = reshape(RegGradientDvfDynamic,[],1)        * scaling_reg;
                         
                        GradientPhiCoeff = GradientPhiCoeff + reshape(regradient_rshp*Psi(dynamic,:),[],1);
                        GradientPsiT(:,dynamic) = GradientPsiT(:,dynamic) + SpatialCoefficients_rshp2.'* regradient_rshp;
                    end

                    
              
                end
                
%                 % Note: regularization can also be done outside parfor loop, this is faster in some cases
%                 % to use this uncomment the if-block below, and comment the if-block above in the parfor loop.
%                 if RegularizationFlag 
% 
%                     % Compute 'effective' spatial coefficients alpha * Psi^T
%                     TotalCoefficients_rshp = reshape(SpatialCoefficients_rshp2*PsiT(:,:),[],NumberOfSpatialDims,NumberOfDynamics);
% 
%                     % Compute the regularization function value R(D_t) and gradient dR/dD_t
%                     [RegValueDynamic,RegGradientDvfDynamic] = fwd_grad_reg(TotalCoefficients_rshp, RegularizationOptions,static_reg_struct,ParallelComputationFlag); 
% 
%                     % Update internal objective function value
%                     ObjfuncValueRegDynamic = RegValueDynamic                                    * scaling_reg;
% 
%                     % Compute dR/d(alpha) and dR/d(Psi^T) from dR/D_t
%                     regradient_rshp = reshape(RegGradientDvfDynamic,[],NumberOfDynamics)        * scaling_reg;
% 
%                     GradientPhiCoeff = GradientPhiCoeff + reshape(regradient_rshp*Psi,[],1);
%                     GradientPsiT = GradientPsiT + SpatialCoefficients_rshp2.'* regradient_rshp;
% 
% 
%                 end

               
            end
            
                          
            %% Project gradients on the cubic B-spline coefficients, and collect everything in one vector
            
            % Phi
            Gradient(1:obj.Ncoefficients_spatial,1)        = Gradient(1:obj.Ncoefficients_spatial,1)+reshape(reshape(obj.SpatialBasis'*reshape(GradientPhi,[],obj.NumberOfSpatialDims*obj.param_struct.NumberOfComponents),[],obj.param_struct.NumberOfComponents).*Phi_scaling,[],1);
            
            % Phi coefficients directly
            Gradient(1:obj.Ncoefficients_spatial,1)        = Gradient(1:obj.Ncoefficients_spatial,1)+reshape(reshape(GradientPhiCoeff,[],obj.param_struct.NumberOfComponents).* Phi_scaling,[],1);
            
            % Psi
            Gradient(obj.Ncoefficients_spatial+1:end,1)    = Gradient(obj.Ncoefficients_spatial+1:end,1) + reshape(obj.TemporalBasis'*GradientPsiT',[],1);

            %% Compute objective function value
            ObjfuncValue = sum(ObjfuncValueDataFid) + sum(ObjfuncValueRegDynamic);
            
            %% Write some feedback to the terminal
            fprintf('f(x) = %14.6e, D(x) = %14.6e , R(x) = %14.6e, ||G_phi||_2 = %14.6e, ||G_psi||_2 = %14.6e , ||G||_2 = %14.6e \n', ...
            ObjfuncValue,...
            sum(ObjfuncValueDataFid),...
            sum(ObjfuncValueRegDynamic),...
            norm(Gradient(1:obj.Ncoefficients_spatial,1),2),...
            norm(Gradient(obj.Ncoefficients_spatial+1:end,1),2), ...
            norm(Gradient(:),2))  ;
            
            
            
            %% Visual feedback
            if obj.param_struct.VisualizationFlag
                for rr=1:obj.param_struct.NumberOfComponents
                    viss = reshape(Phi(:,rr),[ones(1,obj.NumberOfSpatialDims)*obj.ImDims,obj.NumberOfSpatialDims,1]);
                    vis(:,:,:,rr) = ((abs((sqrt(sum(viss.^2,obj.NumberOfSpatialDims+1:numel(size(viss))))))));
                end
                try
                    figure(91),
                    subplot(4,3,1);plot(abs(ForwardSignal(1:300,1)));xlabel('Sample index');ylabel('Magnitude [a.u.]');title('Forward model vs. Kspace data');hold on;
                    plot(abs(KspaceData(1:300)));hold off;legend('forward model','kspace data');
                    subplot(4,3,2),plot(PsiT.');title('Temporal profiles');legend;xlabel('Dynamic index');ylabel('Magnitude [a.u.]');

                    if obj.NumberOfSpatialDims==3
                        subplot(4,3,4);imagesc(rot90(squeeze(vis(:,end/2,:,1)),1));axis image;colormap gray; axis off;colorbar;title('Spatial component 1 - Sagittal')
                        subplot(4,3,5);imagesc(rot90(squeeze(vis(:,:,end/2,1)),0));axis image;colormap gray; axis off;colorbar;title('Spatial component 1 - Axial')
                        subplot(4,3,6);imagesc(rot90(squeeze(vis(end/2,:,:,1)),1));axis image;colormap gray; axis off;colorbar;title('Spatial component 1 - Coronal')
                        subplot(4,3,3);plot(MotionFieldCoefficients{1}(:));title('Spline coefficients');xlabel('Coefficient index');ylabel('Magnitude [a.u.]');
                        if obj.param_struct.NumberOfComponents > 1
                            
                            for i=1:min(2,(obj.param_struct.NumberOfComponents-1))
                                subplot(4,3,7+(i-1)*3);imagesc(rot90(squeeze(vis(:,end/2,:,i+1)),1));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1),' - Sagittal'])
                                subplot(4,3,8+(i-1)*3);imagesc(rot90(squeeze(vis(:,:,end/2,i+1)),0));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1),' - Axial'])
                                subplot(4,3,9+(i-1)*3);imagesc(rot90(squeeze(vis(end/2,:,:,i+1)),1));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1),' - Coronal'])
                            end
                        end
                    else
                        subplot(4,3,3);plot(MotionFieldCoefficients{1}(:));title('Spline coefficients');xlabel('Coefficient index');ylabel('Magnitude [a.u.]');

                        for i=0:min(5,obj.param_struct.NumberOfComponents-1)
                            subplot(4,3,4+i);imagesc(rot90(squeeze(vis(:,:,:,i+1)),0));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1)])
                        end

                    end

                    set_figure_fullscreen;
                    drawnow;

                
                catch
                    warning('Caught an error in the visualization');
                end
            end
            
        toc(lowrankstart)
        end % end forward_and_gradient_lowrank

        
        %% Forward and gradient - multi dynamic & low-rank
        function [ObjfuncValue,Gradient,ForwardSignal,DvfOperator] = forward_and_gradient_affinespline(obj,AffineSplineMotionFieldCoefficients,KspaceData)
            % Main function to perform lowrank MR-MOTUS: returns objective funtion value and gradients
            %
            % Inputs
            %   AffineSplineMotionFieldCoefficients: (d+1)*num_temp_splines x d
            %
            %   KspaceData                  - target kspace data to compute residuals [#readoutsamples x #interleaves]
            %
            % Outputs
            %   ObjfuncValue                - Objective function value 0.5*|| F(MotionFieldCoefficients) - KspaceData||_2^2
            %   Gradient                    - Gradient of the objective function w.r.t. MotionFieldCoefficients (see code for size specifications)
            %   ForwardSignal               - F(MotionFieldCoefficients)
            %   DvfOperator                 - Cell with MotionFieldOperator objects for each dynamic
            

            
            [Phi,PsiT,TimeDepAffineMatrix] = obj.ExpandMotionfieldCoefficientsAffineSpline(AffineSplineMotionFieldCoefficients);
            
            Psi = PsiT.';
           

            lowrankstart=tic;
            
            
                
            
            %% Initialize...
       
            ForwardSignal               = zeros(size(KspaceData));
            % - Gradient of objective function w.r.t. ALL cubic B-spline coefficients, i.e. [dE/d(alpha);dE/d(beta)]
            % - Size: (d*N_phi + N_psi)*R
            Gradient                    = zeros([numel(AffineSplineMotionFieldCoefficients),1]);
            % - Gradient of objective function w.r.t. spatial cubic B-spline coefficients, i.e. dE/d(alpha)
            % - Size: d*N_phi*R x 1
            % Set to zero to check if regularizations produce correct
            % results:
            data_fid_lambda             = 1;
            % - Objective function values: total=data_fid + reg, data fidelity, regularization
            ObjfuncValue                = 0;
            ObjfuncValueDataFid         = zeros(obj.NumberOfDynamics,1);
            ObjfuncValueRegDynamic      = zeros(obj.NumberOfDynamics,1);
                
            
            
            
            % Make sure MATLAB doesn't broadcast the complete object, i.e.
            % don't make a copy for every worker:
            KspaceTrajectory        = obj.KspaceTrajectory;
            PreconditionMatrix      = obj.PreconditionMatrix;
            NumberOfDynamics        = obj.NumberOfDynamics;
            ReferenceGrid           = obj.ReferenceGrid;
            NumberOfComponents      = 3;
            NumberOfSpatialDims     = obj.NumberOfSpatialDims;
            RegularizationOptions   = obj.RegularizationOptions;
            referenceImageMask      = obj.referenceImageMask;
            RegularizationFlag      = obj.RegularizationFlag;
            static_transfer_struct  = obj.static_transfer_struct;
            static_reg_struct       = obj.static_reg_struct;
            fwd_and_gradient        = @ MRMOTUS_Operator.forward_and_gradient_singledynamic;
            fwd_grad_reg            = @ MRMOTUS_Operator.forward_and_gradient_regularization_singledynamic;
            ParallelComputationFlag = obj.param_struct.ParallelComputationFlag;
            
            % objective function scalings - don't touch
            scaling_reg             = 1e4/NumberOfDynamics;
            scaling_datafid = 0;
            for i=1:obj.NumberOfDynamics
                scaling_datafid        = scaling_datafid + norm(obj.PreconditionMatrix{i}*KspaceData(:,i))^2;
            end
            scaling_datafid = 1/sqrt(scaling_datafid)*1e4;
            
            
           
            


            
            %% Main computation block
            
            % computations in serie
            if ~obj.param_struct.ParallelComputationFlag
             
                % loop over all dynamics
                for dynamic = 1:obj.NumberOfDynamics 


                   % compute motion-field for this dynamic
                    MotionField_current                         = reshape(Phi*PsiT(:,:,dynamic),size(ReferenceGrid)) - ReferenceGrid;
                    
                    % compute obj and gradient value for this dynamic
                    % w.r.t. complete motion-field
                    [ObjfuncValueDynamic,GradientDvfDynamic,ForwardSignal(:,dynamic),DvfOperator{dynamic}] = fwd_and_gradient(MotionField_current,squeeze(KspaceTrajectory(:,:,dynamic)),KspaceData(:,dynamic),PreconditionMatrix{dynamic},static_transfer_struct,ParallelComputationFlag);
                    ObjfuncValueDynamic = ObjfuncValueDynamic       * data_fid_lambda * scaling_datafid;
                    GradientDvfDynamic  = GradientDvfDynamic        * data_fid_lambda * scaling_datafid;
                    
                    % seperate complete gradient to phi and psi components
                    GradientPsiT(:,:,dynamic)                     = Phi'*GradientDvfDynamic;
                    
                    ObjfuncValueDataFid(dynamic) = ObjfuncValueDynamic;                             

                    
                    if ~isempty(obj.RegularizationFlag) && obj.RegularizationFlag
                        
                        % Compute 'effective' spatial coefficients alpha * Psi(:,t)^T for this dynamics
                        TotalCoefficients_rshp = reshape(SpatialCoefficients_rshp2*PsiT(:,dynamic),[],NumberOfSpatialDims);
                        
                        % Compute the regularization function value R(D_t) and gradient dR/dD_t
                        [RegValueDynamic,RegGradientDvfDynamic] = fwd_grad_reg(TotalCoefficients_rshp, RegularizationOptions,static_reg_struct,ParallelComputationFlag); 
                        
                        % Update internal objective function value
                        ObjfuncValueRegDynamic(dynamic) = RegValueDynamic                                                           * scaling_reg;
                        
                        % Compute dR/d(alpha) and dR/d(Psi^T) from dR/D_t
                        GradientPhiCoeff = GradientPhiCoeff + kron(PsiT(:,dynamic),RegGradientDvfDynamic(:))                        * scaling_reg;
                        GradientPsiT(:,dynamic) = GradientPsiT(:,dynamic) + SpatialCoefficients_rshp2.'*RegGradientDvfDynamic(:)    * scaling_reg;
                    end
                    
                end           
                
            % computations in parallel
            else
%                 spe = speye(size(PreconditionMatrix{1},1),size(PreconditionMatrix{1},1));
                a=1;
                % loop over all dynamics
                parallel_reg = 1;
                parfor dynamic = 1:obj.NumberOfDynamics % some notes: matlab transfers all variables in the loop, regardless of whether they are accessed or within if statements
                    
                    
                    % compute motion-field for this dynamic
                    MotionField_current                         = reshape(Phi*PsiT(:,:,dynamic),size(ReferenceGrid)) - ReferenceGrid;
                    
                    % compute obj and gradient value for this dynamic
                    % w.r.t. complete motion-field
                    [ObjfuncValueDynamic,GradientDvfDynamic,ForwardSignal(:,dynamic),DvfOperator{dynamic}] = fwd_and_gradient(MotionField_current,squeeze(KspaceTrajectory(:,:,dynamic)),KspaceData(:,dynamic),PreconditionMatrix{dynamic},static_transfer_struct,ParallelComputationFlag);
                    ObjfuncValueDynamic = ObjfuncValueDynamic       * data_fid_lambda * scaling_datafid;
                    GradientDvfDynamic  = GradientDvfDynamic     * data_fid_lambda * scaling_datafid;
                    
                    % seperate complete gradient to phi and psi components
                    GradientPsiT(:,dynamic)                     = Phi'*GradientDvfDynamic;
                    
                    % update obj. function value
                    ObjfuncValueDataFid(dynamic) = ObjfuncValueDynamic;                             
                    
                    if RegularizationFlag && parallel_reg
                        
                        % Compute 'effective' spatial coefficients alpha * Psi(:,t)^T for this dynamics
                        TotalCoefficients_rshp = reshape(SpatialCoefficients_rshp2*PsiT(:,dynamic),[],NumberOfSpatialDims);
                        % Compute the regularization function value R(D_t) and gradient dR/dD_t
                        [RegValueDynamic,RegGradientDvfDynamic] = fwd_grad_reg(TotalCoefficients_rshp, RegularizationOptions,static_reg_struct,ParallelComputationFlag); 
                        
                        % Update internal objective function value
                        ObjfuncValueRegDynamic(dynamic) = RegValueDynamic            * scaling_reg;

                        
                        % Compute dR/d(alpha) and dR/d(Psi^T) from dR/D_t
                        regradient_rshp = reshape(RegGradientDvfDynamic,[],1)        * scaling_reg;
                         
                        GradientPhiCoeff = GradientPhiCoeff + reshape(regradient_rshp*Psi(dynamic,:),[],1);
                        GradientPsiT(:,dynamic) = GradientPsiT(:,dynamic) + SpatialCoefficients_rshp2.'* regradient_rshp;
                    end

                    
              
                end
                
%                 % Note: regularization can also be done outside parfor loop, this is faster in some cases
%                 % to use this uncomment the if-block below, and comment the if-block above in the parfor loop.
%                 if RegularizationFlag 
% 
%                     % Compute 'effective' spatial coefficients alpha * Psi^T
%                     TotalCoefficients_rshp = reshape(SpatialCoefficients_rshp2*PsiT(:,:),[],NumberOfSpatialDims,NumberOfDynamics);
% 
%                     % Compute the regularization function value R(D_t) and gradient dR/dD_t
%                     [RegValueDynamic,RegGradientDvfDynamic] = fwd_grad_reg(TotalCoefficients_rshp, RegularizationOptions,static_reg_struct,ParallelComputationFlag); 
% 
%                     % Update internal objective function value
%                     ObjfuncValueRegDynamic = RegValueDynamic                                    * scaling_reg;
% 
%                     % Compute dR/d(alpha) and dR/d(Psi^T) from dR/D_t
%                     regradient_rshp = reshape(RegGradientDvfDynamic,[],NumberOfDynamics)        * scaling_reg;
% 
%                     GradientPhiCoeff = GradientPhiCoeff + reshape(regradient_rshp*Psi,[],1);
%                     GradientPsiT = GradientPsiT + SpatialCoefficients_rshp2.'* regradient_rshp;
% 
% 
%                 end

               
            end
            
                          
            %% Project gradients on the cubic B-spline coefficients, and collect everything in one vector
            
            % Psi
            Gradient = Gradient + reshape(permute(squeeze(mtimesx(obj.TemporalBasis,'c',permute(GradientPsiT,[3 1 2]))),[2 1 3]),[],1);
            

            %% Compute objective function value
            ObjfuncValue = sum(ObjfuncValueDataFid) + sum(ObjfuncValueRegDynamic);
            
            %% Write some feedback to the terminal
            fprintf('f(x) = %14.6e, D(x) = %14.6e , R(x) = %14.6e, ||G_phi||_infty = %14.6e, ||G_psi||_infty = %14.6e , ||G||_infty = %14.6e \n', ...
            ObjfuncValue,...
            sum(ObjfuncValueDataFid),...
            sum(ObjfuncValueRegDynamic),...
            norm(Gradient(1:obj.Ncoefficients_spatial,1),'inf'),...
            norm(Gradient(obj.Ncoefficients_spatial+1:end,1),'inf'), ...
            norm(Gradient(:),'inf'))  ;
            
            
            
            %% Visual feedback
%             if obj.param_struct.VisualizationFlag
%                 for rr=1:obj.param_struct.NumberOfComponents
%                     viss = reshape(Phi(:,rr),[ones(1,obj.NumberOfSpatialDims)*obj.ImDims,obj.NumberOfSpatialDims,1]);
%                     vis(:,:,:,rr) = ((abs((sqrt(sum(viss.^2,obj.NumberOfSpatialDims+1:numel(size(viss))))))));
%                 end
%                 try
%                     figure(91),
%                     subplot(4,3,1);plot(abs(ForwardSignal(1:300,1)));xlabel('Sample index');ylabel('Magnitude [a.u.]');title('Forward model vs. Kspace data');hold on;
%                     plot(abs(KspaceData(1:300)));hold off;legend('forward model','kspace data');
%                     subplot(4,3,2),plot(PsiT.');title('Temporal profiles');legend;xlabel('Dynamic index');ylabel('Magnitude [a.u.]');
% 
%                     if obj.NumberOfSpatialDims==3
%                         subplot(4,3,4);imagesc(rot90(squeeze(vis(:,end/2,:,1)),1));axis image;colormap gray; axis off;colorbar;title('Spatial component 1 - Sagittal')
%                         subplot(4,3,5);imagesc(rot90(squeeze(vis(:,:,end/2,1)),0));axis image;colormap gray; axis off;colorbar;title('Spatial component 1 - Axial')
%                         subplot(4,3,6);imagesc(rot90(squeeze(vis(end/2,:,:,1)),1));axis image;colormap gray; axis off;colorbar;title('Spatial component 1 - Coronal')
%                         subplot(4,3,3);plot(MotionFieldCoefficients{1}(:));title('Spline coefficients');xlabel('Coefficient index');ylabel('Magnitude [a.u.]');
%                         if obj.param_struct.NumberOfComponents > 1
%                             
%                             for i=1:min(2,(obj.param_struct.NumberOfComponents-1))
%                                 subplot(4,3,7+(i-1)*3);imagesc(rot90(squeeze(vis(:,end/2,:,i+1)),1));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1),' - Sagittal'])
%                                 subplot(4,3,8+(i-1)*3);imagesc(rot90(squeeze(vis(:,:,end/2,i+1)),0));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1),' - Axial'])
%                                 subplot(4,3,9+(i-1)*3);imagesc(rot90(squeeze(vis(end/2,:,:,i+1)),1));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1),' - Coronal'])
%                             end
%                         end
%                     else
%                         subplot(4,3,3);plot(MotionFieldCoefficients{1}(:));title('Spline coefficients');xlabel('Coefficient index');ylabel('Magnitude [a.u.]');
% 
%                         for i=0:min(5,obj.param_struct.NumberOfComponents-1)
%                             subplot(4,3,4+i);imagesc(rot90(squeeze(vis(:,:,:,i+1)),0));axis image;colormap gray; axis off;colorbar;title(['Spatial component ',num2str(i+1)])
%                         end
% 
%                     end
% 
%                     set_figure_fullscreen;
%                     drawnow;
% 
%                 
%                 catch
%                     warning('Caught an error in the visualization');
%                 end
%             end
            
        toc(lowrankstart)
        end % end forward_and_gradient_affinespline
        
        %% Expand motion field coefficients with affine (x) spline basis
        function [Phi, TimeDepAffineCoeff, TimeDepAffineMatrix]=ExpandMotionfieldCoefficientsAffineSpline(obj,AffineSplineMotionFieldCoefficients)
            % Inputs:
            %   - AffineSplineMotionFieldCoefficients: (d+1)*num_temp_splines * d x 1
            % Outputs:
            %   - Phi: spatial motion-field component
            %   - TimeDepAffineCoeff
            %       - 1 x d Cell, one element for each motion direction
            %       - per cell-element: time-dependent affine coefficients, one set of coefficients per dynamic, size: (d+1) x NumberOfDynamics
            %   - TimeDepAffineMatrix, size: NumDynamics x d x (d+1)
            
            
            AffineSplineMotionFieldCoefficients = reshape(AffineSplineMotionFieldCoefficients,[],obj.NumberOfSpatialDims);

            
            for i=1:size(AffineSplineMotionFieldCoefficients,2) % loop over motion direction x,y,z
                TimeDepAffineCoeff(:,i,:)   = reshape(AffineSplineMotionFieldCoefficients(:,i),[],size(obj.TemporalBasis,2)) * obj.TemporalBasis.';
                TimeDepAffineMatrix(i,:,:)  = permute(TimeDepAffineCoeff(:,i,:),[2 1 3]); 
            end
            
            Phi = obj.SpatialBasis;
            
        end % end ExpandMotionfieldCoefficientsAffineSpline
        
        %% Expand motion field coefficients
        function [Phi,Phi_rshp,Psi,PsiT,SpatialCoefficients_rshp2,MotionFieldCoefficients]=ExpandMotionfieldCoefficients(obj,MotionFieldCoefficients)
            % Function to expand motion field spline basis coefficients to a full motion-field.
            
            % Reshaping operator (in case a matrix is provided instaed of
            % cells)
            Matrix2Cell = @(x) {x(1:obj.Ncoefficients_spatial),x(obj.Ncoefficients_spatial+1:end)};
            
            if ~iscell(MotionFieldCoefficients)
                MotionFieldCoefficients = Matrix2Cell(MotionFieldCoefficients);
            end
            
            % Some reshaping for convenience later...
            SpatialCoefficients_rshp    = reshape(MotionFieldCoefficients{1},[],obj.NumberOfSpatialDims*obj.param_struct.NumberOfComponents);
            SpatialCoefficients_rshp2   = reshape(MotionFieldCoefficients{1},[],obj.param_struct.NumberOfComponents);
            
            % build motion-fields
            Phi_rshp                    = obj.SpatialBasis*SpatialCoefficients_rshp;
            Phi                         = reshape(Phi_rshp,[],obj.param_struct.NumberOfComponents);
            Psi                        = (obj.TemporalBasis*reshape(MotionFieldCoefficients{2},[],obj.param_struct.NumberOfComponents));
            PsiT                        = Psi.';
            
        end % end ExpandMotionfieldCoefficients
        
        %% initialization functions
        function obj=set_reference_image(obj,ReferenceImage)
           
            obj.ReferenceImage      = ReferenceImage;
            if size(obj.ReferenceImage,2)>1
                for j=1:size(obj.ReferenceImage,2)
                    obj.DiagReferenceImage{j}  = spdiags(obj.ReferenceImage(:,j),0,numel(obj.ReferenceImage(:,j)),numel(obj.ReferenceImage(:,j)));
                end
            else
                obj.DiagReferenceImage  = spdiags(obj.ReferenceImage,0,numel(obj.ReferenceImage),numel(obj.ReferenceImage));
            end
            obj.referenceImageMask  = abs(obj.ReferenceImage(:,1)*0+1)>0;
        end
        
        function obj=initialize_parallel_computation(obj)
            
            % create some structs to prevent MATLAB to make multiple copies
            % of the same variables in parfor loops
            % for forward model...
            obj.static_transfer_struct=[];
            obj.static_transfer_struct.ReferenceGrid        = obj.ReferenceGrid;
            obj.static_transfer_struct.referenceImageMask   = obj.referenceImageMask;
            obj.static_transfer_struct.DiagReferenceImage   = obj.DiagReferenceImage;
            obj.static_transfer_struct.ReferenceImage       = obj.ReferenceImage;
            obj.static_transfer_struct.NumberOfSpatialDims  = obj.NumberOfSpatialDims;

            % ... and for regularization
            obj.static_reg_struct=[];
            obj.static_reg_struct.ImDims                    = obj.ImDims;
            obj.static_reg_struct.NumberOfSpatialDims       = obj.NumberOfSpatialDims;
            obj.static_reg_struct.referenceImageMask        = obj.referenceImageMask;
            obj.static_reg_struct.referenceImage            = denorm(abs(obj.ReferenceImage(:,1)).*(obj.referenceImageMask(:)));%abs(obj.ReferenceImage)>3.8e-5 & obj.referenceImageMask;%0.12;% ;%0.0025;%ones(numel(obj.ReferenceImage),1);%.0025;;%
  
            no_threads=obj.param_struct.NumberOfThreads;
            parallel_reg = 1;
            pool=gcp('nocreate');
            if obj.param_struct.ParallelComputationFlag 
                
                disp('+Detected parallel computation mode...')
                disp('+Initializing parallel computations...');
               
                if isempty(pool)
                    pool=parpool(min(obj.NumberOfDynamics,no_threads));
                else
                    disp('(Re)starting parallel pool');
                    delete(pool.Cluster.Jobs);
                    delete(pool);
                    pool=parpool(min(obj.NumberOfDynamics,no_threads));
                end
                
                if isfield(obj.RegularizationOptions,'Types')  
                    
                    % determinant
                    if isfield(obj.RegularizationOptions.Types,'Determinant')
                    if obj.RegularizationOptions.Types.Determinant.Lambda>0 && parallel_reg
                    disp('+     ... Transfering Jacobian determinant regularization files to parallel workers')
                        det_fn = fieldnames(obj.RegularizationOptions.Types.Determinant);
                        for ifn=1:numel(det_fn)                             % if parallel_reg == 1 and Determinant regularization active, transfer all required data to the workers (saves memory later)
                            if ~strcmpi(strip(det_fn{ifn}),'Lambda')
                                obj.RegularizationOptions.Types.Determinant.(det_fn{ifn}) = parallel.pool.Constant(obj.RegularizationOptions.Types.Determinant.(det_fn{ifn}));   
                            end
                        end  
                    elseif obj.RegularizationOptions.Types.Determinant.Lambda>0 && ~parallel_reg
                        det_fn = fieldnames(obj.RegularizationOptions.Types.Determinant);
                        for ifn=1:numel(det_fn)                            
                            if ~strcmpi(strip(det_fn{ifn}),'Lambda')
                                obj.RegularizationOptions.Types.Determinant.(det_fn{ifn}).Value = obj.RegularizationOptions.Types.Determinant.(det_fn{ifn});
                            end
                        end  
                    end
                    end
                    
                    % total variation
                    if isfield(obj.RegularizationOptions.Types,'TV')
                    if obj.RegularizationOptions.Types.TV.Lambda>0
                        if parallel_reg
                                disp('+     ... Transfering total variation regularization files to parallel workers')

                            obj.RegularizationOptions.Types.TV.D = parallel.pool.Constant(obj.RegularizationOptions.Types.TV.D);   
                        else
                            buff = obj.RegularizationOptions.Types.TV.D;
                            obj.RegularizationOptions.Types.TV=rmfield(obj.RegularizationOptions.Types.TV,'D');
                            obj.RegularizationOptions.Types.TV.D.Value= buff;
                            clearvars buff
                        end


                    else
                        obj.RegularizationOptions.Types.TV.D = [];
                    end
                    end
                end
                
                obj.static_transfer_struct=parallel.pool.Constant(obj.static_transfer_struct);


            else
                disp('+Detected serial computation mode')

                if ~isempty(pool)
                    delete(pool);
                end

                if isfield(obj.RegularizationOptions,'Types')  
                    if obj.RegularizationOptions.Types.TV.Lambda>0
                        buff = obj.RegularizationOptions.Types.TV.D;
                        obj.RegularizationOptions.Types.TV.D = [];
                        obj.RegularizationOptions.Types.TV.D.Value = buff;
                        clearvars buff
                    end
                
                
                    if isfield(obj.RegularizationOptions.Types,'Determinant')
                        if obj.RegularizationOptions.Types.Determinant.Lambda>0
                            det_fn = fieldnames(obj.RegularizationOptions.Types.Determinant);
                            for ifn=1:numel(det_fn)                            
                                if ~strcmpi(strip(det_fn{ifn}),'Lambda')
                                    buf = obj.RegularizationOptions.Types.Determinant.(det_fn{ifn});
                                    eval(['obj.RegularizationOptions.Types.Determinant.',det_fn{ifn},'=struct()']);
                                    obj.RegularizationOptions.Types.Determinant.(det_fn{ifn}).Value = buf;
                                end
                            end  
                        end
                    end

                end

            end
        end
        
        function obj=initialize_bases(obj)
            
            if (isfield(obj.param_struct,'spatial_affine_basis') && ~obj.param_struct.spatial_affine_basis) || ~isfield(obj.param_struct,'spatial_affine_basis')
                obj = initialize_spatial_spline_bases(obj);
            else
                obj = initialize_spatial_affine_basis(obj);
            end
            
            obj = initialize_temporal_spline_basis(obj);
            
            obj.SpatialBasis = obj.SpatialBasis/norm(obj.SpatialBasis(:));
            obj.TemporalBasis = obj.TemporalBasis/norm(obj.TemporalBasis(:));
            
            
            
            
        end
        
        
        
                
        
        function obj=initialize_spatial_spline_bases(obj)
            
            disp('+Initializing spatial spline basis...');
            
            basis_options_spatial.N                 = obj.ImDims;
            basis_options_spatial.spline_orders     = obj.param_struct.NumberOfSpatialSplines;
            
            if isfield(obj.param_struct,'NumberOfSpatialSplines_Z')
                basis_options_spatial.spline_orders_z = obj.param_struct.NumberOfSpatialSplines_Z;
            end
            
            basis_options_spatial.dimension         = obj.NumberOfSpatialDims;
            basis_options_spatial.spatial_ordering  = [2 1 3 4];



            obj.SpatialBasisHandle                  = ConstructBasisHandle(basis_options_spatial);
            
            if obj.NumberOfSpatialDims==2
                SpatialBasis                            = obj.SpatialBasisHandle(0,0);
            else
                SpatialBasis                            = obj.SpatialBasisHandle(0,0,0);
            end
            SpatialBases                            = {SpatialBasis};  
            obj.SpatialBasis                        = cat(2,SpatialBases{:});
            %NumberOfSpatialBases                    = size(SpatialBases,2);

            
        end
        
        function obj = initialize_spatial_affine_basis(obj)
            
            disp('+Initializing spatial affine basis...');
            
            %n_AffineParameters = (obj.NumberOfSpatialDims + 1)*obj.NumberOfSpatialDims ;
            
            basis_for_1_motionfield_direction = [obj.ReferenceGrid,obj.ReferenceGrid(:,1)*0+(max(abs(obj.ReferenceGrid(:))))/100]; % = [x_coords(:) y_coords(:) (z_coords(:)) 1(:)], to be multiplied with temporal basis via kronecker product
                
            obj.SpatialBasis = basis_for_1_motionfield_direction;
            
        end
        
        function obj = initialize_temporal_spline_basis(obj)
            
%             
%             if ~exist('obj.param_struct.temporal_spline_basis')
%                 obj.param_struct.temporal_spline_basis = 1;
%             end
            

            if obj.param_struct.NumberOfTemporalSplines>0 && obj.param_struct.temporal_spline_basis
                disp('+Initializing temporal spline basis...');
                basis_options_temporal.N                = obj.NumberOfDynamics;
                basis_options_temporal.spline_orders    = obj.param_struct.NumberOfTemporalSplines;
                basis_options_temporal.dimension        = 1;
                basis_options_temporal.spatial_ordering = [1];

                TemporalBasisHandle                     = ConstructBasisHandle(basis_options_temporal);
                TemporalBasis                           = TemporalBasisHandle(0);
                TemporalBases                           = {TemporalBasis}; 
                obj.TemporalBasis                        = cat(2,TemporalBases{:});
                NumberOfTemporalBases                   = size(TemporalBases,2);

            else % don't initialize spline basis in time
                disp('+Initializing identity temporal basis...');
                obj.TemporalBasis                        = speye([obj.NumberOfDynamics obj.NumberOfDynamics]);            
                NumberOfTemporalBases                    = 1;

            end
            
        end
        
        
        function obj = initialize_solution_variables(obj)
            
            obj.Ncoefficients_spatial = size(obj.SpatialBasis,2)*obj.NumberOfSpatialDims*obj.param_struct.NumberOfComponents;
            obj.Ncoefficients_temporal = size(obj.TemporalBasis,2)*obj.param_struct.NumberOfComponents;

            
            if ~isfield(obj.param_struct,'solution_variables_initialization') || isempty(obj.param_struct.solution_variables_initialization)
                disp('+Initializing solution variables as random vectors...');

                PhiCoefficients_0 = reshape(rand([obj.Ncoefficients_spatial,1])-0.5,[],obj.param_struct.NumberOfComponents);
                PsiCoefficients_0 = reshape(rand([obj.Ncoefficients_temporal,1])-0.5,[],obj.param_struct.NumberOfComponents);

                PhiCoefficients_0 = normalize(PhiCoefficients_0,1,'norm',2)*10;
                PsiCoefficients_0 = normalize(PsiCoefficients_0,1,'norm',2);
                
                if ~isfield(obj.param_struct,'init_scaling')
                    obj.param_struct.init_scaling = 1;
                end

                obj.SolutionVariables_init = [PhiCoefficients_0(:);PsiCoefficients_0(:)]*obj.param_struct.init_scaling;
            else
                disp('+Initializing solution variables from param_struct...');

                obj.SolutionVariables_init = obj.param_struct.solution_variables_initialization;
                
                
            end
            
                
            

        end
        
        function obj=initialize_regularization(obj)
            
            disp('+Initializing regularization...');
            
            %% Intialize regularization
            obj.RegularizationOptions.GridSpacing      = [1 1 1]*0.35*1.25;
            obj.RegularizationOptions.Types.Determinant = [];
            obj.RegularizationOptions.Types.Determinant.Lambda   = obj.param_struct.lambda_det;
            if obj.RegularizationOptions.Types.Determinant.Lambda > 0 

                disp('+     ... Jacobian of determinant');
                obj.RegularizationFlag = 1;
                if obj.NumberOfSpatialDims==2
                    obj.RegularizationOptions.Types.Determinant.BasisDxY    = obj.SpatialBasisHandle(1,0);
                    obj.RegularizationOptions.Types.Determinant.BasisXDy    = obj.SpatialBasisHandle(0,1);
                    obj.RegularizationOptions.Types.Determinant.A_1         = [obj.RegularizationOptions.Types.Determinant.BasisDxY.',obj.RegularizationOptions.Types.Determinant.BasisXDy.']; 
%                     obj.RegularizationOptions.Types.Determinant.A_2         = [obj.RegularizationOptions.Types.Determinant.A_1(:,end/2+1:end),obj.RegularizationOptions.Types.Determinant.A_1(:,1:end/2)];
                elseif obj.NumberOfSpatialDims==3
                    obj.RegularizationOptions.Types.Determinant.BasisDxYZ   = obj.SpatialBasisHandle(1,0,0);
                    obj.RegularizationOptions.Types.Determinant.BasisXDyZ   = obj.SpatialBasisHandle(0,1,0);
                    obj.RegularizationOptions.Types.Determinant.BasisXYDz   = obj.SpatialBasisHandle(0,0,1);
                    obj.RegularizationOptions.Types.Determinant.A_1         = [obj.RegularizationOptions.Types.Determinant.BasisDxYZ.',obj.RegularizationOptions.Types.Determinant.BasisXDyZ.',obj.RegularizationOptions.Types.Determinant.BasisXYDz.'];
                end
            end

            obj.RegularizationOptions.Types.TV = [];
            obj.RegularizationOptions.Types.TV.Lambda = obj.param_struct.lambda_TV;

            if obj.RegularizationOptions.Types.TV.Lambda>0
                disp('+     ... Total variation');
                obj.RegularizationFlag = 1;
                if obj.NumberOfSpatialDims==2
                    obj.RegularizationOptions.Types.TV.Eps    = obj.param_struct.eps_TV;
                    D_TV = {(obj.SpatialBasisHandle(1,0)),(obj.SpatialBasisHandle(0,1))};
                    obj.RegularizationOptions.Types.TV.D = D_TV;
                elseif obj.NumberOfSpatialDims==3
                    obj.RegularizationOptions.Types.TV.Eps    = obj.param_struct.eps_TV;
                    D_TV = {(obj.SpatialBasisHandle(1,0,0)),(obj.SpatialBasisHandle(0,1,0)),(obj.SpatialBasisHandle(0,0,1))};
                    obj.RegularizationOptions.Types.TV.D = D_TV;
                end

            end
           
            
        end
        
        
        
        function obj=initialize_preconditioning(obj)
            
            
            
            if obj.param_struct.PreconditionParam==0
                disp('+Initializing identity preconditioning matrix...');
            else
                disp('+Initializing preconditioning matrix...');
            end
            

            obj = obj.initialize_dcf();

            % make precond matrix
            for i=1:obj.NumberOfDynamics
                obj.PreconditionMatrix{i}   = spdiags((reshape(obj.DCF(:,i),[],1)).^(obj.param_struct.PreconditionParam),0,numel(obj.DCF(:,i)),numel(obj.DCF(:,i)));
            end

        end
        
        function obj=initialize_dcf(obj)
            
            disp('+     ... Computing/loading densitiy compensation function');
            
            % get dcf
            if ~isfield(obj.param_struct,'dcf') || isempty(obj.param_struct.dcf)
                dcf_coords = obj.KspaceTrajectory;
                dcf_coords(:,3,:)=zeros(size(dcf_coords(:,1,:)));
                obj.DCF = reshape(iterative_dcf_estimation(demax(dcf_coords)*obj.ImDims(1),6),[],size(obj.KspaceTrajectory,3));
                obj.DCF = obj.DCF/max(abs(obj.DCF(:)));
            else
                obj.DCF = reshape(obj.param_struct.dcf,[],size(obj.KspaceTrajectory,3));
                obj.DCF = obj.DCF/max(abs(obj.DCF(:)));
                disp('+   ... Density compensation function loaded from parameter struct');
            end
            
        end
        

            

    end
    
    methods(Static)
        %% Forward and gradient - single dynamic
        function [ObjfuncValDynamic,GradientDynamic,ForwardSignal,DVF_Operator] = forward_and_gradient_singledynamic(MotionField,KspaceTrajectory,KspaceData,PrecondMatrix,static_transfer_struct,ParallelComputationFlag,varargin) % should be faster paralellized than separately forward and adjoint
            % Function that evaluates the forward model for a single dynamic, this will be called in the larger functions above
 
            
            % Fetch the static_transfer_struct for this worker
            if ParallelComputationFlag
                static_transfer_struct = static_transfer_struct.Value;
            end

            % Initialize the warping operator F(.)
            DVF_Operator                = MotionFieldOperator(KspaceTrajectory,static_transfer_struct.ReferenceGrid(static_transfer_struct.referenceImageMask(:),:),MotionField(static_transfer_struct.referenceImageMask(:),:)); % forward: maps motionfields to data, adjoint: maps data to motionfields
            
            % Compute the forward model signal with current motion-field estimate: F(D_t)
            if numel(varargin)>0
                ForwardSignal               = DVF_Operator*static_transfer_struct.ReferenceImage(static_transfer_struct.referenceImageMask(:),varargin{1});
            else
                ForwardSignal               = DVF_Operator*static_transfer_struct.ReferenceImage(static_transfer_struct.referenceImageMask(:));
            end
            
            % Compute objective function value for this dynamic: 0.5 * || P ( F(D_t) - s_t ) ||_2^2
            ObjfuncValDynamic           = 0.5*norm(PrecondMatrix*(ForwardSignal(:)-KspaceData(:)))^2 ;
            GradientDynamic = zeros(size(static_transfer_struct.ReferenceGrid));
            
            % Compute the gradients for this dynamic w.r.t. motion-field D_t
            if numel(varargin)>0
                GradientDynamic(static_transfer_struct.referenceImageMask(:),:) = reshape(real(conj(static_transfer_struct.DiagReferenceImage{dynamic}(static_transfer_struct.referenceImageMask(:),static_transfer_struct.referenceImageMask(:)))*(DVF_Operator.'*(PrecondMatrix'*PrecondMatrix*(ForwardSignal-KspaceData)))),[],static_transfer_struct.NumberOfSpatialDims);
            else
                GradientDynamic(static_transfer_struct.referenceImageMask(:),:) = reshape(real(conj(static_transfer_struct.DiagReferenceImage(static_transfer_struct.referenceImageMask(:),static_transfer_struct.referenceImageMask(:)))*(DVF_Operator.'*(PrecondMatrix'*PrecondMatrix*(ForwardSignal-KspaceData)))),[],static_transfer_struct.NumberOfSpatialDims);
            end
            
            
        end % end forward_and_gradient_singledynamic
        
        
        %% Regularization - Forward and gradient - single dynamic
        function [ObjfuncValDynamic,GradientDynamic] = forward_and_gradient_regularization_singledynamic(MotionField,RegularizationOptions,static_struct,ParallelComputationFlag)
        % Function that evaluates the regularization and all corresponding gradients
        % MotionFie
        
        reg_types               = fieldnames(RegularizationOptions.Types);
        RegularizationEvaluation    = 0;
        ObjfuncValDynamic       = 0;
        GradientDynamic         = zeros(size(MotionField));
        
        
        for i=1:numel(reg_types)  
        
            switch reg_types{i}
                
                case 'L2_general' % old, not tested recently
                    
                    if RegularizationOptions.Types.(reg_types{i}).Lambda>0
                        RegularizationEvaluation    = RegularizationEvaluation + RegularizationOptions.Types.(reg_types{i}).DiscretizationMatrix(referenceImageMask(:),referenceImageMask(:))*MotionField(referenceImageMask(:),:);
                        ObjfuncValDynamic           = ObjfuncValDynamic + 0.5 * RegularizationOptions.Types.(reg_types{i}).Lambda * norm(RegularizationEvaluation,'fro')^2;
                        GradientDynamic(static_struct.referenceImageMask(:),:) = GradientDynamic(static_struct.referenceImageMask(:),:)  +  RegularizationOptions.Types.(reg_types{i}).Lambda * RegularizationOptions.Types.(reg_types{i}).DiscretizationMatrixAdjoint(static_struct.referenceImageMask(:),static_struct.referenceImageMask(:)) * RegularizationEvaluation;

                    end
                    
                case 'TV'
                    
                    if RegularizationOptions.Types.(reg_types{i}).Lambda>0
                        size_motionfield_coeff = size(MotionField);
                        MotionField = reshape(MotionField,size(MotionField,1),[]);
                        
%                         [regularizationEnergy,regularizationGradient]   = vectorial_tv_operator(MotionField,RegularizationOptions.Types.(reg_types{i}).Eps,RegularizationOptions.Types.(reg_types{i}).D.Value);
                        [regularizationEnergy,regularizationGradient]   = tv_vector_operator_coefficients(MotionField,RegularizationOptions.Types.(reg_types{i}).Eps,static_struct.ImDims*ones(1,static_struct.NumberOfSpatialDims),RegularizationOptions.GridSpacing,RegularizationOptions.Types.(reg_types{i}).D.Value);
                        ObjfuncValDynamic                               = ObjfuncValDynamic   + RegularizationOptions.Types.(reg_types{i}).Lambda * regularizationEnergy;

                        GradientDynamic = GradientDynamic  +  RegularizationOptions.Types.(reg_types{i}).Lambda * reshape(regularizationGradient,size_motionfield_coeff);

                    end
                    
                
                    
                case 'L2_coefficients' % old, not tested recently
                    
                    if RegularizationOptions.Types.(reg_types{i}).Lambda>0
                        RegularizationEvaluation               = RegularizationEvaluation + RegularizationOptions.Types.(reg_types{i}).DiscretizationMatrix(:,:)*reshape(MotionField,size(RegularizationOptions.Types.(reg_types{i}).DiscretizationMatrix(:,:),2),[]);
                        ObjfuncValDynamic           = ObjfuncValDynamic + 0.5 * RegularizationOptions.Types.(reg_types{i}).Lambda * norm(RegularizationEvaluation(:))^2;
                        GradientDynamic = GradientDynamic  +  reshape(RegularizationOptions.Types.(reg_types{i}).Lambda * RegularizationOptions.Types.(reg_types{i}).DiscretizationMatrixAdjoint(:,:) * RegularizationEvaluation,size(MotionField));
                        
                    end
                    
                case 'Determinant'
                    if RegularizationOptions.Types.Determinant.Lambda>0
                        
                        if size(MotionField,2)==3 %3D
                            dxdx = permute(RegularizationOptions.Types.Determinant.BasisDxYZ.Value*squeeze(MotionField(:,1,:))+1,[1 3 2]);
                            dxdy = permute(RegularizationOptions.Types.Determinant.BasisXDyZ.Value*squeeze(MotionField(:,1,:)),[1 3 2]);
                            dxdz = permute(RegularizationOptions.Types.Determinant.BasisXYDz.Value*squeeze(MotionField(:,1,:)),[1 3 2]);
                            dydx = permute(RegularizationOptions.Types.Determinant.BasisDxYZ.Value*squeeze(MotionField(:,2,:)),[1 3 2]);
                            dydy = permute(RegularizationOptions.Types.Determinant.BasisXDyZ.Value*squeeze(MotionField(:,2,:))+1,[1 3 2]);
                            dydz = permute(RegularizationOptions.Types.Determinant.BasisXYDz.Value*squeeze(MotionField(:,2,:)),[1 3 2]);
                            dzdx = permute(RegularizationOptions.Types.Determinant.BasisDxYZ.Value*squeeze(MotionField(:,3,:)),[1 3 2]);
                            dzdy = permute(RegularizationOptions.Types.Determinant.BasisXDyZ*squeeze(MotionField(:,3,:)),[1 3 2]);
                            dzdz = permute(RegularizationOptions.Types.Determinant.BasisXYDz.Value*squeeze(MotionField(:,3,:))+1,[1 3 2]);
                            

                            dtminone =     bsxfun(@times, static_struct.referenceImage,...
                                                                    dxdx    .* (dydy.*dzdz-dydz.*dzdy) - ...
                                                                    dxdy    .* (dydx.*dzdz-dydz.*dzdx) + ...
                                                                    dxdz    .* (dydx.*dzdy-dydy.*dzdx) - 1);
                                         
                            v = size(dtminone,3)  ;                             

                            B = reshape(bsxfun(@times,bsxfun(@times,static_struct.referenceImage,dtminone) ...
                                                              ,cat(2, (dydy.*dzdz-dydz.*dzdy),-(dydx.*dzdz-dydz.*dzdx),(dydx.*dzdy-dydy.*dzdx)  ,...  % derivative to x motion-field | derivative to y motion-field | derivative to z motion-field
                                                                      -dxdy.*dzdz+dxdz.*dzdy ,  dxdx.*dzdz-dxdz.*dzdx ,  -dxdx.*dzdy+dxdy.*dzdx ,...
                                                                       dxdy.*dydz-dxdz.*dydy , -dxdx.*dydz+dxdz.*dydx , dxdx.*dydy - dxdy.*dydx  )) ...
                                      ,[],3*size(dtminone,3)); % 
                                                                   
                            
                            clearvars dxdx dxdy dxdz dydx dydy dydz dzdx dzdy dzdz 
                            
                            ObjfuncValDynamic = ObjfuncValDynamic + RegularizationOptions.Types.Determinant.Lambda*(0.5*norm(dtminone(:)).^2);
                            clearvars dtminone
                            GradientDynamic = GradientDynamic  + reshape(RegularizationOptions.Types.Determinant.Lambda * ([RegularizationOptions.Types.Determinant.BasisDxYZ.Value.',RegularizationOptions.Types.Determinant.BasisXDyZ.Value.',RegularizationOptions.Types.Determinant.BasisXYDz.Value.'] * B),[],3,v);

                            
                        elseif size(MotionField,2)==2 %2D
                            for d=1:size(MotionField,3)

                                dxdx = RegularizationOptions.Types.Determinant.BasisDxY.Value*MotionField(:,1,d);
                                dxdy = RegularizationOptions.Types.Determinant.BasisXDy.Value*MotionField(:,1,d);
                                
                                dydx = RegularizationOptions.Types.Determinant.BasisDxY.Value*MotionField(:,2,d);
                                dydy = RegularizationOptions.Types.Determinant.BasisXDy.Value*MotionField(:,2,d);
                                

                                Determinant = reshape_to_square((dxdx + 1) .* (dydy + 1) - dxdy .* dydx,2);

                                dtminone = Determinant(:)-1;
                                ObjfuncValDynamic = ObjfuncValDynamic + RegularizationOptions.Types.Determinant.Lambda*0.5*norm(static_struct.referenceImage.*dtminone).^2;

                                B_1 = reshape(bsxfun(@times,static_struct.referenceImage.*static_struct.referenceImage.*dtminone,[(dydy(:)+1),-dydx(:)]),[],1); % derivative to x coefficients: to dxdx & dxdy. Bases: DxY, XDy
                                B_2 = reshape(bsxfun(@times,static_struct.referenceImage.*static_struct.referenceImage.*dtminone,[-dxdy(:),(dxdx(:)+1)]),[],1); % derivative to y coefficients: to dydx & dydy. Bases: DxY, XDy

                                GradientDynamic(:,:,d) = GradientDynamic(:,:,d)  + RegularizationOptions.Types.Determinant.Lambda * RegularizationOptions.Types.Determinant.A_1.Value * [B_1 , B_2];
                                
                            
                            end
                        
                        end
                        
                    end
                    
                    
                    

                    

                    
            end
            
        end
        
        end % forward_and_gradient_regularization_singledynamic
        
        
        function RefGrid=MakeReferenceGrid(N,d,affine_flag)
            if nargin<3
                affine_flag = 1;
            end
            
            % return a reference grid with image dimension N and d spatial dimensions

            if d==3 % 3D
                [tmp_x,tmp_y,tmp_z]     = meshgrid(-N/2:N/2-1,-N/2:N/2-1,-N/2:N/2-1);
                RefGrid                 = [tmp_x(:),tmp_y(:),tmp_z(:)];
                clearvars tmp_x tmp_y tmp_z
            elseif d==2 % 2D
                [tmp_x,tmp_y]           = meshgrid(-N/2:N/2-1,-N/2:N/2-1);
                
                
                if affine_flag
                    tmp_x = tmp_x + 0.5*0;
                    tmp_y = tmp_y + 0.5*0;
                end
                
                RefGrid              = [tmp_x(:),tmp_y(:)];
                clearvars tmp_x tmp_y
                    
            end
        end
        
        function SnapshotData=SimulateData(ReferenceImage,DVF,pars)
            % Generate with the MRMOTUS forward model, given all input parameters:
            %   
            % Definitions:     
            %   N:  Image dimensions             [integer]
            %   d:  Number of spatial dimensions [2/3]
            %   R:  Number of ranks in the DVF
            %   T:  Number of dynamics

            % Inputs:
            %   ReferenceImage:         Complex-valued image, size N x N [x N] x 1 x [xT]
            %   DVF:                    Cell with 1 component (no low-rank), or two
            %                           elements (low-rank). 
            %                           - Low-rank mode
            %                               size(DVF{1})= N*d x R 
            %                               size(DVF{2})= T   x R
            %                           - Full DVF mode
            %                               size(DVF{1})= N x d x T
            %   pars:                   struct with all the required parameters.
            %   - SamplesPerReadout:    Number of samples per readout
            %   - ReadoutsPerDynamic:   Number of readouts per dynamic
            %   - KspaceCoords:         Kspace trajectory coordinates, size pars.SamplesPerReadout*pars.ReadoutsPerDynamic x d x T
            %
            %
            % Outputs:
            %   SnapshotData:           Complex-valued snapshot data, size
            %                           SamplesPerReadout*ReadoutsPerDynamic x T


            NumberOfSpatialDims = size(pars.KspaceCoords,2);
            NumberOfDynamics    = size(pars.KspaceCoords,3);
            N               = size(ReferenceImage,1);
            d               = NumberOfSpatialDims;
            ReferenceGrid   = MRMOTUS_Operator.MakeReferenceGrid(N,d);
           
            dyn_ref_flag = size(ReferenceImage,5)>1;
            

            SnapshotData    = zeros(size(pars.KspaceCoords,1),NumberOfDynamics);
            low_rank_flag   = numel(DVF)>1;
            
            if low_rank_flag 
                T = size(DVF{2},1);
            else
                T = NumberOfDynamics;%size(DVF{1},3);
            end
            
            if low_rank_flag
                    error('not implemented')
            end
            
            NC = size(ReferenceImage,4);
            for coils_iii=1:NC
            disp(['Simulating data for coil ',num2str(coils_iii),'/',num2str(NC)])
            for dynamic_iii=1:T
            disp(['Simulating data for dynamic ',num2str(dynamic_iii),'/',num2str(T)])
%                 for readout_iii=1:pars.ReadoutsPerDynamic
%                     indices = [1:pars.SamplesPerReadout]+(readout_iii-1)*pars.SamplesPerReadout;

%                     if low_rank_flag
%                         dvf = DVF{1}*DVF{2}(readout_iii+(dynamic_iii-1)*pars.ReadoutsPerDynamic,:).';
%                     else
%                         dvf = DVF{1}(:,dynamic_iii);
%                     end

                
                
%                 dvf = DVF{1}(:,:,dynamic_iii);

                if ~dyn_ref_flag
                    SnapshotData(:,dynamic_iii,1,coils_iii)=MotionFieldOperator(pars.KspaceCoords(:,:,dynamic_iii),ReferenceGrid,DVF{1}(:,:,dynamic_iii))*reshape(ReferenceImage(:,:,:,coils_iii,:),[],1);
                else
                    SnapshotData(:,dynamic_iii,1,coils_iii)=MotionFieldOperator(pars.KspaceCoords(:,:,dynamic_iii),ReferenceGrid,DVF{1}(:,:,dynamic_iii)*0)*reshape(ReferenceImage(:,:,:,coils_iii,dynamic_iii),[],1);
                end
                
            end
            end

        end
    end
    
    
    
    
end

            