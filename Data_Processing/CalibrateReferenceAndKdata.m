function [ReferenceImage,Kdata,varargout]=CalibrateReferenceAndKdata(ReferenceImage,Kdata,Traj,pars)
    % Function to calibrate reference image and kdata such that the scaling is
    % consistent with the forward model, i.e. F(displacements=0) = ReferenceImage.
    % This is not always guaranteed due different to scalings in (NU)FFT operators
    % 
    % Inputs:
    %   ReferenceImage                      - NxN(xN) Complex reference image
    %   Kdata                               - [#readoutsamples x #readouts] k-space data
    %   Traj                                - [#spatialdims x #readoutsamples x #readouts] trajectory
    %   pars                                - struct with the following fields
    %       .calibration_readouts           - indices of readouts to use for the calibration
    %       .calibration_indices_on_readout - indices on readouts to use for the calibration
    %
    % Outputs:
    %   ReferenceImage                      - Scaled input such that F(displacements=0) = ReferenceImage
    %   Kdata                               - Scaled input such that F(displacements=0) = ReferenceImage
    %   varargout{1}                        - Scaling due to norm of data
    %   varargout{2}                        - Scaling due to fwd model <-> data inconsistency 
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
    
    if size(Traj,1)<=3 
        Traj = permute(Traj,[2 1 3]);
    end

    Traj = demax(Traj)/2;
    if nargin<4
        pars = [];
    end
    
%     scaling_handle = @(x) mean(abs(x),'all');
    scaling_handle = @(x) norm(x,'fro');
    
    pars = set_default(pars,'calibration_readouts',1:min(size(Kdata,2),100));
    
    traj_magn = sqrt(sum(Traj(:,:,pars.calibration_readouts(1)).^2,2));
    [~,ind]=min(traj_magn,[],1);
    pars = set_default(pars,'calibration_indices_on_readout',ind+[-9:9]);
    
    ordering = [2 1 3];
    calibration_coordinates(:,:,:) = Traj(pars.calibration_indices_on_readout,:,pars.calibration_readouts);
    calibration_coordinates(:,:,:) = calibration_coordinates(:,ordering(1:size(Traj,2)),:);
    ReferenceImage = ForceSquareShape(ReferenceImage);
    

    ref_grid = MRMOTUS_Operator.MakeReferenceGrid(size(ReferenceImage,1),numel(size(ReferenceImage)));

    % generate data on the calibration coordinates with 0 motion
    MFO = MotionFieldOperator(calibration_coordinates,ref_grid,repmat(ref_grid*0,[1 1 numel(pars.calibration_readouts)]));
    ref_data = reshape(MFO*ReferenceImage(:),numel(pars.calibration_indices_on_readout),numel(pars.calibration_readouts));
%     ref_data(:,1)=[];
%     pars.calibration_readouts(1) = [];
    magn_ref = scaling_handle(ref_data);
    
    % scaling due to data
    magn_data =  scaling_handle(Kdata(pars.calibration_indices_on_readout,pars.calibration_readouts));
    magn_data_target = 1;
    scaling_data = 1 / magn_data * magn_data_target ;


    % scaling due to inconsistency between motion field operator and kdata
    scaling_operators       = magn_data / magn_ref  ;



    ReferenceImage  = ReferenceImage * scaling_data * scaling_operators;
    Kdata           = Kdata * scaling_data;
    
%     figure;
%     PlotTrajectory(calibration_coordinates);
    
    figure;
    plot(matrix_to_vec(abs(Kdata(pars.calibration_indices_on_readout,pars.calibration_readouts(1:min(numel(pars.calibration_readouts),20)))))); hold on;
    plot(abs(matrix_to_vec(ref_data(:,1:min(numel(pars.calibration_readouts),20))))* scaling_data * scaling_operators );    
    varargout{1} = scaling_data;
    varargout{2} = scaling_operators;


end



% 
%     ReferenceImage = ForceSquareShape(ReferenceImage);
%     
% 
%     ref_grid = MRMOTUS_Operator.MakeReferenceGrid(size(ReferenceImage,1),numel(size(ReferenceImage)));
% 
%     % construct forward model operator for only the k0 index
%     MFO = MotionFieldOperator(zeros(1,size(ref_grid,2)),ref_grid,ref_grid*0);
%     
%     % generate data on k0 with the MR-MOTUS fwd model
%     magn_ref = abs(MFO*ReferenceImage(:));
% 
%     % determine scaling due to data
%     magn_data = double(max(abs(Kdata(:))));
%     magn_data_target = 1;
%     scaling_data = 1 / magn_data * magn_data_target ;
% 
% 
%     % scaling due to inconsistency between motion field operator and kdata
%     scaling_operators       = double(mean(abs(Kdata(k0_index_on_readout:size(Kdata,1):end))) / magn_ref );
% 
% 
%     % resolve the scalings
%     ReferenceImage  = ReferenceImage * scaling_data * scaling_operators;
%     Kdata           = Kdata * scaling_data;