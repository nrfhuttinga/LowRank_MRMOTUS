function [ReferenceImage,Kdata]=CalibrateReferenceAndKdata(ReferenceImage,Kdata,k0_index_on_readout)
    % function to calibrate reference image and kdata such that the scaling is
    % consistent with the forward model, i.e. F(displacements=0) = ReferenceImage.
    % This is not always guaranteed due different to scalings in (NU)FFT operators
    % 
    % Inputs:
    %   ReferenceImage          - NxN(xN) Complex reference image
    %   Kdata                   - [#readoutsamples x #readouts] k-space data
    %   k0_index_on_readout     - index of k0 on the readouts (first dim of 'Kdata')
    %
    % Outputs:
    %   ReferenceImage          - Scaled input such that F(displacements=0) = ReferenceImage
    %   Kdata                   - Scaled input such that F(displacements=0) = ReferenceImage
    %
    % Niek Huttinga, UMC Utrecht, 2020
    
    ReferenceImage = ForceSquareShape(ReferenceImage);
    

    ref_grid = MRMOTUS_Operator.MakeReferenceGrid(size(ReferenceImage,1),numel(size(ReferenceImage)));

    % construct forward model operator for only the k0 index
    MFO = MotionFieldOperator(zeros(1,size(ref_grid,2)),ref_grid,ref_grid*0);
    
    % generate data on k0 with the MR-MOTUS fwd model
    magn_ref = abs(MFO*ReferenceImage(:));

    % determine scaling due to data
    magn_data = double(max(abs(Kdata(:))));
    magn_data_target = 1;
    scaling_data = 1 / magn_data * magn_data_target 


    % scaling due to inconsistency between motion field operator and kdata
    scaling_operators       = double(mean(abs(Kdata(k0_index_on_readout:size(Kdata,1):end))) / magn_ref )


    % resolve the scalings
    ReferenceImage  = ReferenceImage * scaling_data * scaling_operators;
    Kdata           = Kdata * scaling_data;


end