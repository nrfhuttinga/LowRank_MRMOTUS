function ReferenceImage =  ReconstructRefImage(DataStruct,parameters)
    % Function to reconstruct the reference image.
    %
    %
    % Inputs:
    %   DataStruct                              - See main preprocessing codes
    %   parameters                              - Struct with the following fields:
    %       .parallel_reconstruction            - Flag for parallel imaging recon of reference image
    %       .readout_indices_ref                - Which readouts to extract from the complete dataset for ref. image recon 
    %       .readout_downsampling               - Factor to downsample the spatial resolution of the reference image
    %       .recon_overgridding                 - Overgridding in cs recon
    %       .ImDims                             - Image dimensions for final ref image
    %       .bart.iterations                    - Number of iterations               
    %       .bart.regularization_lambda         - L1-Wavelet regularization parameter
    %       .bart.version                       - BART version [required because there is a slight change in commands over the versions]
    %
    %
    % Outputs:
    %   ReferenceImage                          - ImDims x ImDims [x ImDims] complex reference image
    %
    %
    % Niek Huttinga, UMC Utrecht, 2020.

    
    if nargin < 2
        parameters = [];
        
    end

    if ~isfield(parameters,'bart')
        parameters.bart = [];
    end
    
    parameters      = set_default(parameters,'parallel_reconstruction',0);
    parameters      = set_default(parameters,'readout_indices_ref',1:size(DataStruct.RawKspaceData,2));
    parameters      = set_default(parameters,'readout_downsampling',1);
    
    parameters      = set_default(parameters,'recon_overgridding',2);
    parameters      = set_default(parameters,'ImDims',80);
    parameters      = set_default(parameters,'center_out_traj',0);

    
    parameters.bart = set_default(parameters.bart,'iterations',550);
    parameters.bart = set_default(parameters.bart,'version',6);
    parameters.bart = set_default(parameters.bart,'regularization_lambda',0.00033);
    
    NumberOfSpatialDims = size(DataStruct.Coordinates,1);

    
    if parameters.center_out_traj
        indices_on_readouts = 1:numel(Crop1D(size(DataStruct.Coordinates,2),parameters.readout_downsampling));
    else
        indices_on_readouts = Crop1D(size(DataStruct.Coordinates,2),parameters.readout_downsampling);
    end
        
        
    ImDimsOvergridding = make_even(parameters.ImDims * parameters.recon_overgridding);
   
 
    %% extract reference image data

    % extract reference img kspace data from total set
    kdata_refimage          = DataStruct.RawKspaceData(indices_on_readouts,parameters.readout_indices_ref,:,:);

    % extract reference img coordinates from total set
    traj_refimage        = double(demax(DataStruct.Coordinates(:,indices_on_readouts,parameters.readout_indices_ref))/2);
    
    
    %% reconstruction
    

    traj_refimage_bart = demax(traj_refimage)*ImDimsOvergridding/2;

    if NumberOfSpatialDims==2
        traj_refimage_bart=cat(1,traj_refimage_bart,zeros(size(traj_refimage_bart(1,:,:))));
    end


    if ~parameters.parallel_reconstruction
        if parameters.bart.version>4
            xstar = bart(['pics -S -g -U -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage,[3 1 2 4]), ones(ones(1,NumberOfSpatialDims)*ImDimsOvergridding));
        else
            xstar = bart(['pics -S -g -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage,[3 1 2 4]), ones(ones(1,NumberOfSpatialDims)*ImDimsOvergridding));
        end
    else
        pars_csm.overgrid = parameters.recon_overgridding;
        senseEst = EstimateCSM(DataStruct.RawKspaceData , DataStruct.Coordinates, pars_csm);
        if parameters.bart.version>4
            xstar = bart(['pics -S -g -U -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage,[3 1 2 4]), senseEst);
        else
            xstar = bart(['pics -S -g -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage,[3 1 2 4]), senseEst);
        end
    end


    ReferenceImage = reshape_to_square(reshape(crop_boundary(reshape_to_square(double(xstar),NumberOfSpatialDims),ones(1,NumberOfSpatialDims)*parameters.ImDims),[],1),NumberOfSpatialDims);
    disp('+Done')
    slicer5d(ReferenceImage)


end