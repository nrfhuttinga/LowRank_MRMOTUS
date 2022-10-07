function [ReferenceImage, varargout] =  ReconstructRefImage(DataStruct,parameters)
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
    %       .pics_flag
    %       .vis_overgridding
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
    
    parameters      = set_default(parameters,'parallel_reconstruction',double(size(DataStruct.RawKspaceData,4)>1));
    parameters      = set_default(parameters,'readout_indices_ref',1:size(DataStruct.RawKspaceData,2));
    parameters      = set_default(parameters,'readout_downsampling',1);
    
    parameters      = set_default(parameters,'recon_overgridding',2);
    parameters      = set_default(parameters,'ImDims',make_even(116/parameters.readout_downsampling));
    parameters      = set_default(parameters,'centerout_flag',0);

    parameters      = set_default(parameters,'vis_overgridding',1.25);
    parameters      = set_default(parameters,'gpu_flag',1);
    
    parameters.bart = set_default(parameters.bart,'iterations',550);
    parameters.bart = set_default(parameters.bart,'version',6);
    parameters.bart = set_default(parameters.bart,'regularization_lambda',0.00033);
    parameters      = set_default(parameters,'pics_flag',1);

    
    NumberOfSpatialDims = size(DataStruct.Coordinates,1);

    
    
    indices_on_readouts = Crop1D(size(DataStruct.Coordinates,2),parameters.readout_downsampling,parameters.centerout_flag);
    
        
        
    
   
    
 
    %% extract reference image data

    % extract reference img kspace data from total set
    kdata_refimage          = DataStruct.RawKspaceData(indices_on_readouts,parameters.readout_indices_ref,:,:,:);

    % extract reference img coordinates from total set
    traj_refimage = double(DataStruct.Coordinates(:, indices_on_readouts,:));
    max_traj = max(abs(traj_refimage(:)));
    ImDimsOvergridding = make_even(max_traj * parameters.recon_overgridding);
    traj_refimage = traj_refimage(:, :, parameters.readout_indices_ref);
    
    
    %% reconstruction
    

    traj_refimage_bart = traj_refimage/max_traj*ImDimsOvergridding/2;

    if NumberOfSpatialDims==2
        traj_refimage_bart=cat(1,traj_refimage_bart,zeros(size(traj_refimage_bart(1,:,:))));
    end

    if ~parameters.pics_flag
        if ~isfield(DataStruct, 'DCF') || isempty(DataStruct.DCF)
            dcf = iterative_dcf_estimation(traj_refimage_bart);
        else
            dcf = DataStruct.DCF(indices_on_readouts,parameters.readout_indices_ref);
        end
    end
    
    
    if parameters.pics_flag && ~parameters.parallel_reconstruction
        
        
        for ch=1:size(kdata_refimage,4)
            disp(['Channel ',num2str(ch),'/',num2str(size(kdata_refimage,4))])
            if parameters.bart.version>4
%                 if parameters.pics_flag
                    if parameters.gpu_flag
                        xstar(:,:,:,ch) = bart(['pics -S -g -U -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage(:,:,:,ch),[3 1 2 4]), ones(ones(1,NumberOfSpatialDims)*ImDimsOvergridding));
                    else
                        xstar(:,:,:,ch) = bart(['pics -S -U -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage(:,:,:,ch),[3 1 2 4]), ones(ones(1,NumberOfSpatialDims)*ImDimsOvergridding));
                    end
                    
%                 else
%                     xstar(:,:,:,ch) = bart('nufft -i -t -g -l20 ', traj_refimage_bart, permute(bsxfun(@times,dcf*0+1,kdata_refimage(:,:,:,ch)),[3 1 2 4]));

%                 end
            else
                xstar(:,:,:,ch) = bart(['pics -S -g -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage(:,:,:,ch),[3 1 2 4]), ones(ones(1,NumberOfSpatialDims)*ImDimsOvergridding));
            end
        end
    elseif  parameters.pics_flag && parameters.parallel_reconstruction
        
        
        pars_csm.overgrid = parameters.recon_overgridding;
        pars_csm.nreadouts = size(DataStruct.RawKspaceData,2);
        if ~isfield(DataStruct.Coils,'Sensitivities') || isempty(DataStruct.Coils.Sensitivities)
             senseEst = EstimateCSM(DataStruct.RawKspaceData(indices_on_readouts,:,:,:,:) , DataStruct.Coordinates(:,indices_on_readouts,:), pars_csm);
        else
            senseEst = DataStruct.Coils.Sensitivities*0+1;
        end
        if parameters.bart.version>4
            xstar = bart(['pics -S -g -U -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage,[3 1 2 4]), senseEst);
        else
            xstar = bart(['pics -S -g -i ',num2str(parameters.bart.iterations),' -d5 -e -R W:7:0:',num2str(parameters.bart.regularization_lambda),' -t'], traj_refimage_bart, permute(kdata_refimage,[3 1 2 4]), senseEst);
        end
    elseif ~parameters.pics_flag
        
        for ch=1:size(kdata_refimage,4)
            disp(['Channel ',num2str(ch),'/',num2str(size(kdata_refimage,4))])
            xstar(:,:,:,ch) = bart(['nufft -a -t -d',num2str(ImDimsOvergridding),' '], traj_refimage_bart, permute(bsxfun(@times,dcf,kdata_refimage(:,:,:,ch)),[3 1 2 4]));
        end
    end
    

    for ch=1:size(xstar,4)
        ReferenceImage(:,:,:,ch) = reshape_to_square(reshape(crop_boundary(reshape_to_square(double(xstar(:,:,:,ch)),NumberOfSpatialDims),ones(1,NumberOfSpatialDims)*parameters.ImDims*parameters.vis_overgridding),[],1),NumberOfSpatialDims);
    end
    disp('+Done')
    
    varargout{1} = '';


end