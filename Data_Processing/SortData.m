function [DataStruct,varargout] = SortData(DataStruct,param_struct)

    DataStruct.RawKspaceData                    = DataStruct.RawKspaceData(:,param_struct.BeginReadoutIdx:end,:,:);
    DataStruct.Coordinates                      = DataStruct.Coordinates(:,:,param_struct.BeginReadoutIdx:end);
    DataStruct.SelfNavigator.SurrogateSignal    = DataStruct.SelfNavigator.SurrogateSignal(param_struct.BeginReadoutIdx:end);
    
    
     % Check trajectory
%     disp('=== Plotting trajectory ===');
%     figure;PlotTrajectory(DataStruct.Coordinates(:,:,1:100))


    if isfield(param_struct,'binning_pars')
        binning_pars = param_struct.binning_pars;
    else
        binning_pars = [];
    end

    binning_pars = set_default(binning_pars,'surrogate_signal',DataStruct.SelfNavigator.SurrogateSignal);
    binning_pars = set_default(binning_pars,'binning_strategy','phase');
    binning_pars = set_default(binning_pars,'thresh',0.005);
    binning_pars = set_default(binning_pars,'resp_phases',param_struct.NumberOfComponents);
    binning_pars = set_default(binning_pars,'return_extreme_phases',0);

    %% sort the data in case of respiratory-resolved reconstruction

    disp('=== Sorting data ===');

    if param_struct.RespResolvedReconstruction

        % sort the data based on the respiratory phase
        [sorting_indices,phase]                       = RespiratoryBinning(binning_pars);


        sorting_indices = reshape((SelfnavToNormal(DataStruct.SelfNavigator.ReadoutIndices(sorting_indices),DataStruct.Sequence.Self_navigation_interval)),[],1);
        sorting_indices(sorting_indices>size(DataStruct.Coordinates,3))=[];

        DataStruct.RawKspaceData            = DataStruct.RawKspaceData(:,sorting_indices,:,:);
        DataStruct.Coordinates              = DataStruct.Coordinates(:,:,sorting_indices,:,:);

        if isfield(param_struct,'dcf')
            param_struct.dcf = param_struct.dcf(:,sorting_indices);
        end
        % force some other parameters for the resp resolved recon
        param_struct.ReadoutsPerDynamic         = floor( size(DataStruct.RawKspaceData,2) / param_struct.NumberOfDynamics );
        param_struct.NumberOfTemporalSplines    = round(4000/800); 
    end


    %% automatic parameters [don't touch]

    NumberOfSpatialDims             = size(DataStruct.Coordinates,1);                                                
    spatial_ordering                = [2 1 3];
    param_struct.IndicesOnReadout   = 1:size(DataStruct.RawKspaceData,1);

    RefImDims = size(DataStruct.ReferenceImage);
    svrs = structvars(param_struct);for i=1:size(svrs,1);eval(svrs(i,:));end

    param_struct.export_suffix = generate_export_suffix(NumberOfDynamics,ReadoutsPerDynamic,BeginReadoutIdx,RespResolvedReconstruction,lambda_det,lambda_TV,eps_TV,NumberOfComponents,NumberOfSpatialSplines,NumberOfTemporalSplines,RefImDims);

    
    %% Reshape the snapshot data and kspace coordinates according to specified parameters



    ReadoutIndices            = 1+[0:param_struct.ReadoutsPerDynamic*param_struct.NumberOfDynamics-1];

    if isfield(param_struct,'dcf')
        param_struct.dcf = param_struct.dcf(:,ReadoutIndices);
        varargout{1}=param_struct;
    end

    DataStruct.Coordinates    = reshape(DataStruct.Coordinates(spatial_ordering(1:NumberOfSpatialDims),param_struct.IndicesOnReadout,ReadoutIndices),NumberOfSpatialDims,numel(param_struct.IndicesOnReadout)*param_struct.ReadoutsPerDynamic,param_struct.NumberOfDynamics);
    DataStruct.RawKspaceData  = double(reshape(DataStruct.RawKspaceData(param_struct.IndicesOnReadout,ReadoutIndices),numel(param_struct.IndicesOnReadout)*param_struct.ReadoutsPerDynamic,param_struct.NumberOfDynamics));

end

