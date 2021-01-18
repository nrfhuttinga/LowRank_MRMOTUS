function DataStruct=ExtractSurrogateSignal(DataStruct,pars)
    % Perform filtering and extract surrogate signal from k=0 by selecting the princpal 
    % component over the coil dimension with the highest magntidue in a specified
    % frequency range
    %
    % Inputs:
    %   DataStruct                          - structure with at least the following fields
    %       .RawKspaceData                  - [#readoutsamples #readouts 1 #coils]
    %       .Sequence.Repetition_time       - Repetition time to deterimine the sampling frequency of the input k-space data
    %   pars                                - parameter struct with the following fields
    %       .visualize                      - visualization flag (default 1)
    %       .k0_index                       - index of k0 on the readout (default N/2+1)
    %       .lowpass_filtering              - lowpass_filtering flag (default 1)
    %       .pc_freq_range                  - freq. range to select principal component with largest magnitude (default [.1 .5] Hz)
    %       .filter_parameters      
    %           .fcuts                      - Transition band for low pass filter [Hz]
    %           .Fs                         - Sampling frequency of RawKspaceData over second dimension
    %           .PassBandRipple             - Pass band ripple for low pass filter
    %           . StopbandAttenuation       - Stopband attenuation for low pass filter
    %
    % Outputs:
    %   DataStruct                          - structure where the following fields have been added
    %       .SelfNavigator
    %           .RawKspaceData              - all k0 profiles extracted from the input k-space [#readouts #coils]
    %           .LowpassFilterFlag          - flag that lists if filtering is performed or not
    %           .LowpassFilterDelay         - introduced low pass filter delay
    %           .SurrogateSignal            - final extracted surrogate signal
    %
    % Niek Huttinga, UMC Utrecht, 2020

    if nargin<2 
        pars=[];
    end
    
    if ~isfield(pars,'filter_parameters')
        pars.filter_parameters=[];
    end
    
    % set default parameters
    pars = set_default(pars,'visualize',1);
    pars = set_default(pars,'k0_index',round(size(DataStruct.RawKspaceData,1)/2+1));
    pars = set_default(pars,'lowpass_filtering',1);
    pars = set_default(pars,'pc_freq_range',[0.1 0.5]);
    
    pars.filter_parameters = set_default(pars.filter_parameters,'fcuts',[0.7 1.5] );                                    % transition band (start_attenuation -> full attenuation)
    pars.filter_parameters = set_default(pars.filter_parameters,'Fs', 1/(DataStruct.Sequence.Repetition_time));         % Sampling Frequency (Hz) = 1/TR
    pars.filter_parameters = set_default(pars.filter_parameters,'PassBandRipple', 0.1);
    pars.filter_parameters = set_default(pars.filter_parameters,'StopbandAttenuation', 100);
    
    % extract all k0 points over time
    DataStruct.SelfNavigator.RawKspaceData              = squeeze(abs(DataStruct.RawKspaceData(pars.k0_index,:,:,:)));
    kdata_rshped                                        = DataStruct.SelfNavigator.RawKspaceData.';

    % low pass filtering
    disp('+   Filtering k-space data');
    if pars.lowpass_filtering
        [kdata_rshped,filter_delay]                         = kaiserbessel_filter_static(kdata_rshped,pars.filter_parameters);
        kdata_rshped(:,end-filter_delay+1:end)              = [];
        kdata_rshped                                        = kdata_rshped.';
        DataStruct.SelfNavigator.LowpassFilterFlag          = 1;
        DataStruct.SelfNavigator.LowpassFilterDelay         = filter_delay;
    else
        DataStruct.SelfNavigator.LowpassFilterFlag          = 0;
    end
    

    pars_surrogate.fs               = pars.filter_parameters.Fs;
    pars_surrogate.pc_freq_range    = pars.pc_freq_range;
    pars_surrogate.visualize        = pars.visualize;
    
    % perform SVD and select the principal component with maximum content in the breathing
    % frequency range [0.1 0.5]
    DataStruct.SelfNavigator.SurrogateSignal = SelectPC(abs(kdata_rshped.'),pars_surrogate);

    

    




end
