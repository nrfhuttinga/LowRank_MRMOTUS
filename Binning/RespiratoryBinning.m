function [index,phase] = RespiratoryBinning(par)
    % Function to sort 'par.surrogate_signal'   
    % 
    % Inputs:
    %   par                         - structure with all required data
    %       .binning_strategy       - strategy for the sorting: amplitude | phase | hybrid (phase, then amplitude) | amplitude binning in mid position,  (default: phase)
    %       .thresh                 - threshold for the peak detection (default: 0.005)
    %       .resp_phases            - number of respiratory phases to bin in (default: 5)
    %       .return_extreme_phase   - 0 - return all phases, 1 - return end-exhale, 2 - return end-inhale (default: 0)
    %
    % Outputs:
    %   index                       - indices in order to perform sorting, i.e. par.surrogate_signal(index) sorts the surrogate signal
    %   phase                       - calculated phase for each of the input samples
    %
    % Written by Bjorn Stemkens, UMC Utrecht, edited by Niek Huttinga, UMC
    % Copyright UMC Utrecht, 2020. For academic purpose only.
    
    par = set_default(par,'binning_strategy','phase');      % amplitude | phase | hybrid (phase, then amplitude) | amplitude binning in mid position
    par = set_default(par,'thresh',0.005);                  % sensitivity of the peak detection
    par = set_default(par,'resp_phases',5);                 % number of respiratory phases to bin in
    par = set_default(par,'return_extreme_phase',0);        % 0 - return all phases, 1 - return end-exhale, 2 - return end-inhale
    
    
    phase=[];


    if strcmpi('amplitude',par.binning_strategy)
        fprintf('+Performing amplitude binning \n');
        [~,index] = sort(par.surrogate_signal,'descend');
    elseif strcmpi('phase',par.binning_strategy)
        fprintf('+Performing phase binning \n');
        maxtab=[];
        mintab=[];
        thresh=par.thresh;
        while isempty(maxtab)
         [maxtab,mintab] = PeakDetection(par.surrogate_signal,thresh);
         thresh=thresh*.75;
        end

        figure;
        plot(par.surrogate_signal); hold on;
        plot(maxtab(:,1),maxtab(:,2),'.','MarkerSize',20);
        plot(mintab(:,1),mintab(:,2),'.','MarkerSize',20);
        hold off;
        title('Peak detection');
        peaks = zeros(length(par.surrogate_signal),1);
        peaks(maxtab(:,1)) = 1;
        peaks(mintab(:,1)) = -1;
        phase = CalculatePhase(peaks);
        [~,index] = sort(phase,'descend'); 
        
        
        SpokesPerBin = round(numel(par.surrogate_signal)/par.resp_phases);

        if par.return_extreme_phase>0 % only return extreme phases
            index =   [ index([1:round(SpokesPerBin/2),end-round(SpokesPerBin/2)+1:end]),...
                        index([round(end/2)+[-round(SpokesPerBin/2):round(SpokesPerBin/2)-1]]) ];
            index = index(:,par.return_extreme_phase);
        end




    elseif strcmpi('hybrid',par.binning_strategy)
        fprintf('+Performing hybrid binning \n');
        maxtab=[];
        mintab=[];
        thresh=par.thresh;
        while isempty(maxtab)
         [maxtab,mintab] = PeakDetection(par.surrogate_signal,thresh);
         thresh=thresh*.75;
        end

        figure;
        plot(par.surrogate_signal); hold on;
        plot(maxtab(:,1),maxtab(:,2),'.','MarkerSize',20);
        plot(mintab(:,1),mintab(:,2),'.','MarkerSize',20);
        hold off;
        title('Peak detection');
        peaks = zeros(length(par.surrogate_signal),1);
        peaks(maxtab(:,1)) = 1;
        peaks(mintab(:,1)) = -1;
        phase = CalculatePhase(peaks);

        inexhale = zeros(length(phase),1);
        for npa = 1:length(phase)
            if phase(npa) >= pi
                inexhale(npa) = 1;
            else
                inexhale(npa) = -1;
            end
        end

        respiration_inhale = par.surrogate_signal(inexhale == 1);
        nline_in = floor(length(respiration_inhale)/par.resp_phases);
        respiration_exhale = par.surrogate_signal(inexhale == -1);
        nline_ex = floor(length(respiration_exhale)/par.resp_phases);
        if nline_in > nline_ex
            nline_in = nline_ex;
        else
            nline_ex = nline_in;
        end

        [val_inhale,~] = unique(sort(respiration_inhale(1:nline_in*par.resp_phases),'descend'));
        [val_exhale,~] = unique(sort(respiration_exhale(1:nline_ex*par.resp_phases),'ascend'));

        real_in_index=[];
        real_ex_index=[];

        for t=1:numel(val_inhale)
            real_in_index = [real_in_index,find(par.surrogate_signal == val_inhale(t))];
        end

        for t=1:numel(val_exhale)
            real_ex_index = [real_ex_index,find(par.surrogate_signal == val_exhale(t))];
        end

        index = [fliplr(real_in_index),real_ex_index];

    elseif strcmpi('amplitude_midpos',par.binning_strategy)

        n_readouts = 0;
        par.resp_phases = par.resp_phases + 1; 
        while n_readouts < par.target_number_of_readouts && par.resp_phases > 5 

            clearvars  ci_95



            par.resp_phases = par.resp_phases - 1;
            [amplitudes,index] = sort(par.surrogate_signal,'descend');
            bin_oversampling = 100;


            n=1;
            mu_amplitudes = mean(amplitudes(:));
            sigma_amplitudes = std(amplitudes(:));

            z_score = 1.95;
            ci_95 = [mu_amplitudes-z_score*sigma_amplitudes/sqrt(n) , mu_amplitudes+z_score*sigma_amplitudes/sqrt(n)];
            motion_amplitude_range = (ci_95(2)-ci_95(1)); 

            ci_95 = ci_95 + [1 -1]*motion_amplitude_range*(2/5);

            bin_width = motion_amplitude_range/par.resp_phases;

            if bin_width > (ci_95(2)-ci_95(1) ) 
                break;
            else
                clearvars amplitude_bin_boundaries indices_in_bin bin_size
            end

            amplitude_bin_boundaries(:,1) = linspace(ci_95(1),ci_95(2)-bin_width,bin_oversampling*par.resp_phases); % bin starts
            amplitude_bin_boundaries(:,2) = min([amplitude_bin_boundaries(:,1)+bin_width,amplitude_bin_boundaries*0+ci_95(2)],[],2);                                % bin ends

            for i=1:size(amplitude_bin_boundaries,1)
                indices_in_bin{i}=amplitudes(:)>=amplitude_bin_boundaries(i,1) & amplitudes(:)<amplitude_bin_boundaries(i,2);
                bin_size(i)=sum(indices_in_bin{i});
            end


            [n_readouts,largest_bin]=max(bin_size);
        end
        
        amplitude_bin_boundaries = amplitude_bin_boundaries(largest_bin,:);
        index=index(amplitudes(:) >= amplitude_bin_boundaries(1) &  amplitudes(:) <= amplitude_bin_boundaries(2));

        figure,plot(par.surrogate_signal,'b');
        hold on;
        plot(1:length(par.surrogate_signal),repmat(amplitude_bin_boundaries(1),1,length(par.surrogate_signal)),'r');
        hold on;
        plot(1:length(par.surrogate_signal),repmat(amplitude_bin_boundaries(2),1,length(par.surrogate_signal)),'r');
        hold on;
        plot(1:length(par.surrogate_signal),repmat(ci_95(1),1,length(par.surrogate_signal)),'g');
        hold on;
        plot(1:length(par.surrogate_signal),repmat(ci_95(2),1,length(par.surrogate_signal)),'g');

        ylabel('Motion magnitude [a.u.]');
        xlabel('Time [a.u.]');
        legend('Surrogate signal','Selected bin begin','Selected bin end','Search window');
        set_paper_plot_export;
        set_figure_fullscreen;
        
    else
        error('>>Wrong input for par.binning_strategy \n');
    end