function surrogate=SelectPC(data,pars)
    % Perform SVD on Casorati matrix 'data' and selects
    % the principal component with largest magnitude in the frequency range 
    % 'pars.pc_freq_range'.
    %
    % Inputs:
    %   data                - k-space data to perform SVD on [samples x time]
    %   pars                - parameter struct with the following fields
    %       .pc_freq_range  - freq. range to select the surrogate in (default: [0.1 0.5])
    %       .fs             - sampling frequency along the second dimension of 'data' (default: 1/4e-3, for TR=4ms)
    %       .visualize      - visualization flag [0/1] (default: 1)
    %
    % Outputs
    %   surrogate           - extract surrogate signal [time x 1]
    %
    % Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.
    
    if nargin<2
        pars=[];
    end

    pars=set_default(pars,  'pc_freq_range', [0.1 0.5]);
    pars=set_default(pars,  'fs',         1/(4.4*1e-3) );
    pars=set_default(pars,  'visualize', 1);
    
    disp(['+   Performing SVD on Casorati matrix of k-space data'])
    [~,~,v]=svd(abs(data),'econ');
    [spectrum,frequencies]=SingleSidedFFT(v,pars.fs,1);
    
%     [~,max_index]=max(sum(abs(spectrum(frequencies(:)<=pars.pc_freq_range(2) & frequencies(:)>= pars.pc_freq_range(1) , :)),1),[],2);
    [~,max_index]=max(max(abs(spectrum(frequencies(:)<=pars.pc_freq_range(2) & frequencies(:)>= pars.pc_freq_range(1) , 1:5)),[],1),[],2);
    if max_index~=2
        disp(['Maximum index: ',num2str(max_index)]);
        warning('Check if correct PC is selected!');
    end
    
    disp(['+   Selecting principal component #',num2str(max_index), ' in provided frequency range'])
    surrogate = v(:,max_index);
    
    
    if pars.visualize
        figure;
        subplot(1,2,1);
        plot(surrogate);title(['Selected principal component #',num2str(max_index)]);xlim([1,length(surrogate)])
        ylabel('Motion magnitude [a.u.]');
        xlabel('Time [a.u.]');
        set_figure_fullscreen;
        set_paper_plot_export_black;
        drawnow;
        pause(1);
        subplot(1,2,2);
        
        a=plot(frequencies,spectrum(:,max_index));
        title('Corresponding single-sided frequency spectrum');
        xlabel('Hz');ylabel('Fourier coefficient magnitude [a.u.]');
        yticks([]);
        xlim([0 2]);
        
        hold on; 
        b=plot([0.1 0.1],[0 max(spectrum(:,max_index),[],1)],'--r');set_paper_plot_export_black;
        
        hold on; 
        plot([0.5 0.5],[0 max(spectrum(:,max_index),[],1)],'--r');
        ylim([0 max(spectrum(:,max_index),[],1)]);
        
        SkipLegendEntry({a,b});
        
        legend('Frequency range');
        
        set_figure_fullscreen;
        set_paper_plot_export_black;
        
        
        
    end
    
    
    
    
    
    

end


    