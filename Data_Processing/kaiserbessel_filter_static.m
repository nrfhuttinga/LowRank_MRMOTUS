function [filtered_signal,delay]=kaiserbessel_filter_static(raw_data,varargin)
% Filter all rows in raw_data with a kaiser bessel filter with filter
% settings below, or with settings in varargin{1} struct
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

if nargin<2
    % Filter settings
    filter_parameters=[];
    filter_parameters.fcuts = [0.5 0.75];                               % [stopband passband passband stopband]
    filter_parameters.Fs = 1/(4.4e-3);                               % Sampling Frequency (Hz) = 1/(2*TR)
    filter_parameters.PassBandRipple = 0.5;
    filter_parameters.StopbandAttenuation=65;
else
    filter_parameters = varargin{1};
end

filter_parameters.Fn = filter_parameters.Fs/2;                      % Nyquist Frequency (Hz)

lpFilt = designfilt('lowpassfir', 'PassbandFrequency', filter_parameters.fcuts(1)/filter_parameters.Fn,...
           'StopbandFrequency', filter_parameters.fcuts(2)/filter_parameters.Fn, 'PassbandRipple', filter_parameters.PassBandRipple, ...
           'StopbandAttenuation', filter_parameters.StopbandAttenuation, 'DesignMethod', 'kaiserwin');
    

delay=ceil(mean(grpdelay(lpFilt,size(raw_data,2),filter_parameters.Fs)))
% 
% if abs(delay-round(delay))>0
%     error(['Non-integer delay: ',num2str(delay),'. Change filter parameters!']);
% end

filtered_signal = zeros(size(raw_data));
% 
for i=1:size(raw_data,1)
for j=1:size(raw_data,3)
for k=1:size(raw_data,4)

    filter_output = filter(lpFilt,squeeze(raw_data(i,:,j,k)));
    filter_output(1:delay)=[];
    filtered_signal(i,1:size(filter_output,2),j,k)=filter_output;
end
end
end

end