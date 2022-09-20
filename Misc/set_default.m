function struct=set_default(struct,field,value)
% Function to set default value 'value' for 'field' in 'struct'.
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

if ~isfield(struct,field)
    try
        if ischar(value)
            disp(['+Setting default parameter for ',field,': ',value]);
        else
            if numel(value)>1
                disp(['+Setting default parameter for ',field,': [',num2str(value(1)),'...',num2str(value(end)),']']);
            else
                disp(['+Setting default parameter for ',field,': ',num2str(value)]);
            end
        end
    catch
    end
    
    struct.(field)=value;
end

end
