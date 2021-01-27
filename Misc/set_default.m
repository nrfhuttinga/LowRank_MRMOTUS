function struct=set_default(struct,field,value)
% Function to set default value 'value' for 'field' in 'struct'.
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

if ~isfield(struct,field)
    try
        disp(['+Setting default parameter for ',field,': ',num2str(value)]);
    catch
    end
    
    struct.(field)=value;
end

end
