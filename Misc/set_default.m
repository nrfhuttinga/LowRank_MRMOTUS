function struct=set_default(struct,field,value)
% Function to set default value 'value' for 'field' in 'struct'.
%
% Niek Huttinga, UMC Utrecht, 2020

if ~isfield(struct,field)
    try
        disp(['+Setting default parameter for ',field,': ',num2str(value)]);
    catch
    end
    
    struct.(field)=value;
end

end
