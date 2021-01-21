function suffix=generate_export_suffix(varargin)
% Generate an export suffix that specifies all parameters
% 
% example: a=3; b=1; suffix = generate_export_suffix(a,b);
% --> suffix='_a=3_b=1'.
%
% Niek Huttinga, UMC Utrecht, 2020.

suffix=[];
for i=1:nargin
    
    variable_name = inputname(i);
    
    if ~ischar(varargin{i})

        val_str = num2str(varargin{i}(1));
        for j=2:numel(varargin{i})
            val_str = [val_str,',',num2str(varargin{i}(j))];
        end
    else
        val_str = varargin{i};
    end
    
    
    suffix=[suffix,'_',variable_name,'=',val_str];
    
end

    