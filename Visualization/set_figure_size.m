function set_figure_size(figure,varargin)

if nargin<2
    Position = [10 10 1200 600];
else
    Position = varargin{1};
end

set(figure,'Position',Position);


end


