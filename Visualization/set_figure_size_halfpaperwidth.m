function set_figure_size_halfpaperwidth(figure,height,width_scaling)

if nargin<1 || isempty(figure)
    figure=gcf;
end

set(figure,'Units','centimeters');


if nargin<2
    height=29.7/2;
else
    height=29.7/height;
end

if nargin<3
    width=21/2;
else
    width=21/width_scaling;
end


set(figure,'Position',[10 10 width height]);

end

