function set_figure_size_paperwidth(figure,height)

if nargin<1
    figure=gcf;
end


set(figure,'Units','centimeters');

if nargin<2
    height=29.7/3;
else
    height=29.7/height;
end
% ax_old=gca;
% set(gca,'Units','centimeters');
% ax=gca;
% scaling_width = 21/ax.Position(3);
set(figure,'Position',[10 10 21*1.3 height]);
% ax.Position(3) = ax.Position(3)*scaling_width;
end

