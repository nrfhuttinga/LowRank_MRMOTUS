function set_paper_plot_export(FontSize,LineWidth)

if nargin < 2
    LineWidth = [1.5; 2];
end

if nargin < 1
    FontSize = 30;
end

set_background_white;
set(gca,'FontSize',FontSize);
set(gca,'FontWeight','bold');
set(gca,'LineWidth',LineWidth(1))
set(findobj(gca, 'Type', 'Line', 'Linestyle', '-'), 'LineWidth', LineWidth(2));

end
