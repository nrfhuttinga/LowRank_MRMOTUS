function set_paper_plot_export(FontSize,LineWidth)

if nargin < 2
    LineWidth = [1.5; 2];
end

if nargin < 1
    FontSize = 30;
end
a=gca;
set_background_white;
set(gca,'FontSize',FontSize);
set(gca,'FontWeight','bold');
a.XAxis.FontWeight = 'normal';
a.YAxis.FontWeight = 'normal';

a.XAxis.FontSize = round(FontSize*0.8);
a.YAxis.FontSize = round(FontSize*0.8);


a.XLabel.FontWeight = 'bold';
a.YLabel.FontWeight = 'bold';

a.XLabel.FontSize = FontSize;
a.YLabel.FontSize = FontSize;


set(gca,'LineWidth',LineWidth(1))
set(findobj(gca, 'Type', 'Line', 'Linestyle', '-'), 'LineWidth', LineWidth(2));

if isfield(a,'Legend')
    a.Legend.FontSize = FontSize/1.15;
    a.Legend.FontWeight='normal';
end
if isfield(a,'Title')
    a.Title.FontSize = FontSize*1.2;
end

% set(gca,'Legend','FontSize',);
% set(gca,'Title','FontSize',);

end
