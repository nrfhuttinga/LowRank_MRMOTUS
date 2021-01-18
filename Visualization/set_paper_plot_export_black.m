function set_paper_plot_export_black(a,b)

if nargin==0
    set_paper_plot_export();
elseif nargin==1
    set_paper_plot_export(a);
else 
    set_paper_plot_export(a,b);
end


set_background_black;
% set(gca,'FontColor','w')
set(gca,'Color',[0 0 0])
set(gca,'XColor','w')
set(gca,'YColor','w')

potential_title_objects = findall(gca, 'type', 'text', 'handlevisibility', 'off');
set(potential_title_objects,'Color','w');


a=findobj(gcf, 'Type', 'Legend');
if ~isempty(a)

    for i=1:numel(a)
        a(i).TextColor = 'white';
        a(i).Color = 'black';
        a(i).EdgeColor = [1 1 1]*0.4;
    end
end
drawnow;
% a.TextColor = 'black';
%     a.Color = 'white';
%     a.EdgeColor = [1 1 1]*0.4;
end
