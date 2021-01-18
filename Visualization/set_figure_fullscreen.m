function set_figure_fullscreen()

% set(gcf,'units','normalized','outerposition',[0 0 1 1])
warning ('off','all');
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1);
pause(1);
warning ('on','all');

end


