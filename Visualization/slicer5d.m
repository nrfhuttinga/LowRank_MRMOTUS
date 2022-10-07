function f1=slicer5d(data,aspect,varargin)
%% Function for displaying 3D/4D data sets
%
% INPUT:    - data (>=2D, with channels as 4th dimension and time as 5th dimension)
%           - aspect ratio in vector format (eg: [1 1 3] if z is 3 times 
%           larger than x and y. (didn't check 5D support for aspect ratio)
%
% OUTPUT:   Currently no output is defined.
%
% Functionality
%   -   Clicking activates current viewport/subplot
%   -   3rd dimension can be skipped through by scrolling the mouse wheel
%   -   4th and 5th dimensions can be skipped through with vertical 
%       left-mouse/scroll-button dragging
%   -   Right clicking in the image creates a crosshair. The value at this
%       point is reported above the 4th viewport/subplot. In case of >3D data
%       this subplot will also show the signal evolution in the 5th dimension
%       (eg, the time series) of the point at the crosshair.
%   -   The currently displayed data is indicated by an asterisk in the plot 
%       of the timeseries. (>3D data only)
%   -   Arrow buttons scale the image (all viewports simultaneously)


% v1.0 July 2014
% Tim Schakel, Frank Simonis, Bjorn Stemkens.
%
% v1.1 April 2015
% fixed for Matlab2014b & higher
%
% v1.2 September 2018
% Niek Huttinga: changed input to 5D support, fourth dimension = channels,
% fifth dimension = time. Left mouse drag scrolls through time, scroll
% button drag scrolls through channels.


    
%% Initialize
    if ~isreal(data)
%         fprintf('The program will automatically use absolute data.\n');
        data = abs(data);
    end

    egg=0;
    if nargin < 2
        aspect=[1,1,1];
        egg=0;
    elseif ischar(aspect) && strcmp(aspect,[char(80),char(86)])
        aspect=[1,1,1];
        egg=1;
    elseif length(aspect)<3
        aspect=[1,1,1];
        egg=0;
    end

    dim=size(data);
    ndim=length(size(data));

    %Check if datasizes are valid
    if ndim>5
        error('Error: cannot handle this type of data. Only 2d-5d data is allowed');
    end
    
    if dim(end)==1 || ndim == 1
        data = reshape_to_square(data);
    end
    
    
    k=round(size(data)/2); % Default display are the center slices
    
    if ndim < 5
        k = [k,ones(1,5-ndim)];
        dim = [dim,ones(1,5-ndim)];
    end
    
    if numel(varargin)>=2
        range = varargin{2};
    else
        range=[min(data(:)),max(data(:))];
    end
    step=abs(max(data(:))-min(data(:)))/100;
    lastPoint = [];
    whichView = 1;
    crossLoc=[k(1),k(2),k(3)];
    p = [k(2),k(1),k(3)
         k(2),k(1),k(3)];
    offsliceCoord = [3 1 2];
    titleLines = {'Transversal view: ', 'Coronal view: ', 'Sagittal view: '};

%% Make figure
    if numel(varargin)==0 || isempty(varargin{1})
        f1=figure('Position',[50 50 800 800],'Name',inputname(1),'Color',[1 1 1],...
            'WindowScrollWheelFcn', @scroller, ....
            'WindowKeyPressFcn',@scale, ...
            'WindowButtonDownFcn', @click);
    else
        f1 = varargin{1};
    end
    
    
    colormap('gray')
    if ndim==3
        k(4)=1;
        handle{1} = subplot(221);
        imHandles{1} = imagesc(squeeze(data(:,:,k(3))),range);
        set(gca,'DataAspectRatio',[aspect(1) aspect(2) 1])
        axis off
        titleHandles{1} = title(num2str(k(3)),'FontSize',12);

        handle{2} = subplot(222);
        imHandles{2} = imagesc(squeeze(data(k(1),:,:))',range);
        set(gca,'DataAspectRatio',[aspect(3) aspect(1) 1],'YDir','normal')
        axis off
        titleHandles{2} = title(num2str(k(1)),'FontSize',12);

        handle{3} = subplot(223);
        imHandles{3} = imagesc(squeeze(data(:,k(2),:))',range);
        set(gca,'DataAspectRatio',[aspect(3) aspect(2) 1],'YDir','normal')
        axis off
        titleHandles{3} = title(num2str(k(2)),'FontSize',12);

        for i = 1:3
            set(titleHandles{i}, 'String', [sprintf('%s%d', titleLines{i}, k(offsliceCoord(i))),' of ', num2str(dim(offsliceCoord(i)))]);
        end
        
        %Plot a fancy logo in the otherwise empty subplot
        handle{4} = subplot(224);
        patch([0,0,50,50],[0,50,50,0],'white');
        patch([0,50,65,15],[50,50,65,65],'white')
        patch([50,65,65,50],[0,15,65,50],'white')
        for m=[6,9,11,13,13.5,14,14.5,15];
            line([m,50+m],[50+m,50+m],'Color','k');
            line([50+m,50+m],[m,50+m],'Color','k');
        end; 
        axis image off
        if egg
            text(4,25,char([80,114,111,106,101,99,116,32,86]),'FontSize',32,'FontName','Arial')
        else
        text(4,25,'Slicer','FontSize',52,'FontName','Arial')
        end
        value = data(k(1),k(2),k(3)); if value>10^5;value=num2str(value,'%10.4e\n');else value=num2str(value);end
        titleHandles{4} = title({['Crosshair at point [',num2str(k(1)),' ',num2str(k(2)),' ',num2str(k(3)),']. Value: ',value],''}, 'FontSize', 12);
    else
        handle{1} = subplot(221);
        imHandles{1} = imagesc(squeeze(data(:,:,k(3),k(4),k(5))),range);
        set(gca,'DataAspectRatio',[aspect(1) aspect(2) 1])
        axis off
        titleHandles{1} = title(['Slice ',num2str(k(3))],'FontSize',12);

        handle{2} = subplot(222);
        imHandles{2} = imagesc(squeeze(data(k(1),:,:,k(4),k(5)))',range);
        set(gca,'DataAspectRatio',[aspect(3) aspect(1) 1],'YDir','normal')
        axis off
        titleHandles{2} = title(['Slice ',num2str(k(1))],'FontSize',12);

        handle{3} = subplot(223);
        imHandles{3} = imagesc(squeeze(data(:,k(2),:,k(4),k(5)))',range);
        set(gca,'DataAspectRatio',[aspect(3) aspect(2) 1],'YDir','normal')
        axis off
        titleHandles{3} = title(['Slice ',num2str(k(2))],'FontSize',12);

        handle{4} = subplot(224);
        imHandles{4} = plot((1:dim(5)),squeeze(data(k(1),k(2),k(3),k(4),:)));
        hold on
        imHandles{5} = plot(k(5),data(k(1),k(2),k(3),k(4),k(5)),'r*','MarkerSize',8); hold off
        value = data(k(1),k(2),k(3),k(4),k(5)); 
        if value>10^5
            value=num2str(value,'%10.4e\n');
        else
            value=num2str(value);
        end
        titleHandles{4} = title({['Timeseries at point [',num2str(k(1)),' ',num2str(k(2)),' ',num2str(k(3)), ' ', num2str(k(4)), ' ', num2str(k(5)),']. Value: ',value],''}, 'FontSize', 12);
        
        for i = 1:3
            set(titleHandles{i}, 'String', [sprintf('%s%d',titleLines{i}, k(offsliceCoord(i))),' of ', num2str(dim(offsliceCoord(i)))]);
        end
    end
    
    %initialize crosshair and hide it
    x = zeros(2,2);
    y = zeros(2,2);
    subplot(221);
    crossHandle{1} = line(x,y);
    subplot(222);
    crossHandle{2} = line(x,y);
    subplot(223);
    crossHandle{3} = line(x,y);
    for i=1:3
        set(crossHandle{i}, 'Color','red');
        set(crossHandle{i}, 'HandleVisibility', 'off');
    end

%% Click function
    function click(src,evnt)
        p = get(gca, 'CurrentPoint');
        p = round(p);
        s = get(gcf, 'SelectionType');
        lastPoint = [p(1); p(3)];
        
        % Determine in which view was clicked using the lastPoint parameter
        temp=gca;
        if (strfind(temp.Title.String,'Transversal')==1)
            whichView = 1;
        elseif (strfind(temp.Title.String,'Coronal')==1)
            whichView = 2;
        elseif (strfind(temp.Title.String,'Sagittal')==1)
            whichView = 3;
        else
            whichView = 4;
        end
        
        %old cold, still works with older versions, where axes handles
        %and/or cell2mat were different
        %whichView = find(cell2mat(handle) == gca); 
        
        xlim = get(gca, 'XLim');
        ylim = get(gca, 'YLim');
        % If you click outside a figure, nothing will be updated (i.e.
        % the figure you clicked last, will change)
        if lastPoint(1) >= xlim(1) && lastPoint(1) <= xlim(2) && ...
                lastPoint(2) >= ylim(1) && lastPoint(2) <= ylim(2) && whichView ~= 4
          
            %If it's a left click activate window, scaling and scrolling
            if strcmp(s,'normal')
                set(f1, 'WindowButtonMotionFcn', @dragCallback);
                set(f1, 'WindowScrollWheelFcn', @scroller);
                set(f1, 'WindowKeyPressFcn',@scale);
                set(f1, 'WindowButtonUpFcn', @stopdrag);
            
            %If it's a right click activate the crosshair
            elseif strcmp(s,'alt')
%                 set(f1, 'WindowButtonMotionFcn', @dragCallback_2);
                set(f1, 'WindowButtonUpFcn', @stopdrag);
                
                border = [get(gca, 'xlim') get(gca, 'ylim')];
                x = [ border(1) p(1)
                      border(2) p(1) ];
                y = [ p(3)    border(3)
                      p(3)    border(4) ];
                
                set(crossHandle{whichView}(1),'xdata',x(:,1));
                set(crossHandle{whichView}(1),'ydata',y(:,1));
                set(crossHandle{whichView}(2),'xdata',x(:,2));
                set(crossHandle{whichView}(2),'ydata',y(:,2));
                
                if whichView == 1
                    set(crossHandle{1}, 'Visible', 'on');
                    set(crossHandle{2}, 'Visible', 'off');
                    set(crossHandle{3}, 'Visible', 'off');
                    updateTitle(p(3),p(1),k(3),k(4),k(5));
                    updateTimeSeries(p(3),p(1),k(3),k(4),k(5));
                    crossLoc=[p(3),p(1),k(3)];
                elseif whichView == 2
                    set(crossHandle{2}, 'Visible', 'on');
                    set(crossHandle{1}, 'Visible', 'off');
                    set(crossHandle{3}, 'Visible', 'off');
                    updateTitle(k(1),p(1),p(3),k(4),k(5));
                    updateTimeSeries(k(1),p(1),p(3),k(4),k(5));
                    crossLoc=[k(1),p(1),p(3)];
                elseif whichView == 3
                    set(crossHandle{3}, 'Visible', 'on');
                    set(crossHandle{1}, 'Visible', 'off');
                    set(crossHandle{2}, 'Visible', 'off');
                    updateTitle(p(1),k(2),p(3),k(4),k(5));
                    updateTimeSeries(p(1),k(2),p(3),k(4),k(5));
                    crossLoc=[p(1),k(2),p(3)];
                else
                    return %do nothing
                end      
            elseif strcmp(s,'extend')
                set(f1, 'WindowButtonMotionFcn', @dragCallback_2);
                set(f1, 'WindowScrollWheelFcn', @scroller);
                set(f1, 'WindowKeyPressFcn',@scale);
                set(f1, 'WindowButtonUpFcn', @stopdrag);   
            end
        end
    end

%% Drag function
    function dragCallback(varargin)
        p = get(gca, 'CurrentPoint');
        p = round(p);
        if ndim>3
            

            motionV = round([p(1); p(3)] - lastPoint);
            if abs(motionV(2)) < 2
                return
            end
            lastPoint = [p(1); p(3)];
            
            % Update the last dimension of k, the time dimension
            k(5) = k(5)+round(motionV(2)/2);
            if k(5) > size(data,5);
                k(5) =1;
            elseif k(5)<1
                k(5)=size(data,5);
            end;
            
            % REMARK, update the title name with the dynamic name
            set(imHandles{1}, 'CData', squeeze(data(:,:,k(3),k(4),k(5)))); %XY view (transversal)
            set(imHandles{2}, 'CData', squeeze(data(k(1),:,:,k(4),k(5)))'); %YZ view (coronal)
            set(imHandles{3}, 'CData', squeeze(data(:,k(2),:,k(4),k(5)))'); %XZ view (sagittal)
            %Set asterisk
            old_Ydat = get(imHandles{4},'Ydata');
            set(imHandles{5},'YData',old_Ydat(k(5)),'XData',k(5));
            updateTitle(crossLoc(1),crossLoc(2),crossLoc(3),k(4),k(5));
            
            

        else
            return %dragging only does something with 5D datasets
        end;
    end

    function dragCallback_2(varargin)
        p = get(gca, 'CurrentPoint');
        p = round(p);
        if ndim>3
           

            motionV = round([p(1); p(3)] - lastPoint);
            if abs(motionV(2)) < 2
                return
            end
            lastPoint = [p(1); p(3)];

            % Update the fourth dimension of k, the time dimension
            k(4) = k(4)+round(motionV(2)/2);
            if k(4) > size(data,4);
                k(4) =1;
            elseif k(4)<1
                k(4)=size(data,4);
            end;

            % REMARK, update the title name with the dynamic name
            set(imHandles{1}, 'CData', squeeze(data(:,:,k(3),k(4),k(5)))); %XY view (transversal)
            set(imHandles{2}, 'CData', squeeze(data(k(1),:,:,k(4),k(5)))'); %YZ view (coronal)
            set(imHandles{3}, 'CData', squeeze(data(:,k(2),:,k(4),k(5)))'); %XZ view (sagittal)
            %Set asterisk
            old_Ydat = get(imHandles{4},'Ydata');
            set(imHandles{5},'YData',old_Ydat(k(5)),'XData',k(5));
            updateTitle(crossLoc(1),crossLoc(2),crossLoc(3),k(4),k(5));
 


%             updateTitle(k(1),k(2),k(3),k(4),k(5));
%             updateTimeSeries(k(1),k(2),k(3),k(4),k(5));
            
            
        else
            return %dragging only does something with 5D datasets
        end;
    end

%% Stop drag function
% Will execute when mousebutton is released: 'buttonUp'
    function stopdrag(varargin)      
        set(f1,'WindowButtonMotionFcn','');
        
        s = get(gcf, 'SelectionType');
        if strcmp(s,'alt')
            p = get(gca, 'CurrentPoint');
            p=round(p);
            lastPoint = [p(1); p(3)];
            
            temp=gca;
            if (strfind(temp.Title.String,'Transversal')==1)
                whichView = 1;
            elseif (strfind(temp.Title.String,'Coronal')==1)
                whichView = 2;
            elseif (strfind(temp.Title.String,'Sagittal')==1)
                whichView = 3;
            else
                whichView = 4;
            end
            %whichView = find(cell2mat(handle) == gca);
            xlim = get(gca, 'XLim');
            ylim = get(gca, 'YLim');
            % If you click outside a figure, nothing will be updated (i.e.
            % the figure you clicked last, will change)
            if lastPoint(1) >= xlim(1) && lastPoint(1) <= xlim(2) && ...
                    lastPoint(2) >= ylim(1) && lastPoint(2) <= ylim(2)
                if ndim>3
                    if whichView == 1
                        updateTimeSeries(p(3),p(1),k(3),k(4),k(5));
                        updateTitle(p(3),p(1),k(3),k(4),k(5));
                    elseif whichView == 2
                        updateTimeSeries(k(1),p(1),p(3),k(4),k(5));
                        updateTitle(k(1),p(1),p(3),k(4),k(5));
                    elseif whichView == 3
                        updateTimeSeries(p(1),k(2),p(3),k(4),k(5));
                        updateTitle(p(1),k(2),p(3),k(4),k(5));
                    end
                    if whichView == 1 && p(1) == 18 && p(3) == 17
                        figure('Name','Such Wow'); spy
                    end
                else %dim=3
                    if whichView == 1
                        updateTitle(p(3),p(1),k(3),k(4),k(5));
                    elseif whichView == 2
                        updateTitle(k(1),p(1),p(3),k(4),k(5));
                    elseif whichView == 3
                        updateTitle(p(1),k(2),p(3),k(4),k(5));
                    end
                end
            end
            
            border = [get(gca, 'xlim') get(gca, 'ylim')];
            x = [ border(1) p(1)
                  border(2) p(1) ];
            y = [ p(3)    border(3)
                  p(3)    border(4) ];

            set(crossHandle{whichView}(1),'xdata',x(:,1));
            set(crossHandle{whichView}(1),'ydata',y(:,1));
            set(crossHandle{whichView}(2),'xdata',x(:,2));
            set(crossHandle{whichView}(2),'ydata',y(:,2));

            if whichView == 1
                set(crossHandle{1}, 'Visible', 'on');
                set(crossHandle{2}, 'Visible', 'off');
                set(crossHandle{3}, 'Visible', 'off');
                updateTitle(p(3),p(1),k(3),k(4),k(5));
                updateTimeSeries(p(3),p(1),k(3),k(4),k(5));
                crossLoc=[p(3),p(1),k(3)];
            elseif whichView == 2
                set(crossHandle{2}, 'Visible', 'on');
                set(crossHandle{1}, 'Visible', 'off');
                set(crossHandle{3}, 'Visible', 'off');
                updateTitle(k(1),p(1),p(3),k(4),k(5));
                updateTimeSeries(k(1),p(1),p(3),k(4),k(5));
                crossLoc=[k(1),p(1),p(3)];
            elseif whichView == 3
                set(crossHandle{3}, 'Visible', 'on');
                set(crossHandle{1}, 'Visible', 'off');
                set(crossHandle{2}, 'Visible', 'off');
                updateTitle(p(1),k(2),p(3),k(4),k(5));
                updateTimeSeries(p(1),k(2),p(3),k(4),k(5));
                crossLoc=[p(1),k(2),p(3)];
            else
                return %do nothing
            end      
        end
    end   


%% Scroll function
    function scroller(src,evnt)
        % Callback function to scroll through the images, update the
        % coordinates within k
        if evnt.VerticalScrollCount>0;
            newslice = k(offsliceCoord(whichView)) - 1;
        else
            newslice = k(offsliceCoord(whichView)) + 1;
        end
        upDateSlice(newslice)
    end

%% Function
    function upDateSlice(newslice)
        % Update the slice.
        % If you are moving in the xy plane, the update should be in the z
        % direction, etc.        
        if newslice > 0 && newslice <= dim(offsliceCoord(whichView))
            k(offsliceCoord(whichView)) = newslice;
        end
        
        subplot(handle{whichView});
        if ndim==3
            if whichView == 1
                set(imHandles{1}, 'CData', squeeze(data(:,:,k(3),k(4),k(5))));
                updateTitle(p(3),p(1),k(3),1,1);
                
            elseif whichView == 2
                set(imHandles{2}, 'CData', squeeze(data(k(1),:,:,k(4),k(5)))');
                updateTitle(k(1),p(1),p(3),1,1);
            else
                set(imHandles{3}, 'CData', squeeze(data(:,k(2),:,k(4),k(5)))');
                updateTitle(p(1),k(2),p(3),1,1);
            end
        else
            if whichView == 1
                set(imHandles{1}, 'CData', squeeze(data(:,:,k(3),k(4),k(5))));
                updateTitle(crossLoc(1),crossLoc(2),k(3),k(4),k(5));
                updateTimeSeries(crossLoc(1),crossLoc(2),k(3),k(4),k(5));
            elseif whichView == 2
                set(imHandles{2}, 'CData', squeeze(data(k(1),:,:,k(4),k(5)))');
                updateTitle(k(1),p(1),p(3),k(4),k(5));
                updateTimeSeries(k(1),p(1),p(3),k(4),k(5));
            else
                set(imHandles{3}, 'CData', squeeze(data(:,k(2),:,k(4),k(5)))');
                updateTitle(p(1),k(2),p(3),k(4),k(5));
                updateTimeSeries(p(1),k(2),p(3),k(4),k(5));
            end
        end
        set(titleHandles{whichView}, 'String', ...
            [sprintf('%s%d', titleLines{whichView}, k(offsliceCoord(whichView))),' of ', num2str(dim(offsliceCoord(whichView)))]);
    end

    function updateTitle(k1,k2,k3,k4,k5)
        if ndim==3
            value=data(k1,k2,k3); 
            if value>10^5; %set precision for large numbers
                value=num2str(value,'%10.4e\n');
            else
                value=num2str(value);
            end
            set(titleHandles{4}, 'String', {['Crosshair at point [',num2str(k1),' ',num2str(k2),' ',num2str(k3),']. Value: ',value],''});
            
        else
            value=data(k1,k2,k3,k4,k5); 
            if value>10^5;
                value=num2str(value,'%10.4e\n');
            else
                value=num2str(value);
            end
            set(titleHandles{4}, 'String', {['Timeseries at point [',num2str(k1),' ',num2str(k2),' ',num2str(k3), ' ', num2str(k4), ' ', num2str(k5),']. Value: ',value],''});
        end
    end

    function updateTimeSeries(k1,k2,k3,k4,k5)
        if ndim>3
            set(imHandles{4}, 'YData',squeeze(data(k1,k2,k3,k4,:)));
            set(imHandles{5}, 'YData',squeeze(data(k1,k2,k3,k4,k5)),'XData',k5);
        end
    end

%% Keyboard functions
    function scale(src,evnt)
        % Scaling function
        if strcmp(evnt.Key,'downarrow');
            range(1)=range(1)-step;            
        elseif strcmp(evnt.Key,'uparrow');
            range(1)=range(1)+step;
            if range(1)>range(2)
                range(1)=range(2)-1;
            end
        elseif strcmp(evnt.Key,'leftarrow');
            range(2)=range(2)-step;
            if range(2)<range(1)
                range(2)=range(1)+1;
            end
        elseif strcmp(evnt.Key,'rightarrow');
            range(2)=range(2)+step;
        else
            %donothing
        end
            
        % To be able to update CLim, one has to grab the subplots, and
        % update the CLim range
        for c = get(f1,'Children')
            try
                set(c,'CLim',range)
                save('CLimRange.mat','range'); %save the CLim to file for later use
            catch exception
                continue
            end
        end
    end

end