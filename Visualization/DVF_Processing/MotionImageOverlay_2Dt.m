function [imgs_out,cm] = MotionImageOverlay_2Dt(Images,motionFields,varargin)
    % Inputs: 
    %   - Image:                    :   N x N x 1 x 1 x T 2d+t images used to overlay motionFields on
    %   - motionFields{t}(:,:,1)    :   motionFields in horizontal direction at time t; 
    %                                   left is negative, right is positive
    %   - motionFields{t}(:,:,2)    :   motionFields in vertical direction at time t; 
    %                                   top is postive, bottom is negative 
    %   - dimension                 :   index of spatial dimension along which to select motion-fields for visualization [1/2/3]
    %   - slice                     :   slice to select along 'dimension'
    %   - varargin{1}               :   threshold for visualization
    %   - varargin{2}               :   downsampling of motion-fields for visualization
    %   - varargin{3}               :   color of the motion-field vectors
    %   - varargin{4}               :   padding to remove on the boundary of the motion-fields [arrow heads should not go outside 'Images']
    %   - varargin{5}               :   motion-field scaling factor for visualization
    %   - varargin{6}               :   number of rotations counterclockwise for visualization
    %   - varargin{7}               :   binary 2D visualization mask corresponding to the 'slice' selected in 'dimension'
    %
    % Outputs:
    %   - cells with images and colormaps
    %
    %
    % Niek Huttinga, UMC Utrecht, 2020

    if nargin < 10
        im_auto_scale = 1;
    else
        im_auto_scale=varargin{8};
    end

    if nargin < 9 || isempty(varargin{7})
        maskk = ones(size(Images(:,:,:,:,1)));
    else
        maskk = varargin{7};
    end


    if nargin < 8
        rotations = 0;
    else
        rotations = varargin{6};
    end


    if nargin < 7
        scaling = 500;
    else
        scaling = varargin{5};
    end

    if nargin < 6
        pad = 10;
    else
        pad = varargin{4};
    end


    if nargin < 5
        color = 'g';
    else
        color = varargin{3};
    end

    if nargin < 4
        motion_field_display_factor = 5;
    else
        motion_field_display_factor = varargin{2};
    end


    if nargin < 3
        threshold   = 0.25;
    else
        threshold    = varargin{1};
    end

    if nargin < 2
        error('Not enough input arguments')
    end
    set(gcf,'Color','w')

    [~,~,~,~,dynamics]=size(Images);

    if rotations == 1 % counterclockwise
        % x -> y
        % y -> -x

        Images = rot90(Images,rotations);


        for i=1:dynamics
            mf_new(:,:,:,1)=-motionFields{i}(:,:,2);
            mf_new(:,:,:,2)=-motionFields{i}(:,:,1);
            motionFields{i}=rot90(mf_new,rotations);
        end



    elseif rotations == 2
        % x -> y -> -x
        % y -> -x -> -y


        Images = rot90(Images,rotations);


        for i=1:dynamics
            mf_new(:,:,:,2)=-motionFields{i}(:,:,2);
            mf_new(:,:,:,1)=motionFields{i}(:,:,1);
            motionFields{i}=rot90(mf_new,rotations);
        end





    elseif rotations == 3
        % x -> y -> -x -> -y
        % y -> -x -> -y -> x

        Images = rot90(Images,rotations);


        for i=1:dynamics
            mf_new(:,:,:,1)=-motionFields{i}(:,:,2);
            mf_new(:,:,:,2)=motionFields{i}(:,:,1);
            motionFields{i}=rot90(mf_new,rotations);
        end   

    elseif rotations == 0
        for i=1:dynamics
            motionFields{i}(:,:,2) = motionFields{i}(:,:,2);
            motionFields{i}(:,:,1) = motionFields{i}(:,:,1);

        end
    end
    
    


    [dim2,dim1,~]=size(Images);



    quiver(zeros(dim2,dim1),zeros(dim2,dim1),0,'color',color,'linewidth',1,'clipping','off');axis([0 dim2 0 dim1]);
    quiver_mask = zeros(dim2,dim1);
    quiver_mask(pad+1:motion_field_display_factor:end-pad,pad+1:motion_field_display_factor:end-pad)=1;

    axis image;
    hax = gca;
    xlimit = hax.XLim;
    ylimit = hax.YLim;
    close all;
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','black');
    for i=1:dynamics


        motionField_HOR = scaling*squeeze(motionFields{i}(:,:,1)).*quiver_mask;
        motionField_VER = scaling*squeeze(motionFields{i}(:,:,2)).*quiver_mask;
        Image           = squeeze(Images(:,:,:,:,i)).*maskk;

        motionField_HOR(abs(Image)<threshold | ~maskk)=0;
        motionField_VER(abs(Image)<threshold | ~maskk)=0;


        if im_auto_scale == 1
            imagesc(xlimit,ylimit,abs(Image),[0 mean(max(abs(Images),[],[1,2])*0.9)]);
        else
            imagesc(xlimit,ylimit,abs(Image),[0 im_auto_scale]);
        end

        colormap('gray')
        hold on;
        quiver(motionField_HOR, motionField_VER,0,'color',color,'linewidth',2,'clipping','off','AutoScaleFactor',5000000); axis([0 dim2/motion_field_display_factor 0 dim1/motion_field_display_factor]);
        axis image;
        axis off;

        [imgs_out{i},cm{i}] = to_cells(gca);

    end
    
end




