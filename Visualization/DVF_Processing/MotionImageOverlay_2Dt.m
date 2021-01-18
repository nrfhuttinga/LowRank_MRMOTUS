function [imgs_out,cm] = MotionImageOverlay_2Dt(Images,motionFields,varargin)
    % Inputs: 
    %   - Image:                 2d+t images used to overlay motionField
    %   - motionFields{t}(:,:,1) :  motionFields in horizontal direction at time t; 
    %                            left is negative, right is positive
    %   - motionFields{t}(:,:,2):   motionFields in vertical direction at time t; 
    %                            top is postive, bottom is negative 
    %   - For images: 3rd dimension is animated
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
        maskk = ones(size(Images(:,:,1)));
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

    maskk=rot90(maskk);
    if rotations == 1 % counterclockwise
        % x -> y
        % y -> -x

        Images = rot90(Images,rotations);
        [~,~,dynamics]=size(Images);


        for i=1:dynamics
            mf_new(:,:,:,1)=motionFields{i}(:,:,2);
            mf_new(:,:,:,2)=-motionFields{i}(:,:,1);
            motionFields{i}=rot90(mf_new,rotations);
        end



    elseif rotations == 2
        % x -> y -> -x
        % y -> -x -> -y


        Images = rot90(Images,rotations);
        [~,~,dynamics]=size(Images);


        for i=1:dynamics
            mf_new(:,:,:,2)=-motionFields{i}(:,:,2);
            mf_new(:,:,:,1)=motionFields{i}(:,:,1);
            motionFields{i}=rot90(mf_new,rotations);
        end





    elseif rotations == 3
        % x -> y -> -x -> -y
        % y -> -x -> -y -> x

        Images = rot90(Images,rotations);
        [~,~,dynamics]=size(Images);


        for i=1:dynamics
    %         mf_old = motionFields{i};
            mf_new(:,:,:,1)=-motionFields{i}(:,:,2);
            mf_new(:,:,:,2)=motionFields{i}(:,:,1);
            motionFields{i}=rot90(mf_new,rotations);
        end   

    end


    [dim2,dim1,dynamics]=size(Images);



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
        Image           = squeeze(Images(:,:,i)).*maskk;

        motionField_HOR(abs(Image)<threshold | ~maskk)=0;
        motionField_VER(abs(Image)<threshold | ~maskk)=0;


        if im_auto_scale == 1
            imagesc(xlimit,ylimit,abs(Image),[0 max(abs(Image(:)))]);
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




