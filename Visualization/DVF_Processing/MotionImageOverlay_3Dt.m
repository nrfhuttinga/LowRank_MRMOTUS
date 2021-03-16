function [imgs_out,cm] = MotionImageOverlay_3Dt(Images,motionFields,dimension,slice,varargin)
    % Inputs: 
    %   - Images:                 	:   N x N x N x 1 x T images on which to overlay motion-fields
    %   - motionFields              :   motionFields as N^3 x 3 x T matrix
    %   - dimension                 :   index of spatial dimension along which to select motion-fields for visualization [1/2/3]
    %   - slice                     :   slice to select along 'dimension'
    %   - varargin{1}               :   threshold for visualization
    %   - varargin{2}               :   downsampling of motion-fields for visualization
    %   - varargin{3}               :   color of the motion-field vectors
    %   - varargin{4}               :   padding to remove on the boundary of the motion-fields [arrow heads should not go outside 'Images']
    %   - varargin{5}               :   motion-field scaling factor for visualization
    %   - varargin{6}               :   number of rotations counterclockwise for visualization
    %   - varargin{7}               :   binary 2D visualization mask corresponding to the 'slice' selected in 'dimension'
    
    % Outputs:
    %   - cells with images and colormaps
    %
    %
    % Niek Huttinga, UMC Utrecht, 2020


    [dim2,dim1,dim3,~,dynamics]=size(Images);

    motionFields_reshaped = reshape(motionFields,dim2,dim1,dim3,size(motionFields,2),dynamics);

        
    if nargin < 12
        im_auto_scale = 1;
    else
        im_auto_scale = varargin{8};
    end

    if nargin<11
        maskk = [];
    else
        maskk = varargin{7};
    end


    if nargin<10 % counterclockwise
        rotations = 0;
    else
        rotations = varargin{6};
    end


    if nargin < 9
        scaling = 500;
    else
        scaling = varargin{5};
    end

    if nargin < 8
        pad = 10;
    else
        pad = varargin{4};
    end


    if nargin < 7
        color = 'g';
    else
        color = varargin{3};
    end

    if nargin < 6
        motion_field_display_factor = 5;
    else
        motion_field_display_factor = varargin{2};
    end


    if nargin < 5
        threshold   = 0.25;
    else
        threshold    = varargin{1};
    end




    if nargin < 4
        error('Not enough input arguments')
    end
    set(gcf,'Color','w')
    pause(1)
    quiver(zeros(dim2,dim1),zeros(dim2,dim1),0,'color',color,'linewidth',1,'clipping','off');axis([0 dim2 0 dim1]);
    quiver_mask = zeros(dim2,dim1);
    quiver_mask(2:motion_field_display_factor:end-1,2:motion_field_display_factor:end-1)=1;
    axis image;
    hax = gca;
    xlimit = hax.XLim;
    ylimit = hax.YLim;
    close all;
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    for i=1:dynamics

        mf = {-motionFields_reshaped(:,:,:,1,i)*scaling,-motionFields_reshaped(:,:,:,2,i)*scaling,motionFields_reshaped(:,:,:,3,i)*scaling};


        Image = squeeze(Images(:,:,:,i));

        [dvf{i},image(:,:,:,:,i)] = ExtractSlices_3D(Image,mf,dimension,slice);



    end

    [imgs_out,cm] = MotionImageOverlay_2Dt(image,dvf,threshold,motion_field_display_factor,color,pad,scaling,rotations,maskk,im_auto_scale);

end
