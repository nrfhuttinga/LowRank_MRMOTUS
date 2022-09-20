function CellsToGif(ImageCell,ColormapCell,varargin)
    % Function to convert cell of Images and Colormaps to GIF. Both should be
    % made with the rgb2ind function that is applied to a getframe(.) object.
    % See also MotionImageOverlay_2Dt function.
    % If no FileName specified, writes to 'working directory'/images.gif. 
    % If no DelayTime specified uses .1.
    % [DelayTime]=ms;
    %
    % Niek Huttinga - UMC Utrecht - 2020

    warning('off')

    if isempty(ColormapCell{1}) % grayscale
        grayscale = 0;
    else
        grayscale = 1;
    end

    ims = size(ImageCell{1},1);

    if grayscale
        for i=1:numel(ImageCell)
            rsz = imresize(ImageCell{i},[ims,ims]);
            [Images(:,:,i),ColormapCell{i}] = imresize(ImageCell{i},ColormapCell{i},[ims,ims]);
        end
    else
        for i=1:numel(ImageCell)
            [Images(:,:,i),ColormapCell{i}] = imresize(ImageCell{i},ColormapCell{i},[ims,ims]);
        end
    end







    AnimateDimension=3;

    if nargin > 2
        DelayTime = varargin{1};
    else
        DelayTime = .01;
    end

    if nargin > 3
        FileName = varargin{2};
    else
        FileName = 'images.gif';

    end

    if nargin > 4
        text_flag = 1;
        text_cell = varargin{3};

        if numel(text_cell)~=size(Images,AnimateDimension)
            numel(text_cell)
            size(Images,AnimateDimension)
            error('Text cell dynamics not same as image dynamics');
        end
    else
        text_flag = 0;
    end


    if text_flag


            text_img = TextToImages(text_cell,round(size(Images,2)/2)*2);

    %         for i=1:size(text_img,3)

    %         end
            for i=1:size(text_img,3)
                ColormapCell{i}(end+1,:)=[1 1 1];
                text_img(:,:,i) = demax(uint8(text_img(:,:,i)))*size(ColormapCell{i},1);
            end

    %         text_img = bsxfun(@times,max(Images,[],[1 2]),demax(uint8(text_img)));
            pad = ceil(abs((size(text_img,2)-size(Images,2))/2));
            Images = cat(1,[uint8(zeros([size(Images,1),pad,size_ext(Images,[3])])),Images],text_img);
    %         slicer5d(Images)
    end




    for loopIndex=1:numel(ImageCell)

    %     ExtractedImage = squeeze(ImageCell{loopIndex});
        ExtractedImage = squeeze(Images(:,:,loopIndex));
%         ExtractedImage(isnan(ExtractedImage(:)))=0;


        % add padding on top for text
        % ExtractedImage = [zeros(20,size(ExtractedImage,2));ExtractedImage];


        if loopIndex == 1
            imwrite(ExtractedImage,ColormapCell{loopIndex},FileName,'gif','LoopCount',Inf,'DelayTime',DelayTime);
        else
            imwrite(ExtractedImage,ColormapCell{loopIndex},FileName,'gif','WriteMode','append','DelayTime',DelayTime);
        end


    end

    warning('on')
end
