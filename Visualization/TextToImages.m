function im_out = TextToImages(text_cell,varargin)

    dynamics = numel(text_cell);
    
    if nargin>1
        width = varargin{1};
    else
        width = 500;
    end
        
    for ii=1:dynamics
        im_out(:,:,ii)=text2im(text_cell{ii},width);
    end
        
end
