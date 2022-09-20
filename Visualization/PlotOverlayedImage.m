function varargout=PlotOverlayedImage(im_1,im_2,alpha_,axes_object,caxis_scale,varargin)
% Function to overlay 'im_2' on 'im_1' with opacity 'alpha_', in 'axes_object'.
%
% Inputs:
%   im_1            - 2D underlay image
%   im_2            - 2D overlay image
%   alpha_          - opacity of the overlayed image
%   axes_objects    - axes to plot on
%   caxis_scale     - range of the colorbar
%   varargin{1}     - horizontal colobar location offset
%   varargin{2}     - flag to plot colorbar or not
%   varargin{3}     - color bar width
%
%
% Niek Huttinga, UMC Utrecht, 2020

% im_1 = imresize(im_1,2);
% im_2 = imresize(im_2,2);

im_1 = demax(abs(im_1));

if numel(varargin)>0 && ~isempty(varargin{1})
    cb_offset = varargin{1};
else
    cb_offset = 0;
end

if numel(varargin)>1 && ~isempty(varargin{2})
    cb_flag = varargin{2};
else
    cb_flag = 0;
end

if nargin<=3 || isempty(axes_object)
    axes_object = gca;
end

if nargin<=4 || isempty(caxis_scale)
    caxis_scale = [min(im_2(:)) max(im_2(:))];
end

ax=axes_object;
p=get(ax,'pos');


if numel(varargin)>2 && ~isempty(varargin{3})
    cb_width = varargin{3};
else
    cb_width = p(3)/15;
end

if numel(varargin)>3 && ~isempty(varargin{4})
    colormp = varargin{4};
else
    colormp = 'jet';
end


delete(ax);
ax1=axes;
set(ax1,'pos',p);
% cla('reset')
UnderlayImage = imagesc(mat2gray(im_1));axis off; %axis image;
axis tight;
colormap(ax1,'gray');
ax2=axes;ax2.Visible='off';set(ax2,'pos',p);
OverlayImage = imagesc(ax2,im_2);% axis image;
axis off;
ax_keep = ax2;
caxis(ax2,caxis_scale);
colormap( ax2, colormp );

if cb_flag
    cb=colorbar( OverlayImage.Parent );
%     daspect(OverlayImage.Parent ,[1 1 1]);
    cb.FontSize = 14;
    cb.FontWeight = 'bold';
    cb.Color='white';
    cb.Position = p;
%     cb.Position(1) = p(1)+cb_offset;
    cb.Position(3) = cb_width;
    cb.Ticks=0:cb.Ticks(end);
    varargout{1}=cb;
end
axis tight;

if numel(varargin)<5
    alpha = ones(size(im_1))*alpha_;
else
    alpha = varargin{5}*alpha_;
end
    
% alpha(abs(im_1)<0.8)=alpha_*(im_1(abs(im_1)<0.8)/0.8);
% alpha(abs(im_1)>=0.8)=alpha_;%(abs(im_1)==0)=0;

set( OverlayImage, 'AlphaData', alpha );
% hold off;

end
