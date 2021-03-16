function PlotOverlayedImage(im_1,im_2,alpha_,axes_object,caxis_scale,varargin)
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

im_1 = demax(abs(im_1));

if numel(varargin)>0
    cb_offset = varargin{1};
else
    cb_offset = 0;
end

if numel(varargin)>1
    cb_flag = varargin{2};
else
    cb_flag = 0;
end

ax=axes_object;
p=get(ax,'pos');


if numel(varargin)>2
    cb_width = varargin{3};
else
    cb_width = p(3)/15;
end

if numel(varargin)>3
    colormp = varargin{4};
else
    colormp = 'jet';
end


delete(ax);
ax1=axes;
set(ax1,'pos',p);
UnderlayImage = imagesc(mat2gray(im_1));axis off; axis image;
colormap(ax1,'gray');
ax2=axes;ax2.Visible='off';set(ax2,'pos',p);
OverlayImage = imagesc(ax2,im_2);axis off; axis image;
caxis(ax2,caxis_scale);
colormap( ax2, colormp );

if cb_flag
    cb=colorbar( OverlayImage.Parent );
    cb.FontSize = 14;
    cb.FontWeight = 'bold';
    cb.Color='white';
    cb.Position = p;
    cb.Position(1) = p(1)+cb_offset;
    cb.Position(3) = cb_width;
    cb.Ticks=[0 1 2 3];
end


alpha = ones(size(im_1))*alpha_;
alpha(abs(im_1)==0)=0;

set( OverlayImage, 'AlphaData', alpha );


end
