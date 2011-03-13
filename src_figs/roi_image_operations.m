function roi_image_operations(action, varargin)
%ROI_IMAGE_OPERATIONS Region-of-interest processing.
%
%
%   This function has bits of code from ROIDEMO and was modified to work to MIRONE  
%
%   Function Subroutines:
%   InitializeWIN - Set up controls and Images.
%
%   ApplyOperation - Look at the popup to see what operation we are doing
%                    and then do it to the original image only in the
%                    Region-of-Interest specified by the mask.
%
%   UpdateOperation - Updates the ROI operation from the current value of
%                     the operation popup menu.
% 
%   ImageReset - Reset Image to its initial state

if (nargin < 1),   action='InitializeWIN';  end

if (nargin == 2 && ishandle(action))
    varargin{2} = varargin{1};  % "varargin{1}" should contain the polygon ROI vertices
    varargin{1} = action;       % "action" should contain the Axis handle
    action='InitializeWIN';
end    

feval(action,varargin{:})
return

% -------------------------------------------------------------------------------------------
function InitializeWIN(opt1,opt2)
% OPT1, when not empty, contains the Axis handle (gca in Mirone)
% OPT2, when OPT1 is not empty, contains the polygon ROI vertices (MUST be a two column matrix)

% If roi_image_operations is already running, bring it to the foreground
h = findobj(allchild(0), 'tag', 'ROI Processing');
if ~isempty(h)
   figure(h(1))
   return
end

screenD = get(0, 'ScreenDepth');
if (screenD > 8),   grayres=256;
else                grayres=128;    end

RoiImageOpsFig = figure( ...
   'Name','ROI Processing', ...
   'NumberTitle','off', 'HandleVisibility', 'on', ...
   'tag', 'ROI Processing', ...
   'Visible','off', 'Resize', 'off',...
   'BusyAction','Queue','Interruptible','off', ...
   'Color', [.8 .8 .8], 'IntegerHandle', 'off', ......
   'DoubleBuffer', 'on', 'Colormap', gray(grayres),'MenuBar','none');

figpos = get(RoiImageOpsFig, 'position');
% Adjust the size of the figure window
figpos(3:4) = [200 190];
% screenSize = get(0,'ScreenSize');
set(RoiImageOpsFig, 'position', figpos);

% Colors
bgcolor = [0.45 0.45 0.45];  % Background color for frames
% wdcolor = [.8 .8 .8];  % Window color
fgcolor = [1 1 1];  % For text

% hs = 44;        % Horizantal Spacing
% vs = 22;        % Vertical Spacing
ifs = 11;       % Intraframe spacing

%====================================
% Parameters for all buttons and menus

Std.Interruptible = 'off';
Std.BusyAction = 'queue';

Ctl = Std;
Ctl.Units = 'Pixels';
Ctl.Parent = RoiImageOpsFig;

Frame = Ctl;
Frame.Style = 'Frame';

Btn = Ctl;
Btn.Style = 'pushbutton';
Btn.Enable = 'off';

Menu = Ctl;
Menu.Style = 'Popupmenu';

Text = Ctl;
Text.Style = 'text';
Text.HorizontalAlignment = 'left';
Text.BackgroundColor = bgcolor;
Text.ForegroundColor = fgcolor;

btnHt = 26;     txtHt = 18;     menuHt = 26;

%=================================
%  The frame
fleft = 10;     fbot = 10;
fwid = 180;     fht = 170;
ud.h.ControlFrame = uicontrol(Frame,'Position', [fleft fbot fwid fht], 'BackgroundColor', bgcolor);

menuWid = 170;  menuBot = 132;
labelBot = 158;    % For the labels above the three top menus

%================================
% Operations Menu:
opleft = 15;
ud.h.OperationsPop = uicontrol(Menu, ...
   'Position',[opleft menuBot menuWid menuHt], ...
   'String',...
     {'Histogram Equalization', ...
      'Adaptive Histogram Equalization', ...
      'Image Adjust', ...
      'Unsharp Masking', ...
      'Lowpass Filter', ...
      'Median Filter', ...
      'Brighten', 'Darken', ...
      'Increase Contrast', ...
      'Decrease Contrast', ...
      'Inpaint'}, ...
   'Callback','roi_image_operations(''UpdateOperation'')');
%       'Boundary Interpolation', ...
%       'Add Gaussian Noise' }, ...
% Text label for Image Menu Popup
uicontrol( Text, ...
   'Position',[opleft labelBot menuWid txtHt], ...
   'String','Region of Interest Operations:');

%====================================
% Buttons - Apply, Info and Close
ud.h.Close=uicontrol(Btn, ...
   'Position',[opleft fbot+ifs menuWid btnHt], ...
   'String','Close', ...
   'Callback','close(gcbf)');

ud.h.Reset=uicontrol(Btn, ...
   'Position',[opleft fbot+2*ifs+btnHt menuWid btnHt], ...
   'String','Reset', ...
   'Callback','roi_image_operations(''ImageReset'')');

ud.h.Apply=uicontrol(Btn, ...
   'Position',[opleft fbot+3*ifs+2*btnHt menuWid btnHt], ...
   'String','Apply', ...
   'Callback','roi_image_operations(''ApplyOperation'')');

ud.ParentAxis = opt1;
ud.ParentFig  = get(opt1,'Parent');

ud.OriginalImage = get(findobj(opt1,'Type','image'), 'Cdata');
x_lim = get(opt1,'xlim');   y_lim = get(opt1,'ylim');
ud.Mask = img_fun('roipoly_j',x_lim,y_lim,ud.OriginalImage,opt2(:,1),opt2(:,2));

set(RoiImageOpsFig, 'Userdata', ud);
set(RoiImageOpsFig, 'visible','on','HandleVisibility','callback');
set(ud.h.Close, 'Enable', 'on');

% -------------------------------------------------------------------------------------------
function UpdateOperation(ImageOpsFig)

if (nargin < 1),   ImageOpsFig = gcbf;  end

ud = get(ImageOpsFig, 'UserData');
opval = get(ud.h.OperationsPop, 'value');
opstr = get(ud.h.OperationsPop, 'string');
ud.Operation = opstr{opval};
set(ImageOpsFig, 'UserData', ud);
set([ud.h.Reset ud.h.Apply], 'Enable', 'on');

% -------------------------------------------------------------------------------------------
function ApplyOperation(ImageOpsFig)

if (nargin < 1),    ImageOpsFig = gcbf;     end
ud = get(ImageOpsFig, 'UserData');
orig = ud.OriginalImage;
mask = ud.Mask;
result = orig;
threeD = (ndims(orig) == 3);
handlesParent = guidata(ud.ParentFig);

gray_test = getappdata(ImageOpsFig,'gray_test');
if (isempty(gray_test))         % First time call, guess image "grayness"
	% Find if we are dealing with a color (gray == 0) or gray image (gray == 1)
    was_fake_gray = 0;          % Used when the orig image is MxNx3 but it is a gray image anyway
	if (threeD)                 % True color image
        % Pick up 64 random colors and see if they are equal (it should be enough as a test)
        n_to_test = min(64,size(orig,1));       % Play safe that we don't have an image smaller than 64
        tmp_m = round(rand(1,n_to_test)*size(orig,1));
        tmp_n = round(rand(1,n_to_test)*size(orig,2));
        tmp_m(tmp_m == 0) = 1;        tmp_n(tmp_n == 0) = 1;
        df = diff(double(orig(tmp_m,tmp_n,:)),1,3);
        if (any(df(:) ~= 0))
            is_gray = 0;                % It is really a true color image
        else
            is_gray = 1;                % Although is has 3 planes, it is a gray image
            orig(:,:,2:3) = [];         % Remove the non-needed pages
            result = orig;
            was_fake_gray = 1;
        end
	else                        % Indexed image
        cm = get(ud.ParentFig,'Colormap');      dif_cm = diff(cm,1,2);
        if isempty(find(dif_cm > 10e-10))
            is_gray = 1;
        else
            is_gray = 0;
        end
        clear cm dif_cm;
	end
	gray_test.is_gray = is_gray;
	gray_test.was_fake_gray = was_fake_gray;
	setappdata(ImageOpsFig,'gray_test',gray_test)
    
else            % We already know the image "grayness"
	is_gray = gray_test.is_gray;
	was_fake_gray = gray_test.was_fake_gray;
end

% Bellow was an attempt to increase the speed, but failed due to the fact that for indexed
% images we have to convert the whole image to YIQ color space and that's the memory consuming part.
% [i,j] = find(mask == 1);
% il_max = max(i);   il_min = min(i);
% ic_max = max(j);   ic_min = min(j);
% orig = orig(il_min:il_max,ic_min:ic_max);
% mask = mask(il_min:il_max,ic_min:ic_max);

switch ud.Operation
	case 'Histogram Equalization'
        if (threeD && ~is_gray)
            result = cvlib_mex('color',ud.OriginalImage,'rgb2YCrCb');       % Y is on the 3rth plane, bug?
            result_t = result(:,:,3);    result_t(mask) = img_fun('histeq_j',result_t(mask));
            result(:,:,3) = result_t;
            result = cvlib_mex('color',result,'YCrCb2rgb');      clear result_t;     %convert back to RGB
        else
            J = img_fun('histeq_j',orig(mask));
            result(mask) = J;           clear J;
        end
	case 'Adaptive Histogram Equalization'
        if (threeD && ~is_gray)
            result = cvlib_mex('color',ud.OriginalImage,'rgb2YCrCb');
            result_t = result(:,:,3);           tmp = img_fun('adapthisteq',result_t);
            result_t(mask) = tmp(mask);         result(:,:,3) = result_t;
            result = cvlib_mex('color',result,'YCrCb2rgb');      clear result_t tmp;     %convert back to RGB
        else
            tmp = img_fun('adapthisteq',result);
            result(mask) = tmp(mask);       clear tmp;
        end
	case {'Lowpass Filter' 'Unsharp Masking'}
        if (ud.Operation(1) == 'L')
            order = 15;  cutoff = 0.3;
            [f1,f2] = freqspace(order,'meshgrid');
            d = find(f1.^2+f2.^2 < cutoff^2);
            Hd = zeros(order);
            Hd(d) = 1;
            % Use hanning(15) as the window.  The window coefficients are
            % inlined to remove dependency on the Signal Processing Toolbox.
            h = img_fun('fwind1',Hd, ...
                   [0.0381 0.1464 0.3087 0.5000 0.6913 0.8536 ...
                    0.9619 1.0000 0.9619 0.8536 0.6913 0.5000 ...
                    0.3087 0.1464 0.0381]);
        else
            h = img_fun('fspecial','unsharp');
        end
        if (threeD && ~is_gray)
            result = cvlib_mex('color',ud.OriginalImage,'rgb2YCrCb');
            result(:,:,3) = img_fun('roifilt2',h, result(:,:,3), mask);
            result = cvlib_mex('color',result,'YCrCb2rgb');     %convert back to RGB
        else
            result = img_fun('roifilt2',h, orig, mask);
        end
	case 'Median Filter'
        if (threeD && ~is_gray)
            result = cvlib_mex('color',ud.OriginalImage,'rgb2YCrCb');
            result_t = result(:,:,3);           tmp = img_fun('medfilt2',result_t,[5 5]);
            result_t(mask) = tmp(mask);         result(:,:,3) = result_t;
            result = cvlib_mex('color',result,'YCrCb2rgb');      clear result_t tmp;     %convert back to RGB
        else
            tmp = img_fun('medfilt2',orig,[5 5]);
            result(mask) = tmp(mask);       clear tmp;
        end
	case 'Image Adjust'
        if (threeD && ~is_gray)
            result = cvlib_mex('color',ud.OriginalImage,'rgb2YCrCb');
            result_t = result(:,:,3);    result_t(mask) = img_fun('imadjust_j',result_t(mask));
            result(:,:,3) = result_t;
            result = cvlib_mex('color',result,'YCrCb2rgb');      clear result_t;     %convert back to RGB
        else
            result(mask) = img_fun('imadjust_j',result(mask));
        end
	case {'Brighten' 'Darken' 'Increase Contrast' 'Decrease Contrast'}
        if (strcmp(ud.Operation(1:2), 'Br')),           par = {[], [.25 1]};
        elseif (strcmp(ud.Operation(1:2), 'Da')),       par = {[], [0 .75]};
        elseif (strcmp(ud.Operation(1:2), 'In')),       par = {[.25 .75],[]};
        elseif (strcmp(ud.Operation(1:2), 'De')),       par = {[],[.25 .75]};
        end
        if (threeD && ~is_gray)
            result = cvlib_mex('color',ud.OriginalImage,'rgb2YCrCb');
            result_t = result(:,:,3);           result_t(mask) = img_fun('imadjust_j',result_t(mask),par{:});
            result(:,:,3) = result_t;
            result = cvlib_mex('color',result,'YCrCb2rgb');      clear result_t;     %convert back to RGB
        else
            tmp = img_fun('imadjust_j',orig, par{:});
            result(mask) = tmp(mask);       clear tmp;
        end
	case 'Inpaint'
        result = cvlib_mex('inpaint',ud.OriginalImage,mask);

	% case 'Boundary Interpolation'
	%     if ~is_gray
	%         if (threeD)
	%             result = rgb2YIQ(ud.OriginalImage);    %  convert to YIQ image
	%         else
	%             result = ind2rgb8(ud.OriginalImage,get(ud.ParentFig,'Colormap'));
	%             result = rgb2YIQ(result);              %  convert to YIQ image
	%         end
	%         result(:,:,1) = roifill(result(:,:,1), mask);
	%         result(:,:,2) = roifill(result(:,:,2), mask);
	%         result(:,:,3) = roifill(result(:,:,3), mask);
	%         result = YIQ2rgb(result);                  %convert back to RGB
	%     else        % processing a gray image
	%         result = roifill(orig, mask);
	%     end
	% case 'Add Gaussian Noise' 
	%    noisy_orig = imnoise(orig, 'gaussian', 0, .01);
	%    result(mask) = noisy_orig(mask);     clear noisy_orig;
	otherwise 
       warndlg('Invalid ROI operation','Warning');
       return
end

%handlesParent.origFig = result;
guidata(ud.ParentFig,handlesParent);

h_img = findobj(ud.ParentAxis,'Type','image');
if (was_fake_gray)                                  % It was a 3D image, but gray
    set(handlesParent.figure1,'colormap',gray(256))
    set(h_img,'CData', result);                     % Probaly the result is not what you expect
else
    set(h_img,'CData', result);
end
set(ud.h.Apply, 'Enable', 'off');

% -------------------------------------------------------------------------------------------
function ImageReset(ImageOpsFig)
if (nargin < 1),    ImageOpsFig = gcbf;     end
ud = get(ImageOpsFig, 'UserData');
orig = ud.OriginalImage;
h_img = findobj(ud.ParentAxis,'Type','image');
set(h_img,'CData', orig);

% --------------------------------------------------------------------
function yiq = rgb2YIQ(rgb)
%   YIQ = rgb2YIQ(RGB) converts the truecolor image RGB to equivalent YIQ image.
%   No error checking
rgb = double(rgb)/255;      % Image has to be converted to double
T = [1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703].';
[m n k] = size(rgb);
yiq = reshape(reshape(rgb,m*n,k)/T,m,n,k);

% --------------------------------------------------------------------
function rgb = YIQ2rgb(yiq)
%   RGB = YIQ2rgb(YIQ) converts the image YIQ to the equivalent truecolor image RGB.
%   The input image must be of class double. The output is of class uint8.
%   No error checking

T = [1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703];
[m,n,k] = size(yiq);
rgb = reshape(yiq(:),m*n,3)*T';

% Make sure the rgb values are between 0.0 and 1.0
rgb = max(0,rgb);
d = find(any(rgb'>1));
rgb(d,:) = rgb(d,:)./(max(rgb(d,:)')'*ones(1,3));

rgb = reshape(rgb,m,n,3);
rgb = uint8(round(rgb*255));      % Put it back to uint8 (Uff)
