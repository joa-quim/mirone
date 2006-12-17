function img = imcapture( h, opt, dpi )
% IMCAPTURE do screen captures at controllable resolution using the undocumented "hardcopy" built-in function.
%
% USAGE:
%   IMG = IMCAPTURE(H) gets a screen capture from object H, where H is a handle to a figure, axis or an image
%   IMG = IMCAPTURE or IMG = IMCAPTURE([]) operates on the current figure.
%   IMG = IMCAPTURE(H, OPT) selects one of three possible operations depending on the OPT string
%           OPT == 'img' returns an image corresponding to the axes contents and a size determined
%                   by the figure size. If DPI is not provided, the screen capture is done at 150 dpi
%           OPT == 'imgAx' returns an image corresponding to the displyed image plus its labels.
%                   This option tries to fill a enormous gap in Matlab base which is the non availability of
%                   an easy way of geting a raster image that respects the data original aspect ratio.
%           OPT == 'all' returns an image with all the figure's contents and size determined by the figure size.
%                   If DPI is not provided, the screen capture is done at 150 dpi
%   IMG = IMCAPTURE(H,OPT,DPI) do the screen capture at DPI resolution. DPI can be either a string or numeric
%   IMG = IMCAPTURE(H,'img',0) returns an image that has exactly the same width and height of the original CData
%                   You can use this, for example, to convert an indexed image to RGB.
%                   Note, this should be equivalent to a call to GETFRAME, but maybe it fails less in
%                   returning EXACTLY an image of the same size of CData (R13 ,as usual, is better than R14).
%   IMG = IMCAPTURE(H,'img',[MROWS NCOLS]) returns an image of the size specified by [mrows ncols]
%                   This a 10 times faster option to the use of imresize(IMG,[mrows ncols],method)
%                   when method is either 'bilinear' or 'bicubic'.
%
%   NOTE. This function also works with the 'imgAx' option for plots and surfaces but DPI settings
%   in more "fluid" because I cannot exactly determine what is the image.
%   For plots, when DPI is not set it defaults to '0', which probably the best choice.
%
%   EXAMPLE:
%       load clown
%       h = imagesc(X);
%       set(gcf,'colormap',map)
%       x = 27*cos(-pi:.1:pi) + 55;
%       y = 27*sin(-pi:.1:pi) + 120;
%       patch('XData',x,'YData',y,'FaceColor','y');     % Let's change the clown nose color
%       % Notice how the yellow circle as turned into an ellipse.
%       % Now do an image screen capture at 300 dpi
%       rgb=imcapture(h,'img',300);
%       % In order to confirm the correct aspect ratio you need either to call axis image
%       figure; image(rgb); axis image
%       % or save into an raster format an open it with outside Matlab. E.G
%       imwrite(rgb,'lixo.jpg')

%   HOW DOES IT WORKS?
%       A lot of testings revealed that the array size returned by 'hardcopy()' function
%       depends on the 'PaperPosition' property.
%       BTW the help hardcopy says
%           "HARDCOPY(H,'filename','format') saves the figure window with handle H
%            to the designated filename in the specified format."
%       and
%           "Do NOT use this function directly. Use PRINT instead."
%
%       Well, fortunately the first is false and I didn't follow the advise.
%       But let us return to the 'PaperPosition' property. By default figures are created
%       (at least in my system) with this values in centimeters
%       get(hFig,'PaperPosition') = [0.6345    6.3452   20.3046   15.2284]
%       Using the example above (the clown) X = 200x320
%           rgb= hardcopy(get(h,'parent'),'lixo.jpg','-dzbuffer');
%       shows that rgb = 900x1201x3. And how did hardcopy picked this dimensions?
%       We can find that X & Y resolution are different
%           DPIx = round(320 / 20.3046 * 2.54) = 40
%           DPIy = round(200 / 15.2284 * 2.54) = 33
%       Shit, we have different resolutions - what are these guys doing?
%       But see this
%           DPIx = round(1201 / 20.3046 * 2.54) = 150
%       So the capture is carried out at 150 dpi, and using this on the DPIy expression
%           N = 15.22842 * 150 / 2.54 = 899
%       We nearly got it. 899 instead of 900. This may be due to inch<->centimeters rounding
%       errors or to a bug. You can see in the imgOnly function code that I add 0.5 to
%       rows number in order get correct results (sometimes it missed by one).
%       OK, now that it is found how the size is picked, the trick is to play arround with
%       the 'PaperPosition' in order to get what we want.
%
%   AUTHOR
%       Joaquim Luis (jluis@ualg.pt)    12-Dez-2006
%       University of Algarve
%
%       Let it work also with plots     13-Dez-2006

    if (nargin == 0 || isempty(h))
        h = get(0,'CurrentFigure');
    elseif (strcmp(get(h,'Type'),'figure'))
    elseif (strcmp(get(h,'Type'),'axes'))
        h = get(h,'Parent');
    elseif (strcmp(get(h,'Type'),'image'))
        h = get(get(h,'Parent'),'Parent');
    else
        h = [];
    end
    if (~ishandle(h))
        error('First argument is not a valid Fig/Axes/Image handle')
    end
    if (nargin == 1),   opt = [];   end
    
    inputargs{1} = h;
    inputargs{2} = 'lixo.jpg';          % The name doesn't really matter, but we need one.
    renderer = get(h, 'Renderer');
    if strcmp( renderer, 'painters' )
        renderer = 'zbuffer';
    end
    inputargs{3} = ['-d' renderer];
    
    if (nargin == 3)
        if (numel(dpi) == 1)
            dpi = num2str(dpi);
            inputargs{4} = ['-r' dpi];
        elseif (numel(dpi) == 2)            % New image size in [mrows ncols]
            inputargs{4} = dpi;
        else
            error('ERROR: third argument must be a ONE or TWO elements vector')
        end
    end

    msg = [];
	hAxes = findobj(h,'Type','axes');
    if (numel(hAxes) == 1 && strcmp( get(hAxes,'Visible'),'off') && (nargin == 1 || isempty(opt)) )
        % Default to 'imgOnly' when we have only one image with invisible axes
        opt = 'img';
    end

    if (nargin > 1)
        switch opt
            case 'img'      % Capture only the image
                [img,msg] = imgOnly(h,[],inputargs{:});
            case 'imgAx'    % Capture an image including its Labels
                [img,msg] = imgOnly(h,'nik',inputargs{:});
            case 'all'      % Capture everything in Figure
                img = allInFig(h,inputargs{:});
            otherwise
                msg = 'Second argument is not a recognized option';
        end
    else
        img = allInFig(h,inputargs{:});
    end
    
    if (~isempty(msg))      % If we had an error inside imgOnly()
        error(msg);        img = [];
    end

% ------------------------------------------------------------------    
function [img,msg] = imgOnly(h,opt,varargin)
    % Capture the image, and optionaly the frame, mantaining the original image aspect ratio.
    % We do that be messing with the Figure's 'PaperPosition' property
    msg = [];
    hAxes = findobj(h,'Type','axes');
    if (isempty(hAxes) || numel(hAxes) > 1)
        msg = 'ERROR: The figure must contain one, and one ONLY axes';
    end
    im = get(findobj(h,'Type','image'),'CData');
	if (~isempty(im))
        nx = size(im,2);                ny = size(im,1);
    else                    % We have something else. A plot, a surface, etc ...
        axUnit = get(hAxes,'Units');    set(hAxes,'Units','pixels')
        axPos = get(hAxes,'pos');       set(hAxes,'Units',axUnit)
        nx = axPos(3);                  ny = axPos(4);
        if (numel(varargin) == 3)
            varargin{4} = '-r0';        % For the plot cases this is probably the best choice
        end
	end

    PU = get(h,'paperunits');       set(h,'paperunits','inch')
    pp = get(h,'paperposition');    PPM = get(h,'PaperPositionMode');
    dpi = round(nx / pp(3));
    % Here is the kee point of all this manip.
    set(h,'paperposition',[pp(1:3) ny / dpi])
    
    axUnit = get(hAxes,'Units');
    axPos = get(hAxes,'pos');           % Save this because we will have to restore it later
    set(hAxes,'Units','Normalized')     % This is the default, but be sure
    fig_c = get(h,'Color');       set(h,'Color','w')
    
    if (isempty(opt))                   % Pure Image only capture. Even if axes are visible, ignore them
        set(hAxes,'pos',[0 0 1 1],'Visible','off')
    elseif (strcmp(get(hAxes,'Visible'),'on'))      % Try to capture an image that respects the data aspect ratio
		h_Xlabel = get(hAxes,'Xlabel');         h_Ylabel = get(hAxes,'Ylabel');
		units_save = get(h_Xlabel,'units');
		set(h_Xlabel,'units','pixels');         set(h_Ylabel,'units','pixels');
		Xlabel_pos = get(h_Xlabel,'pos');
		Ylabel_pos = get(h_Ylabel,'Extent');
		
		if (abs(Ylabel_pos(1)) < 20)    % Stupid hack, but there is a bug somewhere
            Ylabel_pos(1) = 30;
		end
		
		y_margin = abs(Xlabel_pos(2))+get(h_Xlabel,'Margin');  % To hold the Xlabel height
		x_margin = abs(Ylabel_pos(1))+get(h_Ylabel,'Margin');  % To hold the Ylabel width
		y_margin = min(max(y_margin,20),30);            % Another hack due to the LabelPos non-sense
        
        figUnit = get(h,'Units');        set(h,'Units','pixels')
        figPos = get(h,'pos');           set(h,'Units',figUnit)
        x0 = x_margin / figPos(3);
        y0 = y_margin / figPos(4);
        set(hAxes,'pos',[x0 y0 1-[x0 y0]-1e-2])
        set(h_Xlabel,'units',units_save);     set(h_Ylabel,'units',units_save);
    else            % Dumb choice. Default to Image only
        set(hAxes,'pos',[0 0 1 1],'Visible','off')
    end
    
    confirm = false;
    try
        if (numel(varargin) == 3)
            varargin{4} = '-r150';
        elseif (numel(varargin) == 4 && strcmp(varargin{4},'-r0'))  % One-to-one capture
            varargin{4} = ['-r' num2str(round(dpi))];
            confirm = true;
            mrows = ny;            ncols = nx;      % To use in "confirm"
        elseif (numel(varargin) == 4 && numel(varargin{4}) == 2)    % New size in mrows ncols
            mrows = varargin{4}(1);
            ncols = varargin{4}(2);
            set(h,'paperposition',[pp(1:2) ncols/dpi mrows/dpi])
            varargin{4} = ['-r' num2str(round(dpi))];
            confirm = true;
        end
        img = hardcopy( varargin{:} );      % Capture
        if (confirm)                        % We asked for a pre-determined size. Check that the result is correct
            dy = mrows - size(img,1);       % DX & DY should be zero or one (when it buggs).
            dx = ncols - size(img,2);
            if (dx ~= 0 || dy ~= 0)         % ML failed (probably R14). Repeat to make it obey
                mrows_desBUG = mrows + dy;
                ncols_desBUG = ncols + dx;
                set(h,'paperposition',[pp(1:2) ncols_desBUG/dpi mrows_desBUG/dpi])
                img = hardcopy( varargin{:} );      % Insist
            end
        end
    catch                                   % If it screws, restore original Fig properties anyway
        set(hAxes,'Units',axUnit,'pos',axPos,'Visible','on')
        set(h,'paperposition',pp,'paperunits',PU,'PaperPositionMode',PPM,'Color',fig_c)
        msg = lasterr;      img = [];
    end
    
    % Reset the original fig properties
    set(hAxes,'Units',axUnit,'pos',axPos,'Visible','on')
    set(h,'paperposition',pp,'paperunits',PU,'PaperPositionMode',PPM,'Color',fig_c)
    
% ------------------------------------------------------------------    
function img = allInFig(h,varargin)
    % Get everything in the Figure
    fig_c = get(h,'Color');       set(h,'Color','w')
    if (numel(varargin) == 3)
        varargin{4} = '-r150';
    end
    img = hardcopy( varargin{:} );    
    set(h,'Color',fig_c)
