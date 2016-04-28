function out = zoom_j(varargin)
%   This is a modified version of zoom that checks the axes label format to see
%   if it is in one of geographical notation (e.g. dd mm ss). If it is, new labels
%   resulting from zoom will keep the current format. It does that by inquiring
%   gcf's appdata for 'LabelFormatType' and than calls ChangeAxesLabels with the apropriate args
%
%   For help type help zoom, for this modified version behaves otherwise like the original

%   Copyright 1984-2002 The MathWorks, Inc.

%   Note: zoom uses the figure buttondown and buttonmotion functions
%
%   ZOOM XON zooms x-axis only
%   ZOOM YON zooms y-axis only
%
%   MODIFIED 20-8-04 The calls to draw_funs where replaced by calls to ChangeAxesLabels
%   (it's the same function as in draw_funs). I had to do it to be able to compile m_gmt
%   without having to compile the hole mirone. That happened via a call to zoom by
%   "ecran" (which is used either by m_gmt and some moduls it calls) and than zoom
%   called draw_funs, that also calls mirone and ...
%
%   The pointer now turns into a magnifier when zoom is 'on'    (4-1-06)

%
% PARSE ARGS - set fig and zoomCmd
%

funHand = [];		% It may evantually come to contain a function handle
switch nargin
	case 0      % no arg in
        zoomCmd = 'toggle';
        fig = get(0,'currentfigure');
        if isempty(fig),    return;    end
	case 1      % one arg in
        % If argument is string, it is a zoom command (i.e. (on, off, down, xdown, etc.).
        if ischar(varargin{1})
			zoomCmd=varargin{1};
        else    % Otherwise, the argument is assumed to be a zoom factor.
			scale_factor=varargin{1};
			zoomCmd='scale';
        end
        fig = get(0,'currentfigure');
        if isempty(fig),    return;    end
	case 2      % two arg in
        fig = varargin{1};
        if (~ishandle(fig)),		error('First argument must be a figure handle.'),	end
        if (isnumeric(varargin{2}))
			scale_factor = varargin{2};
			zoomCmd = 'scale';
        else
			zoomCmd = varargin{2};
        end
	case 3      % three arg in. Third, if empty, center zoom on current point, otherwise -- function handle world
		fig = varargin{1};
		if ( isempty(varargin{3}) )				% Do a zoom centered on the current point or the image center
			scale_factor = varargin{2};
			zoomCmd = 'ptscale';
			anchor_pt = [];
		elseif ( numel(varargin{3}) == 2 && isnumeric(varargin{3}(1)) )	% Center on the varargin{3} point
			scale_factor = varargin{2};
			zoomCmd = 'ptscale';
			anchor_pt = varargin{3};
		else
			zoomCmd = varargin{2};
			funHand = varargin{3}{1};			% Tchan tchan -- a function handle
			if ( numel(varargin{3}) > 1 )
				funHvarargins = varargin{3}(2:end);
			else
				funHvarargins = [];
			end
		end
	otherwise   % too many args
        error(nargchk(0, 3, nargin));
end

%------------------------------------------------------------------------------%
zoomCmd = lower(zoomCmd);

if strcmp(zoomCmd,'off')			% turn zoom off
    doZoomOff(fig);
    %scribefiglisten_j(fig,'off');
    state = getappdata(fig,'ZOOMFigureState');
    if ~isempty(state)
        % since we didn't set the pointer, make sure it does not get reset
        % ptr = get(fig,'pointer');     % NOW WE DO (JL)
        % restore figure and non-uicontrol children don't restore uicontrols
        % because they were restored already when zoom was turned on
        uirestore_j(state,'nouicontrols');
        %set(fig,'pointer',ptr)
        if isappdata(fig,'ZOOMFigureState')
            rmappdata(fig,'ZOOMFigureState');
        end
        % get rid of on state appdata if it exists the non-existance
        % of this appdata indicates that zoom is off.
        if isappdata(fig,'ZoomOnState'),    rmappdata(fig,'ZoomOnState');   end
    end
	if (~isempty(funHand))					% If we had a function handle saved, kill it now
		hz = get(get(fig,'currentaxes'),'ZLabel');
		rmappdata(hz,'ExtFunHand')
		rmappdata(hz,'ExtFunHvararg')
	end
    return      % done, go home.
end

%---------------------------------------------------------------------------------%
% set some things we need for other zoomCommands
%
ax = get(fig,'currentaxes');
rbbox_mode = 0;
% initialize unconstrained state
zoomx = 1; zoomy = 1;
% catch 3d zoom
if ~isempty(ax) && any(get(ax,'view')~=[0 90]) ...
        && ~(strcmp(zoomCmd,'scale') || strcmp(zoomCmd,'fill'))
    fZoom3d = 1;
else
    fZoom3d = 0;
end

%----------------------------------------------------------------------------------%
% the zoomCmd is 'toggle'
if strcmp(zoomCmd,'toggle'),
	state = getappdata(fig,'ZOOMFigureState');
	if isempty(state),	zoom_j(fig,'on');
	else				zoom_j(fig,'off');
	end
	return
end

%----------------------------------------------------------------------------------%
% Set/check a few more things before switching to one of remaining zoom actions
%
% Set zoomx,zoomy and zoomCmd for constrained zooms
if strcmp(zoomCmd,'xdown'),
    zoomy = 0; zoomCmd = 'down'; % Constrain y
elseif strcmp(zoomCmd,'ydown')
    zoomx = 0; zoomCmd = 'down'; % Constrain x
end

% Catch bad argin/argout match
if (nargout ~= 0) && ~isequal(zoomCmd,'getmode') && ...
        ~isequal(zoomCmd,'getlimits') && ~isequal(zoomCmd,'getconnect')
    error('ZOOM only returns an output if the command is getmode, getlimits, or getconnect');
end

%----------------------------------------------------------------------------------%
% Switch for rest of zoomCommands
%
switch zoomCmd
case 'down'
    % Activate axis that is clicked in
    allAxes = findobj(fig,'-depth',1,'type','axes');
    ZOOM_found = 0;

    % this test may be causing failures for 3d axes
    for i=1:numel(allAxes)
        ax=allAxes(i);
        ZOOM_Pt1 = get(ax,'CurrentPoint');
        xlim = get(ax,'xlim');
        ylim = get(ax,'ylim');
        if (xlim(1) <= ZOOM_Pt1(1,1) && ZOOM_Pt1(1,1) <= xlim(2) && ...
                ylim(1) <= ZOOM_Pt1(1,2) && ZOOM_Pt1(1,2) <= ylim(2))
            ZOOM_found = 1;
            set(fig,'currentaxes',ax);
            break
        end
    end

    if (ZOOM_found == 0),   return;    end

    % Check for selection type
    selection_type = get(fig,'SelectionType');
    zoomMode = getappdata(fig,'ZOOMFigureMode');

    axz = get(ax,'ZLabel');

	if fZoom3d
		viewData = getappdata(axz,'ZOOMAxesView');
		if isempty(viewData)
			viewProps = {'CameraTarget' 'CameraTargetMode' 'CameraViewAngle' 'CameraViewAngleMode'};
			setappdata(axz,'ZOOMAxesViewProps', viewProps);
			setappdata(axz,'ZOOMAxesView', get(ax,viewProps));
		end
		if isempty(zoomMode) || strcmp(zoomMode,'in');
			zoomLeftFactor = 1.5;
			zoomRightFactor = .75;
		elseif strcmp(zoomMode,'out');
			zoomLeftFactor = .75;
			zoomRightFactor = 1.5;
		end
		switch selection_type
			case 'open'
				set(ax,getappdata(axz,'ZOOMAxesViewProps'), getappdata(axz,'ZOOMAxesView'));
			case 'normal'
				newTarget = mean(get(ax,'CurrentPoint'),1);
				set(ax,'CameraTarget',newTarget);
				camzoom(ax,zoomLeftFactor);
			otherwise
				newTarget = mean(get(ax,'CurrentPoint'),1);
				set(ax,'CameraTarget',newTarget);
				camzoom(ax,zoomRightFactor);
		end
		return
	end

	if isempty(zoomMode) || strcmp(zoomMode,'in');
		switch selection_type
			case 'normal'           % Zoom in
				m = 1;
				scale_factor = 2;   % the default zooming factor
			case 'open'             % Zoom all the way out
				zoom_j(fig,'out');
				return;
			otherwise               % Zoom partially out
				m = -1;
				scale_factor = 2;
		end
	elseif strcmp(zoomMode,'out')
		switch selection_type
			case 'normal'           % Zoom partially out
				m = -1;
				scale_factor = 2;
			case 'open'             % Zoom all the way out
				zoom_j(fig,'out');
				return;
			otherwise               % Zoom in
				m = 1;
				scale_factor = 2; % the default zooming factor
		end
	else % unrecognized zoomMode
		return
	end

    ZOOM_Pt1 = get_currentpoint(ax);
    ZOOM_Pt2 = ZOOM_Pt1;
    center = ZOOM_Pt1;

    if (m == 1)        % Zoom in
        units = get(fig,'units'); set(fig,'units','pixels')
        rbbox([get(fig,'currentpoint') 0 0],get(fig,'currentpoint'));
        ZOOM_Pt2 = get_currentpoint(ax);
        set(fig,'units',units)

        % Note the currentpoint is set by having a non-trivial up function.
        if (min(abs(ZOOM_Pt1-ZOOM_Pt2)) >= min(.01*[diff(get_xlim(ax)) diff(get_ylim(ax))]))
            % determine axis from rbbox
            a = [ZOOM_Pt1;ZOOM_Pt2]; a = [min(a);max(a)];

            % Failed attempt to do a ruberband isotropic (but it nearly works)
            limits = getappdata(get(ax,'ZLabel'),'ZOOMAxesData');
            if isempty(limits)        % Use quick method if possible
                xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2);
                ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2);
            else
                xmin = limits(1); xmax = limits(2);
                ymin = limits(3); ymax = limits(4);
            end
            
            a_dif = diff(a);    % remember diff operates by rows
            rato = a_dif(1) / a_dif(2);
            a_med = sum(a) / 2;
            if (rato > 1)       % rectangle wider then heigh
                a = [a(1,1) a_med(2)-a_dif(1)/2; a(2,1) a_med(2)+a_dif(1)/2];
            elseif (rato < 1)   % the oposite
                a = [a_med(1)-a_dif(2)/2 a(1,2); a_med(1)+a_dif(2)/2 a(2,2)];
            end
            
            % Now we must check if the new 'a' does not exceed 'limits'. If it did the square
            % will become a rectangle, and we must square it again
            a_med = sum(a) / 2;     a_dif = diff(a);
            if (a(1,1) < xmin)          % Left side is out of limits
                a(1,1) = xmin;      a_med = sum(a) / 2;     a_dif = diff(a);
                a = [xmin a_med(2)-a_dif(1)/2; a(2,1) a_med(2)+a_dif(1)/2];
            elseif (a(2,1) > xmax)      % Right side is out of limits
                a(2,1) = xmax;      a_med = sum(a) / 2;     a_dif = diff(a);
                a = [a(1,1) a_med(2)-a_dif(1)/2; xmax a_med(2)+a_dif(1)/2];                
            elseif (a(1,2) < ymin)      % yMin side is out. Min can either be top or bot depending on Ydir
                a(1,2) = ymin;      a_med = sum(a) / 2;     a_dif = diff(a);
                a = [a_med(1)-a_dif(2)/2 ymin; a_med(1)+a_dif(2)/2 a(2,2)];                
            elseif (a(2,2) > ymax)      % yMax side is out
                a(2,2) = ymax;      a_med = sum(a) / 2;     a_dif = diff(a);
                a = [a_med(1)-a_dif(2)/2 a(1,2); a_med(1)+a_dif(2)/2 ymax];
            end
            % And finaly, compute the scale_factor
            a_dif = diff(a);
            scale_factor = min((xmax - xmin) / a_dif(1), (ymax - ymin) / a_dif(2));
            
            % Undo the effect of get_currentpoint for log axes
            if strcmp(get(ax,'XScale'),'log'),
                a(1:2) = 10.^a(1:2);
            end
            if strcmp(get(ax,'YScale'),'log'),
                a(3:4) = 10.^a(3:4);
            end
            rbbox_mode = 1;
        end
    end
    limits = zoom_j(fig,'getlimits');

case 'scale',
    if (all(get(ax,'view')==[0 90]))				% 2D zooming with scale_factor
        % Activate axis that is clicked in
        ZOOM_found = 0;
        ax = gca;
        xlim = get(ax,'xlim');
        ylim = get(ax,'ylim');
        ZOOM_Pt1 = [sum(xlim)/2 sum(ylim)/2];
        center = ZOOM_Pt1;

        if (xlim(1) <= ZOOM_Pt1(1,1) && ZOOM_Pt1(1,1) <= xlim(2) && ...
                ylim(1) <= ZOOM_Pt1(1,2) && ZOOM_Pt1(1,2) <= ylim(2))
            ZOOM_found = 1;
        end

        if ZOOM_found==0, return, end

        if (scale_factor >= 1),     m = 1;
        else                        m = -1;
        end
    else % 3D
        old_CameraViewAngle = get(ax,'CameraViewAngle')*pi/360;
        ncva = atan(tan(old_CameraViewAngle)*(1/scale_factor))*360/pi;
        set(ax,'CameraViewAngle',ncva);
        return
    end
    limits = zoom_j(fig,'getlimits');
case 'getmode'
    state = getappdata(fig,'ZOOMFigureState');
    if isempty(state)
        out = 'off';
    else
        mode = getappdata(fig,'ZOOMFigureMode');
        if isempty(mode),   out = 'on';
        else                out = mode;
        end
    end
    return
case 'on'
    state = getappdata(fig,'ZOOMFigureState');
    if isempty(state),
        % turn off all other interactive modes
%        state = uiclearmode_j(fig,'docontext','zoom_j',fig,'off');
		state = uisuspend_j(fig, true);
        % restore button down functions for uicontrol children of the figure
        uirestore_j(state,'uicontrols');
        setappdata(fig,'ZOOMFigureState',state);
    end
    
    % Here we change the pointer, but only after uiclearmode_j that calls uisuspend
    % This way uirestore will reset the pointer to its previous state
    lupa = getPointer;
    set(fig,'Pointer', 'custom','PointerShapeCData',lupa)

    set(fig,'windowbuttondownfcn','zoom_j(gcbf,''down'')', ...
        'windowbuttonupfcn','ones;', ...
        'windowbuttonmotionfcn','', ...
        'buttondownfcn','', ...
        'interruptible','on');
    set(ax,'interruptible','on');
    % set an appdata so it will always be possible to determine
    % whether or not zoom is on and in what type of 'on' state
    % it is. This appdata will not exist when zoom is off
    setappdata(fig,'ZoomOnState','on');
    %scribefiglisten_j(fig,'on');
	setappdata(fig,'ZOOMFigureMode','in');
	if ( ~isempty(funHand) )					% If we have a function handle to execute later, save it now
		hz = get(ax,'ZLabel');
		setappdata(hz,'ExtFunHand', funHand)
		setappdata(hz,'ExtFunHvararg', funHvarargins)
	end
    return
case 'inmode'
    zoom_j(fig,'on');
	setappdata(fig,'ZOOMFigureMode','in');
    return
case 'outmode'
    zoom_j(fig,'on');
	setappdata(fig,'ZOOMFigureMode','out');
    return
case 'reset',
    axz = get(ax,'ZLabel');
    if isappdata(axz,'ZOOMAxesData')
        rmappdata(axz,'ZOOMAxesData');
    end
    return
case 'xon',
    zoom_j(fig,'on') % Set up userprop
    set(fig,'windowbuttondownfcn','zoom_j(gcbf,''xdown'')', 'windowbuttonupfcn','ones;', ...
        'windowbuttonmotionfcn','','buttondownfcn','', 'interruptible','on');
    set(ax,'interruptible','on')
    % set an appdata so it will always be possible to determine whether or not zoom is on
    % and in what type of 'on' state it is. This appdata will not exist when zoom is off
    setappdata(fig,'ZoomOnState','xon');
    return
case 'yon',
    zoom_j(fig,'on') % Set up userprop
    set(fig,'windowbuttondownfcn','zoom_j(gcbf,''ydown'')', 'windowbuttonupfcn','ones;', ...
        'windowbuttonmotionfcn','','buttondownfcn','', 'interruptible','on');
    set(ax,'interruptible','on')
    % set an appdata so it will always be possible to determine whether or not zoom is on
    % and in what type of 'on' state it is. This appdata will not exist when zoom is off
    setappdata(fig,'ZoomOnState','yon');
    return
case 'out',
    limits = zoom_j(fig,'getlimits');
    center = [sum(get_xlim(ax))/2 sum(get_ylim(ax))/2];
    m = -inf; % Zoom totally out
case 'getlimits', % Get axis limits
    axz = get(ax,'ZLabel');
    limits = getappdata(axz,'ZOOMAxesData');
    % Do simple checking of userdata
    if (size(limits,2) == 4 && size(limits,1) <= 2)
        if all(limits(1,[1 3]) < limits(1,[2 4])),
            out = limits(1,:);
            return   % Quick return
        else
            getlimits = -1; % Don't munge data
        end
	else
		if isempty(limits),	getlimits = 1;
		else				getlimits = -1;
		end
    end
    % If I've made it to here, we need to compute appropriate axis limits.
    if isempty(getappdata(axz,'ZOOMAxesData')),
        % Use quick method if possible
        xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2);
        ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2);
    elseif strcmp(get(ax,'xLimMode'),'auto') && strcmp(get(ax,'yLimMode'),'auto'),
        % Use automatic limits if possible
        xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2);
        ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2);
    else
        % Use slow method only if someone else is using the userdata
        h = get(ax,'Children');
        xmin = inf; xmax = -inf; ymin = inf; ymax = -inf;
        for i=1:length(h),
            t = get(h(i),'Type');
            if ~strcmp(t,'text'),
                if strcmp(t,'image'), % Determine axis limits for image
                    x = get(h(i),'Xdata'); y = get(h(i),'Ydata');
                    x = [min(min(x)) max(max(x))];
                    y = [min(min(y)) max(max(y))];
                    [ma,na] = size(get(h(i),'Cdata'));
					if (na > 1),	dx = diff(x)/(na-1);
					else			dx = 1;
					end
					if (ma > 1),	dy = diff(y)/(ma-1);
					else			dy = 1;
					end
					x = x + [-dx dx]/2; y = y + [-dy dy]/2;
                end
                xmin = min(xmin,min(min(x)));
                xmax = max(xmax,max(max(x)));
                ymin = min(ymin,min(min(y)));
                ymax = max(ymax,max(max(y)));
            end
        end
        % Use automatic limits if in use (override previous calculation)
        if strcmp(get(ax,'xLimMode'),'auto'),
            xlim = get_xlim(ax); xmin = xlim(1); xmax = xlim(2);
        end
        if strcmp(get(ax,'yLimMode'),'auto'),
            ylim = get_ylim(ax); ymin = ylim(1); ymax = ylim(2);
        end
    end
    limits = [xmin xmax ymin ymax];
    if (getlimits ~= -1)   % Don't munge existing data.
        % Store limits ZOOMAxesData store it with the ZLabel, so that it's cleared if
        % the user plots again into this axis.  If that happens, this state is cleared
        axz = get(ax,'ZLabel');
        setappdata(axz,'ZOOMAxesData',limits);
    end
    out = limits;
    return
case 'getconnect', % Get connected axes
    axz = get(ax,'ZLabel');
    limits = getappdata(axz,'ZOOMAxesData');
    if all(size(limits)==[2 4]), % Do simple checking
        out = limits(2,[1 2]);
    else
        out = [ax ax];
    end
    return
case 'fill',
    old_view = get(ax,'view');
    view(45,45);
    set(ax,'CameraViewAngleMode','auto');
    set(ax,'CameraViewAngle',get(ax,'CameraViewAngle'));
    view(old_view);
    return

case 'ptscale'
    xlim = get(ax,'xlim');
    ylim = get(ax,'ylim');
	if (isempty(anchor_pt))
	    ZOOM_Pt1 = get_currentpoint(ax);
	else
		ZOOM_Pt1 = anchor_pt(:)';
	end
    center = ZOOM_Pt1;

    if (xlim(1) <= ZOOM_Pt1(1,1) && ZOOM_Pt1(1,1) <= xlim(2) && ...
            ylim(1) <= ZOOM_Pt1(1,2) && ZOOM_Pt1(1,2) <= ylim(2))
    else
        ZOOM_Pt1 = [sum(xlim)/2 sum(ylim)/2];
        center = ZOOM_Pt1;
    end

    if (scale_factor >= 1),		m = 1;
    else						m = -1;
    end
    limits = zoom_j(fig,'getlimits');
otherwise
    error(['Unknown option: ',zoomCmd,'.']);
end

%---------------- Actual zoom operation -------------------------------------------%

if ~rbbox_mode,
    xmin = limits(1); xmax = limits(2);
    ymin = limits(3); ymax = limits(4);
    if (m == (-inf))
        dx = xmax-xmin;
        dy = ymax-ymin;
    else
        dx = diff(get_xlim(ax))*(scale_factor.^(-m-1)); dx = min(dx,xmax-xmin);
        dy = diff(get_ylim(ax))*(scale_factor.^(-m-1)); dy = min(dy,ymax-ymin);
    end
    % Limit zoom.
    center = max(center,[xmin ymin] + [dx dy]);
    center = min(center,[xmax ymax] - [dx dy]);
    a = [max(xmin,center(1)-dx) min(xmax,center(1)+dx) ...
            max(ymin,center(2)-dy) min(ymax,center(2)+dy)];

    % Check for log axes and return to linear values.
    if strcmp(get(ax,'XScale'),'log'),
        a(1:2) = 10.^a(1:2);
    end
    if strcmp(get(ax,'YScale'),'log'),
        a(3:4) = 10.^a(3:4);
    end
end

% Check for axis equal and update a as necessary
if strcmp(get(ax,'plotboxaspectratiomode'),'manual') && strcmp(get(ax,'dataaspectratiomode'),'manual')
    ratio = get(ax,'plotboxaspectratio') ./ get(ax,'dataaspectratio');
    dx = a(2)-a(1);
    dy = a(4)-a(3);
    [kmax,k] = max([dx dy]./ratio(1:2));
    if (k == 1)
        dy = kmax*ratio(2);
        a(3:4) = mean(a(3:4))+[-dy dy]/2;
    else
        dx = kmax*ratio(1);
        a(1:2) = mean(a(1:2))+[-dx dx]/2;
    end
end

% Update circular list of connected axes
list = zoom_j(fig,'getconnect');		% Circular list of connected axes.
if (zoomx)
    if (a(1) == a(2)),      return,     end % Short circuit if zoom is moot.
    set(ax,'xlim',a(1:2))
    % Save current X labels in appdata for easear access when we want to change them in ChangeAxesLabels
    labelType = getappdata(ax,'LabelFormatType');
    if isempty(labelType)       % Prevent an error in the next switch when labelType = []
        labelType = ' ';
    else                        % Set it to 'auto' (that is numeric) so ChangeAxesLabels knows how to reformat to a geog string
        set(ax, 'XTickLabelMode', 'auto')      % It was implicitly set to 'manual' by ChangeAxesLabels
    end
    setappdata(ax,'XTickOrig',get(ax,'XTickLabel'));
    switch labelType
        case {'DegDec' 'DegMin' 'DegMinDec' 'DegMinSec' 'DegMinSecDec'}
            ChangeAxesLabels(fig,ax,labelType,'X')
    end
    h = list(1);
    while h ~= ax,
        set(h,'xlim',a(1:2))
        % Get next axes in the list
        hz = get(h,'ZLabel');
        next = getappdata(hz,'ZOOMAxesData');
        if all(size(next)==[2 4]), h = next(2,1); else h = ax; end
    end
end
if (zoomy)
    if (a(3) == a(4)),  return,     end % Short circuit if zoom is moot.
    set(ax,'ylim',a(3:4))

    % Save current Y labels in appdata for easear access when we want to change them in ChangeAxesLabels
    labelType = getappdata(ax,'LabelFormatType');
    if isempty(labelType)       % Prevent an error in the next switch when labelType = []
        labelType = ' ';
    else                        % Set it to 'auto' (that is numeric) so ChangeAxesLabels knows how to refrmat to a geog string
        set(ax, 'YTickLabelMode', 'auto')      % It was implicitly set to 'manual' by ChangeAxesLabels
    end
    setappdata(ax,'YTickOrig',get(ax,'YTickLabel'));
    switch labelType
        case {'DegDec' 'DegMin' 'DegMinDec' 'DegMinSec' 'DegMinSecDec'}
            ChangeAxesLabels(fig,ax,labelType,'Y')
    end
    h = list(2);
    while (h ~= ax)
        set(h,'ylim',a(3:4))
        % Get next axes in the list
        hz = get(h,'ZLabel');
        next = getappdata(hz,'ZOOMAxesData');
        if all(size(next)==[2 4]), h = next(2,2); else h = ax; end
    end
end

	% -------- My private add-ons to the zoom function ------------- JL
	% Keep track of the current (accumulated) zooming factor. To be eventualy used by
	% other functions that may be interested in this info.
	% NOTE: The rbbox_mode is not foreseen here
	hz = get(ax,'ZLabel');
	oldZoomFactor = getappdata(hz,'AxesScaleFactor');
	if (~isempty(oldZoomFactor))
        if (m ~= -inf),     oldZoomFactor = oldZoomFactor * scale_factor;
        else                oldZoomFactor = 1;
        end
        setappdata(hz,'AxesScaleFactor',oldZoomFactor);
	else
        % No need to test for m because it should only pass here on the first call to zoom
        setappdata(hz,'AxesScaleFactor',scale_factor);
	end

    % Update amplification info displayed in figure name
    fname = get(fig,'Name');    ind = strfind(fname,'@ ');
    if (~isempty(ind))
		hImg = findobj(ax,'Type','Image');
		imageWidth  = size(get(hImg, 'CData'), 2);
		imageHeight = size(get(hImg, 'CData'), 1);
		axUnits = get(ax, 'Units');     set(ax, 'Units', 'pixels');
		axPos = get(ax,'Pos');          set(ax, 'Units', axUnits);
		imgMagRatio =  ((axPos(3) / imageWidth + axPos(4) / imageHeight) * 0.5 * 1);
        lims = getappdata(ax,'ThisImageLims');
		if (~isempty(lims))
			magRatio =  round( ( diff(lims(1:2)) / diff(a(1:2)) + diff(lims(3:4)) / diff(a(3:4)) ) * 0.5 * imgMagRatio * 100);
			fname = [fname(1:ind) sprintf('  %d%%',magRatio)];
			set(fig,'Name',fname)
		end
    end

	% Update slider value (if the figure has one (or two))
	if ( zoomx || zoomy ),		imscroll_j(ax,'ZoomSetSliders');	end

	% If we have a function handle, execute it
	funHand = getappdata(hz,'ExtFunHand');
	if ( isa(funHand, 'function_handle') && (zoomx || zoomy) )
		funHvarargins = getappdata(hz,'ExtFunHvararg');
		if ( ~isempty(funHvarargins) )				% If the function wants inputs
			feval(funHand, funHvarargins{:});
		else										% Function doesn't want inputs
			feval(funHand);
		end
	end

%----------------------------------------------------------------------------------%
% Helper Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = get_currentpoint(ax)
%GET_CURRENTPOINT Return equivalent linear scale current point
	p = get(ax,'currentpoint'); p = p(1,1:2);
	if strcmp(get(ax,'XScale'),'log'),
        p(1) = log10(p(1));
	end
	if strcmp(get(ax,'YScale'),'log'),
        p(2) = log10(p(2));
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xlim = get_xlim(ax)
%GET_XLIM Return equivalent linear scale xlim
	xlim = get(ax,'xlim');
	if strcmp(get(ax,'XScale'),'log'),
        xlim = log10(xlim);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ylim = get_ylim(ax)
%GET_YLIM Return equivalent linear scale ylim
	ylim = get(ax,'ylim');
	if strcmp(get(ax,'YScale'),'log'),
        ylim = log10(ylim);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doZoomOff(fig)
	if ~isempty(getappdata(fig,'ZoomFigureMode'))
        rmappdata(fig,'ZOOMFigureMode');
	end

% -------------------------------------------------------------------------------------------------------
function lupa = getPointer()
lupa = [...
   NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
   NaN NaN NaN NaN   1   1   1   1 NaN NaN NaN NaN NaN NaN NaN NaN
   NaN NaN   1   1 NaN NaN NaN NaN   1   1 NaN NaN NaN NaN NaN NaN
   NaN   1 NaN NaN NaN NaN NaN NaN NaN NaN   1 NaN NaN NaN NaN NaN
   NaN   1 NaN NaN NaN   1   1 NaN NaN NaN   1 NaN NaN NaN NaN NaN
     1 NaN NaN NaN NaN   1   1 NaN NaN NaN NaN   1 NaN NaN NaN NaN
     1 NaN NaN   1   1   1   1   1   1 NaN NaN   1 NaN NaN NaN NaN
     1 NaN NaN   1   1   1   1   1   1 NaN NaN   1 NaN NaN NaN NaN
     1 NaN NaN NaN NaN   1   1 NaN NaN NaN NaN   1 NaN NaN NaN NaN
   NaN   1 NaN NaN NaN   1   1 NaN NaN NaN   1 NaN NaN NaN NaN NaN
   NaN   1 NaN NaN NaN NaN NaN NaN NaN NaN   1 NaN NaN NaN NaN NaN
   NaN NaN   1   1 NaN NaN NaN NaN   1   1   1   1 NaN NaN NaN NaN
   NaN NaN NaN NaN   1   1   1   1 NaN NaN   1   1   1 NaN NaN NaN
   NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1 NaN NaN
   NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1 NaN
   NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   1 NaN NaN];

% -------------------------------------------------------------------------------------------------------
function ChangeAxesLabels(hFig,hAx,type,eixo)
% This function formats the axes labels strings using a geographical notation
% EIXO is a single character (X or Y) that indicates on which axe (abssissa or
% ordenate) this function will operate with.

tick = getappdata(hAx,[eixo,'TickOrig']);
if (isempty(tick)),     return;     end     % Not a Mirone (type) figure
n_tick = size(tick,1);
sep = ':';

switch type
    case 'DegDec'
        setappdata(hAx,'LabelFormatType','DegDec')                      % Save it so zoom can know the label type
    case 'DegMin'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMM',0,'str');   % e_str is a structure with string fields
        str_e = [e_str.dd repmat(sep, n_tick,1) e_str.mm];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMin')                     % Save it so zoom can know the label type
    case 'DegMinDec'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMM.x',2,'str'); % e_str is a structure with string fields
        str_e = [e_str.dd repmat(sep, n_tick,1) e_str.mm];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMinDec')                  % Save it so zoom can know the label type
    case 'DegMinSec'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMMSS',0,'str'); % e_str is a structure with string fields
        str_e = [e_str.dd repmat(sep, n_tick,1) e_str.mm repmat(sep, n_tick,1) e_str.ss];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMinSec')                  % Save it so zoom can know the label type
    case 'DegMinSecDec'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMMSS.x',1,'str');   % e_str is a structure with string fields
        str_e = [e_str.dd repmat(sep, n_tick,1) e_str.mm repmat(sep, n_tick,1) e_str.ss];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMinSecDec')               % Save it so zoom can know the label type
end
