function [out] = pan(arg1,arg2)
%PAN Interactively pan the view of a plot
%  PAN ON turns on mouse-based panning.
%  PAN XON turns on x-only panning    
%  PAN YON turns on y-only panning
%  PAN OFF turns it off.
%  PAN by itself toggles the state.
%
%  PAN(FIG,...) works on specified figure handle.

% Copyright 2003-2004 The MathWorks, Inc.

% Undocumented syntax
%  PAN(FIG,STYLE); where STYLE = 'x'|'y'|'xy', Note: syntax doesn't turn pan on like 'xon'
%  OUT = PAN(FIG,'getstyle')  'x'|'y'|'xy'
%  OUT = PAN(FIG,'ison')  true/false

% Higly hacked version to work in R13 and only in 2D (that's all that is needed to Mirone)
out = [];

if (nargin == 0) 
    fig = gcf; % caller did not specify handle
    locSetState(fig,'toggle');
elseif nargin==1
    if ishandle(arg1)
        locSetState(arg1,'toggle');
    else
        fig = gcf; % caller did not specify handle
        locSetState(fig,arg1); 
    end 
elseif nargin==2 
    if ~ishandle(arg1), error('Unknown figure.'); end
    switch arg2
        case 'getstyle'   
			pandata = locGetData(arg1);
			out = pandata.style;
        case 'ison'
			pandata = locGetData(arg1);
			out = strcmp(pandata.enable,'on');
        otherwise
           locSetState(arg1,arg2);
    end
end

%-----------------------------------------------%
function locSetState(target,state)
%Enables/disables panning callbacks
% target = figure || axes
% state = 'on' || 'off'

fig = target;
pandata = locGetData(fig);
uistate = pandata.uistate;

if strcmpi(state,'xon')
    state = 'on';
    pandata.style = 'x';
elseif strcmpi(state,'yon')
    state = 'on';
    pandata.style = 'y';
elseif strcmpi(state,'x')
    pandata.style = 'x';
	setappdata(fig,'PanFigureState',pandata);
    return			% All done
elseif strcmpi(state,'y')
    pandata.style = 'y';
    setappdata(fig,'PanFigureState',pandata);
    return			% All done
elseif strcmpi(state,'xy')
    pandata.style = 'xy';
    setappdata(fig,'PanFigureState',pandata);
    return			% All done
elseif strcmpi(state,'toggle')
    if strcmpi(pandata.enable,'on')
        state = 'off';
    else 
        state = 'on';
    end
end

setappdata(fig,'PanFigureState',pandata);
   
if strcmpi(state,'on') 
    % turn off all other interactive modes
    if isempty(uistate)
        % Specify deinstaller callback: pan(fig,'off') 
%        uistate = uiclearmode_j(fig,'pan',target,'off');
		uistate = uisuspend_j(fig, true);
        
        % restore button down functions for uicontrol children of the figure
        uirestore_j(uistate,'uicontrols');
    
        pandata.uistate = uistate;
        pandata.enable = 'on';
        %setappdata(fig,'PanFigureState',pandata);
    end
        
    % Enable pan mode
    set(fig,'WindowButtonDownFcn',{@locWindowButtonDownFcn,fig});
    set(fig,'WindowButtonMotionFcn',{@locWindowButtonMotionFcn,target});
    set(fig,'WindowButtonUpFcn',{@locWindowButtonUpFcn,fig});
    set(fig,'KeyPressFcn',{@locKeyPressFcn,target});
    setptr(fig,'hand');
         
    locDoPanOn(fig);
    setappdata(fig,'PanFigureState',pandata);
else 
    if ~isempty(uistate)
        % restore figure and non-uicontrol children
        % don't restore uicontrols because they were restored
        % already when rotate3d was turned on
%        uirestore_j(uistate,'nouicontrols');
        uirestore_j(uistate);
    end
    pandata.uistate = [];
    pandata.enable = 'off';
    setappdata(fig,'PanFigureState',pandata);
    
    % Turn off UI
    locDoPanOff(fig);
end

%-----------------------------------------------%
function locWindowButtonDownFcn(obj,evd,fig)
% Begin panning

	pandata = locGetData(fig);
	ax = locGetAxes(fig);
	if isempty(ax) || ~ishandle(ax)
		return;
	end

	if (strcmpi(get(fig,'selectiontype'), 'normal'))		% left click
		pandata.axes = ax;
		setptr(fig,'closedhand');
	end

	% Store original axis limits
	pandata.orig_axlim = [get(ax,'XLim') get(ax,'YLim')];

	setappdata(fig,'PanFigureState',pandata);

%-----------------------------------------------%
function locDoPan(ax,origlim,newlim)
% Pan
	set(ax,'XLim',newlim(1:2),'YLim',newlim(3:4),'XLimMode','manual','YLimMode','manual');

%-----------------------------------------------%
function locWindowButtonUpFcn(obj,evd,fig)
% Stop panning

	setptr(fig,'hand');
	pandata = locGetData(fig);

	% Axes may be empty if we are double clicking
	if ishandle(pandata.axes)
		newlim = [get(pandata.axes,'XLim') get(pandata.axes,'YLim')];
		origlim = pandata.orig_axlim;
		ax = pandata.axes;

		if (~isequal(newlim,origlim))
			locDoPan(ax,origlim,newlim);
		end
	end

	% Clear all transient pan state
	pandata.axes = [];
	pandata.last_pixel = [];

	setappdata(fig,'PanFigureState',pandata);

%-----------------------------------------------%
function locWindowButtonMotionFcn(obj,evd,target)
% This gets called everytime we move the mouse in pan mode, 
% regardless of whether any buttons are pressed. 

	fig = target;			% Remember,	target == fig | axes
	pandata = locGetData(fig);
	ax = pandata.axes;

	if isempty(ax) || ~ishandle(ax)
		return
	else
		fig = get(ax,'parent');
	end

	% Only pan if we have a previous pixel point
	ok2pan = ~isempty(pandata.last_pixel); 

	% Get current point in pixels
	curr_pixel = get(fig,'CurrentPoint');
	curr_pixel = curr_pixel(1:2);

	if ok2pan         
		delta_pixel = curr_pixel - pandata.last_pixel;
		orig_units = get(ax,'units');
		set(ax,'units','pixel');
		set(ax,'units',orig_units);
		locDataPan(ax,delta_pixel(1),delta_pixel(2),pandata.style);
		imscroll_j(ax,'PanSetSliders')
	end

	pandata.last_pixel = curr_pixel;
	setappdata(fig,'PanFigureState',pandata);

%-----------------------------------------------%
function locKeyPressFcn(obj,evd,target)
% Pan if the user clicks on arrow keys

	fig = target;
	pandata = locGetData(fig);
	ax = get(fig,'CurrentAxes');

	% Bail out if no handle to axes
	if ~ishandle(ax),	return,		end

	ch = get(fig,'CurrentCharacter');

	% Bail out if empty
	if isempty(ch),		return,		end

	style = pandata.style;

	switch uint16(ch)
	  case 28  % ascii left arrow
		locDataPan(ax,-2,0,style);
	  case 29 % ascii right arrow
		locDataPan(ax,2,0,style);  
	  case 30 % ascii up arrow 
		locDataPan(ax,0,2,style);
	  case 31 % ascii down arrow
		locDataPan(ax,0,-2,style);
	end

	drawnow			% required to prevent "flickering"

%-----------------------------------------------%
function locDataPan(ax, delta_pixel1, delta_pixel2, style)
% This is where the panning computation occurs.

	orig_units = get(ax,'units');
	set(ax,'units','pixel');
	range_pixel = get(ax,'position');

	set(ax,'units',orig_units);

	orig_lim1 = get(ax,'xlim');
	orig_lim2 = get(ax,'ylim');

	curr_lim1 = orig_lim1;
	curr_lim2 = orig_lim2;

	% For log plots, transform to linear scale
	if strcmp(get(ax,'xscale'),'log')
		is_abscissa_log = true;
		curr_lim1 = log10(curr_lim1);
	else
		is_abscissa_log = false;
	end
	if strcmp(get(ax,'yscale'),'log')
		is_ordinate_log = true;
		curr_lim2 = log10(curr_lim2);
	else
		is_ordinate_log = false;
	end

	range_data1 = abs(diff(curr_lim1));
	range_data2 = abs(diff(curr_lim2));

	pixel_width = range_pixel(3);
	pixel_height = range_pixel(4);

	delta_data1 = delta_pixel1 * range_data1 / pixel_width;
	delta_data2 = delta_pixel2 * range_data2 /  pixel_height;

	% Consider direction of axis: [{'normal'|'reverse'}]
	dir1 = get(ax, 'xdir');
	if strcmp(dir1,'reverse')
		new_lim1 = curr_lim1 + delta_data1;
	else
		new_lim1 = curr_lim1 - delta_data1;
	end

	dir2 = get(ax, 'ydir');
	if strcmp(dir2,'reverse')
		new_lim2 = curr_lim2 + delta_data2;
	else
		new_lim2 = curr_lim2 - delta_data2;
	end

	% For log plots, untransform limits
	if is_abscissa_log
		new_lim1 = 10.^new_lim1;
		%curr_lim2 = 10.^curr_lim2;
	end
	if is_ordinate_log
		%curr_lim1 = 10.^curr_lim1;
		new_lim2 = 10.^new_lim2;
	end

	h = findobj(ax,'Type','Image');
	if ~isempty(h) % Determine axis limits for image
		lims = getappdata(ax,'ThisImageLims');      % I think this is a more reliable value
		if (isempty(lims))
			lims = objbounds(ax);
		end
		x = lims(1:2);        y = lims(3:4);

		%If we are within the bounds of the image to begin with. This is to
		%prevent odd behavior if we panned outside the bounds of the image
		if ( x(1) <= orig_lim1(1) && x(2) >= orig_lim1(2) && y(1) <= orig_lim2(1) && y(2) >= orig_lim2(2) )
			dx = new_lim1(2) - new_lim1(1);
			if new_lim1(1) < x(1)
				new_lim1(1) = x(1);
				new_lim1(2) = new_lim1(1) + dx;
			end
			if new_lim1(2) > x(2)
				new_lim1(2) = x(2);
				new_lim1(1) = new_lim1(2) - dx;
			end
			dy = new_lim2(2) - new_lim2(1);
			if new_lim2(1) < y(1)
				new_lim2(1) = y(1);
				new_lim2(2) = new_lim2(1) + dy;
			end
			if new_lim2(2) > y(2)
				new_lim2(2) = y(2);
				new_lim2(1) = new_lim2(2) - dy;
			end
		end
	end

	% Set new limits
	fig = get(ax, 'Parent');
	labelType = getappdata(ax,'LabelFormatType');
	if strcmp(style,'x')
		set(ax,'xlim',new_lim1);
		set(ax,'ylim',orig_lim2);
		if ~isempty(labelType)
			set(ax, 'XTickLabelMode','auto')                    % Needed for ChangeAxesLabels to work
			ChangeAxesLabels(fig,ax,labelType,'X')
		end
	elseif strcmp(style,'y')
		set(ax,'xlim',orig_lim1);
		set(ax,'ylim',new_lim2);
		if ~isempty(labelType)
			set(ax, 'YTickLabelMode','auto')                    % Needed for ChangeAxesLabels to work
			ChangeAxesLabels(fig,ax,labelType,'Y')
		end
	else
		set(ax,'xlim',new_lim1);
		set(ax,'ylim',new_lim2);
		if ~isempty(labelType)
			set(ax, 'XTickLabelMode','auto', 'YTickLabelMode','auto')   % Needed for ChangeAxesLabels to work
			ChangeAxesLabels(fig,ax,labelType,'X')
			ChangeAxesLabels(fig,ax,labelType,'Y')
		end
	end

%-----------------------------------------------%
function [ax] = locGetAxes(fig)
% Return the axes that the mouse is currently over
% Return empty if no axes found (i.e. axes has hidden handle)

% TBD, this code is replicated in zoom.m and rotate3d.m

ax = [];
if ~ishandle(fig),    return;   end

% Determine which axes mouse ptr hits. We can't use the 
% figure's "currentaxes" property since the current
% axes may not be under the mouse ptr.
allAxes = findobj(datachildren(fig),'type','axes');
orig_fig_units = get(fig,'units');
set(fig,'units','pixels');
for (i = 1:numel(allAxes))
   candidate_ax = allAxes(i);
   
   % Check behavior support
   doIgnore = false;
        
	if ~doIgnore
		cp = get(fig,'CurrentPoint');
		orig_ax_units = get(candidate_ax,'units');
		set(candidate_ax,'units','pixels')
		pos = get(candidate_ax,'position');
		set(candidate_ax,'units',orig_ax_units)
		if cp(1) >= pos(1) && cp(1) <= pos(1)+pos(3) && ...
			cp(2) >= pos(2) && cp(2) <= pos(2)+pos(4)
			set(fig,'currentaxes',candidate_ax);
			ax = candidate_ax;
			break
		end
	end
end

% restore state
set(fig,'units',orig_fig_units)

%-----------------------------------------------%
function [pandata] = locGetData(fig)
pandata = getappdata(fig,'PanFigureState');
if isempty(pandata)
    pandata.uistate = [];
    pandata.axes = [];
    pandata.orig_axlim = [];
    pandata.last_pixel = [];
    pandata.style = 'xy';
    pandata.enable = 'off';
end

%-----------------------------------------------%
function locDoPanOn(fig)
	set(findall(fig,'tag','figMenuPan'),'Checked','on');

	% Remove when uitoolfactory is in place
	set(findall(fig,'tag','figToolPan'),'State','on');  

%-----------------------------------------------%
function locDoPanOff(fig)
	set(findall(fig,'tag','figMenuPan'),'Checked','off');

	% Remove when uitoolfactory is in place
	set(findall(fig,'tag','figToolPan'),'State','off');
  
%---------------------------------------------------------------------------------
function ChangeAxesLabels(hFig,hAx,type,eixo)
% This function formats the axes labels strings using a geographical notation
% EIXO is a single character (X or Y) that indicates on which axe (abssissa or
% ordenate) this function will operate with.

tick = get(hAx,[eixo,'TickLabel']);
n_tick = size(tick,1);

switch type
    case 'DegDec'
        % Nothing to do
    case 'DegMin'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMM',0,'str');   % e_str is a structure with string fields
        str_e = [e_str.dd repmat(' ',n_tick,1) e_str.mm];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMin')                     % Save it so zoom can know the label type
    case 'DegMinDec'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMM.x',2,'str'); % e_str is a structure with string fields
        str_e = [e_str.dd repmat(' ',n_tick,1) e_str.mm];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMinDec')                  % Save it so zoom can know the label type
    case 'DegMinSec'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMMSS',0,'str'); % e_str is a structure with string fields
        str_e = [e_str.dd repmat(' ',n_tick,1) e_str.mm repmat(' ',n_tick,1) e_str.ss];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMinSec')                  % Save it so zoom can know the label type
    case 'DegMinSecDec'
        e_str = degree2dms(str2num( ddewhite(tick) ),'DDMMSS.x',1,'str');   % e_str is a structure with string fields
        str_e = [e_str.dd repmat(' ',n_tick,1) e_str.mm repmat(' ',n_tick,1) e_str.ss];
        set(hAx,[eixo,'TickLabel'],str_e);
        setappdata(hAx,'LabelFormatType','DegMinSecDec')               % Save it so zoom can know the label type
end
