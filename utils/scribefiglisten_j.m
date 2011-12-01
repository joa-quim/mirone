function scribefiglisten_j(fig,onoff)
%SCRIBEFIGLISTEN listeners for figures and their axes children.
% SCRIBEFIGLISTEN(fig,onoff) creates child added/removed listeners 
% (if they do not already exist) for fig and its non-legend,
% non-colorbar axes children, and enables (onoff=true) or 
% disables(onoff=false)them. Firing listeners turns off zoom_j for the figure.
% Called by zoom_j
%

%   Glen M. DeLoid 02-01-2001
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.5 $  $Date: 2002/04/08 21:44:36 $

	% create listeners if they don't already exist
	if isempty(findprop(handle(fig),'ScribeFigListeners'))
		% create listeners property
		pl = schema.prop(handle(fig),'ScribeFigListeners','MATLAB array');
		pl.AccessFlags.Serialize = 'off';
		% create listeners array
		l.chadd = handle.listener(fig,'ObjectChildAdded', {@scribeFigChildAdded, fig});
		l.chremove = handle.listener(fig,'ObjectChildRemoved', {@scribeFigChildRemoved, fig});
		% set listeners property to listeners array
		set(handle(fig),'ScribeFigListeners',l);
	end

	% set figure listeners enabled/disabled
	l = get(handle(fig),'ScribeFigListeners');
	set(l.chadd,'Enabled',onoff);
	set(l.chremove,'Enabled',onoff);

	% get list of axes children
	ax = findobj(fig,'type','axes');

	if ~isempty(ax)
		% create listeners for axes if they don't exist
		for (i = 1:numel(ax))
			% don't create them for axes that are legends or colorbars
			if ~any(strcmp(get(get(ax(i),'children'),'tag'),'TMW_COLORBAR')) && ~strcmp(get(ax(i),'tag'),'legend')
				if isempty(findprop(handle(ax(i)),'ScribeFigAxListeners'))
					% create listeners property
					pl = schema.prop(handle(ax(i)),'ScribeFigAxListeners','MATLAB array');
					pl.AccessFlags.Serialize = 'off';
					% create listeners array
					l.chadd = handle.listener(ax(i),'ObjectChildAdded', {@scribeFigAxChildAdded, fig});
					l.chremove = handle.listener(ax(i),'ObjectChildRemoved', {@scribeFigAxChildRemoved, fig});
					% set listeners property to listeners array
					set(handle(ax(i)),'ScribeFigAxListeners',l);
				end
				% set figure axes listeners enabled/disabled
				l = get(handle(ax(i)),'ScribeFigAxListeners');
				set(l.chadd,'Enabled',onoff);
				set(l.chremove,'Enabled',onoff);
			end
		end
	end

%------------------------------------------------------------------------%
function scribeFigChildAdded(src,event,fig)
% figure add child callback if the added child is an axes or is 
% a uicontextmenu or uicontrol turn off zoom_j

	chh = handle(event.child);
	tag = get(chh,'tag');
	type = get(chh,'type');

	% don't do anything if child being added is a temporary text or uicontrol
	% object (e.g. used to get the dimensions of a string) as in legend
	if (strcmpi(tag,'temphackytext') && strcmpi(type,'text')) || ...
			(strcmpi(tag,'temphackyui') && strcmpi(type,'uicontrol'))

		% if obect being added is a legend, and zoom_j is on, turn it
		% off and then on again to make state data saved by uiclearmode current.
	elseif strcmpi(tag,'legend')
		if strcmpi('out',zoom_j(fig,'getmode'))
			zoom_j(fig,'off');
			zoom_j(fig,'outmode');
		elseif isappdata(fig,'ZoomOnState')
			zoomstate = getappdata(fig,'ZoomOnState');
			zoom_j(fig,'off');
			zoom_j(fig,zoomstate);
		end

		% otherwise if object is an axes, uicontextmenu or uicontrol turn zoom_j off.
	elseif  event.child.isa('hg.axes') | event.child.isa('hg.uicontextmenu') | event.child.isa('hg.uicontrol')

		zoom_j(fig,'off');
	end

%------------------------------------------------------------------------%
function scribeFigChildRemoved(src,event,fig)
% figure remove child callback if the child is an axes or is
% a uicontextmenu or uicontrol turn off zoom_j

	chh = handle(event.child);
	tag = get(chh,'tag');
	type = get(chh,'type');

	% don't do anything if child being removed is a temporary text or ui
	% object (e.g. used to get the dimensions of a string) as in legend
	if (strcmpi(tag,'temphackytext') && strcmpi(type,'text')) || ...
			(strcmpi(tag,'temphackyui') && strcmpi(type,'uicontrol'))

		% if object being removed is a legend, and zoom_j is on, turn
		% it off and then on again to make state data saved by uiclearmode current.
	elseif strcmpi(tag,'legend')
		if strcmpi('out',zoom_j(fig,'getmode'))
			zoom_j(fig,'off');
			zoom_j(fig,'outmode');
		elseif isappdata(fig,'ZoomOnState')
			zoomstate = getappdata(fig,'ZoomOnState');
			zoom_j(fig,'off');
			zoom_j(fig,zoomstate);
		end

		% otherwise if object is an axes, uicontextmenu or uicontrol turn zoom_j and rotate3d off.
	elseif  event.child.isa('hg.axes') | event.child.isa('hg.uicontextmenu') | event.child.isa('hg.uicontrol')

		zoom_j(fig,'off');
	end

%------------------------------------------------------------------------%
function scribeFigAxChildAdded(src,event,fig)
% axes add child callback turn off zoom_j
	chh = handle(event.child);
	tag = get(chh,'tag');
	type = get(chh,'type');

	% don't turn zoom_j off if object being added is a temporary text or ui object
	% or a legend delete proxy
	if ~strcmpi(tag,'LegendDeleteProxy') && ~(strcmpi(tag,'temphackytext') && ...
			strcmpi(type,'text')) && ~(strcmpi(tag,'temphackyui') && strcmpi(type,'uicontrol'))
		zoom_j(fig,'off');
	end

%------------------------------------------------------------------------%
function scribeFigAxChildRemoved(src,event,fig)
% axes remove child callback turn off zoom_j
	chh = handle(event.child);
	tag = get(chh,'tag');
	type = get(chh,'type');

	% don't turn zoom_j off if object being added is a temporary text or ui object
	% or a legend delete proxy
	if ~strcmpi(tag,'LegendDeleteProxy') && ~(strcmpi(tag,'temphackytext') && ...
			strcmpi(type,'text')) && ~(strcmpi(tag,'temphackyui') && strcmpi(type,'uicontrol'))
		zoom_j(fig,'off');
	end
