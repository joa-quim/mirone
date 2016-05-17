function magnify(hAx)
%  Creates a magnification box centered on last clicked point.
%
% Press '+'/'-' to increase/decrease magnification.
% Press '>'/'<' to increase/decrease box size.
% Left-click to delete the magnification on figure.

% Modified from 'magnify' of Rick Hindman - FEX-5961 (BSD License)
% Probably it uses too much memory because it duplicates the all image and
% it should be possible to make it access the original. But maybe the lazzy
% copy mechanism saves us from this memory duplication.

% Coffeeright J. Luis 2004-2016

% $Id$

	if (nargin == 0), hAx = gca; end
	oldFigureUIState = uisuspend_j(get(hAx, 'parent'), true);
	setappdata(hAx, 'Magnify_UIState', oldFigureUIState)

	hFig = get(hAx, 'parent');
	set(get(hAx, 'parent'), 'WindowButtonDownFcn', {@wbd_magnifyCB, hAx}, ...
		'WindowButtonMotionFcn', {@wbm_magnifyCB, hAx}, 'KeyPressFcn', {@KeyPress_magnifyCB, hAx});

	hAxMag = copyobj(hAx,hFig);

	set(hFig, 'Pointer','crosshair', 'CurrentAxes', hAxMag);
	setappdata(hAx, 'Magnified', hAxMag);
	set(hAxMag, 'UserData',[2,0.2], 'Color',get(hAx,'Color'), 'Box','on');  % magnification, frame size
	%xlabel(''); ylabel(''); zlabel(''); title('');
	set(hAx, 'Color',get(hAx,'Color')*0.95);
	set(hFig, 'CurrentAxes', hAx);
	wbm_magnifyCB([], [], hAx);

function wbd_magnifyCB(obj, evt, hAx)
	hFig = get(hAx, 'parent');
	hAxMag = getappdata(hAx, 'Magnified');
	set(hAx, 'Color', get(hAxMag,'Color'));
	set(hFig, 'Pointer','arrow', 'CurrentAxes',hAx);
	uirestore_j(getappdata(hAx, 'Magnify_UIState'))
	if ~strcmp(get(hFig,'SelectionType'),'alt')
		delete(hAxMag)
		setappdata(hAx, 'Magnified', '')
	end

function wbm_magnifyCB(obj, evt, hAx)
	hFig = get(hAx, 'parent');
	hAxMag = getappdata(hAx, 'Magnified');
	if ~isempty(hAxMag)
		a2_param = get(hAxMag,'UserData');
		f_pos  = get(hFig,'Position');
		a1_pos = get(hAx,'Position');
		fUnit = get(hFig,'Units');			set(hFig,'Units','pixels')
		f_cp  = get(hFig,'CurrentPoint');
		a1_cp = get(hAx,'CurrentPoint');
		set(hFig, 'Units', fUnit)
		set(hAxMag,'Position',[(f_cp./f_pos(3:4)) 0 0] + a2_param(2)*a1_pos(3)*[-1 -1 2 2]);
		a2_pos = get(hAxMag,'Position');
		set(hAxMag,'XLim',a1_cp(1,1)+(1/a2_param(1))*(a2_pos(3)/a1_pos(3))*diff(get(hAx,'XLim'))*[-0.5 0.5]);
		set(hAxMag,'YLim',a1_cp(1,2)+(1/a2_param(1))*(a2_pos(4)/a1_pos(4))*diff(get(hAx,'YLim'))*[-0.5 0.5]);
	end

function KeyPress_magnifyCB(obj, evt, hAx)
	hFig = get(hAx, 'parent');
	hAxMag = getappdata(hAx, 'Magnified');
	if ~isempty(hAxMag)
		a2_param = get(hAxMag, 'UserData');
		if (strcmp(get(hFig,'CurrentCharacter'),'+') || strcmp(get(hFig,'CurrentCharacter'),'='))
			a2_param(1) = a2_param(1)*1.2;
		elseif (strcmp(get(hFig,'CurrentCharacter'),'-') || strcmp(get(hFig,'CurrentCharacter'),'_'))
			a2_param(1) = a2_param(1)/1.2;
		elseif (strcmp(get(hFig,'CurrentCharacter'),'<') || strcmp(get(hFig,'CurrentCharacter'),','))
			a2_param(2) = a2_param(2)/1.2;
		elseif (strcmp(get(hFig,'CurrentCharacter'),'>') || strcmp(get(hFig,'CurrentCharacter'),'.'))
			a2_param(2) = a2_param(2)*1.2;
		end
		set(hAxMag, 'UserData', a2_param);
		wbm_magnifyCB([], [], hAx);
	end
