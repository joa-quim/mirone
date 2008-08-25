function move_obj(arg)
% callback function for draggable objects
% Any object can be made draggable via
%   set(obj, 'ButtonDownFcn', 'move_obj(1)');
% using deltas allows us to drag big objects
%
% The original version came from the lightspeed toolbox from Tom Minka, but I
% changed it to make it more robust and to respect previous 'windowbutton*' settings
% J. Luis

this_obj = gco;
last_pos = getappdata(this_obj,'pos');

% handle events
switch arg
	case 2    % button motion
		hAx = gca;
		pos = get(hAx, 'CurrentPoint');
		pos = pos(1,:);
	
		switch get(this_obj, 'type')
			case 'text',	obj_pos = get(this_obj, 'Pos')';
			case 'line', 
				obj_pos(1,:) = get(this_obj,'xdata');
				obj_pos(2,:) = get(this_obj,'ydata');
			otherwise error(['cannot handle type ' get(this_obj,'type')])
		end
	
		% if the scale is logarithmic then the delta is a ratio
		if strcmp(get(hAx,'xscale'), 'log')
			new_pos(1,:) = obj_pos(1,:) * (pos(1)/last_pos(1));
		else
			new_pos(1,:) = obj_pos(1,:) + (pos(1)-last_pos(1));
		end
		if strcmp(get(hAx,'yscale'), 'log')
			new_pos(2,:) = obj_pos(2,:) * (pos(2)/last_pos(2));
		else
			new_pos(2,:) = obj_pos(2,:) + (pos(2)-last_pos(2));
		end
		switch get(this_obj, 'type')
			case 'text',	set(this_obj, 'Pos', new_pos);
			case 'line',
				set(this_obj, 'xdata', new_pos(1,:), 'ydata', new_pos(2,:));
		end
		setappdata(this_obj,'pos', pos);
	case 1    % buttondown
		% start moving
		% dragging looks better with double buffering on
		hFig = gcf;
		set(hFig, 'DoubleBuffer', 'on');
		last_pos = get(gca,'CurrentPoint');
		setappdata(this_obj,'pos', last_pos(1,:));
		% set callbacks
		state = uisuspend_fig(hFig);		% Remember initial figure state
		set(hFig, 'pointer', 'fleur');
		set(hFig, 'windowbuttonmotionfcn', 'move_obj(2)');
		set(hFig, 'windowbuttonupfcn', {@wbu_move_obj,this_obj,state});
end

%--------------------------------------------------
function wbu_move_obj(obj,eventdata,h,state)
	uirestore_fig(state);			% Restore the figure's initial state
