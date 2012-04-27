function stopBar =  aguentabar(varargin)
% AGUENTABAR Display a await bar.
%
% Usage:
%   STOP = AGUENTABAR(X) will set the length of the bar in the most recently
%   created aguentabar window to the fractional length X.
%   The handle to the aguentabar figure is returned in STOP.

%   X specifies what fraction (0.0 - 1.0) of the task is complete.
%   Typically, the figure will be updated according to that value. However, if
%   X == 0.0, a new figure is created (an existing figure would be
%   closed first). If X == 1.0, the AGUENTABAR figure will close.
%   The color of the AGUENTABAR is choosen randomly when it is created or
%   reset. Clicking inside the figure will cause a random color change.
%   For best results, call aguentabar(0) (or just aguentabar) before starting
%   a task. This sets the proper starting time to calculate time remaining.
%
%   AGUENTABAR(X,'title','title string') will create|update the title text in
%   the aguentabar figure, in addition to setting the fractional length to X.
%   AGUENTABAR(X,'title','title string') with a negative X updates only the 'title string'
%
%   STOP = AGUENTABAR(..., 'CreateCancelBtn') adds a cancel button to the figure.
%   Hiting this Cancel button will signal the figure to be killed on its next call
%   and the return value of STOP = NaN. So the calling code can use the STOP value
%   test whether or not to interrupt a lengthy calculation.
%
% Example:
%   n = 1000;
%   aguentabar % Create figure and set starting time
%   for i = 1:n
%       pause(0.01) % Do stuff
%       aguentabar(i/n) % Update figure
%   end
%
%	HISTORY:
% 		ORIGINAL VERSION: Steve Hoelzer
% 		WORKING VERSION: m-file modified by Quan Quach on 12/12/07 (quan.quach@gmail.com)
% 		HOWEVER, the Quan-Quach simply does not work as advertized (well, at least I had lots of failurs)
% 		It has a couple of bugs and by no means it can do what it it claims (return a 1 when the figure is killed)
%
% 	So I re-wrote an important part of this code
% 		Now -
%			1) No more persistent variables
%			2) Killing the figure via 'Cancel' button returns '1' which allows calling code to stop
%			3) One can provide a 'Title' and even change that title between calls
%
% 	05-Jan-2008		Joaquim Luis (jluis@ualg.pt)

%	Copyright (c) 2004-2012 by J. Luis
%
% 	This program is part of Mirone and is free software; you can redistribute
% 	it and/or modify it under the terms of the GNU Lesser General Public
% 	License as published by the Free Software Foundation; either
% 	version 2.1 of the License, or any later version.
% 
% 	This program is distributed in the hope that it will be useful,
% 	but WITHOUT ANY WARRANTY; without even the implied warranty of
% 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% 	Lesser General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

	% Parse inputs
	n_argin = nargin;
	fracdone = 0;		haveCancel = 0;			titulo = [];		% Default values
	if (n_argin)
		for ( i = 1:numel(varargin) )
			if ( isnumeric(varargin{i}) )
				fracdone = varargin{i};
			elseif ( strncmpi(varargin{i}, 'title',5) )
				try		titulo = varargin{i+1};				% In case title string was not provided
				catch,	titulo = '';
				end
			elseif ( strncmpi(varargin{i}, 'createcancelbtn',7) )
				haveCancel = 1;
			end
		end
	end

	hFig = [];		hPatch = [];	starttime = [];		lastupdate = [];	killBill = false;
	f = findobj(allchild(0),'flat','Tag','AGUENTAbar');
	if (~isempty(f))
		hFig = f(1);
		if (fracdone ~= 0)
			ud = get(hFig, 'UserData');
			hPatch = ud(1);				starttime = ud(2:7);
			lastupdate = ud(8:13);		killBill = ud(14);
			if (fracdone < 0)			% Change only the title string but keep the percentage
				fracdone = get(hPatch, 'UserData');
			end
		else
			% Progress bar needs to be reset, close figure and set handle to empty
			delete(hFig)		% Close progress bar
			hFig = [];			% Set to empty so a new progress bar is created
		end
	end

	if ( killBill )				% Ah Ha, user hit (and accepted) the Cancel button
		if (nargout),	stopBar = nan;	end
		delete(hFig)
		return
	end

	percentdone = floor(100*fracdone);

	% Create new progress bar if needed
	if ( isempty(hFig) )
		screenSize = get(0,'ScreenSize');
		x0 = screenSize(3)/2-160;			y0 = screenSize(4)/2-20;

		figPos = [x0 y0 320 20];			axPos = [5 4 315 15];
		butaoPos = [132 4 66 21];
		if (~isempty(titulo) && ~haveCancel)		% With title
			figPos = [x0 y0 320 45];		axPos = [5 4 315 15];
		elseif (isempty(titulo) && haveCancel)		% With Cancel button
			figPos = [x0 y0 320 50];		axPos = [5 30 315 15];
		elseif (~isempty(titulo) && haveCancel)		% With title plus Cancel button
			figPos = [x0 y0 320 70];		axPos = [5 30 315 15];
		end
	
        % Initialize progress bar
        hFig = figure('Units', 'pixels', 'Position', figPos, 'NumberTitle','off', 'Resize','off',...
			'MenuBar','none', 'Tag','AGUENTAbar', 'HandleVisibility','callback', 'BackingStore','off');
        hAxes = axes('Parent',hFig, 'Units','pixels', 'Position',axPos, 'XLim',[0 1], 'YLim',[0 1],...
            'Box','on', 'ytick',[], 'xtick',[]);
        hPatch = patch('XData',[0 0 0 0], 'YData',[0 0 1 1], 'Parent',hAxes, 'FaceColor',[.1 1 .1], 'EraseMode','none');
	
		if (haveCancel)
			uicontrol('Parent',hFig, 'Units','pixels', 'String','Cancel', ...
						'Position', butaoPos, 'Callback', {@closeBar, hFig})
		end
	
		if ( ~isempty(titulo) ),	set(get(hAxes,'Title'),'String',titulo,'Interpreter','none'),	end
	
		% enable this code if you want the bar to change colors when the user clicks on the progress bar
		set([hFig hAxes hPatch],'ButtonDownFcn',{@changecolor,hPatch})
		changecolor(0,0,hPatch)			% Select a random color
	
        % Set time of last update to ensure a redraw
        lastupdate = clock - 1;
	
        % Task starting time reference
        if (isempty(starttime) || (fracdone == 0)),		starttime = clock;		end

	elseif ( ~isempty(titulo) )			% Update title on an existing aguentabar figure
		set(get(findobj(hFig, 'type', 'axes'),'Title'),'String',titulo,'Interpreter','none')
	end

	% Enforce a minimum time interval between updates, but allows for
	% the case when the bar reaches 100% so that the user can see it
	if (etime(clock,lastupdate) < 0.1 && (percentdone < 100))
		set(hFig, 'UserData', [hPatch starttime lastupdate killBill])
		if (nargout),	stopBar = hFig;		end
		return
	end

	% Update progress patch
	set(hPatch,'XData',[0 fracdone fracdone 0], 'UserData',fracdone)

	% Update progress figure title bar
	if (fracdone == 0)
		titlebarstr = ' 0%';
	else
		runtime = etime(clock,starttime);
		timeleft = runtime/fracdone - runtime;
		timeleftstr = sec2timestr(timeleft);
		titlebarstr = sprintf('%2d%%    %s remaining', percentdone, timeleftstr);
	end
	set(hFig,'Name',titlebarstr)

	drawnow				% Force redraw to show changes

	% If task completed, close figure and clear vars, then exit
	if (percentdone == 100)		% Task completed
		delete(hFig)			% Close figure
		if (nargout),	stopBar = hFig;		end
		return
	end

	lastupdate = clock;			% Record time of this update

	set(hFig, 'UserData', [hPatch starttime lastupdate killBill])		% Store for the next round
	if (nargout),	stopBar = hFig;		end

% ------------------------------------------------------------------------------
function changecolor(obj,evt,hPatch)
% Change the color of the progress bar patch
	colorlim = 2.8;				% Must be <= 3.0 - This keeps the color from being too light
	thiscolor = rand(1,3);
	while (sum(thiscolor) > colorlim)
		thiscolor = rand(1,3);
	end
	set(hPatch,'FaceColor',thiscolor);

% ------------------------------------------------------------------------------
function timestr = sec2timestr(sec)
% Convert a time measurement from seconds into a human readable string.

% Convert seconds to other units
d = floor(sec/86400);		sec = sec - d*86400;		% Days and remaing seconds
h = floor(sec/3600);		sec = sec - h*3600;			% Hours and remaing seconds
m = floor(sec/60);			sec = sec - m*60;			% Minutes and remaing seconds
s = floor(sec);				% Seconds

% Create time string
if d > 0
	timestr = sprintf('%d day, %.1f hr', d, (h+m/60));
elseif h > 0
	timestr = sprintf('%d hr, %d min',h, m);
elseif m > 0
    if m > 9
		timestr = sprintf('%d min',m);
    else
		timestr = sprintf('%d min, %d sec',m,s);
    end
else
	timestr = sprintf('%d sec',s);
end

% ----------------------------------------------------------------------------------
function closeBar(obj, evt, hFig)
% Signal, using data stored in UserData, that next time 'aguentabar' is called
% it is to return a code to the calling program indicating a STOP request.
	resp = questdlg('Do you want to stop this process (in the next call)?', 'Stop process', 'Yes','No','Yes');
	if strcmp(resp,'Yes')
		ud = get(hFig, 'UserData');
		ud(14) = true;
		set(hFig, 'UserData', ud)
	end
