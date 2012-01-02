function varargout = ground_motion(varargin)
% Helper figure to control ground motion computations
%
% References:
%	David J. Wald and Trevor I. Allen. Topographic Slope as a Proxy for 
%	Seismic Site Conditions and Amplification. 
%	Bulletin of the Seismological Society of America; October 2007; v. 97; no. 5; p. 1379-1395;
%	DOI: 10.1785/0120060267

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

	if isempty(varargin)
		errordlg('Ground Motion: wrong number of arguments.','Error'),		return
	end

	handMir  = varargin{1};
	if (handMir.no_file),		return,		end
	if (~handMir.validGrid)
		errordlg('Ground Motion: This operation is deffined only for images derived from DEM grids.','ERROR')
		return
	end
	if (~handMir.geog)
		errordlg('Ground Motion: Currently this operation is only deffined for geographic grids.','ERROR')
		return
	end
	if (handMir.head(6) <= 0)
		errordlg('Ground Motion: Grid has only negative depths (underwater), nothing to do here. Bye Bye.','ERROR')
		return
	end

	hObject = figure('Vis','off');
	ground_motion_LayoutFcn(hObject);
	handles = guihandles(hObject);
	move2side(handMir.figure1, hObject, 'right')

	Z = getappdata(handMir.figure1,'dem_z');
    
	if (isempty(Z))
		errordlg('Ground Motion: Grid was not saved in memory. Increase "Grid max size" and start over.','ERROR')
		delete(hObject);    return
	end
	handles.mag = 7;
	handles.handMir = handMir;

	%------------ Give a Pro look (3D) to the frame boxes  -------------------------------
	new_frame3D(hObject, [handles.text1 handles.text2])
	%------------- END Pro look (3D) -------------------------------------------------------

	guidata(hObject, handles);

	set(hObject,'Visible','on');
	if (nargout),   varargout{1} = hObject;     end

% --------------------------------------------------------------------
function radio_active_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_stable handles.radio_guess],'Value', 0)

% --------------------------------------------------------------------
function radio_stable_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_active handles.radio_guess],'Value', 0)

% --------------------------------------------------------------------
function radio_guess_CB(hObject, handles)
	if (~get(hObject,'Val')),	set(hObject,'Val',1),	return,		end
	set([handles.radio_active handles.radio_stable],'Value', 0)

% --------------------------------------------------------------------
function edit_mag_CB(hObject, handles)
	xx = str2double(get(hObject,'String'));
	if (isnan(xx) || xx < 2),		set(hObject,'String',handles.mag),	return,	end
	handles.mag = xx;
	guidata(handles.figure1, handles)

% --------------------------------------------------------------------
function push_OK_CB(hObject, handles)
	% See what to compute
	opt_O = '-O';
	if (get(handles.check_INT,'Value')),	opt_O = [opt_O 'i'];	end
	if (get(handles.check_PGA,'Value')),	opt_O = [opt_O 'a'];	end
	if (get(handles.check_PGV,'Value')),	opt_O = [opt_O 'v'];	end
	if (strcmp(opt_O,'-O')),	opt_O = '-Oi';			end			% Default

	opt_F = sprintf('-F%d',get(handles.popup_meca,'Value'));		% Meca type
	opt_M = sprintf('-M%d',handles.mag);							% Magnitude

	h = findobj(handles.handMir.axes1,'type','line','Tag','FaultTrace');
	if (numel(h) > 1)
		warndlg('You have more than one Fault trace in figure. Using only the first one found.','Warning')
	end
	if (isempty(h))			% Hmm no Fault trace. Search for plain lines and use first one found
		h = findobj(handles.handMir.axes1,'type','line');
	end
	if (isempty(h))
		errordlg('You need a ''Fault'' to use this option. Bye Bye.','Error'),		return
	end
	h = h(1);
	x = get(h,'XData');		y = get(h,'YData');		xy = [x(:) y(:)];

	aguentabar(0,'title','Computing VS30 velocities')
	[X,Y,vs30,head] = calc_vs30(handles);

	aguentabar(0.5,'title','Now computing what you asked')
	n_out = numel(opt_O(3:end));
	[varargout{1:n_out}] = shake_mex(vs30, head, xy, opt_O, opt_F, opt_M);
	aguentabar(1,'title','Done'),	pause(0.01)

	ind_0 = [];
	if (handles.handMir.head(5) < 0)		% We have an underwater zone. Make it zero on products
		Z = getappdata(handles.handMir.figure1,'dem_z');
		ind_0 = (Z < 0);
	end

	tmp.X = X;		tmp.Y = Y;
	for (k = 1:n_out)
		if (~isempty(ind_0))				% Sea side
			varargout{k}(ind_0) = NaN;
		end
		zz = grdutils(varargout{k},'-L');		head(5:6) = double(zz(1:2));
		if (opt_O(k+2) == 'i')
			n = ones(1,16);
			pal = ([221*n 208*n 195*n 181*n 168*n 154*n 141*n 128*n 114*n 101*n 87*n 74*n 60*n 47*n 34*n 20*n] / 255)';
			tmp.cmap = [pal pal pal];		% Build a gray color map
			tmp.head = head;		tmp.name = 'Intensity';		mirone(varargout{k},tmp);
			X = tmp.X;		Y = tmp.Y;
			clear tmp;						% We don't want to apply the gray cmap to the grids
			tmp.X = X;		tmp.Y = Y;
		end
		if (opt_O(k+2) == 'a')
			tmp.head = head;		tmp.name = 'PGA';			mirone(varargout{k},tmp);
		end
		if (opt_O(k+2) == 'v')
			tmp.head = head;		tmp.name = 'PGV';			mirone(varargout{k},tmp);
		end
	end

% ---------------------------------------------------------------------------------------------
function [X,Y,vs30,head] = calc_vs30(handles)
% Compute the VS30 (m/sec) array from topography slope according to Wald & Allen, 2007

	[X,Y,slope,head] = GetSlope(handles.handMir);		% Get the slope in m/m

	% Values from the Slope Range of Table 2 in Wald & Allen, 2007
	if (get(handles.radio_stable,'Value'))
		rng_slope = [1e-6 2e-3; 2e-3 4e-3; 4e-3 7.2e-3; 7.2e-3 0.013; 0.013 0.018; 0.018 0.025];
	elseif (get(handles.radio_active,'Value'))
		rng_slope = [3.2e-5 2.2e-3; 2.2e-3 4.3e-3; 6.3e-3 0.018; 0.018 0.05; 0.05 0.1; 0.1 0.138];
	else
		ms = mean(slope(:));		% Find mean slope to apply point 2) of Wald & Allen's recipe
		if (ms < 0.05)
			rng_slope = [1e-6 2e-3; 2e-3 4e-3; 4e-3 7.2e-3; 7.2e-3 0.013; 0.013 0.018; 0.018 0.025];
		else
			rng_slope = [3.2e-5 2.2e-3; 2.2e-3 4.3e-3; 6.3e-3 0.018; 0.018 0.05; 0.05 0.1; 0.1 0.138];
		end
	end

	rng_class = [180 240; 240 300; 300 360; 360 490; 490 620; 620 760];

	vs30 = zeros(size(slope));
	for (m = 1:size(rng_slope,1))
		ind = (slope >= rng_slope(m,1) & slope < rng_slope(m,2));
		ind = find(ind(:));
		delta_class = rng_class(m,2) - rng_class(m,1);
		for (k = 1:numel(ind))
			x = slope(ind(k));
			vs30(ind(k)) = rng_class(m,1) + (x - rng_class(m,1)) / delta_class;
		end
	end

	% Now do the (const) lower and upper limits
	ind = (slope < rng_slope(1));		vs30(ind) = rng_class(1);
	ind = (slope > rng_slope(end));		vs30(ind) = rng_class(end);
	vs30 = single(vs30);

% --------------------------------------------------------------------
function [X,Y,slope,head] = GetSlope(handles)
% This 'handles' is the mirone handles

	[X,Y,Z,head] = load_grd(handles);
	
	D2R = pi/180;
	if (handles.geog)
		near1km = 1 / 120;			% Half a minute arc ~ 1 km
		if (head(9) < near1km / 2)	% If grid has resolution than ~450 m
			winSize = round( min(near1km / head(8), near1km / head(9)) );
			if (rem(winSize, 2) == 0)	winSize = winSize + 1;		end		% Force odd sized window
			opt_W = sprintf('-W%d', winSize);
			opt_A = '-A6';			% Slope
			opt_N = sprintf('-N%d', handles.have_nans);
			slope = double(mirblock(Z, head, opt_A, opt_N, opt_W, '-G'));
		else						% Coarser grid (should had a warning message)
			[X,Y,Z,head] = load_grd(handles,'double');
			slope = gradient_geo(Y,X,Z,'slope');
		end
		slope = tan(slope*D2R);
	else			% Not used, but maybe one day
		nz = getnormals(X,Y,Z);			slope = acos(nz);
		slope = tan(slope);
	end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata)
	if isequal(get(hObject,'CurrentKey'),'escape')
		delete(hObject);
	end

% --- Creates and returns a handle to the GUI figure. 
function ground_motion_LayoutFcn(h1)

set(h1, 'Position',[520 649 427 151],...
'Color',get(0,'factoryUicontrolBackgroundColor'),...
'KeyPressFcn',@figure1_KeyPressFcn,...
'MenuBar','none',...
'Name','Ground motion',...
'NumberTitle','off',...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1');

uicontrol('Parent',h1, 'Position',[180 61 240 79], 'Style','frame');
uicontrol('Parent',h1, 'Position',[7 61 163 79], 'Style','frame');

uicontrol('Parent',h1, 'Position',[14 112 66 15],...
'FontName','Helvetica',...
'String','Intensity',...
'Style','checkbox',...
'Tooltip','Calculate Intensity map',...
'Value',1,...
'Tag','check_INT');

uicontrol('Parent',h1, 'Position',[14 90 150 16],...
'FontName','Helvetica',...
'String','Peak Ground Acceleration',...
'Style','checkbox',...
'Tooltip','Calculate PGA map',...
'Tag','check_PGA');

uicontrol('Parent',h1, 'Position',[14 69 150 16],...
'FontName','Helvetica',...
'String','Peak Ground Velocity',...
'Style','checkbox',...
'Tooltip','Calculate PGV map',...
'Tag','check_PGV');

uicontrol('Parent',h1, 'Position',[32 130 110 18],...
'FontName','Helvetica',...
'FontSize',9,...
'String','What to compute',...
'Style','text',...
'Tag','text1');

uicontrol('Parent',h1, 'Position',[7 5 132 23],...
'BackgroundColor',[1 1 1],...
'String',{'Unknown'; 'Strike-slip'; 'Normal'; 'Thrust'},...
'Style','popupmenu',...
'Value',1,...
'Tag','popup_meca');

uicontrol('Parent',h1, 'Position',[6 30 132 18],...
'FontName','Helvetica',...
'FontSize',9,...
'String','Focal mechanism type',...
'Style','text');

uicontrol('Parent',h1, 'Position',[190 112 225 17],...
'Call',@ground_motion_uiCB,...
'FontName','Helvetica',...
'String','Force "Active Tectonic" in VS30 calculus',...
'Style','radiobutton',...
'Tag','radio_active');

uicontrol('Parent',h1, 'Position',[190 90 225 17],...
'Call',@ground_motion_uiCB,...
'FontName','Helvetica',...
'String','Force "Stable Continent" in VS30 calculus',...
'Style','radiobutton',...
'Tag','radio_stable');

uicontrol('Parent',h1, 'Position',[190 69 110 17],...
'Call',@ground_motion_uiCB,...
'FontName','Helvetica',...
'String','Let me guess',...
'Style','radiobutton',...
'Tooltip','I''ll do the guessing based on the average slope',...
'Value',1,...
'Tag','radio_guess');

uicontrol('Parent',h1, 'Position',[224 130 150 18],...
'FontName','Helvetica',...
'FontSize',9,...
'String','How to compute VS30?',...
'Style','text',...
'Tag','text2');

uicontrol('Parent',h1, 'Position',[243 8 47 21],...
'BackgroundColor',[1 1 1],...
'Call',@ground_motion_uiCB,...
'String','7',...
'Style','edit',...
'Tooltip','Event Magnitude',...
'Tag','edit_mag');

uicontrol('Parent',h1, 'Position',[181 11 60 15],...
'FontName','Helvetica',...
'HorizontalAlignment','right',...
'ListboxTop',0,...
'String','Magnitude:',...
'Style','text');

uicontrol('Parent',h1, 'Position',[341 9 80 21],...
'Call',@ground_motion_uiCB,...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','bold',...
'String','OK',...
'Tag','push_OK');

function ground_motion_uiCB(hObject, eventdata)
% This function is executed by the callback and than the handles is allways updated.
	feval([get(hObject,'Tag') '_CB'],hObject, guidata(hObject));
