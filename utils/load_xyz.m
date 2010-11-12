function varargout = load_xyz(handles, opt, opt2)
% Read a generic ascii file that can be, or not, a multiseg file
%
%	Multi-segs files accept -G, -W & -S GMT type options.
%	It does also deal with the case of ploting the isochrons.dat
%
%	HANDLES	->	Should be the Mirone handles. However, when this function is used with output
%				HANDLES can be just a simple structure with the field 'no_file' = false
%
%	Optional
%		OPT can be either [] in which case the fiename will be asked here or contain the filename
%		OPT2 can take several values 'arrows' to read a x,y,u,v file
%			'AsLine'		plots a "regular" line
%			'AsPoint'		plots the line vertex only using small filled circles 
%			'AsMaregraph'	plots yellow dots used in the Tsunami modeling tools
%			'FaultTrace'	plots lines/polylines used by the elastic deformation tools
%			'Isochron'		plots a isochrons polyline from the internal db.
%							Attention, this option needs that OPT is not empty
%			If not given defaults to 'AsLine'
%
% If first line in file is of the form
%		'>U_N_I_K'	plot a single line NaN separated
%		'>ARROW'	plot an arrow field
%		'>VIMAGE'	tell Fleder to plot a scene with a VIMAGE
%		'>-:'		swap 1st and 2nd columns (assumed as (y,x) -> (x,y))
%		'>CLOSE'	plot patches instead of lines (idependently of pline being closed or not)

%	Copyright (c) 2004-2010 by J. Luis
%
%	This program is free software; you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation; version 2 of the License.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	Contact info: w3.ualg.pt/~jluis/mirone
% --------------------------------------------------------------------

%	EXAMPLE CODE OF HOW TO CREATE A TEMPLATE FOR UICTX WHEN THESE ARE TOO MANY
% 	cmenuHand = get(h, 'UIContextMenu');
% 	setappdata(handles.figure1, 'cmenuHand', cmenuHand)
% 	cmenuHand = uicontextmenu('Parent',handles.figure1);
% 	set(h, 'UIContextMenu', cmenuHand);
% 	%uimenu(cmenuHand, 'Label', 'Set all UIcontexts', 'Call', {@resetUIctx,h,handles.axes1});
% 	uimenu(cmenuHand, 'Label', 'Set all UIcontexts', 'Call', 'hand=guidata(gco); set(gco, ''UIContextMenu'', getappdata(hand.axes1, ''cmenuHand''))' );
% 
% function resetUIctx(obj,evt,h,hAxes)
% 	cmenuHand = getappdata(hAxes, 'cmenuHand');
% 	set(h, 'UIContextMenu', cmenuHand)

	% ------------------- Some defaults -----------------------------------
	tol = 0.5;
	do_project = false;         % We'll estimate below if this holds true
	got_arrow = false;
	got_isoc  = false;
	orig_no_mseg = false;		% Flag to know if original file had multiseg strings to store
	line_type = 'AsLine';
	tag = 'polyline';
	struc_vimage = [];
	is_bin = false;				% To flag binary files
	% ---------------------------------------------------------------------

	% ------------------- Parse inputs ------------------------------------
	if (nargin >= 2 && isempty(opt))            % Read a ascii file
		[FileName, PathName, handles] = put_or_get_file(handles, ...
			{'*.dat;*.DAT', 'Data files (*.dat,*.DAT)';'*.*', 'All Files (*.*)'},'Select File','get');
		if isequal(FileName,0),		return,		end
		fname = [PathName FileName];
	elseif (nargin >= 2)		% Read a ascii file of which we already know the name (drag N'drop)
		fname = opt;
		PathName = fileparts(fname);			% We need the 'PathName' below
	end
	if (nargin == 3)
		if (strcmp(opt2, 'AsArrow'))		got_arrow = true;		% This case does not care about 'line_type'
		else								line_type = opt2;
		end
		if (strncmpi(line_type,'isochron',4))
			tag = 'isochron';		fname = [handles.path_data 'isochrons.dat'];
			got_isoc  = true;		PathName = handles.path_data;
		end
	end
	% ---------------------------------------------------------------------

	[bin, n_column, multi_seg, n_headers] = guess_file(fname);
	if isempty(bin) && isempty(n_column) && isempty(multi_seg) && isempty(n_headers)
		errordlg(['Error reading file ' fname],'Error'),	return
	end
	if (isa(bin,'struct') || bin ~= 0)				% ---****** BINARY FILE *******---
		if (isa(bin,'struct'))
			bin = guess_bin(bin.nCols, bin.type);	% Ask user to confirm/modify guessing
		else
			bin = guess_bin(false);					% Ask user what's in file
		end
		if (isempty(bin))		% User quit
			varargout = {};		return
		end
		n_column = bin.nCols;
		multi_seg = 0;		n_headers = 0;		is_bin = true;
	end

	if (n_column == 1 && multi_seg == 0)			% Take it as a file names list
		fid = fopen(fname);
		c = fread(fid,'*char')';	fclose(fid);
		names = strread(c,'%s','delimiter','\n');   clear c fid;
	else
		names = {fname};
	end

	if (handles.no_file)		% Start empty but below we'll find the true data region
		if (ischar(handles.DefLineColor) && handles.DefLineColor(1) == 'w')
			handles.DefLineColor = 'k';		% To not plot a white line over a white background
		end
		XMin = 1e50;			XMax = -1e50;    YMin = 1e50;            YMax = -1e50;

		for (k = 1:numel(names))
			fname = names{k};
			j = strfind(fname,filesep);
			if (isempty(j)),    fname = sprintf('%s%s',PathName, fname);   end		% Need to add path as well 
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (multi_seg)
				[numeric_data, multi_segs_str] = text_read(fname,NaN,n_headers,'>');
			elseif (~is_bin)
				numeric_data = text_read(fname,NaN,n_headers);
			else				% Try luck with a binary file
				fid = fopen(fname);		numeric_data = fread(fid,['*' bin.type]);		fclose(fid);
				numeric_data = reshape(numeric_data,bin.nCols,numel(numeric_data)/bin.nCols)';
			end

			if (~isa(numeric_data,'cell'))			% File was not multi-segment.
				numeric_data = {numeric_data};
				multi_segs_str = {'> Nikles '};		% Need something (>= 8 chars) to not error further down
				orig_no_mseg = true;
			end
			for i=1:length(numeric_data)
				tmpx = numeric_data{i}(:,1);	tmpy = numeric_data{i}(:,2);
				XMin = min(XMin,min(tmpx));		XMax = max(XMax,max(tmpx));
				YMin = min(YMin,min(tmpy));		YMax = max(YMax,max(tmpy));
			end
		end

		if ( multi_seg && strncmp(multi_segs_str{1},'>-:',3) )		% See if we need to swap x<->y
			tmp = XMin;		XMin = YMin;	YMin = tmp;				% Need to swap min/max
			tmp = XMax;		XMax = YMax;	YMax = tmp;
		end

		dx = XMax - XMin;			dy = YMax - YMin;
		if (dx == 0 || dy == 0)
			errordlg('File is has only one point or all XXs are equal or all YYs are equal','Error')
			return
		end
		
		if (XMin > -179.5 && XMax < 359.5 && YMin > -89.5 && YMax < 89.5)
			XMin = XMin - dx / 200;		XMax = XMax + dx / 200;		% Give an extra 0.5% padding margin
			YMin = YMin - dy / 200;		YMax = YMax + dy / 200;
		end
		xx = [XMin XMax];			yy = [YMin YMax];
		region = [xx yy];
		handles.geog = aux_funs('guessGeog',region);

		if (got_isoc)				% We know it's geog (Global Isochrons)
			xx = [-180 180];		yy = [-90 90];
			if (~handles.geog)				handles.geog = 1;
			elseif (handles.geog == 2)		xx = [0 360];
			end
			region = [xx yy];
		end
		mirone('FileNewBgFrame_CB', handles, [region handles.geog])	% Create a background
		hMirFig = handles.figure1;
		drawnow						% Otherwise it takes much longer to plot and other shits
	else							% Reading over an established region
		XYlim = getappdata(handles.axes1,'ThisImageLims');
		xx = XYlim(1:2);			yy = XYlim(3:4);
		if (handles.is_projected && (got_isoc == 1 || handles.defCoordsIn > 0) )
			do_project = true;
		end
		XMin = XYlim(1);			XMax = XYlim(2);		% In case we need this names below for line trimming
	end

	% --------------------------- Main loop over data files ------------------------------
	for (k = 1:numel(names))
		fname = names{k};
		if (handles.no_file && k == 1)			% Rename figure with draged file name
			[pato,barName] = fileparts(fname);
			old_name = get(hMirFig,'Name');		ind = strfind(old_name, '@');
			set(hMirFig,'Name',[barName old_name(ind-1:end)])
		end

		if (~handles.no_file)					% Otherwise we already read it
			j = strfind(fname,filesep);
			if (isempty(j)),    fname = sprintf('%s%s',PathName, fname);   end		% Need to add path as well 
			if (isempty(n_headers)),    n_headers = NaN;    end
			if (multi_seg)
				[numeric_data, multi_segs_str] = text_read(fname,NaN,n_headers,'>');
			elseif (~is_bin)
				numeric_data = text_read(fname,NaN,n_headers);
			else				% Try luck with a binary file
				fid = fopen(fname);		numeric_data = fread(fid,['*' bin.type]);		fclose(fid);
				numeric_data = reshape(numeric_data,bin.nCols,numel(numeric_data)/bin.nCols)';
			end
		end

		if (~isa(numeric_data,'cell'))			% File was not multi-segment. Now pretend it was but with no info
			numeric_data = {numeric_data};
			multi_segs_str = {'> Nikles '};		% Need something in it to not error below
			orig_no_mseg = true;
		end
		n_isoc = 0;     n_segments = length(numeric_data);
		hLine = ones(n_segments,1)*NaN;			% This is the maximum we can have
		n_clear = false(n_segments,1);
		do_patch = false;						% Default to line object

		% Test if conversion into a single, NaN separated, line its wanted 
		if (strncmp(multi_segs_str{1}, '>U_N_I_K', 8))
			for (i = 1:n_segments-1)
				numeric_data{i} = [numeric_data{i}(:,1:2); nan nan];
			end
			numeric_data{1} = cat(1,numeric_data{:});
			% Rip the U_N_I_K identifier
			if (numel(multi_segs_str{1}) > 8)					% We may have line type specifications
				multi_segs_str{1} = ['> ' multi_segs_str{1}(9:end)];
			else
				multi_segs_str{1} = '> ';
			end
			n_segments = 1;				% Pretend we have only one segment

		elseif (strncmpi(multi_segs_str{1}, '>ARROW', 6) || got_arrow)		% ARROW field (the got_arrow can came via varargin)
			if (~got_arrow),	multi_segs_str{1}(2:6) = [];	end			% Rip the ARROW identifier
			got_arrow = true;
			if (n_column < 4)
				errordlg('Files for arrow plot need 4 columns with the traditial (x,y,u,v)','ERROR'),	return
			end
			UV = cell(n_segments,1);
			for (i = 1:n_segments)		% Split the XY & UV columns to be compatible with the other options
				UV{i} = numeric_data{i}(:,3:4);
				numeric_data{i}(:,3:end) = [];
			end
			struc_arrow = struct('spacingChanged',[], 'hQuiver', [], 'hAx', handles.axes1);

		elseif (strncmp(multi_segs_str{1}, '>VIMAGE', 7))
			[z_Vmin, r] = strtok(multi_segs_str{k}(8:end));		z_Vmin = str2double(z_Vmin);
			[z_Vmax, r] = strtok(r);							z_Vmax = str2double(z_Vmax);
			vimage = strtok(r);
			if (isnan(z_Vmin) || isnan(z_Vmax))
				errordlg('Load VIMAGE error. First 2 fields must contain Z_START & Z_END info.','Error'),	return
			end
			if (~ischar(vimage) || ~exist(vimage,'file'))
				errordlg('Load VIMAGE error. Third field must contain an existing picture file name.','Error'),	return
			end
			struc_vimage = struct('z_min', z_Vmin, 'z_max', z_Vmax, 'vimage', vimage);

		elseif (strncmp(multi_segs_str{1}, '>-:', 3))					% File has y,x instead of x,y
			multi_segs_str{1}(2:3) = [];								% Rip the swap -: identifier
			for (i = 1:n_segments)				% Swapp 1st and 2th columns. Do differently would be very complex
				tmp = numeric_data{i}(:,1);
				numeric_data{i}(:,1) = numeric_data{i}(:,2);
				numeric_data{i}(:,2) = tmp;
			end
			clear tmp

		elseif (strncmp(multi_segs_str{1}, '>CLOSE', 6))				% Closed or not, plot a patch
			multi_segs_str{1}(2:6) = [];								% Rip the swap CLOSE identifier
			do_patch = true;

		end

		% If OUT is requested there is nothing left to be done here  
		if (nargout)
			if (orig_no_mseg),		numeric_data = numeric_data{1};		end
			[varargout{1:nargout}] = numeric_data;		return
		end

		drawnow
		for (i = 1:n_segments)
			tmpz = [];
			if (do_project)         % We need to project
				try
					numeric_data{i} = geog2projected_pts(handles,numeric_data{i});
					if (any( isinf(numeric_data{1}(1:min(20,size(numeric_data{1},1)))) ))
						warndlg('Your data was probably already projected. Right-click on the axes frame and uncheck the ''Load files in Geogs'' ','Warning')
					end
				catch
					errordlg(lasterr,'ERROR');    return
				end
			end

			indx = false;	indy = false;			% Default to no need for map clipping
			difes = [numeric_data{i}(1,1)-numeric_data{i}(end,1) numeric_data{i}(1,2)-numeric_data{i}(end,2)];
			if (any(abs(difes) > 1e-4))				% Not a closed polygon
				if (handles.no_file)
					tmpx = numeric_data{i}(:,1);	tmpy = numeric_data{i}(:,2);
				else
					[tmpx,tmpy,indx,indy] = ...		% Get rid of points that are outside the map limits
						aux_funs('in_map_region',handles,numeric_data{i}(:,1),numeric_data{i}(:,2),tol,[xx yy]);
				end
			else
				tmpx = numeric_data{i}(:,1);		tmpy = numeric_data{i}(:,2);
			end
			if (isempty(tmpx)),     n_clear(i) = true;     continue,		end     % Store indexes for clearing vanished segments info
			if ( numel(numeric_data{i}(1,:)) >= 3 )		% If we have a Z column
				tmpz = numeric_data{i}(:,3);
				if (~isempty(indx) || ~isempty(indy))	tmpz(indx) = [];	tmpz(indy) = [];	end	% If needed, clip outside map data			
			end

			[lThick, cor, multi_segs_str{i}] = parseW(multi_segs_str{i}(min(2,numel(multi_segs_str{i})):end)); % First time, we can chop the '>' char
			if (isempty(lThick)),	lThick = handles.DefLineThick;	end		% IF not provided, use default
			if (isempty(cor)),		cor = handles.DefLineColor;		end		%           "

			if (~do_patch || got_arrow)				% Line plottings
				% See if we need to wrap arround the earth roundness discontinuity. Using 0.5 degrees from border. 
				if (handles.geog == 1 && ~do_project && (XMin < -179.5 || XMax > 179.5) )
						[tmpy, tmpx, tmpz] = map_funs('trimwrap', tmpy, tmpx, [-90 90], [XMin XMax], tmpz, 'wrap');
				elseif (handles.geog == 2 && ~do_project && (XMin < 0.5 || XMax > 359.5) )
						[tmpy, tmpx, tmpz] = map_funs('trimwrap', tmpy, tmpx, [-90 90], [XMin XMax], tmpz, 'wrap');
				end

				n_isoc = n_isoc + 1;
				if (~got_arrow)
					switch line_type
						case {'AsLine' 'Isochrons'}
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1,'Linewidth',lThick,...
									'Color',cor,'Tag',tag,'Userdata',n_isoc);
						case 'AsPoint'
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1, 'LineStyle','none', 'Marker','o',...
								'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4,'Tag','Pointpolyline');
							draw_funs(hLine(i),'DrawSymbol')			% Set marker's uicontextmenu (tag is very important)
						case 'AsMaregraph'
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1, 'LineStyle','none', 'Marker','o',...
								'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10,'Tag','Maregraph');
							draw_funs(hLine(i),'DrawSymbol')			% Set marker's uicontextmenu					
						case 'FaultTrace'
							hLine(i) = line('XData',tmpx,'YData',tmpy,'Parent',handles.axes1,'Color',cor,'LineWidth',lThick,'Tag','FaultTrace');
							draw_funs(hLine(i),'line_uicontext')		% Set lines's uicontextmenu
							% Create empty patches that will contain the surface projection of the fault plane
							hp = zeros(1, numel(tmpx)-1);
							for (j = 1:numel(tmpx)-1),	hp(j) = patch('XData', [], 'YData',[]);    end
							setappdata(hLine(i),'PatchHand',hp);
					end
				else
					struc_arrow.color = cor;
					hQuiver = draw_funs([], 'loc_quiver', struc_arrow, tmpx, tmpy, UV{i}(:,1), UV{i}(:,2));
					set(hQuiver,'Tag','Seta','Userdata',n_isoc)
					setappdata(hQuiver(1),'MyHead',hQuiver(2))		% Store the arrows heads handle
					hLine(i) = hQuiver(1);
				end

				if (~isempty(tmpz))		set(hLine(i),'ZData',tmpz');	end
				if (~orig_no_mseg)
					setappdata(hLine(i),'LineInfo',multi_segs_str{i})  % To work with the sessions and will likely replace old mechansim
				end

				% Finish the Vertical image section (if it exists obviously)
				if (~isempty(struc_vimage))
					vimage = getappdata(handles.axes1,'VIMAGE');
					if (isempty(vimage))			% First one
						struc_vimage.hLine = hLine(i);
						setappdata(handles.axes1, 'VIMAGE', struc_vimage)
					else
						struc_vimage.hLine = hLine(i);
						vimage(end+1) = struc_vimage;
						setappdata(handles.axes1, 'VIMAGE', vimage)
					end
				end

			else							% Closed line (parch)
				Fcor = parseG(multi_segs_str{i});
				if (isempty(Fcor)),      Fcor = 'none';   end
				hPat = patch('XData',tmpx,'YData',tmpy,'Parent',handles.axes1,'Linewidth',lThick,'EdgeColor',cor,'FaceColor',Fcor);
				if (~isempty(tmpz))
					set(hPat,'UserData',tmpz');
				end	
				draw_funs(hPat,'line_uicontext')
				n_clear(i) = true;			% Must delete this header info because it only applyies to lines, not patches
			end
		end
		multi_segs_str(n_clear) = [];		% Clear the unused info

		% In case of Lines (and Isocs) uicontexts have not been set yet. Do it now.
		ind = isnan(hLine);    hLine(ind) = [];      % Clear unused rows in hLine (due to over-dimensioning)
		if ( ~isempty(hLine) && (strcmp(line_type, 'AsLine') || strcmp(line_type, 'Isochrons') || got_arrow) )
			if (orig_no_mseg)
				draw_funs(hLine,'line_uicontext')		% Here hLine is actually only a scalar
			else
				draw_funs(hLine,'isochron',multi_segs_str)
			end
		end
	end
	% --------------------- End main loop over files -----------------------------------------

	set(handles.figure1,'pointer','arrow')

	if (handles.no_file)		% Be very carefull, do not trust on the 'geog' estimate donne in show_image (too soon)
		geog = handles.geog;	handles = guidata(handles.figure1);
		handles.geog = geog;	guidata(handles.figure1,handles)
	end

% --------------------------------------------------------------------
function [cor, str2] = parseG(str)
% Parse the STR string in search of color. If not found or error COR = [].
% STR2 is the STR string less the -Gr/g/b part
	cor = [];   str2 = str;
	ind = strfind(str,'-G');
	if (isempty(ind)),      return;     end     % No -G option
	try									% There are so many ways to have it wrong that I won't bother testing
		[strG, rem] = strtok(str(ind:end));
		str2 = [str(1:ind(1)-1) rem];   % Remove the -G<str> from STR

		strG(1:2) = [];					% Remove the '-G' part from strG
		% OK, now 'strG' must contain the color in the r/g/b form
		ind = strfind(strG,'/');
		if (isempty(ind))           % E.G. -G100 form
			cor = eval(['[' strG ']']);
			cor = [cor cor cor] / 255;
		else
			% This the relevant part in num2str. I think it is enough here
			cor = [eval(['[' strG(1:ind(1)-1) ']']) eval(['[' strG(ind(1)+1:ind(2)-1) ']']) eval(['[' strG(ind(2)+1:end) ']'])];
			cor = cor / 255;
		end
		if (any(isnan(cor))),   cor = [];   end
	end

% --------------------------------------------------------------------
function [thick, cor, str2] = parseW(str)
% Parse the STR string in search for a -Wpen. Valid options are -W1,38/130/255 -W3 or -W100/255/255
% If not found or error THICK = [] &/or COR = [].
% STR2 is the STR string less the -W[thick,][r/g/b] part
	thick = [];     cor = [];   str2 = str;
	ind = strfind(str,'-W');
	if (isempty(ind)),      return;     end     % No -W option
	try                                 % There are so many ways to have it wrong that I won't bother testing
		[strW, rem] = strtok(str(ind:end));
		str2 = [str(1:ind(1)-1) rem];   % Remove the -W<str> from STR

		strW(1:2) = [];                 % Remove the '-W' part from strW
		% OK, now 'strW' must contain the pen in the thick,r/g/b form
		ind = strfind(strW,',');
		if (~isempty(ind))          % First thing before the comma must be the line thickness
			thick = eval(['[' strW(1:ind(1)-1) ']']);
			strW = strW(ind(1)+1:end);  % Remove the line thickness part
		else                        % OK, no comma. So we have either a thickness XOR a color
			ind = strfind(strW,'/');
			if (isempty(ind))       % No color. Take it as a thickness
				thick = eval(['[' strW ']']);
			else                    % A color
				cor = [eval(['[' strW(1:ind(1)-1) ']']) eval(['[' strW(ind(1)+1:ind(2)-1) ']']) eval(['[' strW(ind(2)+1:end) ']'])];
				cor = cor / 255;
				if (any(isnan(cor))),   cor = [];   end
				% We are done here. RETURN
				return
			end
		end
		% Come here when -Wt,r/g/b and '-Wt,' have already been riped
		ind = strfind(strW,'/');
		if (~isempty(ind))
			% This the relevant part in num2str. I think it is enough here
			cor = [eval(['[' strW(1:ind(1)-1) ']']) eval(['[' strW(ind(1)+1:ind(2)-1) ']']) eval(['[' strW(ind(2)+1:end) ']'])];
			cor = cor / 255;
		end
		% Notice that we cannot have -W100 represent a color because it would have been interpret above as a line thickness
		if (any(isnan(cor))),   cor = [];   end
	end
