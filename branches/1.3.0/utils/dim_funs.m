function  dim_funs(opt, varargin)
% Functions to deal with -R -I (Ivan the Terrible) stuff of Windows that need it
% There are no outputs. All changes are saved by guidata.

% handles must have these fields
% handles.x_min  handles.x_max handles.y_min  handles.y_max
% handles.x_min_or handles.x_max_or handles.y_min_or  handles.y_max_or
% handles.dms_xinc handles.dms_yinc	handles.one_or_zero
	
switch opt
    case 'xMin'
        edit_x_min(varargin{:})
    case 'xMax'
        edit_x_max(varargin{:})
    case 'yMin'
        edit_y_min(varargin{:})
    case 'yMax'
        edit_y_max(varargin{:})
    case 'xInc'
        edit_x_inc(varargin{:})
    case 'yInc'
        edit_y_inc(varargin{:})
    case 'nRows'
        edit_Nrows(varargin{:})
    case 'nCols'
        edit_Ncols(varargin{:})
end

% --------------------------------------------------------------------------------------
function edit_x_min(hObject, handles)
	x_min_or = handles.x_min_or;
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)            % when dd:mm or dd:mm:ss was given
		% Calculate printing format
		nDigit = round( log10(abs(x_min_or)) );		% Number of digits of the integer part
		frmt = sprintf('%%.%dg',nDigit+8);			% it will be of the type '%.Ng'
        x_min = 0;
        if str2double(val{1}) > 0
            for i = 1:numel(val)   x_min = x_min + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:numel(val)   x_min = x_min - abs(str2double(val{i})) / (60^(i-1));   end
        end
        if (x_min < x_min_or);  set(hObject,'String',sprintf(frmt,x_min_or));    return;     end; 
        handles.x_min = x_min;
        if ~isempty(handles.x_max) && x_min >= handles.x_max
            errordlg('West Longitude >= East Longitude ','Error in Longitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nc = get(handles.edit_Ncols,'String');
        if ~isempty(handles.x_max) && ~isempty(nc)       % x_max and ncols boxes are filled
            % Compute Ncols, but first must recompute x_inc
            x_inc = ivan_the_terrible((handles.x_max - x_min),round(abs(str2double(nc))),1, handles.one_or_zero);
            %xx = floor((handles.x_max - str2double(xx)) / (str2double(get(handles.edit_x_inc,'String')))+0.5) + handles.one_or_zero;
            set(handles.edit_x_inc,'String',sprintf(frmt,x_inc))
            guidata(hObject, handles);
        elseif ~isempty(handles.x_max)      % x_max box is filled but ncol is not, so put to the default (100)
            x_inc = ivan_the_terrible((handles.x_max - x_min),100,1, handles.one_or_zero);
            set(handles.edit_x_inc,'String',sprintf(frmt,x_inc))
            set(handles.edit_Ncols,'String','100')
            guidata(hObject, handles);
        end
	else                % box is empty, so clear also x_inc and ncols
        set(handles.edit_x_inc,'String','');     set(handles.edit_Ncols,'String','');
        set(hObject,'String','');   guidata(hObject, handles);
	end

% --------------------------------------------------------------------------------------
function edit_x_max(hObject, handles)
	x_max_or = handles.x_max_or;
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)
		% Calculate printing format
		nDigit = round( log10(abs(x_max_or)) );		% Number of digits of the integer part
		frmt = sprintf('%%.%dg',nDigit+8);			% it will be of the type '%.Ng'
        x_max = 0;
        if str2double(val{1}) > 0
            for i = 1:numel(val)   x_max = x_max + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:numel(val)   x_max = x_max - abs(str2double(val{i})) / (60^(i-1));   end
        end
        if (x_max > x_max_or);  set(hObject,'String',sprintf(frmt,x_max_or));    return;     end; 
        handles.x_max = x_max;
        if ~isempty(handles.x_min) && x_max <= handles.x_min 
            errordlg('East Longitude <= West Longitude','Error in Longitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nc = get(handles.edit_Ncols,'String');
        if ~isempty(handles.x_min) && ~isempty(nc)       % x_max and ncols boxes are filled
            % Compute Ncols, but first must recompute x_inc
            x_inc = ivan_the_terrible((x_max - handles.x_min),round(abs(str2double(nc))),1, handles.one_or_zero);
            %xx = floor((handles.x_min - str2double(xx)) / (str2double(get(handles.edit_x_inc,'String')))+0.5) + handles.one_or_zero;
            set(handles.edit_x_inc,'String',sprintf(frmt,x_inc))
            guidata(hObject, handles);    
        elseif ~isempty(handles.x_min)      % x_min box is filled but ncol is not, so put to the default (100)
            x_inc = ivan_the_terrible((x_max - handles.x_min),100,1, handles.one_or_zero);
            set(handles.edit_x_inc,'String',sprintf(frmt,x_inc))
            set(handles.edit_Ncols,'String','100')
            guidata(hObject, handles);
        end
	else                % box is empty, so clear also x_inc and ncols
        set(handles.edit_x_inc,'String','');     set(handles.edit_Ncols,'String','');
        set(hObject,'String','');   guidata(hObject, handles);
	end

% --------------------------------------------------------------------
function edit_y_min(hObject, handles)
	% Read value either in decimal or in the dd:mm or dd_mm:ss formats and do some tests
	y_min_or = handles.y_min_or;
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)
		% Calculate printing format
		nDigit = round( log10(abs(y_min_or)) );		% Number of digits of the integer part
		frmt = sprintf('%%.%dg',nDigit+8);			% it will be of the type '%.Ng'
        y_min = 0;
        if str2double(val{1}) > 0
            for i = 1:numel(val)   y_min = y_min + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:numel(val)   y_min = y_min - abs(str2double(val{i})) / (60^(i-1));   end
        end
        if (y_min < y_min_or);  set(hObject,'String',sprintf(frmt,y_min_or));    return;     end; 
        handles.y_min = y_min;
        if ~isempty(handles.y_max) && y_min >= handles.y_max
            errordlg('South Latitude >= North Latitude','Error in Latitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nr = get(handles.edit_Nrows,'String');
        if ~isempty(handles.y_max) && ~isempty(nr)       % y_max and nrows boxes are filled
            % Compute Nrows, but first must recompute y_inc
            y_inc = ivan_the_terrible((handles.y_max - y_min),round(abs(str2double(nr))),1, handles.one_or_zero);
            %xx = floor((handles.y_max - str2double(xx)) / (str2double(get(handles.edit_y_inc,'String')))+0.5) + handles.one_or_zero;
            set(handles.edit_y_inc,'String',sprintf(frmt,y_inc))
            guidata(hObject, handles);
        elseif ~isempty(handles.y_max)      % y_max box is filled but nrows is not, so put to the default (100)
            y_inc = ivan_the_terrible((handles.y_max - y_min),100,1, handles.one_or_zero);
            set(handles.edit_y_inc,'String',sprintf(frmt,y_inc))
            set(handles.edit_Nrows,'String','100')
            guidata(hObject, handles);
        end
	else                % box is empty, so clear also y_inc and nrows
        set(handles.edit_y_inc,'String','');     set(handles.edit_Nrows,'String','');
        set(hObject,'String','');   guidata(hObject, handles);
	end

% --------------------------------------------------------------------
function edit_y_max(hObject, handles)
	y_max_or = handles.y_max_or;
	xx = get(hObject,'String');     val = test_dms(xx);
	if ~isempty(val)
		% Calculate printing format
		nDigit = round( log10(abs(y_max_or)) );		% Number of digits of the integer part
		frmt = sprintf('%%.%dg',nDigit+8);			% it will be of the type '%.Ng'
        y_max = 0;
        if str2double(val{1}) > 0
            for i = 1:numel(val)   y_max = y_max + str2double(val{i}) / (60^(i-1));    end
        else
            for i = 1:numel(val)   y_max = y_max - abs(str2double(val{i})) / (60^(i-1));   end
        end
        if (y_max > y_max_or);  set(hObject,'String',sprintf(frmt,y_max_or));    return;     end; 
        handles.y_max = y_max;
        if ~isempty(handles.y_min) && y_max <= handles.y_min 
            errordlg('North Latitude <= South Latitude','Error in Latitude limits')
            set(hObject,'String','');   guidata(hObject, handles);  return
        end
        nr = get(handles.edit_Nrows,'String');
        if ~isempty(handles.y_min) && ~isempty(nr)       % y_min and nrows boxes are filled
            % Compute Nrows, but first must recompute y_inc
            y_inc = ivan_the_terrible((y_max - handles.y_min),round(abs(str2double(nr))),1, handles.one_or_zero);
            %xx = floor((handles.y_min - str2double(xx)) / (str2double(get(handles.edit_y_inc,'String')))+0.5) + handles.one_or_zero;
            set(handles.edit_y_inc,'String',sprintf(frmt,y_inc))
            guidata(hObject, handles);
        elseif ~isempty(handles.y_min)      % y_min box is filled but nrows is not, so put to the default (100)
            y_inc = ivan_the_terrible((y_max - handles.y_min),100,1, handles.one_or_zero);
            set(handles.edit_y_inc,'String',sprintf(frmt,y_inc))
            set(handles.edit_Nrows,'String','100')
            guidata(hObject, handles);
        end
	else                % This box is empty, so clear also y_inc and nrows
        set(handles.edit_y_inc,'String','');     set(handles.edit_Nrows,'String','');
        set(hObject,'String','');   guidata(hObject, handles);
	end

% --------------------------------------------------------------------
function edit_x_inc(hObject, handles)
	dms = 0;    x_inc = 0;
	xx = get(hObject,'String');     val = test_dms(xx);
	if isempty(val)
        set(hObject, 'String', '');    return
	end

	if numel(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
	for i = 1:numel(val)   x_inc = x_inc + str2double(val{i}) / (60^(i-1));    end
	if ~isempty(handles.x_min) && ~isempty(handles.x_max)
        % Make whatever x_inc given compatible with GMT_grd_RI_verify
        x_inc = ivan_the_terrible((handles.x_max - handles.x_min), x_inc,2, handles.one_or_zero);
        if ~dms         % case of decimal unities
			% Calculate printing format
			nDigit = round( log10(abs(x_inc)) );		% Number of digits of the integer part
			frmt = sprintf('%%.%dg',nDigit+10);			% it will be of the type '%.Ng'
            set(hObject,'String',sprintf(frmt,x_inc))
        else            % inc was in dd:mm or dd:mm:ss format
            ddmm = dec2deg(x_inc);
            set(hObject,'String',ddmm)
        end
        ncol = round((handles.x_max - handles.x_min) / x_inc) + handles.one_or_zero;
        set(handles.edit_Ncols,'String',sprintf('%d',ncol))
	end
	handles.dms_xinc = dms;
	guidata(hObject, handles);
	if isempty(get(handles.edit_y_inc,'String'))     set(handles.edit_y_inc,'String',xx);    end

% --------------------------------------------------------------------
function edit_Ncols(hObject, handles)
	xx = get(hObject,'String');
	if ~isempty(get(handles.edit_x_min,'String')) && ~isempty(get(handles.edit_x_max,'String')) && ...
            ~isempty(get(handles.edit_x_inc,'String')) && ~isempty(xx)
        x_inc = ivan_the_terrible((handles.x_max - handles.x_min),round(abs(str2double(xx))),1, handles.one_or_zero);
        if handles.dms_xinc        % x_inc was given in dd:mm:ss format
            ddmm = dec2deg(x_inc);
            set(handles.edit_x_inc,'String',ddmm)
        else                    % x_inc was given in decimal format
			% Calculate printing format
			nDigit = round( log10(abs(x_inc)) );		% Number of digits of the integer part
			frmt = sprintf('%%.%dg',nDigit+10);			% it will be of the type '%.Ng'
            set(handles.edit_x_inc,'String',sprintf(frmt,x_inc));
        end
        guidata(hObject, handles);
	end

% --------------------------------------------------------------------
function edit_y_inc(hObject, handles)
	dms = 0;    y_inc = 0;
	xx = get(hObject,'String');     val = test_dms(xx);
	if isempty(val)
        set(hObject, 'String', '');    return
	end

	if numel(val) > 1    dms = 1;      end         % inc given in dd:mm or dd:mm:ss format
	for i = 1:numel(val)   y_inc = y_inc + str2double(val{i}) / (60^(i-1));    end
	if ~isempty(handles.y_min) && ~isempty(handles.y_max)
        % Make whatever y_inc given compatible with GMT_grd_RI_verify
        y_inc = ivan_the_terrible((handles.y_max - handles.y_min), y_inc, 2, handles.one_or_zero);
        if ~dms         % case of decimal unities
			% Calculate printing format
			nDigit = round( log10(abs(y_inc)) );		% Number of digits of the integer part
			frmt = sprintf('%%.%dg',nDigit+10);			% it will be of the type '%.Ng'
            set(hObject,'String',sprintf(frmt,y_inc))
        else            % inc was in dd:mm or dd:mm:ss format
            ddmm = dec2deg(y_inc);
            set(hObject,'String',ddmm)
        end
        nrow = round((handles.y_max - handles.y_min) / y_inc) + handles.one_or_zero;
        set(handles.edit_Nrows,'String',sprintf('%d',nrow))
	end
	handles.dms_yinc = dms;
	guidata(hObject, handles);

% --------------------------------------------------------------------
function edit_Nrows(hObject, handles)
	xx = get(hObject,'String');
	if ~isempty(get(handles.edit_y_min,'String')) && ~isempty(get(handles.edit_y_max,'String')) && ...
            ~isempty(get(handles.edit_y_inc,'String')) && ~isempty(xx)
        y_inc = ivan_the_terrible((handles.y_max - handles.y_min),round(abs(str2double(xx))),1, handles.one_or_zero);
        if handles.dms_yinc        % y_inc was given in dd:mm:ss format
            ddmm = dec2deg(y_inc);
            set(handles.edit_y_inc,'String',ddmm)
        else                    % y_inc was given in decimal format
			% Calculate printing format
			nDigit = round( log10(abs(y_inc)) );		% Number of digits of the integer part
			frmt = sprintf('%%.%dg',nDigit+10);			% it will be of the type '%.Ng'
            set(handles.edit_y_inc,'String',sprintf(frmt,y_inc));
        end
        guidata(hObject, handles);
	end
