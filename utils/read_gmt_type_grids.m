function [handles, X, Y, Z, head, misc] = read_gmt_type_grids(handles,fullname,opt)
    % OPT indicates that only the grid info is outputed.
    % If it is OPT = 'hdr' outputs info in the struct format, else outputs in the head format
    
    infoOnly = 0;
    if (nargin == 3),   infoOnly = 1;    end
    
	X = [];		Y = [];		Z = [];		head = [];		misc = [];
	[fid, msg] = fopen(fullname, 'r');
	if (fid < 0),   errordlg([fullname ': ' msg],'ERROR');
		return
	end

	% Because GMT and Surfer share the .grd extension, find out which grid kind we are dealing with
	ID = fread(fid,4,'*char');      ID = ID';      fclose(fid);

	switch upper(ID)
		case 'CDF',		tipo = 'CDF';
		case {'DSBB' 'DSRB'},	tipo = 'SRF_BIN';
			ID = 'CDF';			% Just a dirty trick
		case 'DSAA',    tipo = 'SRF_ASCII';
		case 'MODE',    tipo = 'ENCOM';
		case 'NLIG',    tipo = 'MAN_ASCII';
	end

	% See if the grid is on one of the OTHER (non netCDF & non Surfer) formats that GMT recognizes
	if (~strcmp(ID(1:3),'CDF') && ~(tipo(1) == 'S' || tipo(1) == 'E' || tipo(1) == 'M') )
		str = ['grdinfo ' fullname];
		[PATH,FNAME,EXT] = fileparts(fullname);
		[s,att] = mat_lyies(str,[handles.path_tmp FNAME '.' EXT '.info']);
		if ~(isequal(s,0))          % File could not be read
			errordlg([fullname ' : Is not a grid that GMT can read!'],'ERROR');
			return
		end
	end

	if (~infoOnly)
		[handles, X, Y, Z, head, misc] = read_grid(handles,fullname,tipo);
	elseif ( strmatch(tipo,{'GMT' 'SRF_BIN'}) )
		if (opt(1) == 's')          % Get the info on the struct form
			X = grdinfo_m(fullname,'hdr_struct');       % Output goes in the second arg
		else                        % Get the info on the vector form
			X = grdinfo_m(fullname,'silent');
		end
	else
		errordlg([fullname ' : Is not a GMT or binary Surfer grid!'],'ERROR');
		return
	end

% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function [handles, X, Y, Z, head, misc] = read_grid(handles,fullname,tipo)

if (isfield(handles,'ForceInsitu'))        % Other GUI windows may not know about 'ForceInsitu'
	if (handles.ForceInsitu),   opt_I = 'insitu';    % Use only in desperate cases.
	else                        opt_I = ' ';
	end
else
    opt_I = ' ';
end
	X = [];     Y = [];     Z = [];     head = [];		misc = [];		% MISC is used only by nc_io

	if (~strcmp(tipo,'GMT'))        % GMT files are open by the GMT machinerie
		[fid, msg] = fopen(fullname, 'r');
		if (fid < 0),   errordlg([fullname ': ' msg],'ERROR');  return,		end
	end

if (strcmp(tipo,'CDF'))
	try				% Use the new nc_io()
		[X, Y, Z, head, misc] = nc_io(fullname, 'r');
		if (isa(Z,'int16')),		handles.was_int16 = 1;
		elseif (isa(Z,'single')),	handles.have_nans = grdutils(Z,'-N');
		elseif (isa(Z,'double')),	Z = single(Z);		% The HORRRRRRRRROOOOOOORRRRR
		end
	catch			% If it have failed try GMT
		str = sprintf('First attempt to load netCDF file failed because ... \n\n %s\n\n Trying now with GMT mex', lasterr);
		warndlg(str,'Info')
    	[X, Y, Z, head] = grdread_m(fullname,'single',opt_I);
    	handles.have_nans = grdutils(Z,'-N');
    	if (head(10) == 2 || head(10) == 8 || head(10) == 16),   handles.was_int16 = 1;  end     % New output from grdread_m
		head(10) = [];
	end
    if (head(7))            % Convert to grid registration
        head(1) = head(1) + head(8) / 2;        head(2) = head(2) - head(8) / 2;
        head(3) = head(3) + head(9) / 2;        head(4) = head(4) - head(9) / 2;
        head(7) = 0;
    end
    handles.grdname = fullname;		handles.image_type = 1;		handles.computed_grid = 0;
elseif (strcmp(tipo,'SRF_BIN'))
	[X, Y, Z, head] = grdread_m(fullname,'single',opt_I);
	handles.have_nans = grdutils(Z,'-N');
    handles.grdname = fullname;		handles.image_type = 1;		handles.computed_grid = 0;
elseif (strcmp(tipo,'SRF_ASCII'))	% Pretend that its a internaly computed grid (no reload)
    s = fgetl(fid);
    n_col_row = fscanf(fid,'%f',2);     x_min_max = fscanf(fid,'%f',2);
    y_min_max = fscanf(fid,'%f',2);     z_min_max = fscanf(fid,'%f',2);
    X = linspace(x_min_max(1),x_min_max(2),n_col_row(1));
    Y = linspace(y_min_max(1),y_min_max(2),n_col_row(2));
    Z = single(fscanf(fid,'%f',inf));   fclose(fid);
    Z = reshape(Z,n_col_row')';
    Z(Z >= 1e38) = NaN;
    handles.have_nans = grdutils(Z,'-N');
    dx = diff(x_min_max) / (n_col_row(1) - 1);
    dy = diff(y_min_max) / (n_col_row(2) - 1);
    head = [x_min_max' y_min_max' z_min_max' 0 dx dy];
    handles.image_type = 1;     handles.computed_grid = 1;    handles.grdname = [];
elseif (strcmp(tipo,'ENCOM'))       % Pretend that its a GMT grid
    ID = fread(fid,180,'*char');        % We don't use this header info, so strip it
    no_val = fread(fid,1,'float32');
    ID = fread(fid,4,'*char');      n_rows = fread(fid,1,'float32');    % ROWS FLAG
    ID = fread(fid,4,'*char');      n_cols = fread(fid,1,'float32');    % COLS FLAG
    ID = fread(fid,4,'*char');      x_min = fread(fid,1,'float32');     % XORIG FLAG
    ID = fread(fid,4,'*char');      y_min = fread(fid,1,'float32');     % YORIG FLAG
    ID = fread(fid,4,'*char');      dx = fread(fid,1,'float32');        % DX FLAG
    ID = fread(fid,4,'*char');      dy = fread(fid,1,'float32');        % DY FLAG
    ID = fread(fid,4,'*char');      rot = fread(fid,1,'float32');       % DEGR FLAG (I'll use it one day)
    Z = single(fread(fid,n_rows*n_cols,'float32'));    fclose(fid);
    Z = reshape(Z, n_cols, n_rows)';
    Z(Z == no_val) = NaN;
    [zzz] = grdutils(Z,'-L+');  z_min = zzz(1);     z_max = zzz(2);     handles.have_nans = zzz(3); clear zzz;
    x_max = x_min + (n_cols-1) * dx;        y_max = y_min + (n_rows-1) * dy;
    X = linspace(x_min,x_max,n_cols);       Y = linspace(y_min,y_max,n_rows);
    head = [x_min x_max y_min y_max z_min z_max 0 dx dy];
    handles.image_type = 1;     handles.computed_grid = 1;    handles.grdname = [];
elseif (strcmp(tipo,'MAN_ASCII'))
    h1 = fgetl(fid);    h2 = fgetl(fid);    h3 = fgetl(fid);    h4 = fgetl(fid);    h5 = fgetl(fid);
    n_rows = str2double(h1(6:10));          n_cols = str2double(h1(16:20));
    y_inc = str2double(h2(6:10));           x_inc = str2double(h2(16:20));
    no_val = str2double(h2(25:31));         azim = str2double(h2(35:40));
    x_min = str2double(h2(46:60));          y_max = str2double(h2(66:80));
    Z = single(fscanf(fid,'%f',inf));       fclose(fid);
    Z = flipud(reshape(Z,n_rows,n_cols));
    Z(Z == no_val) = NaN;
    x_max = x_min + (n_cols-1) * x_inc;     y_min = y_max - (n_rows-1) * y_inc;
    if (azim ~= 0)
        Z = transform_fun('imrotate',Z,azim,'bilinear','loose');
        n_cols = size(Z,2);      n_rows = size(Z,1);
        azim_rad = azim * pi / 180;
        rot = [cos(azim_rad) sin(azim_rad); ...
                -sin(azim_rad) cos(azim_rad)];
        % Compute the des-rotated grid limits
        UL = rot * [x_min; y_max];                      % Upper Left corner
        UR = rot * [x_max; y_max];                      % Upper Right corner
        LR = rot * [x_max; y_min];                      % Lower Right  corner
        LL = rot * [x_min; y_min];                      % Lower Left  corner
        x_min = min([UL(1) UR(1) LR(1) LL(1)]);      x_max = max([UL(1) UR(1) LR(1) LL(1)]);
        y_min = min([UL(2) UR(2) LR(2) LL(2)]);      y_max = max([UL(2) UR(2) LR(2) LL(2)]);
        x_inc = (x_max - x_min) / (n_cols - 1);         % We need to recompute those
        y_inc = (y_max - y_min) / (n_rows - 1);
    end
    [zzz] = grdutils(Z,'-L+');  z_min = zzz(1);     z_max = zzz(2);     handles.have_nans = zzz(3); clear zzz;
    X = linspace(x_min,x_max,n_cols);       Y = linspace(y_min,y_max,n_rows);
    head = [x_min x_max y_min y_max z_min z_max 0 x_inc y_inc];
    handles.image_type = 1;     handles.computed_grid = 1;    handles.grdname = [];
end
