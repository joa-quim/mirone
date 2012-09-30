function [handles, X, Y, Z, head, misc] = read_gmt_type_grids(handles,fullname,opt)
% OPT indicates that only the grid info is outputed.
% MISC - which exists only when nc_io was used - is a struct with:
%		'desc', 'title', 'history', 'srsWKT', 'strPROJ4' fields
% If OPT == 'hdr' outputs info in the struct format, else outputs in the head format
%
% The HANDLES fields 'grdname', 'image_type', 'have_nans' and 'computed_grid' are reset
% and 'was_int16' may or not
%
% When used to read netCDF grids HANDLES can be []. Useful to use this function as a standalone

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

% $Id$

    infoOnly = 0;
    if (nargin == 3),   infoOnly = 1;    end
    
	X = [];		Y = [];		Z = [];		head = [];		misc = [];
	[fid, msg] = fopen(fullname, 'r');
	if (fid < 0),   errordlg([fullname ': ' msg],'ERROR');
		return
	end

	% Because GMT and Surfer share the .grd extension, find out which grid kind we are dealing with
	ID = fread(fid,3,'*char');      ID = ID';      fclose(fid);
	tipo = 'U';

	switch upper(ID)
		case 'CDF',		tipo = 'CDF';
		case 'DSB',		tipo = 'SRF6';		% DSBB
		case 'DSR',		tipo = 'SRF7';		% DSRB
			ID = 'CDF';			% Just a dirty trick
		case 'DSA',    tipo = 'SRF_ASCII';		% DSAA
		case 'MOD',    tipo = 'ENCOM';			% MODE
		case 'NLI',    tipo = 'MAN_ASCII';		% NLIG
		case char([2 0 0])		% Risky wild guess that grid is GEOSOFT 2-bytes
			ID = 'CDF';			% Again the dirty trick
			tipo = 'GSOFT';
	end
	if (strcmp(ID(2:3),'HD'))
		tipo = 'CDF';	ID = 'CDF';		% And the trick again
	end

	% See if the grid is on one of the OTHER (non netCDF & non Surfer) formats that GMT recognizes
	if (~strcmp(ID(1:3),'CDF') && ~(tipo(1) == 'S' || tipo(1) == 'E' || tipo(1) == 'M') )
		str = ['grdinfo ' fullname];
		[PATH,FNAME,EXT] = fileparts(fullname);
		s = mat_lyies(str,[handles.path_tmp FNAME EXT '.info']);
		if ~(isequal(s,0))          % File could not be read
			errordlg([fullname ' : Is not a grid that GMT can read!'],'ERROR');
			return
		end
	end

	if (~infoOnly)
		[handles, X, Y, Z, head, misc] = read_grid(handles,fullname,tipo);
	elseif ( any(strcmp(tipo,{'CDF' 'SRF6' 'SRF7'})) )
		if (opt(1) == 's')          % Get the info on the struct form
			X = grdinfo_m(fullname,'hdr_struct');       % Output goes in the second arg
		else                        % Get the info on the vector form
			X = grdinfo_m(fullname,'silent');
		end
	else
		errordlg([fullname ' : Is not a GMT or binary Surfer grid!'],'ERROR');
		return
	end

	if ( (abs(head(5)) < 1e-12) && (abs(head(6)) < 1e-12) )		% Badly behaved grid. No min/max info
		if (handles.have_nans),		zz = grdutils(Z,'-L');
		else						zz = [min(Z(:)) max(Z(:))];
		end
		head(5:6) = double(zz(1:2));
	end
% _________________________________________________________________________________________________	
% -*-*-*-*-*-*-$-$-$-$-$-$-#-#-#-#-#-#-%-%-%-%-%-%-@-@-@-@-@-@-(-)-(-)-(-)-&-&-&-&-&-&-{-}-{-}-{-}-
function [handles, X, Y, Z, head, misc] = read_grid(handles,fullname,tipo)

	X = [];     Y = [];     Z = [];     head = [];		opt_I = ' ';	misc = [];		% MISC is used only by nc_io
	if (~isempty(handles) && isfield(handles,'ForceInsitu'))		% Other GUI windows may not know about 'ForceInsitu'
		if (handles.ForceInsitu),   opt_I = 'insitu';		end		% Use only in desperate cases.
	end

	if (~strcmp(tipo,'CDF'))        % GMT files are open by the GMT machinerie
		[fid, msg] = fopen(fullname, 'r');
		if (fid < 0),   errordlg([fullname ': ' msg],'ERROR');  return,		end
	end

if (strcmp(tipo,'CDF'))
	try				% Use the new nc_io()
		[X, Y, Z, head, misc] = nc_io(fullname, 'r');
		if (isa(Z,'int16')),		handles.was_int16 = 1;
		elseif (isa(Z,'single')),	handles.have_nans = grdutils(Z,'-N');
		elseif (isa(Z,'double')),	Z = single(Z);		% The HORRRRRRRRROOOOOOOOOORRRRR
		end
	catch			% If it have failed try GMT
		str = sprintf(['First attempt to load netCDF file failed because ... \n\n\n %s\n\n\n       Trying now with GMT mex ...' ...
		'\n\nBTW. Please inform me about this error so that I can try to correct it.\nThanks.'], lasterr);
		warndlg(str,'Info')
		if ( ~isempty(strfind(lasterr, 'Out of memory')) ) % If its a memory problem, no use to insist
			error(lasterr)
		end
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
elseif (strcmp(tipo,'SRF6'))
	fread(fid,4,'*char');
	n_cols = fread(fid,1,'int16');			n_rows = fread(fid,1,'int16');
	head = (fread(fid,6,'double'))';
	Z = fread(fid,n_rows*n_cols,'*float32');	fclose(fid);
    Z = reshape(Z, n_cols, n_rows)';
	handles.have_nans = grdutils(Z,'-N');		% Check for the degenerated case where a Surfer grid has native NaNs
	ind = (Z >= 1e38);
	if (any(ind(:)))
    	Z(ind) = NaN;    handles.have_nans = 1;
	end
	head(7:9) = [0 diff(head(1:2))/(n_cols - 1) diff(head(3:4))/(n_rows - 1)];
    X = linspace(head(1),head(2),n_cols);    Y = linspace(head(3),head(4),n_rows);
elseif ( strcmp(tipo,'SRF7') || (tipo(1) == 'U') )
	[X, Y, Z, head] = grdread_m(fullname,'single',opt_I);
	handles.have_nans = grdutils(Z,'-N');
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
    Z = fread(fid,n_rows*n_cols,'*float32');    fclose(fid);
    Z = reshape(Z, n_cols, n_rows)';
    Z(Z == no_val) = NaN;
    [zzz] = grdutils(Z,'-L+');  z_min = zzz(1);     z_max = zzz(2);     handles.have_nans = zzz(3); clear zzz;
    x_max = x_min + (n_cols-1) * dx;        y_max = y_min + (n_rows-1) * dy;
    X = linspace(x_min,x_max,n_cols);       Y = linspace(y_min,y_max,n_rows);
    head = [x_min x_max y_min y_max z_min z_max 0 dx dy];
elseif (strcmp(tipo,'MAN_ASCII'))
    h1 = fgetl(fid);    h2 = fgetl(fid);    fgetl(fid);    fgetl(fid);    fgetl(fid);
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
elseif (strcmp(tipo,'GSOFT'))
	fclose(fid);
	[X,Y,Z,head, azim] = read_geosoft(fullname);
	if (azim ~= 0)
		warndlg('This GEOSOFT grids is rotated but currently I do not take that into account, sorry', 'Warning')
	end
    zzz = grdutils(Z,'-L+');     head(5:6) = zzz(1:2);	handles.have_nans = zzz(3);
end

handles.grdname = fullname;		handles.image_type = 1;		handles.computed_grid = 0;

% ------------------------------------------------------------------------------------
function [X,Y,Z,head,rot] = read_geosoft(filename)
% Read a Geosoft 2-byte binary grid
%
% Based on function getgrd2 of Adam O'Neill
% http://www.mathworks.com/matlabcentral/fileexchange/1527-getgrd2/content/getgrd2.m

	grdhead = grdheadset;
	grdelem = grdelemset;
	% Open data file and get necessary header info
	fid = fopen(filename,'r','n');

	fseek(fid,grdhead.ne,'bof');		n_col = fread(fid,1,grdelem.ne);		% No. of columns
	fseek(fid,grdhead.nv,'bof');		n_row = fread(fid,1,grdelem.nv);		% No. of rows
	fseek(fid,grdhead.kx,'bof');		kx = fread(fid,1,grdelem.kx);			% Orientation sense
	% kx = 1, vectors run left-right
	% kx = -1, vectors run up-down

	% Geographic information
	fseek(fid,grdhead.de,'bof');		dx = fread(fid,1,grdelem.de);
	fseek(fid,grdhead.dv,'bof');		dy = fread(fid,1,grdelem.dv);
	fseek(fid,grdhead.x0,'bof');		x0 = fread(fid,1,grdelem.x0);
	fseek(fid,grdhead.y0,'bof');		y0 = fread(fid,1,grdelem.y0);

	fseek(fid,grdhead.rot,'bof');		rot = fread(fid,1,grdelem.rot);

	% Data (Z) scaling
	fseek(fid,grdhead.zbase,'bof');		zbase = fread(fid,1,grdelem.zbase);
	fseek(fid,grdhead.zmult,'bof');		zmult = fread(fid,1,grdelem.zmult);

	% Go to start of data stream and read into array
	fseek(fid,grdhead.totallength,'bof');
	Z = fread(fid,[n_col n_row], grdelem.datasize); 
	fclose(fid);

	% Make the output array
	if (kx == 1),		Z = Z';
	elseif (kx ~= -1)	error('GEOSOFT Grid: kx not valid')
	end

	ind = (Z == -32767);		% Not sure is they do it in Geosoft
	if (~any(ind)),		ind = false;	end		% If no novalues free memory right away
	Z = single(Z);
	if (ind),	Z(ind) = NaN;	end
	cvlib_mex('CvtScale',Z, 1./zmult, zbase)

	X = linspace(x0, x0 + dx * (n_col - 1), n_col);
	Y = linspace(y0, y0 + dy * (n_row - 1), n_row);
	head = [x0 X(end) y0 Y(end) 0 0 0 dx dy];		% z_min & z_max will be computed by calling f

% ------------------------------------------------------------------------------------
function grdelem = grdelemset()
grdelem = struct(...
	'es','int32',... % Data storage
	'sf','int32',...
	'ne','int32',...
	'nv','int32',...
	'kx','int32',...
	'de','float64',... % Geographic information
	'dv','float64',...
	'x0','float64',...
	'y0','float64',...
	'rot','float64',...
	'zbase','float64',... % Data (Z) scaling
	'zmult','float64',...
	'label','48*char',...     %char*48
	'mapno','16*char',...     %char*16
	'proj','int32',... % Optional parameters
	'unitx','int32',...
	'unity','int32',...
	'unitz','int32',...
	'nvpts','int32',...
	'izmin','int32',...
	'izmax','int32',...
	'izmed','int32',...
	'izmea','int32',...
	'zvar','float64',...
	'prcs','int32',...
	'datasize','int16'); % Geosoft 2 byte grids ONLY

function [grdhead] = grdheadset()
grdhead = struct(...
	'es',0,...
	'sf',4,...
	'ne',8,...
	'nv',12,...
	'kx',16,...
	'de',20,...
	'dv',28,...
	'x0',36,...
	'y0',44,...
	'rot',52,...
	'zbase',60,...
	'zmult',68,...
	'label',76,...
	'mapno',124,...
	'proj',140,...
	'unitx',144,...
	'unity',148,...
	'unitz',152,...
	'nvpts',156,...
	'izmin',160,...
	'izmax',164,...
	'izmed',168,...
	'izmea',172,...
	'zvar',176,...
	'prcs',184,...
	'totallength',512); % Geosoft 2 byte grids ONLY
