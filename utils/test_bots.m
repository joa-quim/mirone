function  test_bots(opt,varargin)

	global gmt_ver;		gmt_ver = 5;

	if (nargin)
		switch opt
			case 'Illum'
				test_Illum
			case 'xyz'
				test_xyz
			case 'writeascii'
				writeascii
			case 'loadPall'
				loadPall
			case 'grdfilter'
				grdfilter
			case 'grdinfo'
				grdinfo
			case 'grdlandmask'
				grdlandmask
			case 'grdsample'
				grdsample
			case 'grdtrend'
				grdtrend
			case 'grdread'
				grdread
			case 'coasts'
				coasts
			case 'gmtlist'
				gmtlist
			case 'implant'
				implant
			case 'fill_poly_hole'
				fill_poly_hole
		end
	else
		test_Illum
		test_xyz
		writeascii
		loadPall
		grdfilter
		grdinfo
		grdlandmask
		grdsample
		grdtrend
		grdread
		coasts
		gmtlist
		implant
	end

% -----------------
function test_Illum
	h = mirone('c:\SVN\mironeWC\swath_grid.grd');
	handles = guidata(h);
	luz = struct('azim',0,'elev',30,'ambient',0.55,'diffuse',0.6,'specular',0.4,'shine',10);
	try		mirone('ImageIllum',luz, handles, 'grdgrad_class'),		pause(1)
	catch,	disp('FAIL: a ilum com o grdgrad')
	end
	try		mirone('ImageIllum',luz, handles, 'grdgrad_lamb'),		pause(1)
	catch,	disp('FAIL: a ilum com o grad lambert')
	end
	try		mirone('ImageIllum',luz, handles, 'grdgrad_peuck'),		pause(1)
	catch,	disp('FAIL: a ilum com o peucker')
	end
	try		mirone('ImageIllum',luz, handles, 'lambertian'),		pause(1)
	catch,	disp('FAIL: a ilum com o lambert')
	end
	try		mirone('ImageIllum',luz, handles, 'manip'),				pause(1)
	catch,	disp('FAIL: a ilum com o manip'),	disp(lasterr)
	end
	try		mirone('ImageIllum',luz, handles, 'hill'),				pause(1)
	catch,	disp(['FAIL: a ilum com o hill -> ' lasterr])
	end
	luz.azim = [0 120 240];		luz.mercedes_type = 1;
	try		mirone('ImageIllumFalseColor',luz, handles)
	catch,	disp(['FAIL: a ilum com o false color -> ' lasterr])
	end
	delete(h)

% ----------------------------
function test_xyz
% Test import dat via load_xyz
	x = rand(2, 30)*100;		fname = 'lixo_test.dat';
	fid = fopen(fname,'w');
	fprintf(fid,'%f %f\n', x);
	fclose(fid);
	h = mirone(fname);		pause(1)
	delete(findobj(h,'type','line'));
	handMir = guidata(h);
	load_xyz(handMir, fname, 'AsPoint'),	pause(1)
	delete(h)
	
	% Arrows
	fid = fopen(fname,'w');
	fprintf(fid,'>ARROW\n');
	fprintf(fid,'%f %f %f %f\n', [-2 -1 -0.06229 -0.02873; -2 -0.85 -0.0822 -0.0299; -2 -0.7 -0.1037 -0.03094;
	-2 -0.55 -0.1251 -0.02925; -2 -0.40 -0.1443 -0.02447; -2 -0.25 -0.1591 -0.01684; -2 -0.10 -0.1676 -0.00709]);
	fclose(fid);
	h = mirone(fname);		pause(1)
	delete(findobj(h,'type','line'));
	handMir = guidata(h);
	load_xyz(handMir, fname),	pause(1)
	delete(h)

	% The girl
	mir_dirs = getappdata(0,'MIRONE_DIRS');
	girl = [mir_dirs.home_dir filesep 'data' filesep 'gp_girl.dat'];
	h = mirone(girl);		pause(1)
	delete(h)

	% A coloured patch
	fid = fopen(fname,'w');
	fprintf(fid,'> -G200/134/34\n');
	fprintf(fid,'%f %f\n', [0 0; 0 1; 1 1; 1 0; 0 0]');
	fclose(fid);
	h = mirone(fname);		pause(1)
	delete(h)

	builtin('delete',fname);

% ----------------------------
function writeascii
% Test save ascii file via double2ascii. Use integers so that we can do numeric comparisons

	fname = 'lixo_test.dat';
	xyz = [(1:4)' (41:44)' (61:64)'];
	double2ascii(fname, xyz);
	zzz = load_xyz([], fname);
	difa = xyz - zzz;
	if (any(difa(:))),		disp('FAIL: 4x3 array without format'),		end
	double2ascii(fname, xyz, '%d %d %f');
	zzz = load_xyz([], fname);
	difa = xyz - zzz;
	if (any(difa(:))),		disp('FAIL: 4x3 array with variable format'),		end
	% CELLS
	double2ascii(fname, {xyz});
	zzz = load_xyz([], fname);
	difa = xyz - zzz;
	if (any(difa(:))),		disp('FAIL: Cell with one 4x3 without format and non-multisegment'),	end
	double2ascii(fname, {xyz}, '%d');
	zzz = load_xyz([], fname);
	difa = xyz - zzz;
	if (any(difa(:))),		disp('FAIL: Cell with one 4x3 with unique format and non-multisegment'),	end
	% CELLS MULTISEGS
	multistr = {'> a'; '> b'};
	double2ascii(fname, {xyz xyz}, '%f', multistr);
	[zzz, sss] = load_xyz([], fname);
	if (~isequal(zzz{1}, xyz)),		disp('FAIL: Cell with 2 4x3 with unique format and multiseg sent in'),	end
	if (~isequal(sss, multistr)),	disp('FAIL: Multiseg string of a 2 cell array and multiseg sent in'),	end
	% NaNs
	double2ascii(fname, [xyz; ones(1,size(xyz,2))*NaN; xyz]);
	zzz = load_xyz([], fname);
	if (numel(zzz(~isnan(zzz))) ~= numel(xyz)*2)
		disp('FAIL: With NaNs and not multiseg')
	end
	double2ascii(fname, [xyz; ones(1,size(xyz,2))*NaN; xyz],'%f','multi');
	zzz = load_xyz([], fname);
	if (~isequal(zzz{1}, xyz) && ~isequal(zzz{2}, xyz)),		disp('FAIL: With NaNs and multiseg'),	end

	builtin('delete',fname);

% ----------------------------
function loadPall
% Load a GMT cpt file
	% Note, if use mola.cpt in 5.2 we get a crash. INVESTIGATE IT
	color_palettes('C:\progs_cygw\GMTdev\gmt4\share\cpt\GMT_ocean.cpt');	% Should be a file in a test dir
	pause(0.2)
	delete(findobj('type','figure','Name','Color Palettes'))

% ----------------------------
function grdfilter
% ...
	h = mirone('swath_grid.grd');
	handles = guidata(h);
	[X,Y,Z,head] = load_grd(handles);
	hf = grdfilter_mir(handles);
	Z = c_grdfilter(Z, head, '-Fb20k', '-D1');	%#ok
	pause(0.3);		delete(h);		delete(hf)

% ----------------------------
function grdinfo
% ...
	fname = 'swath_grid.grd';
	h = mirone(fname);
	grid_info(guidata(h))
	showBak = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on');
	hFig = findobj(get(0,'Children'),'flat', 'tag','Wdmsgfig');
	pause(0.3);
	delete([hFig h])
	set(0,'ShowHiddenHandles',showBak);

% ----------------------------
function grdlandmask
% ...
	h = grdlandmask_win;
	[mask,head,X,Y] = c_grdlandmask('-R-180/180/-90/90', '-I0.5/0.5');
	pause(0.3);		delete(h);

% ----------------------------
function grdsample
% ...
	h = mirone('swath_grid.grd');
	handles = guidata(h);
	[X,Y,Z,head] = load_grd(handles);
	hf = grdsample_mir(handles);
	Z = c_grdsample(Z, head, '-N191/131');	%#ok
	pause(0.3);		delete(h);		delete(hf)

% ----------------------------
function grdtrend
% ...
	h = mirone('swath_grid.grd');
	handles = guidata(h);
	[X,Y,Z,head] = load_grd(handles);
	hf = grdtrend_mir(handles);
	Z = c_grdtrend(Z, head, '-T', '-N3');	%#ok
	pause(0.3);		delete(h);		delete(hf)

% ----------------------------
function grdread
% ...
	fname = 'swath_grid.grd';
	[X, Y, Z, head] = c_grdread(fname,'single');
	h = mirone(Z);
	pause(0.3);		delete(h);

% ----------------------------
function coasts
% ...
	opt_res = '-Di';	opt_N = '-Na';		opt_I = '-Ia';
	opt_R = '-R-10/10/30/50';
	coast = c_shoredump(opt_R,opt_res,'-A1/1/1');
	boundaries = c_shoredump(opt_R,opt_N,opt_res);
	rivers = c_shoredump(opt_R,opt_I,opt_res);
	h = figure; hold on
	plot(coast(1,:), coast(2,:))
	plot(boundaries(1,:), boundaries(2,:))
	plot(rivers(1,:), rivers(2,:))
	pause(0.5);		delete(h);

% ----------------------------
function gmtlist
% ...
	gmtedit('C:\j\cd120\revisao\cd120.gmt');

% ----------------------------
function implant
% ...

	[Z, hdrStruct] = gen_UMF2d(1.8, 0.05, 0.9, 1000);		% Pretend it's geogs
	Zf = c_grdfilter(Z,hdrStruct.head,'-Fb20','-D1');		% Filter boxcar 20 km
	hand1.Z = Zf(1:5:end, 1:5:end);

 	X = 1:5:1000; 	Y = 1:5:1000;
	hand1.head = [1 X(end) 1 Y(end) hdrStruct.head(5:6) 0 5 5];
 	hand2.X = 350:650;
 	hand2.Y = 350:650;
	hand2.Z = Z(350:650, 350:650);
	hand2.head = double([hand2.X(1) hand2.X(end) hand2.Y(1) hand2.Y(end) min(hand2.Z(:)) max(hand2.Z(:)) 0 1 1]);

	Z = transplants([], 'grid', hand1, hand2);
% 	[Z,aa] = transplants([], 'grid', hand1, hand2);
% 	mirone(Z)
	hand1.X = X;	hand1.Y = Y;
	mirone(Z, hand1)

% ----------------------------
function fill_poly_hole
% ...

	[Z, hdrStruct] = gen_UMF2d(1.8, 0.05, 0.9, 1024);		% Pretend it's geogs (Remark, piking 1001 returns Z = NaNs everywhere)
	x = [-10 -7 -5 -10]';		y = [38 40 37.5 38]';
	x_lim = [hdrStruct.X(1) hdrStruct.X(end)];		y_lim = [hdrStruct.Y(1) hdrStruct.Y(end)];
	mask = img_fun('roipoly_j',x_lim,y_lim,Z,x,y);
	Z(mask) = NaN;
	hand1.X = hdrStruct.X;		hand1.Y = hdrStruct.Y;
	hand1.head = hdrStruct.head;
	hFig = mirone(Z, hand1);
	[Z2, hdrStruct2] = gen_UMF2d(1.8, 0.05, 0.95, 1024, [-10.5 -4.5 37 40.5]);
	hdrStruct2.Z = Z2;
	handles = guidata(hFig);
	transplants([], 'one_sharp', handles, hdrStruct2)
	
	
