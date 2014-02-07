function  test_bots(opt,varargin)

if (nargin)
	switch opt
		case 'Illum'
			test_Illum
		case 'xyz'
			test_xyz
		case 'writeascii'
			writeascii
	end
else
	test_Illum
	test_xyz
	writeascii
end

% -----------------
function test_Illum
	h = mirone('swath_grid.grd');
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
	try		mirone('ImageIllum',luz, handles, 'lambertian'),			pause(1)
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

