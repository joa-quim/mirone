function make_mexs64(opt, varargin)
%  make_mexs -- Make the Mirone-mex supplement directly in Matlab
%  
%	OPT  -- one of 'all', 'gmt', 'gdal', 'mexnc', 'simple' or a specifif module name
%	varargin -- '64' or '32' to compile for 64 r 32 bits [default '32']
%				and/or compiler options like '-g', '-v', '-c'
%
%	A major drawback of this file is that the diverse libraries are very hardcoded
%	for the names on my computer, and to private names of some .lib files.
%
%  Author:	Joaquim Luis
%  Date:	02-April-2010

if (nargin == 0)
	opt = 'usage';	% Quite poor message though
end

mexe = 'mex ';
bit = '32';

if (nargin > 1)			% Used to provide compiler options like '-g' or '-v'
	for (i = 1:numel(varargin))
		if ( isnan(str2double(varargin{i})) )		% compiler options like '-g' or '-v'
			mexe = [mexe ' ' varargin{i} ' '];%#ok
		else
			bit = varargin{i};
		end
	end
end

% ------------- Adjust for your own path -----------------------------------------------
% Include path for GMT. Directory where the several *.h GMT files reside 
patoINC_GMT = 'c:\progs_cygw\GMTdev\GMT\';

% Lib path for GMT. Must contain the GMT *.lib library files (A bit of mess, since 64 libs are under ...\GMT_win64
patoLIB_GMT = 'c:\progs_cygw\GMTdev\GMT_win\lib\';
if (strcmp(bit,'64')),		patoLIB_GMT = 'c:\progs_cygw\GMTdev\GMT_win64\lib\';	end

% path for NETCDF base dir. Sub-dir 'lib' and 'include' must exist with, respectively, libnetcdf_w64.lib and header files
pato_NETCDF = ['C:\progs_cygw\netcdf-3.6.3\compileds\VC10_' bit '\'];

% path for the HDF4 base dir. Here I use my own compiled libs that have an extra '_mir' in their
% names to avoid name clashing with the often incompatible versions shiped by Matlab
pato_HDF = ['C:\programs\compa_libs\HDF-4.2.5\compileds\vc10_' bit '\'];

% path for GDAL. Sub-directories 'lib' and 'include' must exist with, respectively, the gdal_i.lib and header files
pato_GDAL = ['c:\programs\GDALtrunk\gdal\compileds\VC10_' bit '\'];

% path for OpenCV. Base directory where the OpenCV library was installed
% INCLUDE & LIB[64] MUST BE DIRECT DESCENDENTS OF "pato_OCV" DIR
pato_OCV = ['C:\programs\OpenCV_SVN\compileds\VC10_' bit '\'];

if (ispc),	COPT = '-DWIN32';	end

INCLUDE_NETCDF = [pato_NETCDF 'include\'];		LIB_NETCDF = [pato_NETCDF 'lib\libnetcdf.lib'];
INCLUDE_HDF = [pato_HDF 'include\'];			LIB_HDF = [pato_HDF 'lib\hm425m.lib ' pato_HDF 'lib\hd425m.lib'];

INCLUDE_GMT = [patoINC_GMT 'src\'];				LIB_GMT = [patoLIB_GMT 'gmt.lib'];
INCLUDE_GMT_MGG = [patoINC_GMT 'src\mgg'];		LIB_GMT_MGG = [patoLIB_GMT 'gmt_mgg.lib'];

INCLUDE_GDAL = [pato_GDAL 'include'];			LIB_GDAL = [pato_GDAL 'lib\gdal_i.lib'];

INCLUDE_CV = [pato_OCV 'include\opencv'];		LIB_CV = [pato_OCV 'lib\cv210.lib'];
INCLUDE_HG = INCLUDE_CV;						LIB_HG = [pato_OCV 'lib\highgui210.lib'];
INCLUDE_CXCORE = INCLUDE_CV;					LIB_CXCORE = [pato_OCV 'lib\cxcore210.lib'];
LIB_CV_HAAR = [pato_OCV 'lib\cvhaartraining.lib'];
% -------------------------- Stop editing -------------------------------------------------------

% GMT mexs
str_gmt = {'grdinfo_m' 'grdproject_m' 'grdread_m' 'grdsample_m' 'grdtrend_m' ...
        'grdwrite_m' 'mapproject_m' 'shoredump' 'surface_m' 'nearneighbor_m' ...
        'grdfilter_m' 'cpt2cmap' 'grdlandmask_m' 'grdppa_m' 'dimfilter_m' 'shake_mex' 'testgmt5'}';

% GMT MGG supplements mexs (currently only one)
str_gmt_mgg = {'gmtlist_m'};

% Gdal mexs
str_gdal = {'gdalread' 'gdalwrite' 'mex_shape'}';

% Gdal c++ mexs
str_gdal_cpp = {'gdalwarp_mex' 'gdaltransform_mex' 'ogrproj'}';

% OpenCV mexs (currently only one)
str_cv = {'cvlib_mex'}';

% netCDF mexes (other than GMT ones)
str_withCDF = {'swan'; 'swan_sem_wbar'};

% HDF mexes (L3 binned HDF SeaWiffs files)
str_withHDF = {'swreadl3b_m'};

% Non LIB dependent mexs (besides matlab libs, of course)
str_simple = {'test_gmt' 'igrf_m' 'scaleto8' 'tsun2' 'wave_travel_time' 'mansinha_m' ...
	'telha_m' 'range_change' 'country_select' 'mex_illuminate' 'grdutils' 'read_isf' ...
	'alloc_mex' 'susan' 'set_gmt' 'mxgridtrimesh' 'trend1d_m', 'gmtmbgrid_m' ...
	'grdgradient_m' 'grdtrack_m' 'spa_mex' 'ind2rgb8' 'mirblock' 'applylutc' 'cq' ...
	'bwlabel1' 'bwlabel2' 'imhistc' 'intlutc' 'inv__lwm' 'grayto8' 'grayto16' 'ordf' ...
	'parityscan' 'ditherc' 'PolygonClip' 'write_mex' 'xyzokb_m'}';

% Non LIB dependent c++ mexs
str_simple_cpp = {'houghmex' 'clipbd_mex' 'akimaspline' 'bwlabelnmex' 'bwboundariesmex' ...
	'imreconstructmex' 'morphmex' 'grayxform' 'resampsep' 'iptcheckinput'}';

% -----------------------------------------------------------------------------------------
include_gmt = ['-I' INCLUDE_GMT ' ' '-I' INCLUDE_NETCDF];
include_gmt_mgg = ['-I' INCLUDE_GMT_MGG];
library_gmt = [LIB_GMT ' ' LIB_NETCDF];
library_gmt_mgg = [LIB_GMT ' ' LIB_GMT_MGG];
library_gdal = LIB_GDAL;
include_gdal = ['-I' INCLUDE_GDAL];
include_cv = ['-I' INCLUDE_CV ' -I' INCLUDE_CXCORE ' -I' INCLUDE_HG];
library_cv = [LIB_CV ' ' LIB_CXCORE ' ' LIB_HG ' ' LIB_CV_HAAR];

opt_mexnc = [' -I' INCLUDE_NETCDF ' ' LIB_NETCDF ' -DDLL_NETCDF'];

opt_gmt = COPT;
opt_gmt_mgg = COPT;
if (ispc)
	opt_gmt = [COPT ' -DDLL_GMT -DDLL_NETCDF'];		% Are the -D... really needed?
	opt_gmt_mgg = [COPT ' -DDLL_GMT -DGMT_MGG'];
end

if (strcmp(opt,'all'))			% Compile the whole family
	for (i=1:numel(str_gmt))		% Compile GMT mexs
		cmd = [mexe [str_gmt{i} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
		eval(cmd)
	end
	for (i=1:numel(str_gmt_mgg))	% Compile GMT MGG mexs
		cmd = [mexe [str_gmt_mgg{i} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
		eval(cmd)
	end
	for (i=1:numel(str_gdal))		% Compile GDAL mexs
		cmd = [mexe [str_gdal{i} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];
		eval(cmd)
	end
	for (i=1:numel(str_gdal_cpp))	% Compile GDAL C++ mexs
		cmd = [mexe [str_gdal_cpp{i} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' COPT];
		eval(cmd)
	end
	for (i=1:numel(str_cv))			% Compile OpenCV mexs
		cmd = [mexe [str_cv{i} '.c']  ' sift\sift.c sift\imgfeatures.c sift\kdtree.c sift\minpq.c ' include_cv ' ' library_cv ' ' COPT];
		eval(cmd)
	end
	for (i=1:numel(str_withCDF))	% Compile netCDF (simple) mexs
		cmd = [mexe [str_withCDF{i} '.c']  ' -I' INCLUDE_NETCDF ' ' LIB_NETCDF ' ' COPT];
		eval(cmd)
	end
	for (i=1:numel(str_withHDF))	% Compile HDF (simple) mexs
		cmd = [mexe [str_withHDF{i} '.c']  ' -I' INCLUDE_HDF ' ' LIB_HDF ' ' COPT];
		eval(cmd)
	end

	make_simple(str_simple, mexe, str_simple_cpp, COPT)

	% Compile the MEXNC mexs
	cmd = [mexe 'mexnc\mexgateway.c mexnc\netcdf2.c mexnc\netcdf3.c mexnc\common.c -output mexnc ' opt_mexnc ' ' COPT];
	eval(cmd)

	% Compile the edison_wraper
	try		cd edison;		compile_edison_wrapper;		cd ..
	catch,	cd ..
	end

elseif (strcmpi(opt,'gmt'))			% Compile only the GMT mexs (and supplements)
	for (i=1:numel(str_gmt))
        cmd = [mexe [str_gmt{i} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
        eval(cmd)
	end
	for (i=1:numel(str_gmt_mgg))	% Compile GMT MGG mexs
        cmd = [mexe [str_gmt_mgg{i} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
        eval(cmd)
	end
elseif (strcmpi(opt,'gdal'))		% Compile only the GDAL mexs
	for (i=1:numel(str_gdal))		% Compile GDAL C mexs
        cmd = [mexe [str_gdal{i} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];
        eval(cmd)
	end
	for (i=1:numel(str_gdal_cpp))	% Compile GDAL C++ mexs
        cmd = [mexe [str_gdal_cpp{i} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' COPT];
        eval(cmd)
	end
elseif (strcmpi(opt,'mexnc'))	% Compile only the MEXNC mexs
	cmd = [mexe 'mexnc\mexgateway.c mexnc\netcdf2.c mexnc\netcdf3.c mexnc\common.c -output mexnc ' opt_mexnc ' ' COPT];
	eval(cmd)
elseif (strcmpi(opt,'simple'))		% Compile the 'simple' mexs
	make_simple(str_simple, mexe, str_simple_cpp, COPT)
elseif (strcmp(opt,'usage'))
	disp('Example usage: make_mexs(''mapproject_m'')')
	disp('	OR: make_mexs(''ALL'') -- Compile all familly')
	disp('	OR: make_mexs(''GMT'') -- Compile GMT mexs')
else							% Compile only one mex
    idx = strmatch(opt,[str_gmt; str_gmt_mgg; str_gdal; str_gdal_cpp; str_simple; str_cv; str_simple_cpp; ...
			str_withCDF; str_withHDF]);
    if (isempty(idx))
		disp('Example usage: make_mexs(''mapproject_m'')');
		error('Bad use, or my fault');
    end
	idx1   = strmatch(opt, str_gmt, 'exact');
	idx2   = strmatch(opt, str_gdal, 'exact');
	idx2pp = strmatch(opt, str_gdal_cpp, 'exact');
	idx3   = strmatch(opt, str_simple, 'exact');
	idx4   = strmatch(opt, str_gmt_mgg, 'exact');
	idx5   = strmatch(opt, str_cv, 'exact');
	idx6   = strmatch(opt, str_simple_cpp, 'exact');
	idx7   = strmatch(opt, str_withCDF, 'exact');
	idx8   = strmatch(opt, str_withHDF, 'exact');
	idx9   = strcmpi(opt, 'polygonclip');
	idx10  = strcmpi(opt, 'ditherc');
	if (idx9 && idx3),		idx3 = [];		end		% Need  a extra arg , so treat sep
	if (idx10 && idx3),		idx3 = [];		end		% 			"
	if (~isempty(idx1))         % Compile GMT mexs
		cmd = [mexe [str_gmt{idx1} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];

    elseif (~isempty(idx4))     % Compile GMT MGG mexs
		cmd = [mexe [str_gmt_mgg{idx4} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];

	elseif (~isempty(idx2))     % Compile GDAL mexs
		cmd = [mexe [str_gdal{idx2} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];

	elseif (~isempty(idx2pp))     % Compile GDAL c++ mexs
		cmd = [mexe [str_gdal_cpp{idx2pp} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' COPT];

    elseif (~isempty(idx3))     % Compile Simple c mexs
        cmd = [mexe [str_simple{idx3} '.c'] ' ' COPT];

	elseif (~isempty(idx5))     % Compile OpenCV mexs
		cmd = [mexe [str_cv{idx5} '.c'] ' sift\sift.c sift\imgfeatures.c sift\kdtree.c sift\minpq.c ' include_cv ' ' library_cv ' ' COPT];

	elseif (~isempty(idx6))     % Compile Simple c++ mexs
		cmd = [mexe [str_simple_cpp{idx6} '.cpp'] ' ' COPT];

	elseif (~isempty(idx7))     % Compile netCDF dependent mexs
		cmd = [mexe str_withCDF{idx7} '.c' ' -I' INCLUDE_NETCDF ' ' LIB_NETCDF ' ' COPT];

	elseif (~isempty(idx8))     % Compile netHDF dependent mexs
		cmd = [mexe str_withHDF{idx8} '.c' ' -I' INCLUDE_HDF ' ' LIB_HDF ' ' COPT];

	elseif (~isempty(idx9) && idx9)
		cmd = [mexe 'PolygonClip.c gpc.c ' COPT];

	elseif (~isempty(idx10))
		cmd = [mexe 'ditherc.c invcmap.c ' COPT];

	else						% Compile Other (simple) mexs
		cmd = [mexe [str_simple{idx3} '.c'] ' ' COPT];

	end
    eval(cmd)
end

% ----------------------------------------------------------------------------
function make_simple(str_simple, mexe, str_simple_cpp, COPT)
% Compile the so called 'simple' mexs. That is, mexs that don't depend on external libs
	for (i=1:numel(str_simple))
		if (strcmp(str_simple{i}, 'ditherc'))
			cmd = [mexe [str_simple{i} '.c '] 'invcmap.c ' COPT];
		elseif (strcmp(str_simple{i}, 'PolygonClip'))
			cmd = [mexe [str_simple{i} '.c  gpc.c '] COPT];
		else
			cmd = [mexe [str_simple{i} '.c '] COPT];
		end
        eval(cmd)
	end
	for (i=1:numel(str_simple_cpp))
		if (strcmp(str_simple_cpp{i}, 'bwlabelnmex'))
	        cmd = [mexe [str_simple_cpp{i} '.cpp '] 'neighborhood.cpp unionfind.c ' COPT];
		elseif (strcmp(str_simple_cpp{i}, 'bwboundariesmex'))
	        cmd = [mexe [str_simple_cpp{i} '.cpp '] 'boundaries.cpp ' COPT];
		elseif (strcmp(str_simple_cpp{i}, 'morphmex'))
	        cmd = [mexe [str_simple_cpp{i} '.cpp '] 'dilate_erode_gray_nonflat.cpp dilate_erode_packed.cpp dilate_erode_binary.cpp neighborhood.cpp vectors.cpp ' COPT];
		elseif (strcmp(str_simple_cpp{i}, 'imreconstructmex'))
	        cmd = [mexe [str_simple_cpp{i} '.cpp '] 'neighborhood.cpp ' COPT];
		else
	        cmd = [mexe [str_simple_cpp{i} '.cpp '] COPT];
		end
        eval(cmd)
	end

