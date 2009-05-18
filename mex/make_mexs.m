function make_mexs(opt)
%  make_gmt_mex -- Make the GMT-mex supplement directly in Matlab
%  This make works on Windows using the VC6 compiler. It should also work
%  for other configurations, but it was not tested.
%
%  Author:	Joaquim Luis (based on make_mexcdf53)
%  Date:	29-April-2005

if (nargin == 0),	opt = 'usage';	end			% Quite poor message though

% ------------- Adjust for your own path -----------------------------------------------

% Include path for GMT. Directory where the several *.h GMT files reside 
patoINC_GMT = 'c:\progs_cygw\GMTdev\GMT\';

% Lib path for GMT - Libs compiled with 'MEX condition'. Must contain the GMT *.lib library files
patoLIB_GMT = 'c:\progs_cygw\GMTdev\GMT_win\libMEX\';
%patoLIB_GMT = 'c:\progs_cygw\GMTdev\GMT_win\lib\';	% Lib path for GMT

% path for NETCDF bae dir. Sub-directories 'lib' and 'include' must exist with, respectively, libnetcdf.lib and header files
% I use ftp://ftp.unidata.ucar.edu/pub/netcdf/contrib/win32/netcdf-3.6.2-beta5_pgi_w32bin.zip
pato_NETCDF = 'c:\progs_interix\netcdf-3.6.2b5_win\';

% path for GDAL. Sub-directories 'lib' and 'include' must exist with, respectively, the gdal_i.lib and header files
pato_GDAL = 'c:\programs\GDALtrunk\gdal\';

% path for OpenCV. Base directory where the OpenCV library was installed (still using v1.0)
pato_OCV = 'C:\programs\OpenCV\';

% path for shapelib. Directory where shapelib.lib shapefil.h files reside
pato_SHAPELIB = 'c:\lixo\shapelib\';

% path for MSVC library dir
pato_VC98LIB = 'C:\programs\VisualStudio\VC98\Lib\';
% -------------------------- Stop editing (at least on Windows) ---------------------------

if (ispc),	COPT = '-DWIN32 -O';
else		COPT = '-O';
end

INCLUDE_GMT = [patoINC_GMT 'src\'];               % Core gmt programs
INCLUDE_GMT_MGG = [patoINC_GMT 'src\mgg'];        % MGG supplements

LIB_GMT = [patoLIB_GMT 'gmt.lib'];
LIB_GMT_MGG = [patoLIB_GMT 'gmt_mgg.lib'];

LIB_NETCDF = [pato_NETCDF 'lib\libnetcdf.lib'];
INCLUDE_NETCDF = [pato_NETCDF 'include\'];

LIB_GDAL = [pato_GDAL 'lib\gdal_i.lib'];
INCLUDE_GDAL = [pato_GDAL 'include'];

INCLUDE_SHAPE = pato_SHAPELIB;
LIB_SHAPE = [pato_SHAPELIB 'shapelib.lib'];

INCLUDE_CV = [pato_OCV 'cv\include'];
INCLUDE_HG = [pato_OCV 'otherlibs\highgui'];
INCLUDE_CXCORE = [pato_OCV 'cxcore\include'];
LIB_CV = [pato_OCV 'lib\cv.lib'];
LIB_CV_HAAR = [pato_OCV 'lib\cvhaartraining.lib'];
LIB_CXCORE = [pato_OCV 'lib\cxcore.lib'];
LIB_HG = [pato_OCV 'lib\highgui.lib'];

% GMT mexs
str_gmt = {'grdinfo_m' 'grdproject_m' 'grdread_m' 'grdsample_m' ...
        'grdtrend_m' 'grdwrite_m' 'mapproject_m' 'shoredump' 'surface_m' ...
        'nearneighbor_m' 'grdfilter_m' 'cpt2cmap' 'grdlandmask_m' 'grdppa_m' 'dimfilter_m'}';

% GMT MGG supplements mexs (currently only one)
str_gmt_mgg = {'gmtlist_m'};

% Gdal mexs
str_gdal = {'gdalread' 'gdalwrite'}';

% Gdal c++ mexs
str_gdal_cpp = {'gdalwarp_mex' 'gdaltransform_mex' 'ogrproj'}';

% Shape mexs (currently only one)
str_shape = {'mex_shape'}';

% OpenCV mexs (currently only one)
str_cv = {'cvlib_mex'}';

% netCDF mexes (other than GMT ones)
str_withCDF = {'swan'; 'swan_sem_wbar'};

% Non LIB dependent mexs (besides matlab libs, of course)
str_simple = {'test_gmt' 'igrf_m' 'scaleto8' 'tsun2' 'wave_travel_time' 'mansinha_m' ...
	'telha_m' 'range_change' 'country_select' 'mex_illuminate' 'grdutils' ...
	'read_isf' 'ind2rgb8' 'alloc_mex' 'susan' 'set_gmt' 'mxgridtrimesh' ...
	'intlutc' 'trend1d_m', 'gmtmbgrid_m' 'grdgradient_m' 'grdtrack_m' 'spa_mex' 'cm4field_m' ...
	'PolygonClip' }';

% Non LIB dependent c++ mexs
str_simple_cpp = {'houghmex' 'clipbd_mex' 'akimaspline'}';
LIB_USER32 = [pato_VC98LIB 'USER32.LIB'];
LIB_GDI32 = [pato_VC98LIB 'GDI32.LIB'];

% -----------------------------------------------------------------------------------------
include_gmt = ['-I' INCLUDE_GMT ' ' '-I' INCLUDE_NETCDF];
include_gmt_mgg = ['-I' INCLUDE_GMT_MGG];
library_gmt = [LIB_GMT ' ' LIB_NETCDF];
library_gmt_mgg = [LIB_GMT ' ' LIB_GMT_MGG];
library_gdal = LIB_GDAL;
include_gdal = ['-I' INCLUDE_GDAL];
include_shape = ['-I' INCLUDE_SHAPE];
library_shape = LIB_SHAPE;
include_cv = ['-I' INCLUDE_CV ' -I' INCLUDE_CXCORE ' -I' INCLUDE_HG];
library_cv = [LIB_CV ' ' LIB_CXCORE ' ' LIB_HG ' ' LIB_CV_HAAR];
library_vc6 = [LIB_USER32 ' ' LIB_GDI32];

opt_gmt = COPT;
opt_gmt_mgg = COPT;
if (ispc)
    opt_gmt = [COPT ' -DDLL_GMT -DDLL_NETCDF'];		% Are the -D... really needed?
    opt_gmt_mgg = [COPT ' -DDLL_GMT -DGMT_MGG'];
end

if (strcmp(opt,'all'))			% Compile the whole family
	for (i=1:numel(str_gmt))		% Compile GMT mexs
        cmd = ['mex ' [str_gmt{i} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
        eval(cmd)
	end
	for (i=1:numel(str_gmt_mgg))	% Compile GMT MGG mexs
        cmd = ['mex ' [str_gmt_mgg{i} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
        eval(cmd)
	end
	for (i=1:numel(str_gdal))		% Compile GDAL mexs
        cmd = ['mex ' [str_gdal{i} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];
        eval(cmd)
	end
	for (i=1:numel(str_gdal_cpp))	% Compile GDAL C++ mexs
        cmd = ['mex ' [str_gdal_cpp{i} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' COPT];
        eval(cmd)
	end
	for (i=1:numel(str_cv))		% Compile OpenCV mexs
        cmd = ['mex ' [str_cv{i} '.c']  ' sift\sift.c sift\imgfeatures.c sift\kdtree.c sift\minpq.c ' include_cv ' ' library_cv ' ' COPT];
        eval(cmd)
	end
	for (i=1:numel(str_simple))	% Compile Other (simple) mexs
        cmd = ['mex ' [str_simple{i} '.c'] ' ' COPT];
        eval(cmd)
	end
	for (i=1:numel(str_simple_cpp))	% Compile Other (simple) c++ mexs
		cmd = ['mex ' [str_simple_cpp{i} '.cpp'] ' ' COPT];
		eval(cmd)
	end
	cmd = ['mex PolygonClip.c gpc.c ' COPT];
	eval(cmd)
elseif (strcmpi(opt,'gmt'))	% Compile only the GMT mexs (and supplements)
    for (i=1:numel(str_gmt))
        cmd = ['mex ' [str_gmt{i} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
        eval(cmd)
    end
    for (i=1:numel(str_gmt_mgg))	% Compile GMT MGG mexs
        cmd = ['mex ' [str_gmt_mgg{i} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
        eval(cmd)
    end
elseif (strcmpi(opt,'gdal'))	% Compile only the GDAL mexs
    for (i=1:numel(str_gdal))		% Compile GDAL C mexs
        cmd = ['mex ' [str_gdal{i} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];
        eval(cmd)
    end
    for (i=1:numel(str_gdal_cpp))	% Compile GDAL C++ mexs
        cmd = ['mex ' [str_gdal_cpp{i} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' COPT];
        eval(cmd)
    end
else                                % Compile only one mex
    idx = strmatch(opt,[str_gmt; str_gmt_mgg; str_gdal; str_gdal_cpp; str_simple; str_shape; str_cv; str_simple_cpp; str_withCDF]);
    if (isempty(idx))
        disp('Example usage: make_mexs(''mapproject_m'')');
        error('Bad use, or my fault');
    end
    idx1   = strmatch(opt, str_gmt, 'exact');
    idx2   = strmatch(opt, str_gdal, 'exact');
    idx2pp = strmatch(opt, str_gdal_cpp, 'exact');
    idx22  = strmatch(opt, str_shape, 'exact');
    idx3   = strmatch(opt, str_simple, 'exact');
    idx4   = strmatch(opt, str_gmt_mgg, 'exact');
    idx5   = strmatch(opt, str_cv, 'exact');
    idx6   = strmatch(opt, str_simple_cpp, 'exact');
    idx7   = strmatch(opt, str_withCDF, 'exact');
    idx8   = strcmpi(opt, 'polygonclip');
    if (~isempty(idx1))         % Compile GMT mexs
        cmd = ['mex ' [str_gmt{idx1} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
    elseif (~isempty(idx4))     % Compile GMT MGG mexs
        cmd = ['mex ' [str_gmt_mgg{idx4} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
    elseif (~isempty(idx2))     % Compile GDAL mexs
        cmd = ['mex ' [str_gdal{idx2} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];
    elseif (~isempty(idx2pp))     % Compile GDAL c++ mexs
        cmd = ['mex ' [str_gdal_cpp{idx2pp} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' COPT];
    elseif (~isempty(idx22))    % Compile Shape mexs
        cmd = ['mex ' [str_shape{idx22} '.c'] ' ' include_shape ' ' library_shape ' ' COPT];
    elseif (~isempty(idx5))     % Compile OpenCV mexs
        cmd = ['mex ' [str_cv{idx5} '.c'] ' sift\sift.c sift\imgfeatures.c sift\kdtree.c sift\minpq.c ' include_cv ' ' library_cv ' ' COPT];
    elseif (~isempty(idx6))     % Compile Simple c++ mexs
        cmd = ['mex ' [str_simple_cpp{idx6} '.cpp'] ' ' library_vc6 ' ' COPT];
    elseif (~isempty(idx7))     % Compile netCDF dependent mexs
        cmd = ['mex ' str_withCDF{idx7} '.c' ' -I' INCLUDE_NETCDF ' ' LIB_NETCDF ' ' COPT];
	elseif (~isempty(idx8))
		cmd = ['mex PolygonClip.c gpc.c ' COPT];
    else                        % Compile Other (simple) mexs
        cmd = ['mex ' [str_simple{idx3} '.c'] ' ' COPT];
    end
    eval(cmd)
end
