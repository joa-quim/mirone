function make_mexs(opt)
%  make_gmt_mex -- Make the GMT-mex supplement directly in Matlab
%  This make works on Windows using the VC6 compiler. It should also work
%  for other configurations, but it was not tested.
%
%  Author:	Joaquim Luis (based on make_mexcdf53)
%  Date:	29-April-2005

if (nargin == 0)	opt = 'usage';	end

% Adjust for your own path
patoINC_GMT = 'd:\progs_cygw\GMTdev\GMT\';				% Include path for GMT
patoLIB_GMT = 'd:\progs_cygw\GMTdev\GMT_win\libMEX\';	% Lib path for GMT - Libs compiled with 'MEX condition'
%patoLIB_GMT = 'd:\progs_interix\GMTdev\GMT_win\lib\';	% Lib path for GMT
%patoLIB_GMT = 'c:\programs\gmt4\lib\';					%
pato_NETCDF = 'D:\progs_interix\netcdf-3.6.2b5_win\';	% path for NETCDF
%pato_GDAL = 'D:\programas\GDALB143\';					% path for GDAL
pato_GDAL = 'D:\programas\GDALtrunk\';					% path for GDAL
pato_OCV = 'C:\programas\OpenCV\';					% path for OpenCV
pato_SHAPELIB = 'D:\lixo\shapelib\';					% path for shapelib
pato_VC98LIB = 'C:\programas\VisualStudio\VC98\Lib\';	% path for MSVC library dir

if (ispc)
	COPT = '-DWIN32 -O';
else
	COPT = '-O';
end
% -------------------------- Stop editing (at least on Windows) ---------------------------

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
LIB_CXCORE = [pato_OCV 'lib\cxcore.lib'];
LIB_HG = [pato_OCV 'lib\highgui.lib'];

% GMT mexs
str_gmt = {'grdgradient_m' 'grdinfo_m' 'grdproject_m' 'grdread_m' 'grdsample_m' ...
        'grdtrack_m' 'grdtrend_m' 'grdwrite_m' 'mapproject_m' 'mapproject_m421' 'shoredump' 'surface_m' ...
        'nearneighbor_m' 'grdfilter_m' 'cpt2cmap' 'grdlandmask_m' 'grdppa_m'}';

% GMT MGG supplements mexs (currently only one)
str_gmt_mgg = {'gmtlist_m'};

% Gdal mexs
str_gdal = {'gdalread' 'gdalwrite'}';

% Gdal c++ mexs
str_gdal_cpp = {'gdalwarp_mex' 'ogrproj'}';

% Shape mexs (currently only one)
str_shape = {'mex_shape'}';

% OpenCV mexs (currently only one)
str_cv = {'cvcolor_mex' 'cvfill_mex' 'cvgetcorners_mex' 'cvresize_mex' 'cvlib_mex'}';

% netCDF mexes (other than GMT ones)
str_withCDF = {'swan'; 'swan_sem_wbar'};

% Non LIB dependent mexs (besides matlab libs, of course)
str_simple = {'test_gmt' 'igrf_m' 'scaleto8' 'tsun2' 'wave_travel_time' 'mansinha_m' ...
        'telha_m' 'range_change' 'country_select' 'mex_illuminate' 'grdutils' ...
        'read_isf' 'ind2rgb8' 'alloc_mex' 'susan' 'set_gmt' 'existmex' 'mxgridtrimesh' 'intlutc'}';

% Non LIB dependent c++ mexs
str_simple_cpp = {'houghmex' 'cimgmatlab_cannyderiche'}';
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
library_cv = [LIB_CV ' ' LIB_CXCORE ' ' LIB_HG];
library_vc6 = [LIB_USER32 ' ' LIB_GDI32];

opt_gmt = COPT;
opt_gmt_mgg = COPT;
if (ispc)
    opt_gmt = [COPT ' -DDLL_GMT -DDLL_NETCDF'];		% Are the -D... really needed?
    opt_gmt_mgg = [COPT ' -DDLL_GMT -DGMT_MGG'];
end

if (strcmp(opt,'all'))              % Compile the whole family
    for (i=1:length(str_gmt))       % Compile GMT mexs
        cmd = ['mex ' [str_gmt{i} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
        eval(cmd)
    end
    for (i=1:length(str_gmt_mgg))    % Compile GMT MGG mexs
        cmd = ['mex ' [str_gmt_mgg{i} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
        eval(cmd)
    end
    for (i=1:length(str_gdal))      % Compile Gdal mexs
        cmd = ['mex ' [str_gdal{i} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];
        eval(cmd)
    end
    for (i=1:length(str_cv))      % Compile OpenCV mexs
        cmd = ['mex ' [str_cv{i} '.c'] ' ' include_cv ' ' library_cv ' ' COPT];
        eval(cmd)
    end
    for (i=1:length(str_simple))    % Compile Other (simple) mexs
        cmd = ['mex ' [str_simple{i} '.c'] ' ' COPT];
        eval(cmd)
    end
    for (i=1:length(str_simple_cpp))    % Compile Other (simple) c++ mexs
        cmd = ['mex ' [str_simple_cpp{i} '.cpp'] ' ' COPT];
        eval(cmd)
    end
elseif (strcmp(lower(opt),'gmt'))   % Compile only the GMT mexs (and supplements)
    for (i=1:length(str_gmt))
        cmd = ['mex ' [str_gmt{i} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
        eval(cmd)
    end
    for (i=1:length(str_gmt_mgg))    % Compile GMT MGG mexs
        cmd = ['mex ' [str_gmt_mgg{i} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
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
    if (~isempty(idx1))         % Compile GMT mexs
        cmd = ['mex ' [str_gmt{idx1} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
    elseif (~isempty(idx4))     % Compile GMT MGG mexs
        cmd = ['mex ' [str_gmt_mgg{idx4} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
    elseif (~isempty(idx2))     % Compile Gdal mexs
        cmd = ['mex ' [str_gdal{idx2} '.c'] ' ' include_gdal ' ' library_gdal ' ' COPT];
    elseif (~isempty(idx2pp))     % Compile Gdal c++ mexs
        cmd = ['mex ' [str_gdal_cpp{idx2pp} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' COPT];
    elseif (~isempty(idx22))    % Compile Shape mexs
        cmd = ['mex ' [str_shape{idx22} '.c'] ' ' include_shape ' ' library_shape ' ' COPT];
    elseif (~isempty(idx5))     % Compile OpenCV mexs
        cmd = ['mex ' [str_cv{idx5} '.c'] ' ' include_cv ' ' library_cv ' ' COPT];
    elseif (~isempty(idx6))     % Compile Simple c++ mexs
        cmd = ['mex ' [str_simple_cpp{idx6} '.cpp'] ' ' library_vc6 ' ' COPT];
    elseif (~isempty(idx7))     % Compile netCDF dependent mexs
        cmd = ['mex ' str_withCDF{idx7} '.c' ' -I' INCLUDE_NETCDF ' ' LIB_NETCDF ' ' COPT];
    else                        % Compile Other (simple) mexs
        cmd = ['mex ' [str_simple{idx3} '.c'] ' ' COPT];
    end
    eval(cmd)
end
