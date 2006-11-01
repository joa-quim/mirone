function make_mexs(opt)
%  make_gmt_mex -- Make the GMT-mex supplement directly in Matlab
%  This make works on Windows using the VC6 compiler. It should also work
%  for other configurations, but it was not tested.
%
%  Author:	Joaquim Luis (based on make_mexcdf53)
%  Date:	29-April-2005

if (nargin == 0)	opt = 'usage';	end

% Adjust for your own path
INCLUDE_GMT = 'd:\progs_interix\GMTdev\GMT\src\';               % Core gmt programs
INCLUDE_GMT_MGG = 'd:\progs_interix\GMTdev\GMT\src\mgg';        % MGG supplements
INCLUDE_NETCDF = 'd:\progs_interix\GMTdev\netcdf_win\include\';
%INCLUDE_GDAL = 'D:\programas\GDALBuild\include';
INCLUDE_GDAL = 'D:\programas\GDALB132\include';
%LIB_GMT = 'c:\programs\gmt4\lib\gmt.lib';

LIB_GMT = 'd:\progs_interix\GMTdev\GMT_win\lib\gmt.lib';
%LIB_GMT_MGG = 'c:\programs\gmt4\lib\gmt_mgg.lib';
LIB_GMT_MGG = 'd:\progs_interix\GMTdev\GMT_win\lib\gmt_mgg.lib';
LIB_NETCDF = 'd:\progs_interix\GMTdev\netcdf_win\lib\netcdf.lib';
%LIB_GDAL = 'D:\programas\GDALBuild\lib\gdal_i.lib';
LIB_GDAL = 'D:\programas\GDALB132\lib\gdal_i.lib';

INCLUDE_SHAPE = 'D:\lixo\shapelib';
LIB_SHAPE = 'D:\lixo\shapelib\shapelib.lib';

% OpenCV
INCLUDE_CV = 'C:\programas\OpenCV\cv\include';
INCLUDE_HG = 'C:\programas\OpenCV\otherlibs\highgui';
INCLUDE_CXCORE = 'C:\programas\OpenCV\cxcore\include';
LIB_CV = 'C:\programas\OpenCV\lib\cv.lib';
LIB_CXCORE = 'C:\programas\OpenCV\lib\cxcore.lib';
LIB_HG = 'C:\programas\OpenCV\lib\highgui.lib';

% GMT mexs
str_gmt = {'grdgradient_m' 'grdinfo_m' 'grdproject_m' 'grdread_m' 'grdsample_m' ...
        'grdtrack_m' 'grdtrend_m' 'grdwrite_m' 'mapproject_m' 'shoredump' 'surface_m' ...
        'nearneighbor_m' 'grdfilter_m' 'cpt2cmap' 'grdppa_m'}';

% GMT MGG supplements mexs (currently only one)
str_gmt_mgg = {'gmtlist_m'};

% Gdal mexs
str_gdal = {'gdalread' 'gdalwrite'}';

% Gdal c++ mexs
str_gdal_cpp = {'gdalvirtual' 'ogrproj' 'leca' 'importwkt'}';

% Shape mexs (currently only one)
str_shape = {'mex_shape'}';

% OpenCV mexs (currently only one)
str_cv = {'cvcolor_mex' 'cvfill_mex' 'cvgetcorners_mex' 'cvresize_mex' 'cvlib_mex'}';

% Non LIB dependent mexs (besides matlab libs, of course)
str_simple = {'test_gmt' 'igrf_m' 'scaleto8' 'swan' 'tsun2' 'wave_travel_time' 'mansinha_m' ...
        'telha_m' 'range_change' 'country_select' 'mex_illuminate' 'grdutils' ...
        'read_isf' 'ind2rgb8' 'alloc_mex' 'susan'}';

% Non LIB dependent c++ mexs
str_simple_cpp = {'houghmex' 'cimgmatlab_cannyderiche'}';
LIB_USER32 = 'C:\programas\VisualStudio\VC98\Lib\USER32.LIB';
LIB_GDI32 = 'C:\programas\VisualStudio\VC98\Lib\GDI32.LIB';

% -------------------------- Stop editing ---------------------------
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

if (ispc)
    %opt_gmt = '-O -DWIN32 -DDLL_GMT -DLL_NETCDF -DDLL_PSL';
    opt_gmt = '-O -DWIN32 -DDLL_GMT -DDLL_NETCDF';
    opt_gmt_mgg = '-O -DWIN32 -DDLL_GMT -DGMT_MGG';
    opt_cv = '-O -DWIN32 -DDLL_CV097 -DDLL_CXCORE097';
    opt_simple = '-O -DWIN32';
else
    opt_gmt = '-O';
    opt_simple = '-O';
    opt_cv = '-O';
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
        cmd = ['mex ' [str_gdal{i} '.c'] ' ' include_gdal ' ' library_gdal ' ' opt_simple];
        eval(cmd)
    end
    for (i=1:length(str_cv))      % Compile OpenCV mexs
        cmd = ['mex ' [str_cv{i} '.c'] ' ' include_cv ' ' library_cv ' ' opt_cv];
        eval(cmd)
    end
    for (i=1:length(str_simple))    % Compile Other (simple) mexs
        cmd = ['mex ' [str_simple{i} '.c'] ' ' opt_simple];
        eval(cmd)
    end
    for (i=1:length(str_simple_cpp))    % Compile Other (simple) c++ mexs
        cmd = ['mex ' [str_simple_cpp{i} '.cpp'] ' ' opt_simple];
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
    idx = strmatch(opt,[str_gmt; str_gmt_mgg; str_gdal; str_gdal_cpp; str_simple; str_shape; str_cv; str_simple_cpp]);
    if (isempty(idx))
        disp('Example usage: make_mexs(''mapproject_m'')');
        error('Bad use, or my fault');
    end;
    idx1 = strmatch(opt,str_gmt,'exact');
    idx2 = strmatch(opt,str_gdal,'exact');
    idx2pp = strmatch(opt,str_gdal_cpp,'exact');
    idx22 = strmatch(opt,str_shape,'exact');
    idx3 = strmatch(opt,str_simple,'exact');
    idx4 = strmatch(opt,str_gmt_mgg,'exact');
    idx5 = strmatch(opt,str_cv,'exact');
    idx6 = strmatch(opt,str_simple_cpp,'exact');
    if (~isempty(idx1))         % Compile GMT mexs
        cmd = ['mex ' [str_gmt{idx1} '.c'] ' ' include_gmt ' ' library_gmt ' ' opt_gmt];
    elseif (~isempty(idx4))     % Compile GMT MGG mexs
        cmd = ['mex ' [str_gmt_mgg{idx4} '.c'] ' ' include_gmt ' ' include_gmt_mgg ' ' library_gmt_mgg ' ' opt_gmt_mgg];
    elseif (~isempty(idx2))     % Compile Gdal mexs
        cmd = ['mex ' [str_gdal{idx2} '.c'] ' ' include_gdal ' ' library_gdal ' ' opt_simple];
    elseif (~isempty(idx2pp))     % Compile Gdal c++ mexs
        cmd = ['mex ' [str_gdal_cpp{idx2pp} '.cpp'] ' ' include_gdal ' ' library_gdal ' ' opt_simple];
    elseif (~isempty(idx22))    % Compile Shape mexs
        cmd = ['mex ' [str_shape{idx22} '.c'] ' ' include_shape ' ' library_shape ' ' opt_simple];
    elseif (~isempty(idx5))     % Compile OpenCV mexs
        cmd = ['mex ' [str_cv{idx5} '.c'] ' ' include_cv ' ' library_cv ' ' opt_cv];
    elseif (~isempty(idx6))     % Compile Simple c++ mexs
        cmd = ['mex ' [str_simple_cpp{idx6} '.cpp'] ' ' library_vc6 ' ' opt_simple];
    else                        % Compile Other (simple) mexs
        cmd = ['mex ' [str_simple{idx3} '.c'] ' ' opt_simple];
    end
    eval(cmd)
end
