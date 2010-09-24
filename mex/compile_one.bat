@echo off
REM --------------------------------------------------------------------------------------
REM
REM	$Id:
REM
REM This is a compile batch that builds one of the MEXs whose source code is distributed in Mirone
REM A major difficulty in using this comes from the fact that several external libraries are needed.
REM They are:
REM		GMT, netCDF, GDAL and OpenCV
REM If a WIN64 version is targeted than the above Libs must have been build in 64-bits as well.
REM
REM Notes: To compile gmtlist_m the file x_system.h must be copyed from 
REM		GMTROOT\src\x_system to %GMT_INC% (see below)
REM        To compile mex_shape the file shapefil.h must be copyed from 
REM		GDALROOT\ogr\ogrsf_frmts\shape to %GDAL_INC%
REM
REM Usage: open the command window set up by the compiler of interest (were all vars are already set)
REM        and run this from there. Example:
REM	   compile_one surface_m
REM
REM
REM Author: Joaquim Luis, 09-MAY-2010
REM --------------------------------------------------------------------------------------

REM ------------- Set the compiler (set to 'icl' to use the Intel compiler) --------------
SET CC=icl
REM --------------------------------------------------------------------------------------

REM If set to "yes", linkage is done againsts ML6.5 Libs (needed in compiled version)
SET R13="no"

REM Set it to "yes" or "no" to build under 64-bits or 32-bits respectively.
SET WIN64="no"

REM The MSVC version. I use this var to select libraries also compiled with this compiler
SET MSVC_VER="1600"

REM Options are "dll", "mexw32" (recent ML version scream when they find .dll) or "mexw64" (when WIN64="yes")
SET MEX_EXT="mexw32"

REM If set some MEXs will print the execution time (in CPU ticks)
REM SET TIMEIT=-DMIR_TIMEIT
SET TIMEIT=

REM
REM Set to "yes" if you want to build a debug version
SET DEBUG="no"
REM
SET LDEBUG=
IF %DEBUG%=="yes" SET LDEBUG=/debug

REM --- Next allows compiling with the compiler you want against the ML6.5 libs (needed in stand-alone version)
IF %R13%=="yes" (
SET MATLIB=C:\SVN\MAT65_pracompa\extern\lib\win32\microsoft
SET MATINC=C:\SVN\MAT65_pracompa\extern\include
SET _MX_COMPAT=
SET MEX_EXT="dll"

) ELSE (

IF %WIN64%=="yes" (
SET MATLIB=C:\PROGRAMS\MATLAB\R2010A\extern\lib\win64\microsoft 
SET MATINC=C:\PROGRAMS\MATLAB\R2010A\extern\include
SET _MX_COMPAT=-DMX_COMPAT_32
SET MEX_EXT="mexw64"

) ELSE (

SET MATLIB=C:\PROGRAMS\MATLAB\R2009B\extern\lib\win32\microsoft 
SET MATINC=C:\PROGRAMS\MATLAB\R2009B\extern\include
SET _MX_COMPAT=-DMX_COMPAT_32
) )

 
REM -------------- Set up libraries here -------------------------------------------------
IF %WIN64%=="yes" (

SET NETCDF_LIB=C:\progs_cygw\netcdf-3.6.3\compileds\VC10_64\lib\libnetcdf.lib
rem SET GMT_LIB=c:\progs_cygw\GMTdev\GMT_win64\lib\gmt.lib
rem SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\GMT_win64\lib\gmt_mgg.lib
SET GMT_LIB=c:\progs_cygw\GMTdev\GMT5\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\GMT5\lib\gmt_mgg.lib
SET GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_64\lib\gdal_i.lib
SET CV_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\cv210.lib
SET CXCORE_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\cxcore210.lib

) ELSE (

IF %MSVC_VER%=="1600" (
SET NETCDF_LIB=C:\progs_cygw\netcdf-3.6.3\compileds\VC10_32\lib\libnetcdf.lib
SET GMT_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt_mgg.lib
SET GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_32\lib\gdal_i.lib
SET CV_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\cv210.lib
SET CXCORE_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\cxcore210.lib
SET LAS_LIB=C:\programs\compa_libs\liblas-src-1.2.1\bin\lib\liblas_i.lib

) ELSE (

SET NETCDF_LIB=C:\progs_cygw\netcdf-3.6.3\lib\libnetcdf_w32.lib
SET GMT_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt_mgg.lib
SET GDAL_LIB=c:\programs\GDALtrunk\gdal\lib\gdal_i.lib
REM I haven't build yet (and maybe I won't) 2.1 libs with VC7.1
SET CV_LIB=C:\programs\OpenCV_SVN\lib\cv200.lib
SET CXCORE_LIB=C:\programs\OpenCV_SVN\lib\cxcore200.lib
SET LAS_LIB=C:\programs\compa_libs\liblas-src-1.2.1\bin\lib\liblas_i.lib
) )

SET NETCDF_INC=C:\progs_cygw\netcdf-3.6.3\include
rem SET GMT_INC=c:\progs_cygw\GMTdev\GMT\include
SET GMT_INC=c:\progs_cygw\GMTdev\GMT5\include
SET GMT_INC2=C:\progs_cygw\GMTdev\GMT5\src\mex
SET GDAL_INC=c:\programs\GDALtrunk\gdal\compileds\VC10_32\include
SET OCV_INC=C:\programs\OpenCV_SVN\compileds\VC10_32\include\opencv
SET LAS_INC=-IC:\programs\compa_libs\liblas-src-1.2.1\bin\include\liblas\capi -IC:\programs\compa_libs\liblas-src-1.2.1\bin\include\liblas
REM ----------------------------------------------------------------------------

REM ____________________________________________________________________________
REM ___________________ STOP EDITING HERE ______________________________________


SET COMPFLAGS=/c /Zp8 /GR /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD 
IF %DEBUG%=="no" SET OPTIMFLAGS=/Ox /Oy- /DNDEBUG
IF %DEBUG%=="yes" SET OPTIMFLAGS=/Z7

IF %WIN64%=="yes" SET arc=X64
IF %WIN64%=="no" SET arc=X86
SET LINKFLAGS=/dll /export:mexFunction /LIBPATH:%MATLIB% libmx.lib libmex.lib libmat.lib /MACHINE:%arc% kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /incremental:NO %LDEBUG% 

%CC% -DWIN32 %COMPFLAGS% -I%GMT_INC% -I%GDAL_INC% -I%NETCDF_INC% -I%OCV_INC% -I%MATINC% %LAS_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %1.c
link  /out:"%1.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% %GDAL_LIB% %LAS_LIB% /implib:templib.x %1.obj

del *.obj *.exp templib.x
