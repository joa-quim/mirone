@echo off
REM --------------------------------------------------------------------------------------
REM
REM	$Id :
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

REM --- List of those who run faster when compiled with VS2010 than Intel ----------------
REM grdutils grdgradient_m
REM --------------------------------------------------------------------------------------


REM ------------- Set the compiler (set to 'icl' to use the Intel compiler) --------------
SET CC=cl
REM --------------------------------------------------------------------------------------

REM If set to "yes", linkage is done againsts ML6.5 Libs (needed in compiled version)
SET R13="no"

REM Set it to 32 or 64 to build under 64-bits or 32-bits respectively.
SET BITS=32

IF %R13%=="yes" SET BITS=32

REM The MSVC version. I use this var to select libraries also compiled with this compiler
SET MSVC_VER="1600"

REM Options are "dll", "mexw32" (recent ML version scream when they find .dll) or "mexw64" (when BITS=64)
SET MEX_EXT="mexw32"

REM If set some MEXs will print the execution time (in CPU ticks)
SET TIMEIT=-DMIR_TIMEIT
SET TIMEIT=

REM To buils with OpenMP support (very few)
SET OMP=
SET OMP=-DHAVE_OPENMP 

REM
REM Set to "yes" if you want to build a debug version
SET DEBUG="no"
REM
SET LDEBUG=
IF %DEBUG%=="yes" SET LDEBUG=/debug

REM --- Next allows compiling with the compiler you want against the ML6.5 libs (needed in stand-alone version)
IF %R13%=="yes" (
SET MATLIB=C:\SVN\pracompila\MAT65\lib\win32\microsoft
SET MATINC=C:\SVN\pracompila\MAT65\include
SET _MX_COMPAT=
SET MEX_EXT="dll"

) ELSE (

IF %BITS%==64 (
SET MATLIB=C:\SVN\pracompila\ML2010a_w64\lib\win64\microsoft 
SET MATINC=C:\SVN\pracompila\ML2010a_w64\include
SET _MX_COMPAT=-DMX_COMPAT_32
SET MEX_EXT="mexw64"

) ELSE (

SET MATLIB=C:\SVN\pracompila\ML2009b_w32\lib\win32\microsoft 
SET MATINC=C:\SVN\pracompila\ML2009b_w32\include
SET _MX_COMPAT=-DMX_COMPAT_32
) )

 
REM -------------- Set up libraries here -------------------------------------------------
IF %BITS%==64 (

SET  NETCDF_LIB=C:\programs\compa_libs\netcdf\compileds\VC10_64\lib\netcdf.lib
SET     GMT_LIB=c:\progs_cygw\GMTdev\gmt4\WIN64\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\gmt4\WIN64\lib\gmt_mgg.lib
SET    GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_64\lib\gdal_i.lib
SET  CXCORE_LIB=C:\programs\compa_libs\opencv\compileds\VC10_64\lib\opencv_core.lib
SET   CVIMG_LIB=C:\programs\compa_libs\opencv\compileds\VC10_64\lib\opencv_imgproc.lib
SET CVCALIB_LIB=C:\programs\compa_libs\opencv\compileds\VC10_64\lib\opencv_calib3d.lib
SET   CVOBJ_LIB=C:\programs\compa_libs\opencv\compileds\VC10_64\lib\opencv_objdetect.lib
SET CVVIDEO_LIB=C:\programs\compa_libs\opencv\compileds\VC10_64\lib\opencv_video.lib
SET CVPHOTO_LIB=C:\programs\compa_libs\opencv\compileds\VC10_64\lib\opencv_photo.lib
SET     LAS_LIB=C:\programs\compa_libs\liblas-src-1.2.1\lib\VC10_64\liblas_i.lib
SET  GEOLIB_LIB=C:\programs\compa_libs\GeographicLib-1.16\compileds\VC10_64\lib\Geographic.lib
SET LASZLIB_LIB=C:\programs\compa_libs\lastools\compileds\VC10_64\lib\laslib_i.lib 

) ELSE (

SET  NETCDF_LIB=C:\programs\compa_libs\netcdf\compileds\VC10_32\lib\netcdf.lib
SET     GMT_LIB=c:\progs_cygw\GMTdev\gmt4\WIN32\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\gmt4\WIN32\lib\gmt_mgg.lib
SET    GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_32\lib\gdal_i.lib
SET  CXCORE_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_core.lib
SET   CVIMG_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_imgproc.lib
SET CVCALIB_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_calib3d.lib
SET   CVOBJ_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_objdetect.lib
SET CVVIDEO_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_video.lib
SET CVPHOTO_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_photo.lib
SET     LAS_LIB=C:\programs\compa_libs\liblas-src-1.2.1\lib\Intel11_32\liblas_i.lib
SET  GEOLIB_LIB=C:\programs\compa_libs\GeographicLib-1.16\compileds\VC10_32\lib\Geographic.lib
SET LASZLIB_LIB=C:\programs\compa_libs\lastools\compileds\VC10_32\lib\laslib_i.lib 

)

SET  NETCDF_INC=C:\programs\compa_libs\netcdf\compileds\VC10_32\include
SET     GMT_INC=c:\progs_cygw\GMTdev\GMT4\include
REM SET GMT_INC=c:\progs_cygw\GMTdev\GMT5\include
SET    GMT_INC2=C:\progs_cygw\GMTdev\GMT5\src\mex
SET    GDAL_INC=c:\programs\GDALtrunk\gdal\compileds\VC10_32\include
SET      CV_INC=C:\programs\compa_libs\opencv\compileds\VC10_32\include
SET       CVInc=C:\programs\compa_libs\opencv\compileds\VC10_32\include\opencv
SET  GEOLIB_INC=C:\programs\compa_libs\GeographicLib-1.16\compileds\VC10_64\include
SET LASZLIB_INC=C:\programs\compa_libs\lastools\compileds\VC10_32\include
SET       INCAS=%INCLUDE%
rem SET     LAS_INC=-IC:\programs\compa_libs\liblas-src-1.2.1\bin\include\liblas\capi -IC:\programs\compa_libs\liblas-src-1.2.1\bin\include\liblas
REM ----------------------------------------------------------------------------

REM ____________________________________________________________________________
REM ___________________ STOP EDITING HERE ______________________________________


SET COMPFLAGS=/c /Zp8 /GR /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD /Qopenmp
SET OPTIM2=/QxSSE4.2 /Qparallel /arch:SSE2 /fp:fast 
IF %DEBUG%=="no" SET OPTIMFLAGS=/Ox /Oy- /DNDEBUG /arch:SSE2 /fp:fast
IF %DEBUG%=="yes" SET OPTIMFLAGS=/Z7

IF %BITS%==64 SET arc=X64
IF %BITS%==32 SET arc=X86
SET LINKFLAGS=/dll /export:mexFunction /LIBPATH:%MATLIB% libmx.lib libmex.lib libmat.lib /MACHINE:%arc% kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib Vfw32.lib /nologo /incremental:NO %LDEBUG% 

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% -I%GDAL_INC% -I%CV_INC% -I%CVInc% %LAS_INC% -I%GEOLIB_INC% -I%LASZLIB_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% -DDLL_GMT %OMP% %1
link  /out:"%~n1.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% %GDAL_LIB% %LAS_LIB% %GEOLIB_LIB% %LASZLIB_LIB% %CXCORE_LIB% %CVIMG_LIB% %CVCALIB_LIB% %CVOBJ_LIB% %CVVIDEO_LIB% %CVPHOTO_LIB% /implib:templib.x %~n1.obj
rem link  /out:"%1.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% %GDAL_LIB% %LAS_LIB% %GEOLIB_LIB% /implib:templib.x %1.obj

SET INCLUDE=%INCAS%
del *.obj *.exp templib.x
