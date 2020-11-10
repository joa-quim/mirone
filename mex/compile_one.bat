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
SET R13="yes"

REM Set it to 32 or 64 to build under 64-bits or 32-bits respectively.
SET BITS=32

REM Set to "yes" if you want to build a debug version
SET DEBUG="no"

REM The compiler version. Used to set up the lin/include paths
SET VC=14

REM --------- Most of times changes STOP here (but you may need to setup things bellow) --
REM --------------------------------------------------------------------------------------

REM If set some MEXs will print the execution time (in CPU ticks)
SET TIMEIT=-DMIR_TIMEIT
SET TIMEIT=

REM To buils with OpenMP support (very few)
SET OMP=
SET OMP=-DHAVE_OPENMP 

IF %R13%=="yes" SET BITS=32

REM Set to "yes" if when you want to build cvlib_mex that needs extra links
SET PROCV="no"

REM Set to "yes" if when you want to build imregionalmaxmex that needs extra links
SET IMREGMAX="no"

REM --- Next allows compiling with the compiler you want against the ML6.5 libs (needed in stand-alone version)
IF %R13%=="yes" (
SET MATLIB=C:\SVN\pracompila\MAT65\lib\win32\microsoft
SET MATINC=C:\SVN\pracompila\MAT65\include
SET _MX_COMPAT=
SET MEX_EXT="dll"

) ELSE (

IF %BITS%==64 (
SET MATLIB=C:\SVN\pracompila\ML2011b_w64\lib\win64\microsoft 
SET MATINC=C:\SVN\pracompila\ML2011b_w64\include
SET _MX_COMPAT=-DMX_COMPAT_32
SET MEX_EXT="mexw64"

) ELSE (

SET MATLIB=C:\SVN\pracompila\ML2012a_w32\lib\win32\microsoft 
SET MATINC=C:\SVN\pracompila\ML2012a_w32\include
SET _MX_COMPAT=-DMX_COMPAT_32
SET MEX_EXT="mexw32"
) )

 
REM -------------- Set up libraries here -------------------------------------------------

SET  NETCDF_LIB=C:\programs\compa_libs\netcdf_GIT\compileds\VC%VC%_%BITS%\lib\netcdf.lib
SET     GMT_LIB=c:\progs_cygw\GMTdev\gmt5\compileds\gmt6\VC%VC%_%BITS%\lib\gmt.lib
IF %VC% == 14 (
SET    GDAL_LIB=c:\programs\compa_libs\gdal_GIT\compileds\VC%VC%_%BITS%\lib\gdal_i.lib
) ELSE (
SET    GDAL_LIB=c:\programs\compa_libs\gdal\compileds\VC%VC%_%BITS%\lib\gdal_i.lib
)
SET  CXCORE_LIB=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\lib\opencv_core.lib
SET   CVIMG_LIB=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\lib\opencv_imgproc.lib
SET CVCALIB_LIB=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\lib\opencv_calib3d.lib
SET   CVOBJ_LIB=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\lib\opencv_objdetect.lib
SET CVVIDEO_LIB=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\lib\opencv_video.lib
SET CVPHOTO_LIB=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\lib\opencv_photo.lib
SET     LAS_LIB=C:\programs\compa_libs\liblas-src-1.2.1\lib\VC10_%BITS%\liblas_i.lib
SET  GEOLIB_LIB=C:\programs\compa_libs\GeographicLib-1.49\compileds\VC%VC%_%BITS%\lib\Geographic_i.lib
SET LASZLIB_LIB=C:\programs\compa_libs\lastools\compileds\VC%VC%_%BITS%\lib\laslib_i.lib 
SET   JULIA_LIB=V:\julia\usr\bin\julia.lib 

SET   JULIA_LIB=

SET  NETCDF_INC=C:\programs\compa_libs\netcdf_GIT\compileds\VC%VC%_%BITS%\include
SET     GMT_INC=c:\progs_cygw\GMTdev\gmt5\compileds\gmt6\VC%VC%_%BITS%\include\gmt
IF %VC% == 14 (
SET    GDAL_INC=c:\programs\compa_libs\gdal_GIT\compileds\VC%VC%_%BITS%\include
) ELSE (
SET    GDAL_INC=c:\programs\compa_libs\gdal\compileds\VC%VC%_%BITS%\include
)
SET      CV_INC=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\include
SET       CVInc=C:\programs\compa_libs\opencv\compileds\VC%VC%_%BITS%\include\opencv
SET  GEOLIB_INC=C:\programs\compa_libs\GeographicLib-1.49\compileds\VC%VC%_%BITS%\include
SET LASZLIB_INC=C:\programs\compa_libs\lastools\compileds\VC%VC%_%BITS%\include

REM SET LASZLIB_INC=C:\programs\compa_libs\lastools_GIT\compileds\VC12_64\include
REM SET LASZLIB_INC2=C:\programs\compa_libs\lastools_GIT\LASlib\inc
SET   JULIA_INC=-IC:\programs\julia64\src -IC:\programs\julia64\src\support -IC:\programs\julia64\usr\include
SET       INCAS=%INCLUDE%
REM ----------------------------------------------------------------------------

REM ____________________________________________________________________________
REM ___________________ STOP EDITING HERE ______________________________________

SET LDEBUG=
IF %DEBUG%=="yes" SET LDEBUG=/debug

IF %PROCV%=="yes" (
SET extra_cv_c=sift\sift.c sift\imgfeatures.c sift\kdtree.c sift\minpq.c 
SET extra_cv_o=sift.obj imgfeatures.obj kdtree.obj minpq.obj 
) ELSE (
SET extra_cv_c=
SET extra_cv_o=
)

SET extra_IMREGMAX=
IF %IMREGMAX%=="yes" (
SET extra_IMREGMAX_src=neighborhood.cpp
SET extra_IMREGMAX_obj=neighborhood.obj
)

SET OMP_F=/openmp
IF %CC%==icl SET OMP_F=/Qopenmp
SET COMPFLAGS=/c /Zp8 /GR /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD %OMP_F% /DHAVE_NETCDF 
SET OPTIM2=/QxSSE4.2 /Qparallel
IF %DEBUG%=="no" SET OPTIMFLAGS=/Ox /fp:precise /DNDEBUG /Ox /Qprec-div-
IF %DEBUG%=="yes" SET OPTIMFLAGS=/Z7

IF %BITS%==64 SET arc=X64
IF %BITS%==32 SET arc=X86
SET LINKFLAGS=/dll /export:mexFunction /LIBPATH:%MATLIB% libmx.lib libmex.lib libmat.lib /MACHINE:%arc% kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib Vfw32.lib /nologo /incremental:NO %LDEBUG% 

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% -I%GDAL_INC% -I%CV_INC% -I%CVInc% -I%GEOLIB_INC% -I%LASZLIB_INC2% -I%LASZLIB_INC% %JULIA_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% -DDLL_GMT %OMP% %extra_cv_c% %extra_IMREGMAX_src% %1

link  /out:"%~n1.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% %GDAL_LIB% %LAS_LIB% %GEOLIB_LIB% %LASZLIB_LIB% %CXCORE_LIB% %CVIMG_LIB% %CVCALIB_LIB% %CVOBJ_LIB% %CVVIDEO_LIB% %CVPHOTO_LIB% %JULIA_LIB% /implib:templib.x %~n1.obj %extra_cv_o% %extra_IMREGMAX_obj%

SET INCLUDE=%INCAS%
del *.obj *.exp templib.x
