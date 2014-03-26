@echo off
REM --------------------------------------------------------------------------------------
REM
REM	$Id: $
REM
REM This is a compile batch that builds all MEXs whose source code is distributed in Mirone
REM A major difficulty in using this comes from the fact that several external libraries are needed.
REM They are:
REM		GMT, netCDF, GDAL and OpenCV
REM If a WIN64 version is targeted than the above Libs must have been build in 64-bits as well.
REM
REM Notes: To compile gmtlist_m the file x_system.h must be copyed from GMTROOT\src\x_system to %GMT_INC% (see below)
REM        To compile mex_shape the file shapefil.h must be copyed from GDALROOT\ogr\ogrsf_frmts\shape to %GDAL_INC%
REM
REM Usage: open the command window set up by the compiler of interest (were all vars are already set)
REM	   and run this from there.
REM	   You cannot build one program individualy but you can build one of the following groups:
REM		simple, swan, edison, GMT, GDAL, OCV, MEXNC, laszreader, mpgwrite
REM	   To do it, give the group name as argument to this batch. E.G. compile_mex GMT
REM
REM
REM Author: Joaquim Luis, 09-MAY-2010
REM --------------------------------------------------------------------------------------

REM ------------- Set the compiler (set to 'icl' to use the Intel compiler) --------------
SET CC=cl

REM If set to "yes", linkage is done againsts ML6.5 Libs (needed in compiled version)
SET R13="no"

REM Set it to 32 or 64 to build under 64-bits or 32-bits respectively.
SET BITS=64

REM Set to "yes" if you want to build a debug version
SET DEBUG="0"

REM --------- Most of times changes STOP here (but you may need to setup things bellow) --
REM --------------------------------------------------------------------------------------

REM If set some MEXs will print the execution time (in CPU ticks)
REM SET TIMEIT=-DMIR_TIMEIT
SET TIMEIT=

REM To buils with OpenMP support (very few)
SET OMP=
SET OMP=-DHAVE_OPENMP 

IF %R13%=="yes" SET BITS=32

REM --- Next allows compiling with the compiler you want against the ML6.5 libs (needed in stand-alone version)
IF %R13%=="yes" (
SET MATLIB=C:\SVN\pracompila\MAT65\lib\win32\microsoft
SET MATINC=C:\SVN\pracompila\MAT65\include
rem SET MATLIB=C:\SVN\pracompila\ML2007b_w32\lib\win32\microsoft
rem SET MATINC=C:\SVN\pracompila\ML2007b_w32\include
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
SET MEX_EXT="mexw32"
) )

 
REM -------------- Set up libraries here -------------------------------------------------
IF %BITS%==64 (

SET  NETCDF_LIB=C:\programs\compa_libs\netcdf\compileds\VC10_64\lib\netcdf.lib
SET     GMT_LIB=c:\progs_cygw\GMTdev\gmt4\WIN64\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\gmt4\WIN64\lib\gmt_mgg.lib
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
SET    GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_32\lib4mex\gdal_i.lib
SET  CXCORE_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_core.lib
SET   CVIMG_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_imgproc.lib
SET CVCALIB_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_calib3d.lib
SET   CVOBJ_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_objdetect.lib
SET CVVIDEO_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_video.lib
SET CVPHOTO_LIB=C:\programs\compa_libs\opencv\compileds\VC10_32\lib\opencv_photo.lib
SET     LAS_LIB=C:\programs\compa_libs\liblas-src-1.2.1\lib\VC10_32\liblas_i.lib
SET  GEOLIB_LIB=C:\programs\compa_libs\GeographicLib-1.16\compileds\VC10_32\lib\Geographic.lib
SET LASZLIB_LIB=C:\programs\compa_libs\lastools\compileds\VC10_32\lib\laslib_i.lib 

)

SET  NETCDF_INC=C:\programs\compa_libs\netcdf\compileds\VC10_32\include
SET     GMT_INC=c:\progs_cygw\GMTdev\gmt4\include
SET    GDAL_INC=c:\programs\GDALtrunk\gdal\compileds\VC10_32\include
SET      CV_INC=C:\programs\compa_libs\opencv\compileds\VC10_32\include
SET       CVInc=C:\programs\compa_libs\opencv\compileds\VC10_32\include\opencv
SET  GEOLIB_INC=C:\programs\compa_libs\GeographicLib-1.16\compileds\VC10_64\include
SET LASZLIB_INC=C:\programs\compa_libs\lastools\compileds\VC10_32\include
SET     LAS_INC=-IC:\programs\compa_libs\liblas-src-1.2.1\bin\include\liblas\capi -IC:\programs\compa_libs\liblas-src-1.2.1\bin\include\liblas
REM ----------------------------------------------------------------------------

REM ____________________________________________________________________________
REM ___________________ STOP EDITING HERE ______________________________________

SET LDEBUG=
IF %DEBUG%=="yes" SET LDEBUG=/debug

SET COMPFLAGS=/c /Zp8 /GR /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD /Qopenmp
IF %DEBUG%=="no" SET OPTIMFLAGS=/Ox /DNDEBUG
IF %DEBUG%=="yes" SET OPTIMFLAGS=/Z7

IF %BITS%==64 SET arc=X64
IF %BITS%==32 SET arc=X86
SET LINKFLAGS=/dll /export:mexFunction /LIBPATH:%MATLIB% libmx.lib libmex.lib libmat.lib /MACHINE:%arc% kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /incremental:NO %LDEBUG%

IF "%1" == ""  GOTO todos
IF %1==GMT  GOTO GMT
IF %1==GDAL GOTO GDAL
IF %1==OCV  GOTO OCV
IF %1==MEXNC  GOTO MEXNC
IF %1==MEXNC4  GOTO MEXNC4
IF %1==swan   GOTO swan
IF %1==edison GOTO edison
IF %1==lasreader GOTO lasreader
IF %1==laszreader GOTO laszreader
IF %1==mpgwrite GOTO mpgwrite

:todos
REM ------------------ "simple" (no external Libs dependency) ------------------
:simple
for %%G in (test_gmt igrf_m scaleto8 tsun2 wave_travel_time mansinha_m telha_m range_change country_select 
	mex_illuminate grdutils read_isf alloc_mex susan set_gmt mxgridtrimesh trend1d_m gmtmbgrid_m 
	grdgradient_m grdtrack_m spa_mex mirblock write_mex xyzokb_m distmin CalcMD5 WindowAPI DateStr2Num) do ( 

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %OMP% %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj
)

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %OMP% PolygonClip.c gpc.c  
link  /out:"PolygonClip.%MEX_EXT%" %LINKFLAGS% /implib:templib.x PolygonClip.obj gpc.obj

REM --------------------- CPPs --------------------------------------------------
for %%G in (clipbd_mex akimaspline) do (

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %OMP% %%G.cpp
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj 
)

IF "%1"=="simple" GOTO END
REM ---------------------- END "simple" --------------------------------------------

REM ---------------------- GMTs ----------------------------------------------------
:GMT
for %%G in (grdinfo_m grdproject_m grdread_m grdsample_m grdtrend_m grdwrite_m mapproject_m shoredump 
	surface_m nearneighbor_m grdfilter_m cpt2cmap grdlandmask_m grdppa_m shake_mex) do (

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %OMP% -DDLL_GMT %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% /implib:templib.x %%G.obj 
)

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %OMP% -DDLL_GMT gmtlist_m.c
link  /out:"gmtlist_m.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% %GMT_MGG_LIB% /implib:templib.x gmtlist_m.obj 
IF "%1"=="GMT" GOTO END
REM ---------------------- END "GMTs" ----------------------------------------------


REM ---------------------- GDALs ---------------------------------------------------
:GDAL
for %%G in (gdalread gdalwrite mex_shape ogrread) do (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%GDAL_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %GDAL_LIB% /implib:templib.x %%G.obj 
)

for %%G in (gdalwarp_mex gdaltransform_mex ogrproj) do (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%GDAL_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.cpp
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %GDAL_LIB% /implib:templib.x %%G.obj 
)
IF "%1"=="GDAL" GOTO END
REM ---------------------- END "GDALs" ----------------------------------------------

REM ---------------------- OpenCV ---------------------------------------------------
:OCV
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%CV_INC% -I%CVInc% %OPTIMFLAGS% %_MX_COMPAT% cvlib_mex.cpp sift\sift.c sift\imgfeatures.c sift\kdtree.c sift\minpq.c 
link  /out:"cvlib_mex.%MEX_EXT%" %LINKFLAGS% %CXCORE_LIB% %CVIMG_LIB% %CVCALIB_LIB% %CVOBJ_LIB% %CVVIDEO_LIB% %CVPHOTO_LIB% /implib:templib.x cvlib_mex.obj sift.obj imgfeatures.obj kdtree.obj minpq.obj 
IF "%1"=="OCV" GOTO END


REM ---------------------- with netCDF ----------------------------------------------
:swan
for %%G in (swan swan_sem_wbar) do (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %OMP% %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% /implib:templib.x %%G.obj 
)
IF "%1"=="swan" GOTO END

REM ---------------------- MEXNC4 ---------------------------------------------------
:MEXNC4
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% -DDLL_NETCDF -DNC4_V2_COMPAT mexnc4\mexgateway.c mexnc4\netcdf3.c mexnc4\netcdf4.c mexnc4\common.c
link  /out:"mexnc.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% /implib:templib.x mexgateway.obj netcdf3.obj netcdf4.obj common.obj
IF "%1"=="MEXNC4" GOTO END

REM ---------------------- MEXNC ----------------------------------------------------
:MEXNC
%CC% -DWIN32 %COMPFLAGS% -DUSE_API_2 -I%MATINC% -I%NETCDF_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% -DDLL_NETCDF -DNC4_V2_COMPAT mexnc\mexgateway.c mexnc\netcdf2.c mexnc\netcdf3.c mexnc\common.c
link  /out:"mexnc.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% /implib:templib.x mexgateway.obj netcdf2.obj netcdf3.obj common.obj
IF "%1"=="MEXNC" GOTO END


REM ---------------------- Edison_wrapper -------------------------------------------
:edison
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% %OPTIMFLAGS% %_MX_COMPAT% edison/edison_wrapper_mex.cpp edison/segm/ms.cpp edison/segm/msImageProcessor.cpp edison/segm/msSysPrompt.cpp edison/segm/RAList.cpp edison/segm/rlist.cpp edison/edge/BgEdge.cpp edison/edge/BgImage.cpp edison/edge/BgGlobalFc.cpp edison/edge/BgEdgeList.cpp edison/edge/BgEdgeDetect.cpp
link  /out:"edison_wrapper_mex.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% /implib:templib.x edison_wrapper_mex.obj Bg*.obj ms*.obj rlist.obj RAList.obj
IF "%1"=="edison" GOTO END


REM ---------------------- LASlib (laszreader) ----------------------------------------
:laszreader
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%LASZLIB_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% laszreader_mex.cpp
link  /out:"laszreader_mex.%MEX_EXT%" %LINKFLAGS% %LASZLIB_LIB% /implib:templib.x laszreader_mex.obj
IF "%1"=="laszreader" GOTO END


REM ---------------------- libLAS (lasreader) ----------------------------------------
:lasreader
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %LAS_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% lasreader_mex.c
link  /out:"lasreader_mex.%MEX_EXT%" %LINKFLAGS% %LAS_LIB% /implib:templib.x lasreader_mex.obj
IF "%1"=="lasreader" GOTO END

REM ---------------------- MPEG (mpgwrite) -------------------------------------------
:mpgwrite
set P=mpgwrite\
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %P%mpgwrite.c %P%mfwddct.c %P%postdct.c %P%huff.c %P%bitio.c %P%mheaders.c %P%iframe.c %P%pframe.c %P%bframe.c %P%psearch.c %P%bsearch.c %P%block.c %P%mpeg.c %P%subsampl.c %P%jrevdct.c %P%frame.c %P%fsize.c
link  /out:"mpgwrite.%MEX_EXT%" %LINKFLAGS% /implib:templib.x mpgwrite.obj mfwddct.obj postdct.obj huff.obj bitio.obj mheaders.obj iframe.obj pframe.obj bframe.obj psearch.obj bsearch.obj block.obj mpeg.obj subsampl.obj jrevdct.obj frame.obj fsize.obj 
IF "%1"=="mpgwrite" GOTO END


:END
del *.obj *.exp templib.x
