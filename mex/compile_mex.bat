@echo off
REM --------------------------------------------------------------------------------------
REM
REM	$Id:
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
REM Usage: open the command window set up by the compiler of interest (were all vars are already set) and run this from there
REM	   There is no mechanism to build one program individualy but you can build one of the following four groups:
REM		simple, GMT, GDAL, OCV
REM	   To do it, give the group name as argument to this batch. E.G. compile_mex GMT
REM
REM
REM Author: Joaquim Luis, 09-MAY-2010
REM --------------------------------------------------------------------------------------

REM ------------- Set the compiler (set to 'icl' to use the Intel compiler) --------------
SET CC=cl
REM --------------------------------------------------------------------------------------

REM If set to "yes", linkage is done againsts ML6.5 Libs (needed in compiled version)
SET R13="yes"

REM Set it to "yes" or "no" to build under 64-bits or 32-bits respectively.
SET WIN64="no"

REM The MSVC version. I use this var to select libraries also compiled with this compiler
SET MSVC_VER="1600"

REM Options are "dll", "mexw32" (recent ML version scream when they find .dll) or "mexw64" (when WIN64="yes")
SET MEX_EXT="mexw32"

REM If set some MEXs will print the execution time (in CPU ticks)
REM SET TIMEIT=-DMIR_TIMEIT
SET TIMEIT=

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
SET MEX_EXT="mexw32"
) )

 
REM -------------- Set up libraries here -------------------------------------------------
IF %WIN64%=="yes" (
SET NETCDF_LIB=C:\progs_cygw\netcdf-3.6.3\compileds\VC10_64\lib\libnetcdf.lib
SET GMT_LIB=c:\progs_cygw\GMTdev\GMT_win64\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\GMT_win64\lib\gmt_mgg.lib
SET GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_64\lib\gdal_i.lib
SET CV_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\cv.lib
SET CXCORE_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\opencv_core211.lib
SET CVIMG_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\opencv_imgproc211.lib
SET CVCALIB_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\opencv_calib3d211.lib
SET CVOBJ_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\opencv_objdetect211.lib
SET CVVIDEO_LIB=C:\programs\OpenCV_SVN\compileds\VC10_64\lib\opencv_video211.lib

) ELSE (

IF %MSVC_VER%=="1600" (
SET NETCDF_LIB=C:\progs_cygw\netcdf-3.6.3\compileds\VC10_32\lib\libnetcdf.lib
SET GMT_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt_mgg.lib
SET GDAL_LIB=c:\programs\GDALtrunk\gdal\compileds\VC10_32\lib\gdal_i.lib
SET CV_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\cv.lib
SET CXCORE_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\opencv_core211.lib
SET CVIMG_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\opencv_imgproc211.lib
SET CVCALIB_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\opencv_calib3d211.lib
SET CVOBJ_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\opencv_objdetect211.lib
SET CVVIDEO_LIB=C:\programs\OpenCV_SVN\compileds\VC10_32\lib\opencv_video211.lib

) ELSE (

SET NETCDF_LIB=C:\progs_cygw\netcdf-3.6.3\lib\libnetcdf_w32.lib
SET GMT_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt.lib
SET GMT_MGG_LIB=c:\progs_cygw\GMTdev\GMT_win\lib\gmt_mgg.lib
SET GDAL_LIB=c:\programs\GDALtrunk\gdal\lib\gdal_i.lib
REM I haven't build yet (and mayne I won't) 2.1 libs with VC7.1
SET CV_LIB=C:\programs\OpenCV_SVN\lib\cv200.lib
SET CXCORE_LIB=C:\programs\OpenCV_SVN\lib\cxcore200.lib
) )

SET NETCDF_INC=C:\progs_cygw\netcdf-3.6.3\include
SET GMT_INC=c:\progs_cygw\GMTdev\GMT\include
SET GDAL_INC=c:\programs\GDALtrunk\gdal\compileds\VC10_32\include
SET CV_INC=C:\programs\OpenCV_SVN\include\opencv
SET INCLUDE=%INCLUDE%;C:\programs\OpenCV_SVN\modules\core\include;C:\programs\OpenCV_SVN\modules\imgproc\include;C:\programs\OpenCV_SVN\modules\features2d\include;C:\programs\OpenCV_SVN\modules\calib3d\include;C:\programs\OpenCV_SVN\modules\objdetect\include;C:\programs\OpenCV_SVN\modules\video\include
REM ----------------------------------------------------------------------------

REM ____________________________________________________________________________
REM ___________________ STOP EDITING HERE ______________________________________


REM link /out:"test_gmt.mexw32" /dll /export:mexFunction /LIBPATH:"C:\PROGRAMS\MATLAB\R2009B\extern\lib\win32\microsoft" libmx.lib libmex.lib libmat.lib /MACHINE:X86 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /incremental:NO /implib:"C:\TMP\MEX_IF~1\templib.x" /MAP:"test_gmt.mexw32.map"  @C:\TMP\MEX_IF~1\MEX_TMP.RSP   

REM cl  -DWIN32 /c /Zp8 /GR /W3 /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD /FoC:\TMP\MEX_IF~1\test_gmt.obj -IC:\PROGRAMS\MATLAB\R2009B\extern\include /O2 /Oy- /DNDEBUG -DMX_COMPAT_32 test_gmt.c 

SET COMPFLAGS=/c /Zp8 /GR /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD 
SET OPTIMFLAGS=/Ox /Oy- /DNDEBUG 
SET DEBUGFLAGS=/Z7 

IF %WIN64%=="yes" SET arc=X64
IF %WIN64%=="no" SET arc=X86
SET LINKFLAGS=/dll /export:mexFunction /LIBPATH:%MATLIB% libmx.lib libmex.lib libmat.lib /MACHINE:%arc% kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /incremental:NO

IF %1==GMT  GOTO GMT
IF %1==GDAL GOTO GDAL
IF %1==OCV  GOTO OCV

REM ------------------ "simple" (no external Libs dependency) ------------------
:simple
for %%G in (test_gmt igrf_m scaleto8 tsun2 wave_travel_time mansinha_m telha_m range_change country_select 
	mex_illuminate grdutils read_isf alloc_mex susan set_gmt mxgridtrimesh trend1d_m gmtmbgrid_m 
	grdgradient_m grdtrack_m spa_mex mirblock write_mex xyzokb_m ind2rgb8 applylutc cq bwlabel1 
	bwlabel2 imhistc intlutc inv__lwm grayto8 grayto16 ordf parityscan) do ( 

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj
)

for %%G in (ditherc PolygonClip) do ( 

IF %%G==ditherc (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.c invcmap.c  
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj invcmap.obj )

IF %%G==PolygonClip (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.c gpc.c  
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj gpc.obj )

)

REM --------------------- CPPs --------------------------------------------------
for %%G in (houghmex clipbd_mex akimaspline grayxform resampsep iptcheckinput) do (

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.cpp
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj 
)

for %%G in (bwlabelnmex bwboundariesmex imreconstructmex morphmex) do (

IF %%G==bwlabelnmex (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %%G.cpp neighborhood.cpp unionfind.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj neighborhood.obj unionfind.obj )

IF %%G==bwboundariesmex (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.cpp boundaries.cpp 
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj boundaries.obj )

IF %%G==imreconstructmex (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.cpp neighborhood.cpp 
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj neighborhood.obj )

IF %%G==morphmex (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.cpp dilate_erode_gray_nonflat.cpp dilate_erode_packed.cpp dilate_erode_binary.cpp neighborhood.cpp vectors.cpp 
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% /implib:templib.x %%G.obj dilate_erode_gray_nonflat.obj dilate_erode_packed.obj dilate_erode_binary.obj neighborhood.obj vectors.obj )
)
IF %1==simple GOTO END
REM ---------------------- END "simple" --------------------------------------------

REM ---------------------- GMTs ----------------------------------------------------
:GMT
for %%G in (grdinfo_m grdproject_m grdread_m grdsample_m grdtrend_m grdwrite_m mapproject_m shoredump 
	surface_m nearneighbor_m grdfilter_m cpt2cmap grdlandmask_m grdppa_m dimfilter_m shake_mex) do (

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% -DDLL_GMT %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% /implib:templib.x %%G.obj 
)

%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% -I%GMT_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% -DDLL_GMT gmtlist_m.c
link  /out:"gmtlist_m.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% %GMT_LIB% %GMT_MGG_LIB% /implib:templib.x gmtlist_m.obj 
IF %1==GMT GOTO END
REM ---------------------- END "GMTs" ----------------------------------------------


REM ---------------------- GDALs ---------------------------------------------------
:GDAL
for %%G in (gdalread gdalwrite mex_shape) do (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%GDAL_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %GDAL_LIB% /implib:templib.x %%G.obj 
)

for %%G in (gdalwarp_mex gdaltransform_mex ogrproj) do (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%GDAL_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.cpp
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %GDAL_LIB% /implib:templib.x %%G.obj 
)
IF %1==GDAL GOTO END
REM ---------------------- END "GDALs" ----------------------------------------------

REM ---------------------- OpenCV ---------------------------------------------------
:OCV
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%CV_INC% %OPTIMFLAGS% %_MX_COMPAT% cvlib_mex.c sift\sift.c sift\imgfeatures.c sift\kdtree.c sift\minpq.c 
link  /out:"cvlib_mex.%MEX_EXT%" %LINKFLAGS% %CV_LIB% %CXCORE_LIB% %CVIMG_LIB% %CVCALIB_LIB% %CVOBJ_LIB% %CVVIDEO_LIB% /implib:templib.x cvlib_mex.obj sift.obj imgfeatures.obj kdtree.obj minpq.obj 
IF %1==OCV GOTO END


REM ---------------------- with netCDF ----------------------------------------------
for %%G in (swan swan_sem_wbar) do (
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% %%G.c
link  /out:"%%G.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% /implib:templib.x %%G.obj 
)

REM ---------------------- MEXNC ----------------------------------------------------
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% -DDLL_NETCDF mexnc\mexgateway.c mexnc\netcdf2.c mexnc\netcdf3.c mexnc\common.c
link  /out:"mexnc.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% /implib:templib.x mexgateway.obj netcdf2.obj netcdf3.obj common.obj


REM ---------------------- Edison_wrapper -------------------------------------------
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% -I%NETCDF_INC% %OPTIMFLAGS% %_MX_COMPAT% edison/edison_wrapper_mex.cpp edison/segm/ms.cpp edison/segm/msImageProcessor.cpp edison/segm/msSysPrompt.cpp edison/segm/RAList.cpp edison/segm/rlist.cpp edison/edge/BgEdge.cpp edison/edge/BgImage.cpp edison/edge/BgGlobalFc.cpp edison/edge/BgEdgeList.cpp edison/edge/BgEdgeDetect.cpp
link  /out:"edison_wrapper_mex.%MEX_EXT%" %LINKFLAGS% %NETCDF_LIB% /implib:templib.x edison_wrapper_mex.obj Bg*.obj ms*.obj rlist.obj RAList.obj


REM ---------------------- Reduce patches --------------------------------------------
%CC% -DWIN32 %COMPFLAGS% -I%MATINC% %OPTIMFLAGS% %_MX_COMPAT% %TIMEIT% reducep/reducep_s.cpp reducep/3D.cpp reducep/AdjModel.cpp reducep/AdjPrims.cpp reducep/avars.cpp reducep/decimate.cpp reducep/heap.cpp reducep/Mat4.cpp reducep/ProxGrid.cpp reducep/quadrics.cpp reducep/smf.cpp
link  /out:"reducep_s.%MEX_EXT%" %LINKFLAGS% /implib:templib.x reducep_s.obj 3D.obj AdjModel.obj AdjPrims.obj avars.obj decimate.obj heap.obj Mat4.obj ProxGrid.obj quadrics.obj smf.obj

:END
del *.obj *.exp templib.x

