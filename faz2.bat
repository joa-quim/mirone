@echo off
REM
REM Copy the various DLLs from their original location 

REM -----------------------------------------------------------------------------------------------
REM If set to "yes", fish the MEXs with the old .dll extension
SET R13="no"

REM Set to 32 or 64 to the corresponding bitage version
SET BITS=32

IF %R13%=="yes" SET BITS=32

REM My directory with simlinks to the original DLLs
SET P=c:\programs\compa_libs\DLLs_%BITS%\
SET Pj=c:\j\bin\

REM The GDAL location 
SET Pgd=C:\programs\GDALtrunk\gdal\compileds\VC10_%BITS%\bin\
REM ------------------------------------------------------------------------------------------------

IF %BITS%==64 (
REM The VS CRT runtime DLLs
SET PVC="C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\redist\x64\Microsoft.VC100.CRT\"
SET P2=C:\SVN\mironeWC64\
SET MEX_EXT="mexw64"
) ELSE (
SET PVC="C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\redist\x86\Microsoft.VC100.CRT\"
SET P2=C:\SVN\mironeWC\
SET MEX_EXT="mexw32"
)
IF %R13%=="yes" SET MEX_EXT="dll"

md prozip  &&  cd prozip

copy /Y %P%gdal_w%BITS%.dll .      && copy /Y %P%netcdf4_w%BITS%.dll .  && copy /Y %P%tbb_w%BITS%.dll .
copy /Y %P%hdf4_w%BITS%.dll .      && copy /Y %P%mfhdf4_w%BITS%.dll .   && copy /Y %P%xdr_w%BITS%.dll .
copy /Y %P%hdf5_hl_w%BITS%.dll .   && copy /Y %P%hdf5_w%BITS%.dll .     && copy /Y %P%xdrbsd_w%BITS%.dll .
copy /Y %P%libjpeg_w%BITS%.dll .   && copy /Y %P%openjp2_w%BITS%.dll .  && copy /Y %P%laslib_w%BITS%.dll .
copy /Y %P%geos_c_w%BITS%.dll .    && copy /Y %P%zlib1_w%BITS%.dll .    && copy /Y %P%xerces-c_3_1_w%BITS%.dll .
copy /Y %P%libcurl_w%BITS%.dll .   && copy /Y %P%proj_w%BITS%.dll .     && copy /Y %P%expat_w%BITS%.dll .
copy /Y %P%libecwj2_w%BITS%.dll .  && copy /Y %P%geos_w%BITS%.dll .

copy /Y %P%opencv_calib3d_w%BITS%.dll .    && copy /Y %P%opencv_core_w%BITS%.dll .
copy /Y %P%opencv_highgui_w%BITS%.dll .    && copy /Y %P%opencv_imgproc_w%BITS%.dll .
copy /Y %P%opencv_objdetect_w%BITS%.dll .  && copy /Y %P%opencv_video_w%BITS%.dll .
copy /Y %P%opencv_flann_w%BITS%.dll .      && copy /Y %P%opencv_features2d_w%BITS%.dll .
copy /Y %P%opencv_photo_w%BITS%.dll .      && copy /Y %P%lti_dsdk.dll .
copy /Y %P%cfitsio_w%BITS%.dll .

copy /Y %P2%libiomp5md.dll .  && copy /Y %P2%libmmd.dll .
copy /Y %P2%country_extract.exe .
copy /Y %P2%gmt.dll .         && copy /Y %P2%gmt_mgg.dll .   && copy /Y %P2%psl.dll .
copy /Y %Pj%gunzip.exe .      && copy /Y %Pj%unzip.exe .     && copy /Y %Pj%wget.exe .
copy /Y %Pgd%gdalinfo.exe .   && copy /Y %Pgd%gdal_translate.exe .
copy /Y %PVC%msvcp100.dll .   && copy /Y %PVC%msvcr100.dll .

REM ------------- Copy the MEXs and the P codes --------------------------------------
md utils  && md lib_mex  &&  cd utils

IF %R13%=="yes" (
copy /Y %P2%lib_mex_dll\*.p .
cd ..\lib_mex
copy /Y %P2%lib_mex_dll\*.dll .
SET R13="R13"

) ELSE (

copy /Y %P2%utils\*.p .
cd ..\lib_mex
copy /Y %P2%lib_mex\*.%MEX_EXT% .
IF %BITS%==32 copy /Y %P2%lib_mex\*.dll .
SET R13=
)
REM ----------------------------------------------------------------------------------

cd ..

REM Create a simple readme
echo unzip the contents of this file to the Mirone's root directory > READEME.txt

zip -q -r ..\mir%R13%_binaries_X%BITS% .

cd ..

rd /S /Q prozip
pause
