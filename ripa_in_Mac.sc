#! /bin/sh
#
# Script to rip the crazy dynamic linking dependencies created by default on MacOS 
# This is in the hope that it could do any good to the mystery of why MacMirone has
# so many problems in running on other machines than the one it was compiled.
#
# Not all dependencies of OpenCV are treated here as I think some they are never called.
# If I'm wrong, complains will raise up.
# Regarding GDAL, references to libexpat, libcurl & libsqlite3 are not fixed too. I have them
# under /usr/lib but don't know if they are still installed when Xcode is not.

install_name_tool -change /Users/j/programs/gmt5/lib/libgmt.5.dylib libgmt.5.dylib gmtmex.mexmaci64

install_name_tool -id libgmt.dylib libgmt.5.dylib
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib	libnetcdf.dylib	libgmt.dylib
install_name_tool -change /usr/local/lib/libgdal.1.dylib	libgdal.dylib	libgmt.dylib

#install_name_tool -id libpsl.dylib libpsl.4.dylib
install_name_tool -id libpostscriptlight.dylib libpostscriptlight.5.dylib
install_name_tool -id libnetcdf.dylib libnetcdf.7.dylib

#install_name_tool -id libNCSCnet.dylib	libNCSCnet.0.0.0.dylib
#install_name_tool -id libNCSEcw.dylib	libNCSEcw.0.0.0.dylib
#install_name_tool -id libNCSEcwC.dylib	libNCSEcwC.0.0.0.dylib
#install_name_tool -id libNCSUtil.dylib	libNCSUtil.0.0.0.dylib

install_name_tool -id libgdal.dylib libgdal.1.dylib 
#install_name_tool -change /usr/local/lib/libNCSEcw.0.dylib	../libNCSEcw.dylib	libgdal.dylib
#install_name_tool -change /usr/local/lib/libNCSEcwC.0.dylib	../libNCSEcwC.dylib	libgdal.dylib
#install_name_tool -change /usr/local/lib/libNCSCnet.0.dylib	../libNCSCnet.dylib	libgdal.dylib
#install_name_tool -change /usr/local/lib/libNCSUtil.0.dylib	../libNCSUtil.dylib	libgdal.dylib

install_name_tool -change /usr/local/lib/libnetcdf.7.dylib	libnetcdf.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libhdf5.10.dylib	libhdf5.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libhdf5_hl.10.dylib	libhdf5_hl.dylib	libgdal.dylib
#install_name_tool -change /usr/local/lib/libmfhdf.0.dylib	../libmfhdf.dylib	libgdal.dylib
#install_name_tool -change /usr/local/lib/libdf.0.dylib		../libdf.dylib		libgdal.dylib

# ------------------ jpeg????
install_name_tool -id libjpeg.dylib	libjpeg.8.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib	libjpeg.dylib		libgdal.dylib
install_name_tool -id libtiff.dylib	libtiff.5.dylib
install_name_tool -change /usr/local/lib/libtiff.5.dylib	libtiff.dylib		libgdal.dylib

# ------------------ OpenCVs
install_name_tool -id libopencv_core.dylib	libopencv_core.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_core.dylib

install_name_tool -id libopencv_imgproc.dylib	libopencv_imgproc.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib	../libopencv_imgproc.dylib	libopencv_imgproc.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_imgproc.dylib

install_name_tool -id libopencv_objdetect.dylib	libopencv_objdetect.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_objdetect.2.4.dylib	../libopencv_objdetect.dylib	libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib	../libopencv_imgproc.dylib	libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib	../libopencv_highgui.dylib	libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib	../libopencv_features2d.dylib	libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib	../libopencv_calib3d.dylib	libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib	../libopencv_flann.dylib	libopencv_objdetect.dylib

install_name_tool -id libopencv_calib3d.dylib	libopencv_calib3d.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib	../libopencv_calib3d.dylib	libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib	../libopencv_imgproc.dylib	libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib	../libopencv_highgui.dylib	libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib	../libopencv_features2d.dylib	libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib	../libopencv_flann.dylib	libopencv_calib3d.dylib

install_name_tool -id libopencv_video.dylib		libopencv_video.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_video.2.4.dylib	../libopencv_video.dylib	libopencv_video.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_video.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib	../libopencv_imgproc.dylib	libopencv_video.dylib

install_name_tool -id libopencv_highgui.dylib	libopencv_highgui.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib	../libopencv_highgui.dylib	libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib	../libopencv_imgproc.dylib	libopencv_highgui.dylib

install_name_tool -id libopencv_features2d.dylib	libopencv_features2d.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib	../libopencv_features2d.dylib	libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib	../libopencv_imgproc.dylib	libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib	../libopencv_highgui.dylib	libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib	../libopencv_flann.dylib	libopencv_features2d.dylib

install_name_tool -id libopencv_flann.dylib			libopencv_flann.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib	../libopencv_flann.dylib	libopencv_flann.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		libopencv_flann.dylib

cd lib_mex

install_name_tool -change /usr/local/lib/libnetcdf.7.dylib ../libnetcdf.dylib mexnc.mexmaci64
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib ../libnetcdf.dylib swan.mexmaci64

install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib	../libopencv_core.dylib		cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib	../libopencv_imgproc.dylib	cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib	../libopencv_calib3d.dylib	cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_objdetect.2.4.dylib	../libopencv_objdetect.dylib	cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_video.2.4.dylib	../libopencv_video.dylib	cvlib_mex.mexmaci64

#for i in grdinfo_m grdproject_m grdread_m grdsample_m grdtrend_m grdwrite_m mapproject_m shoredump surface_m nearneighbor_m grdfilter_m cpt2cmap grdlandmask_m grdppa_m gmtlist_m shake_mex
#do
	#prg=$i.mexmaci64
	#install_name_tool -change /Users/j/programs/GMT4/lib/libgmt.4.dylib ../libgmt.dylib $prg
	#install_name_tool -change /Users/j/programs/GMT4/lib/libpsl.4.dylib ../libpsl.dylib $prg
	#install_name_tool -change /usr/local/lib/libnetcdf.7.dylib ../libnetcdf.dylib $prg
#done

for i in gdalread gdalwrite gdalwarp_mex ogrproj gdaltransform_mex mex_shape
do
	prg=$i.mexmaci64
	install_name_tool -change /usr/local/lib/libgdal.1.dylib ../libgdal.dylib $prg
done

cd ..
