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
install_name_tool -change /usr/local/lib/libfftw3f.3.dylib	libfftw3f.dylib		libgmt.dylib
install_name_tool -change /usr/local/lib/libfftw3f_threads.3.dylib	libfftw3f_threads.dylib		libgmt.dylib
install_name_tool -change /usr/local/lib/libpcre.1.dylib	libpcre.dylib		libgmt.dylib
install_name_tool -change /Users/j/programs/gmt5/lib/libpostscriptlight.5.dylib	libpostscriptlight.dylib	libgmt.dylib

install_name_tool -id libjpeg.dylib	libjpeg.8.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib	libjpeg.dylib		libgdal.dylib
install_name_tool -id libtiff.dylib	libtiff.5.dylib
install_name_tool -change /usr/local/lib/libtiff.5.dylib	libtiff.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libproj.9.dylib	libproj.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libjson-c.2.dylib	libjson-c.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libfreexl.1.dylib	libfreexl.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libgeos_c.1.dylib	libgeos_c.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libwebp.5.dylib	libwebp.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libepsilon.1.dylib	libepsilon.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libodbc.2.dylib	libodbc.dylib		libgdal.dylib
install_name_tool -change /usr/local/lib/libodbcinst.2.dylib	libodbcinst.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libjasper.1.dylib	libjasper.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libgif.4.dylib		libgif.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libgeotiff.2.dylib	libgeotiff.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libpng16.16.dylib	libpng16.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libcfitsio.2.dylib	libcfitsio.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/liblzma.5.dylib	liblzma.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libdap.11.dylib	libdap.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libdapserver.7.dylib	libdapserver.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libdapclient.3.dylib	libdapclient.dylib	libgdal.dylib
install_name_tool -change /usr/local/lib/libspatialite.7.dylib	libspatialite.dylib	libgdal.dylib
install_name_tool -change /usr/local/opt/sqlite/lib/libsqlite3.0.dylib libsqlite3.dylib libgdal.dylib

install_name_tool -id libcfitsio.dylib	libcfitsio.2.dylib
install_name_tool -id libproj.dylib	libproj.9.dylib
install_name_tool -id libjson-c.dylib	libjson-c.2.dylib
install_name_tool -id libfreexl.dylib	libfreexl.1.dylib
install_name_tool -id libgeos_c.dylib	libgeos_c.1.dylib
install_name_tool -id libwebp.dylib	libwebp.5.dylib
install_name_tool -id libjasper.dylib	libjasper.1.dylib
install_name_tool -id libgeotiff.dylib	libgeotiff.2.dylib
install_name_tool -id liblzma.dylib	liblzma.5.dylib
install_name_tool -id libdap.dylib	libdap.11.dylib
install_name_tool -id libdapserver.dylib	libdapserver.7.dylib
install_name_tool -id libdapclient.dylib	libdapclient.3.dylib
install_name_tool -id libspatialite.dylib	libspatialite.7.dylib
install_name_tool -id libsqlite3.dylib		libsqlite3.0.dylib
install_name_tool -id libpcre.dylib	libpcre.1.dylib
install_name_tool -id libxml2.dylib	libxml2.2.dylib
install_name_tool -id libgeos-3.dylib	libgeos-3.4.2.dylib
install_name_tool -id liblwgeom-2.dylib	liblwgeom-2.1.5.dylib
install_name_tool -id libfftw3f.dylib	libfftw3f.3.dylib
install_name_tool -id libfftw3f_threads.dylib	libfftw3f_threads.3.dylib
install_name_tool -id libpopt.dylib	libpopt.0.dylib
install_name_tool -id libsz.dylib	libsz.2.dylib


install_name_tool -id libhdf5.dylib	libhdf5.10.dylib
install_name_tool -id libhdf5_hl.dylib	libhdf5_hl.10.dylib
install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5.dylib
install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5_hl.dylib
install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5.10.dylib libhdf5.dylib libhdf5_hl.dylib

install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libnetcdf.dylib
install_name_tool -change /usr/local/lib/libhdf5_hl.10.dylib libhdf5_hl.dylib libnetcdf.dylib
install_name_tool -change /usr/local/lib/libhdf5.10.dylib libhdf5.dylib libnetcdf.dylib

install_name_tool -change /usr/local/Cellar/geos/3.4.2/lib/libgeos-3.4.2.dylib libgeos-3.dylib libgeos_c.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib libjpeg.dylib libjasper.dylib

install_name_tool -change /usr/local/lib/libproj.9.dylib libproj.dylib libgeotiff.dylib
install_name_tool -change /usr/local/lib/libtiff.5.dylib libtiff.dylib libgeotiff.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib libjpeg.dylib libgeotiff.dylib

install_name_tool -change /usr/local/lib/libjpeg.8.dylib libjpeg.dylib libtiff.dylib
install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libdap.dylib

install_name_tool -change /usr/local/Cellar/libdap/3.12.1_1/lib/libdap.11.dylib libdap.dylib libdapserver.dylib
install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libdapserver.dylib

install_name_tool -change /usr/local/Cellar/libdap/3.12.1_1/lib/libdap.11.dylib libdap.dylib libdapclient.dylib
install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libdapclient.dylib

install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libspatialite.dylib
install_name_tool -change /usr/local/lib/libfreexl.1.dylib libfreexl.dylib libspatialite.dylib
install_name_tool -change /usr/local/lib/libproj.9.dylib libproj.dylib libspatialite.dylib
install_name_tool -change /usr/local/opt/sqlite/lib/libsqlite3.0.dylib libsqlite3.dylib libspatialite.dylib
install_name_tool -change /usr/local/opt/liblwgeom/lib/liblwgeom-2.1.5.dylib liblwgeom-2.dylib libspatialite.dylib
install_name_tool -change /usr/local/lib/libgeos_c.1.dylib libgeos_c.dylib libspatialite.dylib

install_name_tool -change /usr/local/lib/libgeos_c.1.dylib libgeos_c.dylib liblwgeom-2.dylib
install_name_tool -change /usr/local/lib/libproj.9.dylib libproj.dylib liblwgeom-2.dylib
install_name_tool -change /usr/local/lib/libjson-c.2.dylib libjson-c.dylib liblwgeom-2.dylib

install_name_tool -change /Users/j/programs/gmt5/lib/libgmt.5.dylib libgmt.dylib supplements.dylib
install_name_tool -change /Users/j/programs/gmt5/lib/libpostscriptlight.5.dylib libpostscriptlight.dylib supplements.dylib
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib libnetcdf.dylib supplements.dylib 
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib supplements.dylib
install_name_tool -change /usr/local/lib/libpcre.1.dylib libpcre.dylib supplements.dylib 
install_name_tool -change /usr/local/lib/libfftw3f.3.dylib libfftw3f.dylib supplements.dylib
install_name_tool -change /usr/local/lib/libfftw3f_threads.3.dylib libfftw3f_threads.dylib supplements.dylib 

install_name_tool -change /usr/local/Cellar/fftw/3.3.4_1/lib/libfftw3f.3.dylib libfftw3f.dylib libfftw3f_threads.dylib 

install_name_tool -change /usr/local/lib/libpopt.0.dylib libpopt.dylib libepsilon.dylib 


pref=

# ------------------ OpenCVs
install_name_tool -id libopencv_core.dylib     libopencv_core.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib  ${pref}libopencv_core.dylib libopencv_core.dylib

install_name_tool -id libopencv_imgproc.dylib  libopencv_imgproc.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib ${pref}ibopencv_imgproc.dylib  libopencv_imgproc.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib ${pref}libopencv_core.dylib       libopencv_imgproc.dylib

install_name_tool -id libopencv_objdetect.dylib libopencv_objdetect.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_objdetect.2.4.dylib ${pref}libopencv_objdetect.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib ${pref}libopencv_core.dylib       libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib ${pref}libopencv_imgproc.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib ${pref}libopencv_highgui.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib ${pref}libopencv_features2d.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib ${pref}libopencv_calib3d.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib ${pref}libopencv_flann.dylib     libopencv_objdetect.dylib

install_name_tool -id libopencv_calib3d.dylib   libopencv_calib3d.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib ${pref}libopencv_calib3d.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib ${pref}libopencv_core.dylib       libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib ${pref}libopencv_imgproc.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib ${pref}libopencv_highgui.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib ${pref}libopencv_features2d.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib ${pref}libopencv_flann.dylib     libopencv_calib3d.dylib

install_name_tool -id libopencv_video.dylib     libopencv_video.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_video.2.4.dylib   ${pref}libopencv_video.dylib   libopencv_video.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib    ${pref}libopencv_core.dylib    libopencv_video.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib ${pref}libopencv_imgproc.dylib libopencv_video.dylib

install_name_tool -id libopencv_highgui.dylib    libopencv_highgui.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib  ${pref}libopencv_highgui.dylib    libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib     ${pref}libopencv_core.dylib       libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib  ${pref}libopencv_imgproc.dylib    libopencv_highgui.dylib

install_name_tool -id libopencv_features2d.dylib  libopencv_features2d.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib  ${pref}libopencv_features2d.dylib  libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib        ${pref}libopencv_core.dylib        libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib     ${pref}libopencv_imgproc.dylib     libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib     ${pref}libopencv_highgui.dylib     libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib       ${pref}libopencv_flann.dylib       libopencv_features2d.dylib

install_name_tool -id libopencv_flann.dylib  libopencv_flann.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib  ${pref}libopencv_flann.dylib  libopencv_flann.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib   ${pref}libopencv_core.dylib   libopencv_flann.dylib

install_name_tool -id libopencv_photo.dylib  libopencv_photo.2.4.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib libjpeg.dylib libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libpng16.16.dylib libpng16.16.dylib libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libtiff.5.dylib libtiff.5.dylib libopencv_highgui.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_highgui.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_features2d.2.4.dylib libopencv_features2d.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_flann.2.4.dylib libopencv_flann.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_highgui.2.4.dylib libopencv_highgui.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_calib3d.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_flann.2.4.dylib libopencv_flann.dylib libopencv_features2d.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_highgui.2.4.dylib libopencv_highgui.dylib libopencv_features2d.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_features2d.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_features2d.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_flann.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_imgproc.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_highgui.2.4.dylib libopencv_highgui.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_objdetect.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_photo.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_photo.dylib

install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_video.dylib
install_name_tool -change /usr/local/Cellar/opencv/2.4.12/lib/libopencv_core.2.4.dylib libopencv_core.dylib libopencv_video.dylib


cd lib_mex


install_name_tool -change /usr/local/lib/libnetcdf.7.dylib ${pref}libnetcdf.dylib mexnc.mexmaci64
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib ${pref}libnetcdf.dylib swan.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib    ${pref}libopencv_core.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib ${pref}libopencv_imgproc.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib ${pref}libopencv_calib3d.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_objdetect.2.4.dylib ${pref}libopencv_objdetect.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_video.2.4.dylib  ${pref}libopencv_video.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_photo.2.4.dylib  ${pref}libopencv_photo.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib set_gmt.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib akimaspline.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib alloc_mex.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib country_select.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib gmtmbgrid_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib grdgradient_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib grdtrack_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib grdutils.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib igrf_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib mansinha_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib mex_illuminate.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib mirblock.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib mxgridtrimesh.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib ogrread.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib range_change.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib read_isf.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib scaleto8.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib susan.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib telha_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib test_gmt.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib trend1d_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib tsun2.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib wave_travel_time.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib distmin.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib PolygonClip.mexmaci64


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
	install_name_tool -change /usr/local/lib/libgdal.1.dylib ${pref}libgdal.dylib $prg
done

cd ..
