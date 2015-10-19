#! /bin/sh
#
# This script if the first of a triology:
# copy_patch_simple.sh, copy_patch_with_hide.sh, patch_mexs.sh   
# Run them in that order from a sub-directory of the mirone root dir
#
# What they do is to copy several dylibs from /usr/local/lib (built with Homebrew) and strip the
# hard-coded paths. Besides that, several of those shared libs are assigned a different name
# (by appending a '_hide') so that Matlab does not managed to f. in the middle with its own outdated
# versions. For example libhdf5.dylib will become libhdf5_hide.dylib

# Copy and patch the dylibs that don't need the replaced names. Those are dealt by "copy_patch_with_hide.sh"

# $Id :

lib_loc=/usr/local/lib

# ------------------------- OpenCV libs -------------------------
cp -f ${lib_loc}/libopencv_core.2.4.dylib .
ln -s -f libopencv_core.2.4.dylib libopencv_core.dylib

cp -f ${lib_loc}/libopencv_imgproc.2.4.dylib .
ln -s -f libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib 

cp -f ${lib_loc}/libopencv_calib3d.2.4.dylib .
ln -s -f libopencv_calib3d.2.4.dylib libopencv_calib3d.dylib

cp -f ${lib_loc}/libopencv_objdetect.2.4.dylib .
ln -s -f libopencv_objdetect.2.4.dylib libopencv_objdetect.dylib 

cp -f ${lib_loc}/libopencv_video.2.4.dylib .
ln -s -f libopencv_video.2.4.dylib libopencv_video.dylib

cp -f ${lib_loc}/libopencv_photo.2.4.dylib .
ln -s -f libopencv_photo.2.4.dylib libopencv_photo.dylib

cp -f ${lib_loc}/libopencv_highgui.2.4.dylib .
ln -s -f libopencv_highgui.2.4.dylib libopencv_highgui.dylib

cp -f ${lib_loc}/libopencv_features2d.2.4.dylib .
ln -s -f libopencv_features2d.2.4.dylib libopencv_features2d.dylib

cp -f ${lib_loc}/libopencv_flann.2.4.dylib .
ln -s -f libopencv_flann.2.4.dylib libopencv_flann.dylib


# ------------------------- GDAL libs ---------------------------
# WE don't copy those that are dealt by "copy_patch_with_hide.sh"

cp -f ${lib_loc}/libjson-c.2.dylib .
ln -s -f libjson-c.2.dylib libjson-c.dylib

cp -f ${lib_loc}/libfreexl.1.dylib .
ln -s -f libfreexl.1.dylib libfreexl.dylib

cp -f ${lib_loc}/libgeos_c.1.dylib .
ln -s -f libgeos_c.1.dylib libgeos_c.dylib

cp -f ${lib_loc}/libwebp.5.dylib .
ln -s -f libwebp.5.dylib libwebp.dylib

cp -f ${lib_loc}/libepsilon.1.dylib .
ln -s -f libepsilon.1.dylib libepsilon.dylib

cp -f ${lib_loc}/libodbc.2.dylib .
ln -s -f libodbc.2.dylib libodbc.dylib

cp -f ${lib_loc}/libodbcinst.2.dylib .
ln -s -f libodbcinst.2.dylib libodbcinst.dylib

cp -f ${lib_loc}/libjasper.1.dylib .
ln -s -f libjasper.1.dylib libjasper.dylib

cp -f ${lib_loc}/libgif.4.dylib .
ln -s -f libgif.4.dylib libgif.dylib

cp -f ${lib_loc}/libjpeg.8.dylib .
ln -s -f libjpeg.8.dylib libjpeg.dylib

cp -f ${lib_loc}/libtiff.5.dylib .
ln -s -f libtiff.5.dylib libtiff.dylib

cp -f ${lib_loc}/libpng16.16.dylib .
ln -s -f libpng16.16.dylib libpng16.dylib

cp -f ${lib_loc}/liblzma.5.dylib .
ln -s -f liblzma.5.dylib liblzma.dylib

cp -f ${lib_loc}/libdap.11.dylib .
ln -s -f libdap.11.dylib libdap.dylib

cp -f ${lib_loc}/libdapserver.7.dylib .
ln -s -f libdapserver.7.dylib libdapserver.dylib

cp -f ${lib_loc}/libdapclient.3.dylib .
ln -s -f libdapclient.3.dylib libdapclient.dylib

cp -f ${lib_loc}/libspatialite.7.dylib .
ln -s -f libspatialite.7.dylib libspatialite.dylib

cp -f /usr/local/opt/libxml2/lib/libxml2.2.dylib .
ln -s -f libxml2.2.dylib libxml2.dylib

cp -f /usr/local/opt/sqlite/lib/libsqlite3.0.dylib .
ln -s -f libsqlite3.0.dylib libsqlite3.dylib

cp -f ${lib_loc}/libsz.2.dylib .
ln -s -f libsz.2.dylib libsz.dylib

cp /usr/local/Cellar/geos/3.4.2/lib/libgeos-3.4.2.dylib .
ln -s -f libgeos-3.4.2.dylib libgeos-3.dylib 

cp /usr/local/opt/liblwgeom/lib/liblwgeom-2.1.5.dylib .
ln -s -f liblwgeom-2.1.5.dylib liblwgeom-2.dylib

cp /usr/local/lib/libpopt.0.dylib .
ln -s -f libpopt.0.dylib libpopt.dylib

cp -f ${lib_loc}/libpcre.1.dylib .
ln -s -f libpcre.1.dylib libpcre.dylib

cp -f ${lib_loc}/libfftw3f.3.dylib .
ln -s -f libfftw3f.3.dylib libfftw3f.dylib

cp -f ${lib_loc}/libfftw3f_threads.3.dylib .
ln -s -f libfftw3f_threads.3.dylib libfftw3f_threads.dylib

cp /usr/local/lib/libImath-2_2.12.dylib .
ln -s -f libImath-2_2.12.dylib libImath-2_2.dylib

cp /usr/local/lib/libIlmImf-2_2.22.dylib .
ln -s -f libIlmImf-2_2.22.dylib libIlmImf-2_2.dylib

cp /usr/local/lib/libIex-2_2.12.dylib .
ln -s -f libIex-2_2.12.dylib libIex-2_2.dylib

cp /usr/local/lib/libHalf.12.dylib .
ln -s -f libHalf.12.dylib libHalf.dylib

cp /usr/local/lib/libIlmThread-2_2.12.dylib .
ln -s -f libIlmThread-2_2.12.dylib libIlmThread-2_2.dylib

chmod +w *.dylib

# GDAL exes (Shit, these carry again the whole trunk of dependencies to patch)
#cp -f /usr/local/bin/gdalinfo .
#cp -f /usr/local/bin/gdal_translate .

# ------------------------------------------------------------------
#               NOW THE PATCHING PART
# ------------------------------------------------------------------


install_name_tool -id libjpeg.dylib        libjpeg.8.dylib
install_name_tool -id libtiff.dylib        libtiff.5.dylib
install_name_tool -id libpng16.dylib       libpng16.16.dylib
install_name_tool -id libjson-c.dylib      libjson-c.2.dylib
install_name_tool -id libfreexl.dylib      libfreexl.1.dylib
install_name_tool -id libgeos_c.dylib      libgeos_c.1.dylib
install_name_tool -id libwebp.dylib        libwebp.5.dylib
install_name_tool -id libjasper.dylib      libjasper.1.dylib
install_name_tool -id liblzma.dylib        liblzma.5.dylib
install_name_tool -id libdap.dylib         libdap.11.dylib
install_name_tool -id libdapserver.dylib   libdapserver.7.dylib
install_name_tool -id libdapclient.dylib   libdapclient.3.dylib
install_name_tool -id libspatialite.dylib  libspatialite.7.dylib
install_name_tool -id libsqlite3.dylib     libsqlite3.0.dylib
install_name_tool -id libpcre.dylib        libpcre.1.dylib
install_name_tool -id libxml2.dylib        libxml2.2.dylib
install_name_tool -id libgeos-3.dylib      libgeos-3.4.2.dylib
install_name_tool -id liblwgeom-2.dylib    liblwgeom-2.1.5.dylib
install_name_tool -id libfftw3f.dylib      libfftw3f.3.dylib
install_name_tool -id libfftw3f_threads.dylib libfftw3f_threads.3.dylib
install_name_tool -id libpopt.dylib        libpopt.0.dylib
install_name_tool -id libsz.dylib          libsz.2.dylib
install_name_tool -id libepsilon.dylib     libepsilon.1.dylib
install_name_tool -id libgif.dylib         libgif.4.dylib
install_name_tool -id libodbc.dylib        libodbc.2.dylib
install_name_tool -id libodbcinst.dylib    libodbcinst.2.dylib
install_name_tool -id libImath-2_2.dylib   libImath-2_2.12.dylib
install_name_tool -id libIlmImf-2_2.dylib  libIlmImf-2_2.22.dylib
install_name_tool -id libIex-2_2.dylib     libIex-2_2.12.dylib
install_name_tool -id libHalf.dylib        libHalf.12.dylib
install_name_tool -id libIlmThread-2_2.dylib libIlmThread-2_2.12.dylib

install_name_tool -change /usr/local/Cellar/geos/3.4.2/lib/libgeos-3.4.2.dylib libgeos-3.dylib libgeos_c.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib libjpeg.dylib libjasper.dylib

install_name_tool -change /usr/local/lib/libjpeg.8.dylib libjpeg.dylib libtiff.dylib
install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libdap.dylib

install_name_tool -change /usr/local/Cellar/libdap/3.12.1_1/lib/libdap.11.dylib libdap.dylib libdapserver.dylib
install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libdapserver.dylib

install_name_tool -change /usr/local/Cellar/libdap/3.12.1_1/lib/libdap.11.dylib libdap.dylib libdapclient.dylib
install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libdapclient.dylib

install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib libspatialite.dylib
install_name_tool -change /usr/local/lib/libfreexl.1.dylib libfreexl.dylib libspatialite.dylib
install_name_tool -change /usr/local/lib/libproj.9.dylib libproj_hide.dylib libspatialite.dylib
install_name_tool -change /usr/local/opt/sqlite/lib/libsqlite3.0.dylib libsqlite3.dylib libspatialite.dylib
install_name_tool -change /usr/local/opt/liblwgeom/lib/liblwgeom-2.1.5.dylib liblwgeom-2.dylib libspatialite.dylib
install_name_tool -change /usr/local/lib/libgeos_c.1.dylib libgeos_c.dylib libspatialite.dylib

install_name_tool -change /usr/local/lib/libgeos_c.1.dylib libgeos_c.dylib   liblwgeom-2.dylib
install_name_tool -change /usr/local/lib/libproj.9.dylib libproj_hide.dylib  liblwgeom-2.dylib
install_name_tool -change /usr/local/lib/libjson-c.2.dylib libjson-c.dylib   liblwgeom-2.dylib

install_name_tool -change /usr/local/Cellar/fftw/3.3.4_1/lib/libfftw3f.3.dylib libfftw3f.dylib libfftw3f_threads.dylib 

install_name_tool -change /usr/local/lib/libpopt.0.dylib libpopt.dylib libepsilon.dylib 

# ------------------ OpenCVs
install_name_tool -id libopencv_core.dylib     libopencv_core.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib  libopencv_core.dylib libopencv_core.dylib

install_name_tool -id libopencv_imgproc.dylib  libopencv_imgproc.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib    ibopencv_imgproc.dylib     libopencv_imgproc.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib       libopencv_core.dylib       libopencv_imgproc.dylib

install_name_tool -id libopencv_objdetect.dylib libopencv_objdetect.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_objdetect.2.4.dylib  libopencv_objdetect.dylib  libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib       libopencv_core.dylib       libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib    libopencv_imgproc.dylib    libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib    libopencv_highgui.dylib    libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib libopencv_features2d.dylib libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib    libopencv_calib3d.dylib    libopencv_objdetect.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib      libopencv_flann.dylib      libopencv_objdetect.dylib

install_name_tool -id libopencv_calib3d.dylib   libopencv_calib3d.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib libopencv_calib3d.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib libopencv_core.dylib       libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib libopencv_highgui.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib libopencv_features2d.dylib libopencv_calib3d.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib libopencv_flann.dylib     libopencv_calib3d.dylib

install_name_tool -id libopencv_video.dylib     libopencv_video.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_video.2.4.dylib   libopencv_video.dylib   libopencv_video.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib    libopencv_core.dylib    libopencv_video.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib libopencv_video.dylib

install_name_tool -id libopencv_highgui.dylib    libopencv_highgui.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib  libopencv_highgui.dylib    libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib     libopencv_core.dylib       libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib  libopencv_imgproc.dylib    libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libImath-2_2.12.dylib        libImath-2_2.dylib         libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libIlmImf-2_2.22.dylib       libIlmImf-2_2.dylib        libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libIex-2_2.12.dylib          libIex-2_2.dylib           libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libHalf.12.dylib             libHalf-2_2.dylib          libopencv_highgui.dylib
install_name_tool -change /usr/local/lib/libIlmThread-2_2.12.dylib    libIlmThread-2_2.dylib     libopencv_highgui.dylib

install_name_tool -id libopencv_features2d.dylib  libopencv_features2d.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_features2d.2.4.dylib  libopencv_features2d.dylib  libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib        libopencv_core.dylib        libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib     libopencv_imgproc.dylib     libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_highgui.2.4.dylib     libopencv_highgui.dylib     libopencv_features2d.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib       libopencv_flann.dylib       libopencv_features2d.dylib

install_name_tool -id libopencv_flann.dylib  libopencv_flann.2.4.dylib
install_name_tool -change /usr/local/lib/libopencv_flann.2.4.dylib  libopencv_flann.dylib  libopencv_flann.dylib
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib   libopencv_core.dylib   libopencv_flann.dylib

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

install_name_tool -change /usr/local/Cellar/ilmbase/2.2.0/lib/libIex-2_2.12.dylib libIex-2_2.dylib libImath-2_2.dylib

install_name_tool -change /usr/local/lib/libImath-2_2.12.dylib     libImath-2_2.dylib     libIlmImf-2_2.dylib
install_name_tool -change /usr/local/lib/libHalf.12.dylib          libHalf.dylib          libIlmImf-2_2.dylib
install_name_tool -change /usr/local/lib/libIex-2_2.12.dylib       libIex-2_2.dylib       libIlmImf-2_2.dylib
install_name_tool -change /usr/local/lib/libIexMath-2_2.12.dylib   libIexMath-2_2.dylib   libIlmImf-2_2.dylib
install_name_tool -change /usr/local/lib/libIlmThread-2_2.12.dylib libIlmThread-2_2.dylib libIlmImf-2_2.dylib

install_name_tool -change /usr/local/Cellar/ilmbase/2.2.0/lib/libIex-2_2.12.dylib libIex-2_2.dylib libIlmThread-2_2.dylib





