#! /bin/sh
#
# Script to fish dynamic libraries and copy them to the Mirone root dir
# and create symbolic links.
# Note that this is crude and very system dependent
# It should only be needed if we want to create a distributable package

lib_loc=/usr/local/lib
gmt_loc=/Users/j/programs/gmt5/lib

# On Mac set it like this
mac_sufix=.dylib
linus_sufix=
# But on Linus, this way
#mac_sufix=
#linus_sufix=.so


# ------------------------- OpenCV libs -------------------------
cp -f ${lib_loc}/libopencv_core${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_core${linus_sufix}.2.4${mac_sufix} libopencv_core${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_imgproc${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_imgproc${linus_sufix}.2.4${mac_sufix} libopencv_imgproc${linus_sufix}${mac_sufix} 

cp -f ${lib_loc}/libopencv_calib3d${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_calib3d${linus_sufix}.2.4${mac_sufix} libopencv_calib3d${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_objdetect${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_objdetect${linus_sufix}.2.4${mac_sufix} libopencv_objdetect${linus_sufix}${mac_sufix} 

cp -f ${lib_loc}/libopencv_video${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_video${linus_sufix}.2.4${mac_sufix} libopencv_video${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_photo${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_photo${linus_sufix}.2.4${mac_sufix} libopencv_photo${linus_sufix}${mac_sufix}


cp -f ${lib_loc}/libopencv_highgui${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_highgui${linus_sufix}.2.4${mac_sufix} libopencv_highgui${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_features2d${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_features2d${linus_sufix}.2.4${mac_sufix} libopencv_features2d${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_flann${linus_sufix}.2.4${mac_sufix} .
ln -s -f libopencv_flann${linus_sufix}.2.4${mac_sufix} libopencv_flann${linus_sufix}${mac_sufix}



# ----------------------- GMT libs ------------------------------
cp -f ${gmt_loc}/libgmt${linus_sufix}.5${mac_sufix} .
ln -s -f libgmt${linus_sufix}.5${mac_sufix} libgmt${linus_sufix}${mac_sufix}

cp -f ${gmt_loc}/libpostscriptlight${linus_sufix}.5${mac_sufix} .
ln -s -f libpostscriptlight${linus_sufix}.5${mac_sufix} libpostscriptlight${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libpcre${linus_sufix}.1${mac_sufix} .
ln -s -f libpcre${linus_sufix}.1${mac_sufix} libpcre${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libfftw3f${linus_sufix}.3${mac_sufix} .
ln -s -f libfftw3f${linus_sufix}.3${mac_sufix} libfftw3f${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libfftw3f_threads${linus_sufix}.3${mac_sufix} .
ln -s -f libfftw3f_threads${linus_sufix}.3${mac_sufix} libfftw3f_threads${linus_sufix}${mac_sufix}


# ------------------------- GDAL libs ---------------------------
cp -f ${lib_loc}/libgdal${linus_sufix}.1${mac_sufix} .
ln -s -f libgdal${linus_sufix}.1${mac_sufix} libgdal${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libproj${linus_sufix}.9${mac_sufix} .
ln -s -f libproj${linus_sufix}.9${mac_sufix} libproj${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libjson-c${linus_sufix}.2${mac_sufix} .
ln -s -f libjson-c${linus_sufix}.2${mac_sufix} libjson-c${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libfreexl${linus_sufix}.1${mac_sufix} .
ln -s -f libfreexl${linus_sufix}.1${mac_sufix} libfreexl${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libgeos_c${linus_sufix}.1${mac_sufix} .
ln -s -f libgeos_c${linus_sufix}.1${mac_sufix} libgeos_c${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libwebp${linus_sufix}.5${mac_sufix} .
ln -s -f libwebp${linus_sufix}.5${mac_sufix} libwebp${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libepsilon${linus_sufix}.1${mac_sufix} .
ln -s -f libepsilon${linus_sufix}.1${mac_sufix} libepsilon${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libodbc${linus_sufix}.2${mac_sufix} .
ln -s -f libodbc${linus_sufix}.2${mac_sufix} libodbc${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libodbcinst${linus_sufix}.2${mac_sufix} .
ln -s -f libodbcinst${linus_sufix}.2${mac_sufix} libodbcinst${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libjasper${linus_sufix}.1${mac_sufix} .
ln -s -f libjasper${linus_sufix}.1${mac_sufix} libjasper${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libnetcdf${linus_sufix}.7${mac_sufix} .
ln -s -f libnetcdf${linus_sufix}.7${mac_sufix} libnetcdf${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libhdf5${linus_sufix}.10${mac_sufix} .
ln -s -f libhdf5${linus_sufix}.10${mac_sufix} libhdf5${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libgif${linus_sufix}.4${mac_sufix} .
ln -s -f libgif${linus_sufix}.4${mac_sufix} libgif${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libjpeg${linus_sufix}.8${mac_sufix} .
ln -s -f libjpeg${linus_sufix}.8${mac_sufix} libjpeg${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libgeotiff${linus_sufix}.2${mac_sufix} .
ln -s -f libgeotiff${linus_sufix}.2${mac_sufix} libgeotiff${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libtiff${linus_sufix}.5${mac_sufix} .
ln -s -f libtiff${linus_sufix}.5${mac_sufix} libtiff${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libpng16${linus_sufix}.16${mac_sufix} .
ln -s -f libpng16${linus_sufix}.16${mac_sufix} libpng16${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libcfitsio${linus_sufix}.2${mac_sufix} .
ln -s -f libcfitsio${linus_sufix}.2${mac_sufix} libcfitsio${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/liblzma${linus_sufix}.5${mac_sufix} .
ln -s -f liblzma${linus_sufix}.5${mac_sufix} liblzma${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libdap${linus_sufix}.11${mac_sufix} .
ln -s -f libdap${linus_sufix}.11${mac_sufix} libdap${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libdapserver${linus_sufix}.7${mac_sufix} .
ln -s -f libdapserver${linus_sufix}.7${mac_sufix} libdapserver${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libdapclient${linus_sufix}.3${mac_sufix} .
ln -s -f libdapclient${linus_sufix}.3${mac_sufix} libdapclient${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libspatialite${linus_sufix}.7${mac_sufix} .
ln -s -f libspatialite${linus_sufix}.7${mac_sufix} libspatialite${linus_sufix}${mac_sufix}

cp -f /usr/local/opt/libxml2/lib/libxml2${linus_sufix}${mac_sufix} .
ln -s -f libxml2${linus_sufix}.2${mac_sufix} libxml2${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libsz${linus_sufix}.2${mac_sufix} .
ln -s -f libsz${linus_sufix}.2${mac_sufix} libsz${linus_sufix}${mac_sufix}

# ---------- ECW ???
#cp -f ${lib_loc}/libNCSCnet${linus_sufix}.0.0.0${mac_sufix} .
#ln -s -f libNCSCnet${linus_sufix}.0.0.0${mac_sufix} libNCSCnet${linus_sufix}${mac_sufix}
#
#cp -f ${lib_loc}/libNCSEcw${linus_sufix}.0.0.0${mac_sufix} .
#ln -s -f libNCSEcw${linus_sufix}.0.0.0${mac_sufix} libNCSEcw${linus_sufix}${mac_sufix}
#
#cp -f ${lib_loc}/libNCSEcwC${linus_sufix}.0.0.0${mac_sufix} .
#ln -s -f libNCSEcwC${linus_sufix}.0.0.0${mac_sufix} libNCSEcwC${linus_sufix}${mac_sufix}
#
#cp -f ${lib_loc}/libNCSUtil${linus_sufix}.0.0.0${mac_sufix} .
#ln -s -f libNCSUtil${linus_sufix}.0.0.0${mac_sufix} libNCSUtil${linus_sufix}${mac_sufix}

cp -f /usr/local/bin/gdalinfo .
cp -f /usr/local/bin/gdal_translate .


# ------------------------- netcdf libs ---------------------------
cp -f ${lib_loc}/libnetcdf${linus_sufix}.7${mac_sufix} .
ln -s -f libnetcdf${linus_sufix}.7${mac_sufix} libnetcdf${linus_sufix}${mac_sufix}


# ---------- HDF5 libs ---------------------------
cp -f ${lib_loc}/libhdf5${linus_sufix}.10${mac_sufix} .
ln -s -f libhdf5${linus_sufix}.10${mac_sufix} libhdf5${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libhdf5_hl${linus_sufix}.10${mac_sufix} .
ln -s -f libhdf5_hl${linus_sufix}.10${mac_sufix} libhdf5_hl${linus_sufix}${mac_sufix}

# --- HDF4 ???
#cp -f ${lib_loc}/libmfhdf${linus_sufix}.0${mac_sufix} .
#ln -s -f libmfhdf${linus_sufix}.0${mac_sufix} libmfhdf${linus_sufix}${mac_sufix}
#
#cp -f ${lib_loc}/libdf${linus_sufix}.0${mac_sufix} .
#ln -s -f libdf${linus_sufix}.0${mac_sufix} libdf${linus_sufix}${mac_sufix}


# ---------- Shapefile libs ------------------------
#cp -f ${lib_loc}/libshp${linus_sufix}.1.2.10${mac_sufix} .
#ln -s -f libshp${linus_sufix}.1.2.10${mac_sufix} libshp${linus_sufix}${mac_sufix}

# ---------- GMT4 libs ------------------------------
#cp -f ${gmt_loc}/libgmt${linus_sufix}.4${mac_sufix} .
#ln -s -f libgmt${linus_sufix}.4${mac_sufix} libgmt${linus_sufix}${mac_sufix}
#
#cp -f ${gmt_loc}/libpsl${linus_sufix}.4${mac_sufix} .
#ln -s -f libpsl${linus_sufix}.4${mac_sufix} libpsl${linus_sufix}${mac_sufix}

