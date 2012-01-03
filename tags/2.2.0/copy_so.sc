#! /bin/sh
#
# Script to fish dynamic libraries and copy them to the Mirone root dir
# and create symbolic links.
# Note that this is crude and very system dependent
# It should only be needed if we want to create a distributable package

lib_loc=/usr/local/lib
gmt_loc=/Users/j/programs/GMTdev/GMT4/lib

# On Mac set it like this
mac_sufix=.dylib
linus_sufix=
# But on Linus, this way
#mac_sufix=
#linus_sufix=.so

# ---------- OpenCV libs -------------------------
cp -f ${lib_loc}/libopencv_calib3d${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_calib3d${linus_sufix}.2.2.0${mac_sufix} libopencv_calib3d${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_core${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_core${linus_sufix}.2.2.0${mac_sufix} libopencv_core${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_imgproc${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_imgproc${linus_sufix}.2.2.0${mac_sufix} libopencv_imgproc${linus_sufix}${mac_sufix} 

cp -f ${lib_loc}/libopencv_objdetect${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_objdetect${linus_sufix}.2.2.0${mac_sufix} libopencv_objdetect${linus_sufix}${mac_sufix} 

cp -f ${lib_loc}/libopencv_video${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_video${linus_sufix}.2.2.0${mac_sufix} libopencv_video${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_highgui${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_highgui${linus_sufix}.2.2.0${mac_sufix} libopencv_highgui${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_features2d${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_features2d${linus_sufix}.2.2.0${mac_sufix} libopencv_features2d${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libopencv_flann${linus_sufix}.2.2.0${mac_sufix} .
ln -s -f libopencv_flann${linus_sufix}.2.2.0${mac_sufix} libopencv_flann${linus_sufix}${mac_sufix}

# ---------- ?? libs ------------------------------
#cp -f /usr/lib/libexpat${linus_sufix}.1.5.2${mac_sufix} .
#ln -s -f libexpat${linus_sufix}.1.0.0${mac_sufix} libexpat${linus_sufix}.1${mac_sufix}

# ---------- GDAL libs ---------------------------
cp -f ${lib_loc}/libgdal${linus_sufix}.1${mac_sufix} .
ln -s -f libgdal${linus_sufix}.1${mac_sufix} libgdal${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libNCSCnet${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSCnet${linus_sufix}.0.0.0${mac_sufix} libNCSCnet${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libNCSEcw${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSEcw${linus_sufix}.0.0.0${mac_sufix} libNCSEcw${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libNCSEcwC${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSEcwC${linus_sufix}.0.0.0${mac_sufix} libNCSEcwC${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libNCSUtil${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSUtil${linus_sufix}.0.0.0${mac_sufix} libNCSUtil${linus_sufix}${mac_sufix}

cp -f /usr/local/bin/gdalinfo .
cp -f /usr/local/bin/gdal_translate .

# ---------- netcdf libs ---------------------------
cp -f ${lib_loc}/libnetcdf${linus_sufix}.7${mac_sufix} .
ln -s -f libnetcdf${linus_sufix}.7${mac_sufix} libnetcdf${linus_sufix}${mac_sufix}

# ---------- HDF5 libs ---------------------------
cp -f ${lib_loc}/libhdf5${linus_sufix}.7${mac_sufix} .
ln -s -f libhdf5${linus_sufix}.7${mac_sufix} libhdf5${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libhdf5_hl${linus_sufix}.7${mac_sufix} .
ln -s -f libhdf5_hl${linus_sufix}.7${mac_sufix} libhdf5_hl${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libmfhdf${linus_sufix}.0${mac_sufix} .
ln -s -f libmfhdf${linus_sufix}.0${mac_sufix} libmfhdf${linus_sufix}${mac_sufix}

cp -f ${lib_loc}/libdf${linus_sufix}.0${mac_sufix} .
ln -s -f libdf${linus_sufix}.0${mac_sufix} libdf${linus_sufix}${mac_sufix}


# ---------- Don't know why this. Should be an internal dependency to gdal
cp -f /usr/local/Cellar/jpeg/8c/lib/libjpeg.8.dylib .

# ---------- Proj4 libs ----------------------------
#cp -f ${lib_loc}/libproj${linus_sufix}.0.6.6${mac_sufix} .
#ln -s -f libproj${linus_sufix}.0.6.6${mac_sufix} libproj${linus_sufix}${mac_sufix}

# ---------- Shapefile libs ------------------------
#cp -f ${lib_loc}/libshp${linus_sufix}.1.2.10${mac_sufix} .
#ln -s -f libshp${linus_sufix}.1.2.10${mac_sufix} libshp${linus_sufix}${mac_sufix}

# ---------- GMT libs ------------------------------
cp -f ${gmt_loc}/libgmt${linus_sufix}.4${mac_sufix} .
ln -s -f libgmt${linus_sufix}.4${mac_sufix} libgmt${linus_sufix}${mac_sufix}

cp -f ${gmt_loc}/libpsl${linus_sufix}.4${mac_sufix} .
ln -s -f libpsl${linus_sufix}.4${mac_sufix} libpsl${linus_sufix}${mac_sufix}
	
