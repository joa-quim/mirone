#! /bin/sh
#
# Script to fish dynamic libraries and copy them to the Mirone root dir
# and create symbolic links.
# Note that this is crude and very system dependent
# It should only be needed if we want to create a destributable package

lib_loc=/usr/local/lib
gmt_loc=/Users/j/programs/GMT/lib

# On Mac set it like this
mac_sufix=.dylib
linus_sufix=
# But on Linus, this way
#mac_sufix=
#linus_sufix=.so

# ---------- OpenCV libs -------------------------
cp -f ${lib_loc}/libcv${linus_sufix}.1.0.0${mac_sufix} .
ln -s -f libcv${linus_sufix}.1.0.0${mac_sufix} libcv${linus_sufix}.1${mac_sufix}

cp -f ${lib_loc}/libcxcore${linus_sufix}.1.0.0${mac_sufix} .
ln -s -f libcxcore${linus_sufix}.1.0.0 libcxcore${linus_sufix}.1${mac_sufix}

# ---------- ?? libs ------------------------------
#cp -f /usr/lib/libexpat${linus_sufix}.1.5.2${mac_sufix} .
#ln -s -f libexpat${linus_sufix}.1.0.0 libexpat${linus_sufix}.1${mac_sufix}

# ---------- GDAL libs ---------------------------
cp -f ${lib_loc}/libgdal${linus_sufix}.1${mac_sufix} .
ln -s -f libgdal${linus_sufix}.1 libgdal${linus_sufix}.1${mac_sufix}

cp -f ${lib_loc}/libNCSCnet${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSCnet${linus_sufix}.0.0.0 libNCSCnet${linus_sufix}.0${mac_sufix}

cp -f ${lib_loc}/libNCSEcw${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSEcw${linus_sufix}.0.0.0 libNCSEcw${linus_sufix}.0${mac_sufix}

cp -f ${lib_loc}/libNCSEcwC${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSEcwC${linus_sufix}.0.0.0 libNCSEcwC${linus_sufix}.0${mac_sufix}

cp -f ${lib_loc}/libNCSUtil${linus_sufix}.0.0.0${mac_sufix} .
ln -s -f libNCSUtil${linus_sufix}.0.0.0 libNCSUtil${linus_sufix}.0${mac_sufix}

# ---------- netcdf libs ---------------------------
#cp -f ${lib_loc}/libnetcdf${linus_sufix}.4.0.0${mac_sufix} .
#ln -s -f libnetcdf${linus_sufix}.4.0.0 libnetcdf${linus_sufix}.4${mac_sufix}

# ---------- Proj4 libs ----------------------------
cp -f ${lib_loc}/libproj${linus_sufix}.0.5.4${mac_sufix} .
ln -s -f libproj${linus_sufix}.0.5.4 libproj${linus_sufix}.0${mac_sufix}

# ---------- Shapefile libs ------------------------
cp -f ${lib_loc}/libshp${linus_sufix}.1.2.10${mac_sufix} .
ln -s -f libshp${linus_sufix}.1.2.10 libshp${linus_sufix}.1${mac_sufix}

# ---------- GMT libs ------------------------------
cp -f ${gmt_loc}/libgmt${linus_sufix}.4${mac_sufix} .
ln -s -f libgmt${linus_sufix}.4 libgmt${linus_sufix}${mac_sufix}

cp -f ${gmt_loc}/libpsl${linus_sufix}.4${mac_sufix} .
ln -s -f libpsl${linus_sufix}.4 libpsl${linus_sufix}${mac_sufix}
	