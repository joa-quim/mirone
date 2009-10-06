#! /bin/sh
#
# Script to fish dynamic libraries and copy them to the Mirone root dir
# and create symbolic links.
# Note that this is crude and very system dependent
# It should only be needed if we want to create a destributable package

lib_loc=/usr/local/lib

# ---------- OpenCV libs -------------------------
cp -f ${lib_loc}/libcv.so.1.0.0 .
ln -s -f libcv.so.1.0.0 libcv.so.1

cp -f ${lib_loc}/libcxcore.so.1.0.0 .
ln -s -f libcxcore.so.1.0.0 libcxcore.so.1

# ---------- ?? libs ------------------------------
cp -f /usr/lib/libexpat.so.1.5.2 .
ln -s -f libexpat.so.1.0.0 libexpat.so.1

# ---------- GDAL libs ---------------------------
cp -f ${lib_loc}/libgdal.so.1.13.0 .
ln -s -f libgdal.so.1.13.0 libgdal.so.1

cp -f ${lib_loc}/libNCSCnet.so.0.0.0 .
ln -s -f libNCSCnet.so.0.0.0 libNCSCnet.so.0

cp -f ${lib_loc}/libNCSEcw.so.0.0.0 .
ln -s -f libNCSEcw.so.0.0.0 libNCSEcw.so.0

cp -f ${lib_loc}/libNCSEcwC.so.0.0.0 .
ln -s -f libNCSEcwC.so.0.0.0 libNCSEcwC.so.0

cp -f ${lib_loc}/libNCSUtil.so.0.0.0 .
ln -s -f libNCSUtil.so.0.0.0 libNCSUtil.so.0

# ---------- netcdf libs ---------------------------
cp -f ${lib_loc}/libnetcdf.so.4.0.0 .
ln -s -f libnetcdf.so.4.0.0 libnetcdf.so.4

# ---------- Proj4 libs ----------------------------
cp -f ${lib_loc}/libproj.so.0.5.4 .
ln -s -f libproj.so.0.5.4 libproj.so.0

# ---------- Shapefile libs ------------------------
cp -f ${lib_loc}/libshp.so.1.0.1 .
ln -s -f libshp.so.1.0.1 libshp.so.1

# ---------- GMT libs ------------------------------
# Hmm, My GMT is not under /usr/local
cp -f /home/bagside/programs/GMT/lib/libgmt.so.4 .
ln -s -f libgmt.so.4 libgmt.so

cp -f /home/bagside/programs/GMT/lib/libpsl.so.4 .
ln -s -f libpsl.so.4 libpsl.so

