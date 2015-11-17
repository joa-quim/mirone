#! /bin/sh
#
# This script is the second of a triology:
# copy_patch_simple_OSX.sh, copy_patch_with_hide_OSX.sh, patch_mexs_OSX.sh   
# Run them in that order from a sub-directory of the mirone root dir
#
# What they do is to copy several dylibs from /usr/local/lib (built with Homebrew) and strip the
# hard-coded paths. Besides that, several of those shared libs are assigned a different name
# (by appending a '_hide') so that Matlab does not managed to f. in the middle with its own outdated
# versions. For example libhdf5.dylib will become libhdf5_hide.dylib
#
# Copy libs built with homebrew and patch them to new names, including its dependencies
#
# WARNING: The GMT lib is here copyied from a local path, NOT from /usr/local/lib

# $Id :

gmt_loc=/Users/j/programs/gmt5/lib

# HDF5
cp /usr/local/Cellar/hdf5/1.8.15/lib/*.10.dylib .
chmod +w *.dylib

install_name_tool -id libhdf5_hide.dylib libhdf5.10.dylib
install_name_tool -id libhdf5_cpp_hide.dylib libhdf5_cpp.10.dylib
install_name_tool -id libhdf5_hl_hide.dylib libhdf5_hl.10.dylib
install_name_tool -id libhdf5_hl_cpp_hide.dylib libhdf5_hl_cpp.10.dylib

# The order is important
install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5.10.dylib libhdf5_hide.dylib libhdf5_cpp.10.dylib
mv libhdf5_cpp.10.dylib libhdf5_cpp_hide.dylib
install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5_cpp_hide.dylib

mv libhdf5.10.dylib libhdf5_hide.dylib
install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5_hide.dylib

install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5.10.dylib libhdf5_hide.dylib libhdf5_hl.10.dylib
mv libhdf5_hl.10.dylib libhdf5_hl_hide.dylib
install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5_hl_hide.dylib

install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5.10.dylib libhdf5_hide.dylib         libhdf5_hl_cpp.10.dylib
install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5_cpp.10.dylib libhdf5_cpp_hide.dylib libhdf5_hl_cpp.10.dylib
install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5_hl.10.dylib libhdf5_hl_hide.dylib   libhdf5_hl_cpp.10.dylib
mv libhdf5_hl_cpp.10.dylib libhdf5_hl_cpp_hide.dylib
install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5_hl_cpp_hide.dylib


# NetCDF
cp /usr/local/lib/libnetcdf.7.dylib .
chmod +w libnetcdf.7.dylib
install_name_tool -id libnetcdf_hide.dylib libnetcdf.7.dylib
install_name_tool -change /usr/local/lib/libhdf5.10.dylib libhdf5_hide.dylib         libnetcdf.7.dylib
install_name_tool -change /usr/local/lib/libhdf5_hl.10.dylib libhdf5_hl_hide.dylib   libnetcdf.7.dylib
mv libnetcdf.7.dylib libnetcdf_hide.dylib
install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libnetcdf_hide.dylib

# FITSIO
cp /usr/local/lib/libcfitsio.2.3.37.dylib .
chmod +w libcfitsio.2.3.37.dylib
install_name_tool -id libcfitsio_hide.dylib libcfitsio.2.3.37.dylib
mv libcfitsio.2.3.37.dylib libcfitsio_hide.dylib

# GEOTIFF
cp /usr/local/lib/libgeotiff.2.dylib .
chmod +w libgeotiff.2.dylib
install_name_tool -id libgeotiff_hide.dylib libgeotiff.2.dylib
install_name_tool -change /usr/local/lib/libproj.9.dylib libproj_hide.dylib  libgeotiff.2.dylib
install_name_tool -change /usr/local/lib/libtiff.5.dylib libtiff.dylib       libgeotiff.2.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib libjpeg.dylib       libgeotiff.2.dylib
mv libgeotiff.2.dylib libgeotiff_hide.dylib

# PROJ4
cp /usr/local/lib/libproj.9.dylib .
chmod +w libproj.9.dylib
install_name_tool -id libproj_hide.dylib libproj.9.dylib
mv libproj.9.dylib libproj_hide.dylib

# LIBLWGEOM
# Tried here to see if it solved the mystery of the compatibility versions, but doesn't. It still says
# libproj_hide.dylib (compatibility version 10.0.0, current version 10.0.0)
# But it should be 11 like the "libproj_hide.dylib" itself. So let this be done in copy_patch_simple.sh
#cp /usr/local/opt/liblwgeom/lib/liblwgeom-2.1.5.dylib .
#chmod +w liblwgeom-2.1.5.dylib
#install_name_tool -id liblwgeom-2.dylib liblwgeom-2.1.5.dylib
#install_name_tool -change /usr/local/lib/libgeos_c.1.dylib libgeos_c.dylib   liblwgeom-2.1.5.dylib
#install_name_tool -change /usr/local/lib/libproj.9.dylib libproj_hide.dylib  liblwgeom-2.1.5.dylib
#install_name_tool -change /usr/local/lib/libjson-c.2.dylib libjson-c.dylib   liblwgeom-2.1.5.dylib
#mv liblwgeom-2.1.5.dylib liblwgeom-2.dylib


# GDAL
cp /usr/local/lib/libgdal.1.dylib .
chmod +w libgdal.1.dylib
install_name_tool -id libgdal.dylib libgdal.1.dylib
install_name_tool -change /usr/local/lib/libjson-c.2.dylib     libjson-c.dylib         libgdal.1.dylib
install_name_tool -change /usr/local/lib/libproj.9.dylib       libproj_hide.dylib      libgdal.1.dylib
install_name_tool -change /usr/local/lib/libfreexl.1.dylib     libfreexl.dylib         libgdal.1.dylib
install_name_tool -change /usr/local/lib/libgeos_c.1.dylib     libgeos_c.dylib         libgdal.1.dylib
install_name_tool -change /usr/local/lib/libwebp.5.dylib       libwebp.dylib           libgdal.1.dylib
install_name_tool -change /usr/local/lib/libepsilon.1.dylib    libepsilon.dylib        libgdal.1.dylib
install_name_tool -change /usr/local/lib/libodbc.2.dylib       libodbc.dylib           libgdal.1.dylib
install_name_tool -change /usr/local/lib/libodbcinst.2.dylib   libodbcinst.dylib       libgdal.1.dylib
install_name_tool -change /usr/local/lib/libjasper.1.dylib     libjasper.dylib         libgdal.1.dylib
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib     libnetcdf_hide.dylib    libgdal.1.dylib
install_name_tool -change /usr/local/lib/libhdf5.10.dylib      libhdf5_hide.dylib      libgdal.1.dylib
install_name_tool -change /usr/local/lib/libgif.4.dylib        libgif.dylib            libgdal.1.dylib
install_name_tool -change /usr/local/lib/libjpeg.8.dylib       libjpeg.dylib           libgdal.1.dylib
install_name_tool -change /usr/local/lib/libgeotiff.2.dylib    libgeotiff_hide.dylib   libgdal.1.dylib
install_name_tool -change /usr/local/lib/libtiff.5.dylib       libtiff.dylib           libgdal.1.dylib
install_name_tool -change /usr/local/lib/libpng16.16.dylib     libpng16.dylib          libgdal.1.dylib
install_name_tool -change /usr/local/lib/libcfitsio.2.dylib    libcfitsio_hide.dylib   libgdal.1.dylib
install_name_tool -change /usr/local/lib/liblzma.5.dylib       liblzma.dylib           libgdal.1.dylib
install_name_tool -change /usr/local/lib/libdap.11.dylib       libdap.dylib            libgdal.1.dylib
install_name_tool -change /usr/local/lib/libdapserver.7.dylib  libdapserver.dylib      libgdal.1.dylib
install_name_tool -change /usr/local/lib/libdapclient.3.dylib  libdapclient.dylib      libgdal.1.dylib
install_name_tool -change /usr/local/lib/libspatialite.7.dylib libspatialite.dylib     libgdal.1.dylib
install_name_tool -change /usr/local/opt/sqlite/lib/libsqlite3.0.dylib libsqlite3.dylib libgdal.1.dylib
install_name_tool -change /usr/local/lib/libpcre.1.dylib       libpcre.dylib           libgdal.1.dylib
install_name_tool -change /usr/local/opt/libxml2/lib/libxml2.2.dylib libxml2.dylib     libgdal.1.dylib
mv libgdal.1.dylib libgdal.dylib

# GMT
cp ${gmt_loc}/libgmt.5.dylib .
chmod +w libgmt.5.dylib
install_name_tool -id libgmt.dylib libgmt.5.dylib
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib      libnetcdf_hide.dylib     libgmt.5.dylib
install_name_tool -change /usr/local/lib/libpcre.1.dylib        libpcre.dylib            libgmt.5.dylib
install_name_tool -change /usr/local/lib/libgdal.1.dylib        libgdal.dylib            libgmt.5.dylib
install_name_tool -change /usr/local/lib/libfftw3f.3.dylib      libfftw3f.dylib          libgmt.5.dylib
install_name_tool -change /usr/local/lib/libfftw3f_threads.3.dylib  libfftw3f_threads.dylib  libgmt.5.dylib
install_name_tool -change ${gmt_loc}/libpostscriptlight.5.dylib libpostscriptlight.dylib libgmt.5.dylib
mv libgmt.5.dylib libgmt.dylib

cp -f ${gmt_loc}/libpostscriptlight.5.dylib .
chmod +w libpostscriptlight.5.dylib
ln -s -f libpostscriptlight.5.dylib libpostscriptlight.dylib

cp ${gmt_loc}/gmt/plugins/supplements.so supplements.dylib
chmod +w supplements.dylib
install_name_tool -change ${gmt_loc}/libgmt.5.dylib libgmt.dylib supplements.dylib
install_name_tool -change ${gmt_loc}/libpostscriptlight.5.dylib    libpostscriptlight.dylib supplements.dylib
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib         libnetcdf_hide.dylib     supplements.dylib 
install_name_tool -change /usr/local/lib/libgdal.1.dylib           libgdal.dylib            supplements.dylib
install_name_tool -change /usr/local/lib/libpcre.1.dylib           libpcre.dylib            supplements.dylib 
install_name_tool -change /usr/local/lib/libfftw3f.3.dylib         libfftw3f.dylib          supplements.dylib
install_name_tool -change /usr/local/lib/libfftw3f_threads.3.dylib libfftw3f_threads.dylib  supplements.dylib 

# And the GMTMEX
cp ${gmt_loc}/../gmt-mex/src/gmtmex.mexmaci64 .
cp ${gmt_loc}/../gmt-mex/src/gmt.m .
chmod +w *.mexmaci64
install_name_tool -change ${gmt_loc}/libgmt.5.dylib libgmt.dylib gmtmex.mexmaci64
install_name_tool -change ${gmt_loc}/libgmt.5.dylib libgmt.dylib gmt.mexmaci64
