#! /bin/sh
#
# This script is the second of a triology:
# copy_patch_simple_OSX.sh, copy_patch_with_hide_OSX.sh, patch_mexs_OSX.sh   
# Run them in that order from a sub-directory of the mirone root dir
#
# However, don't forget to first: update GMT, rebuild gmtmex, build the mexs in the mex dir
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
cp /usr/local/opt/hdf5/lib/*.100.dylib .
chmod +w *.dylib

# Commented the _cpp ones because they are aparently not needed and are now 11 instead of 10

install_name_tool -id libhdf5_hide.dylib libhdf5.100.dylib
#install_name_tool -id libhdf5_cpp_hide.dylib libhdf5_cpp.10.dylib
install_name_tool -id libhdf5_hl_hide.dylib libhdf5_hl.100.dylib
#install_name_tool -id libhdf5_hl_cpp_hide.dylib libhdf5_hl_cpp.10.dylib

# The order is important
#install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5.10.dylib libhdf5_hide.dylib libhdf5_cpp.10.dylib
#mv libhdf5_cpp.10.dylib libhdf5_cpp_hide.dylib
#install_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5_cpp_hide.dylib

var_shit=`otool -L libhdf5.100.dylib | tail -n +3 | awk '{print $1}' | grep libsz`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libsz.dylib  libhdf5.100.dylib
fi
mv libhdf5.100.dylib libhdf5_hide.dylib

var_shit=`otool -L libhdf5_hl.100.dylib | tail -n +3 | awk '{print $1}' | grep libhdf5`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hide.dylib libhdf5_hl.100.dylib
fi

var_shit=`otool -L libhdf5_hl.100.dylib | tail -n +3 | awk '{print $1}' | grep libsz`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libsz.dylib libhdf5_hl.100.dylib
fi
mv libhdf5_hl.100.dylib libhdf5_hl_hide.dylib

#install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5.10.dylib libhdf5_hide.dylib         libhdf5_hl_cpp.10.dylib
#install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5_cpp.10.dylib libhdf5_cpp_hide.dylib libhdf5_hl_cpp.10.dylib
#install_name_tool -change /usr/local/Cellar/hdf5/1.8.15/lib/libhdf5_hl.10.dylib libhdf5_hl_hide.dylib   libhdf5_hl_cpp.10.dylib
#mv libhdf5_hl_cpp.10.dylib libhdf5_hl_cpp_hide.dylib
#ßinstall_name_tool -change /usr/local/lib/libsz.2.dylib libsz.dylib libhdf5_hl_cpp_hide.dylib


# NetCDF
cp /usr/local/lib/libnetcdf.11.dylib .
chmod +w libnetcdf.11.dylib
install_name_tool -id libnetcdf_hide.dylib libnetcdf.11.dylib

var_shit=`otool -L libnetcdf.11.dylib | tail -n +3 | awk '{print $1}' | grep libhdf5.1`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hide.dylib  libnetcdf.11.dylib
fi

var_shit=`otool -L libnetcdf.11.dylib | tail -n +3 | awk '{print $1}' | grep libhdf5_hl`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hl_hide.dylib  libnetcdf.11.dylib
fi

var_shit=`otool -L libnetcdf.11.dylib | tail -n +3 | awk '{print $1}' | grep libsz`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libsz.dylib libnetcdf.11.dylib
fi

mv libnetcdf.11.dylib libnetcdf_hide.dylib

# FITSIO
cp /usr/local/lib/libcfitsio.5.3.39.dylib .
chmod +w libcfitsio.5.3.39.dylib
install_name_tool -id libcfitsio_hide.dylib libcfitsio.5.3.39.dylib
mv libcfitsio.5.3.39.dylib libcfitsio_hide.dylib

# GEOTIFF
cp /usr/local/lib/libgeotiff.2.dylib .
chmod +w libgeotiff.2.dylib
install_name_tool -id libgeotiff_hide.dylib libgeotiff.2.dylib

var_shit=`otool -L libgeotiff.2.dylib | tail -n +3 | awk '{print $1}' | grep libproj`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libproj_hide.dylib  libgeotiff.2.dylib
fi

var_shit=`otool -L libgeotiff.2.dylib | tail -n +3 | awk '{print $1}' | grep libtiff`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libtiff.dylib  libgeotiff.2.dylib
fi

var_shit=`otool -L libgeotiff.2.dylib | tail -n +3 | awk '{print $1}' | grep libjpeg`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libjpeg.dylib  libgeotiff.2.dylib
fi

mv libgeotiff.2.dylib libgeotiff_hide.dylib

#TIFF
cp -f /usr/local/lib/libtiff.5.dylib .
chmod +w libtiff.5.dylib
install_name_tool -id libtiff.dylib libtiff.5.dylib

var_shit=`otool -L libtiff.5.dylib | tail -n +3 | awk '{print $1}' | grep libjpeg`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libjpeg.8.dylib libtiff.5.dylib 
fi

mv libtiff.5.dylib libtiff.dylib

# PROJ4
cp -f /usr/local/lib/libproj.12.dylib .
chmod +w libproj.12.dylib
install_name_tool -id libproj_hide.dylib libproj.12.dylib
mv libproj.12.dylib libproj_hide.dylib

# GEOS
cp -f /usr/local/lib/libgeos_c.1.dylib .
chmod +w libgeos_c.1.dylib
install_name_tool -id libgeos_c_hide.dylib libgeos_c.1.dylib

var_shit=`otool -L libgeos_c.1.dylib | tail -n +3 | awk '{print $1}' | grep libgeos`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgeos-3.dylib libgeos_c.1.dylib
fi
mv libgeos_c.1.dylib libgeos_c_hide.dylib

#XERCES
cp /usr/local/lib/libxerces-c-3.1.dylib .
chmod +w libxerces-c-3.1.dylib
install_name_tool -id libxerces-c_hide.dylib libxerces-c-3.1.dylib

var_shit=`otool -L libxerces-c-3.1.dylib | tail -n +3 | awk '{print $1}' | grep libxerces-c`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libxerces-c-3.1.dylib libxerces-c.dylib 
fi
mv libxerces-c-3.1.dylib libxerces-c_hide.dylib


# GDAL
cp /usr/local/lib/libgdal.1.dylib .
chmod +w libgdal.1.dylib
install_name_tool -id libgdal.dylib libgdal.1.dylib

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libjson`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libjson-c.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libproj`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libproj_hide.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libfreexl`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libfreexl.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libgeos_c`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgeos_c_hide.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libwebp`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libwebp.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libepsilon`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit}  libepsilon.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libodbc`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit}  libodbc.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libodbcinst`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit}  libodbcinst.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libjasper`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libjasper.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libnetcdf`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libnetcdf_hide.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libhdf5`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hide.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libgif`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgif.dylib   libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libjpeg`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libjpeg.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libgeotiff`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgeotiff_hide.dylib   libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libtiff`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libtiff.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libpng16`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libpng16.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libcfitsio`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libcfitsio_hide.dylib   libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep liblzma`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} liblzma.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libdap`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libdap.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libdapserver`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libdapserver.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libdapclient`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libdapclient.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libspatialite`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libspatialite.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libsqlite3`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libsqlite3.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libpcre`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libpcre.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libxml2`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libxml2.dylib libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libxerces`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libxerces-c_hide.dylib   libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libodbc`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libodbc.dylib  libgdal.1.dylib
fi

var_shit=`otool -L libgdal.1.dylib | tail -n +3 | awk '{print $1}' | grep libodbcinst`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libodbcinst.dylib  libgdal.1.dylib
fi

mv libgdal.1.dylib libgdal.dylib

# GMT
cp ${gmt_loc}/libgmt.5.dylib .
chmod +w libgmt.5.dylib
install_name_tool -id libgmt.dylib libgmt.5.dylib

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libnetcdf`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libnetcdf_hide.dylib     libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libpcre`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit}  libpcre.dylib            libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libgdal`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgdal.dylib  libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libfftw3f.3`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit}  libfftw3f.dylib  libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libfftw3f_threads`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit}  libfftw3f_threads.dylib  libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libpostscriptlight`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libpostscriptlight.dylib libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libhdf5.1`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hide.dylib  libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libhdf5_hl`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hl_hide.dylib  libgmt.5.dylib
fi

var_shit=`otool -L libgmt.5.dylib | tail -n +3 | awk '{print $1}' | grep libsz`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libsz.dylib libgmt.5.dylib
fi

mv libgmt.5.dylib libgmt.dylib

cp -f ${gmt_loc}/libpostscriptlight.5.dylib .
chmod +w libpostscriptlight.5.dylib
install_name_tool -id libpostscriptlight.dylib libpostscriptlight.5.dylib

var_shit=`otool -L libpostscriptlight.5.dylib | tail -n +3 | awk '{print $1}' | grep libpostscriptlight`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libpostscriptlight.dylib libpostscriptlight.5.dylib
fi

mv libpostscriptlight.5.dylib libpostscriptlight.dylib

cp ${gmt_loc}/gmt/plugins/supplements.so .
chmod +w supplements.so

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libgmt`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgmt.dylib supplements.so
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libpostscriptlight`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libpostscriptlight.dylib supplements.so
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libnetcdf`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libnetcdf_hide.dylib   supplements.so 
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libgdal`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgdal.dylib   supplements.so
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libhdf5.1`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hide.dylib  supplements.so
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libhdf5_hl`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libhdf5_hl_hide.dylib  supplements.so
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libsz`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libsz.dylib supplements.so
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libpcre`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libpcre.dylib  supplements.so 
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libfftw3f.3`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit}  libfftw3f.dylib  supplements.so
fi

var_shit=`otool -L supplements.so | tail -n +2 | awk '{print $1}' | grep libfftw3f_threads`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libfftw3f_threads.dylib  supplements.so 
fi

# And the GMTMEX
cp ${gmt_loc}/../gmt-mex/src/gmtmex.mexmaci64 .
cp ${gmt_loc}/../gmt-mex/src/gmt.m .
chmod +w *.mexmaci64

var_shit=`otool -L supplements.dylib | tail -n +2 | awk '{print $1}' | grep libgmt`
if [[ ${var_shit:0:1} == *"/"* ]]; then
install_name_tool -change ${var_shit} libgmt.dylib gmtmex.mexmaci64
fi
