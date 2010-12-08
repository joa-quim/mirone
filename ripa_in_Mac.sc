#! /bin/sh
#
# Script to rip the crazy dynamic linking dependencies created by default on MacOS 
# This is in the hope that it could do any good to the mystery of why MacMirone has
# so many problems in running on other machines that the one it was compiled.

#install_name_tool -id libgmt.4.dylib libgmt.4.dylib
#install_name_tool -change /usr/local/lib/libnetcdf.4.dylib libnetcdf.4.dylib libgmt.4.dylib
#install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.1.dylib libgmt.4.dylib

install_name_tool -id libpsl.4.dylib libpsl.4.dylib

install_name_tool -id libgdal.1.dylib libgdal.1.dylib 
install_name_tool -change /usr/local/lib/libNCSEcw.0.dylib libNCSEcw.0.dylib libgdal.1.dylib
install_name_tool -change /usr/local/lib/libNCSCnet.0.dylib libNCSCnet.0.dylib libgdal.1.dylib
install_name_tool -change /usr/local/lib/libNCSUtil.0.dylib libNCSUtil.0.dylib libgdal.1.dylib
install_name_tool -change /usr/local/lib/libnetcdf.4.dylib libnetcdf.4.dylib libgdal.1.dylib

install_name_tool -id libNCSCnet.0.dylib libNCSCnet.0.0.0.dylib
install_name_tool -id libNCSEcw.0.dylib libNCSEcw.0.0.0.dylib
install_name_tool -id /usr/local/lib/libNCSEcwC.0.dylib libNCSEcwC.0.0.0.dylib
install_name_tool -id /usr/local/lib/libNCSUtil.0.dylib libNCSUtil.0.0.0.dylib

install_name_tool -id libcv.4.dylib libcv.4.dylib
install_name_tool -change /usr/local/lib/libcxcore.4.dylib libcxcore.4.dylib libcv.4.dylib
install_name_tool -change /usr/lib/libz.1.dylib libz.1.dylib libcv.4.dylib

install_name_tool -id libcxcore.4.dylib libcxcore.4.dylib
install_name_tool -change /usr/lib/libz.1.dylib libz.1.dylib libcxcore.4.dylib

install_name_tool -id libnetcdf.4.dylib libnetcdf.4.dylib

install_name_tool -id libsz.2.dylib libsz.2.0.0.dylib

install_name_tool -id libz.1.dylib libz.1.dylib

cd lib_mex

install_name_tool -change /usr/local/lib/libnetcdf.4.dylib libnetcdf.4.dylib mexnc.mexmaci64
install_name_tool -change /usr/local/lib/libnetcdf.4.dylib libnetcdf.4.dylib swan.mexmaci64
install_name_tool -change /usr/local/lib/libcv.4.dylib libcv.4.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libcxcore.4.dylib libcxcore.4.dylib cvlib_mex.mexmaci64

for i in grdinfo_m grdproject_m grdread_m grdsample_m grdtrend_m grdwrite_m mapproject_m shoredump surface_m nearneighbor_m grdfilter_m cpt2cmap grdlandmask_m grdppa_m gmtlist_m shake_mex
do
	prg=$i.mexmaci64
	install_name_tool -change /Users/j/programs/GMT/lib/libgmt.4.dylib libgmt.4.dylib $prg
	install_name_tool -change /Users/j/programs/GMT/lib/libpsl.4.dylib libpsl.4.dylib $prg
	install_name_tool -change /usr/local/lib/libnetcdf.4.dylib libnetcdf.4.dylib $prg
done

for i in gdalread gdalwrite gdalwarp_mex ogrproj gdaltransform_mex
do
	prg=$i.mexmaci64
	install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.1.dylib $prg
done

cd ..
