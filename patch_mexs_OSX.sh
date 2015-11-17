#! /bin/sh
#
# This script is the third of a triology:
# copy_patch_simple_OSX.sh, copy_patch_with_hide_OSX.sh, patch_mexs_OSX.sh   
# Run them in that order from a sub-directory of the mirone root dir
#
# What they do is to copy several dylibs from /usr/local/lib (built with Homebrew) and strip the
# hard-coded paths. Besides that, several of those shared libs are assigned a different name
# (by appending a '_hide') so that Matlab does not managed to f. in the middle with its own outdated
# versions. For example libhdf5.dylib will become libhdf5_hide.dylib
#
# Patch the MEXs in lib_mex

# $Id :

cd ../lib_mex

install_name_tool -change /usr/local/lib/libnetcdf.7.dylib libnetcdf_hide.dylib mexnc.mexmaci64
install_name_tool -change /usr/local/lib/libnetcdf.7.dylib libnetcdf_hide.dylib swan.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_core.2.4.dylib    libopencv_core.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_imgproc.2.4.dylib libopencv_imgproc.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_calib3d.2.4.dylib libopencv_calib3d.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_objdetect.2.4.dylib libopencv_objdetect.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_video.2.4.dylib  libopencv_video.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libopencv_photo.2.4.dylib  libopencv_photo.dylib cvlib_mex.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib set_gmt.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib akimaspline.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib alloc_mex.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib country_select.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib gmtmbgrid_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib grdgradient_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib grdtrack_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib grdutils.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib igrf_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib mansinha_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib mex_illuminate.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib mirblock.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib mxgridtrimesh.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib ogrread.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib range_change.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib read_isf.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib scaleto8.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib susan.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib telha_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib test_gmt.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib trend1d_m.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib tsun2.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib wave_travel_time.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib distmin.mexmaci64
install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib PolygonClip.mexmaci64


for i in gdalread gdalwrite gdalwarp_mex ogrproj gdaltransform_mex mex_shape
do
	prg=$i.mexmaci64
	install_name_tool -change /usr/local/lib/libgdal.1.dylib libgdal.dylib $prg
done

cd ../tmp
