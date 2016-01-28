VPATH = ../../mex

HOSTTYPE = $(shell octave-config -p CANONICAL_HOST_TYPE)

SOURCES = set_gmt.c igrf_m.c scaleto8.c tsun2.c wave_travel_time.c mansinha_m.c telha_m.c range_change.c country_select.c \
	mex_illuminate.c grdutils.c read_isf.c alloc_mex.c susan.c mxgridtrimesh.c trend1d_m.c gmtmbgrid_m.c \
	grdgradient_m.c grdtrack_m.c spa_mex.c mirblock.c write_mex.c xyzokb_m.c ind2rgb8.c applylutc.c cq.c bwlabel1.c \
	bwlabel2.c imhistc.c intlutc.c inv__lwm.c grayto8.c grayto16.c ordf.c parityscan.c distmin.c CalcMD5.c
MEX = $(patsubst %.c,%.mex,$(SOURCES))
TARGET_MEX = $(patsubst %,$(HOSTTYPE)/%,$(MEX))

all: $(HOSTTYPE) $(TARGET_MEX)
	echo $(TARGET_MEX)

$(HOSTTYPE):
	mkdir -p $(HOSTTYPE)

$(HOSTTYPE)/%.mex: %.c
	mkoctfile --mex -o -V $@ $<

clean:
	rm -rf $(HOSTTYPE)
	rm *.o


# mkoctfile --mex ../../mex/gdalread.c -I/c/programs/GDALtrunk/gdal/compileds/VC10_32/include -L/c/programs/GD ALtrunk/gdal/compileds/VC10_32/lib -lgdal_i -o i686-pc-mingw32/gdalread.mex
