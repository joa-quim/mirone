VPATH = ../../mex

HOSTTYPE = $(shell octave-config -p CANONICAL_HOST_TYPE)

SOURCES = set_gmt.c
MEX = $(patsubst %.c,%.mex,$(SOURCES))
TARGET_MEX = $(patsubst %,$(HOSTTYPE)/%,$(MEX))

all: $(HOSTTYPE) $(TARGET_MEX)
	echo $(TARGET_MEX)

$(HOSTTYPE):
	mkdir -p $(HOSTTYPE)

$(HOSTTYPE)/%.mex: %.c
	mkoctfile --mex -o $@ $<

clean:
	rm -rf $(HOSTTYPE)
	rm *.o
