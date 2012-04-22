/*--------------------------------------------------------------------
 *	$Id:$
 *
 *	Copyright (c) 2004-2012 by J. Luis
 *
 * 	This program is part of Mirone and is free software; you can redistribute
 * 	it and/or modify it under the terms of the GNU Lesser General Public
 * 	License as published by the Free Software Foundation; either
 * 	version 2.1 of the License, or any later version.
 * 
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * 	Lesser General Public License for more details.
 *
 *	Contact info: w3.ualg.pt/~jluis/mirone
 *--------------------------------------------------------------------*/

/* Program:	gdalread.c
 * Purpose:	matlab callable routine to read files supported by gdal
 * 		and dumping all band data of that dataset.
 *
 * Revision 22 12/11/2011 Added option -s to force output as float. Force error when fail to open file
 * Revision 21 23/02/2010 nedCDF bug is perhaps fixed. Limit previous solution to to pre 1.7 version
 * Revision 20 17/02/2009 Added -L option to deal with MODIS L2 left-right flipping
 * Revision 19 15/01/2009 Added the "Name" field to the attributes struct
 * Revision 18 18/03/2008 Another attempt to patch the broken netCDF driver
 * Revision 17 26/01/2008 Moved the y_min > y_max test to before the correct_bounds case.  
 *                        opt_r (size by pixels) was extracting a grid one row & column too big 
 * Revision 16 29/10/2007 The netCDF driver is still highly broken. Apply patch to make it a bit less bad
 * Revision 15 21/07/2007 Added the -r option. Like -R but uses pixel units
 * Revision 14 24/04/2007 Was crashing when called with a coards NETCDF file and -M
 * 			  Patch to deal with NETCDF driver adfGeoTransform bug
 * Revision 13 09/03/2007 Added MinMax field to the Band struct field
 * Revision 12 22/01/2007 Left-right flip ENVISAT GCPs because they have X positive to left (can we kill the guy?)
 * Revision 11 21/11/2006 Removed globals
 *                        Fixed a bug with pointer name confusion(dptr & dptr2)
 * Revision 10 04/11/2006 The 'ProjectionRef' metadata attempts to report in PrettyWkt
 * Revision  9 15/09/2006 Fixed (or at least reduced) the memory leak
 *                        Added the att = gdalread('','-M'); form. only att.Drivers is non-empty 
 * Revision  8 12/09/2006 Added the "Subdatasets" & "ImageStructure" to the attributes struct
 * Revision  7 11/09/2006 Requested bands could be larger than actual number of bands
 *                        Ignore option -C (scale) if data is already unsigned int
 * Revision  6 05/09/2006 Added the GCPs and Metadata to the attributes structure
 * Revision  5 08/05/2006 Added Color Table to the attributes structure
 * Revision  4 01/02/2006 Corrected uggly bug in index computations
 *			  Added -B option
 *			  Acelareted the indeces computation
 * Revision  3 07/06/2005 Added the -R (GMT style) option.
 * Revision  2 28/02/2005 Improved memory usage.
 *			  Passing options via -M<echanism>
 *			  Added a "insitu" option
 *			  Outputs an attribute structure with Dataset metadata
 *			  Taken (and modified) from mexgdal of John Evans (johnevans@acm.org)
 * Revision  1 25/6/2003  Joaquim Luis
 *
 */

#define CNULL	((char *)NULL)
#define i32_swap(x, y) {int tmp; tmp = x, x = y, y = tmp;}
#define ui32_swap(x, y) {unsigned int tmp; tmp = x, x = y, y = tmp;}
#define f_swap(x, y) {float tmp; tmp = x, x = y, y = tmp;}
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

/*	This copysign thing is a good mess. I don't know if it exists Windows or not.
	cpl_config.h defines it to _copysign, but then the compiler doesn't know it.
	So lets call it something else */
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#include "mex.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"

int record_geotransform ( char *gdal_filename, GDALDatasetH hDataset, double *adfGeoTransform );
mxArray * populate_metadata_struct (char * ,int , int, int, int, int, double, double, double, double, double, double);
int ReportCorner(GDALDatasetH hDataset, double x, double y, double *xy_c, double *xy_geo);
void grd_FLIPUD_I32(int data[], int nx, int ny);
void grd_FLIPUD_UI32(unsigned int data[], int nx, int ny);
void grd_FLIPUD_F32(float data[], int nx, int ny);
void troca_insituI32(int *a, int n, int m);
void troca_insituUI32(unsigned int *a, int n, int m);
void troca_insituF32(float *a, int n, int m);
void to_col_majorI32(int *in, int *out, int n_col, int n_row, int flipud, int insitu, int *nVector, int *mVector, int offset);
void to_col_majorUI32(unsigned int *in, unsigned int *out, int n_col, int n_row, int flipud, int insitu, int *nVector, int *mVector, int offset);
void to_col_majorF32(float *in, float *out, int n_col, int n_row, int flipud, int insitu, int *nVector, int *mVector, int offset);
void ComputeRasterMinMax(char *tmp, GDALRasterBandH hBand, double adfMinMax[2], int nXSize, int nYSize, double, double);
int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (double w, double e, double s, double n);
int decode_columns (char *txt, int *whichBands, int n_col);
int GMT_strtok (const char *string, const char *sep, int *start, char *token);
double ddmmss_to_degree (char *text);


/* ------------------- Matlab Gateway routine ---------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	char	**argv, *gdal_filename;
	const char	*format;
	int	ndims, bGotNodata, flipud = FALSE, metadata_only = FALSE, got_R = FALSE;
	int	nPixelSize, nBands, nXYSize, i, m, n, nn, nReqBands = 0, got_r = FALSE;
	int	dims[]={0,0,0,0,0,0,0}, argc = 0, n_arg_no_char = 0;
	int	error = FALSE, gdal_dump = FALSE, insitu = FALSE, scale_range = FALSE;
	int	pixel_reg = FALSE, correct_bounds = FALSE, fliplr = FALSE, forceSingle = FALSE;
	int	anSrcWin[4], xOrigin = 0, yOrigin = 0, i_x_nXYSize;
	int	nBufXSize, nBufYSize, jump = 0, *whichBands = NULL, *nVector, *mVector;
	int	n_commas, n_dash, nX, nY;
	int dataType;
	int	nXSize = 0, nYSize = 0;
	int	bGotMin, bGotMax;	/* To know if driver transmited Min/Max */
	char	*tmp, *outByte, *p;
	static int runed_once = FALSE;	/* It will be set to true if reaches end of main */
	float	*tmpF32, *outF32;
	double	dfNoData, range, aux, adfMinMax[2];
	double	dfULX = 0.0, dfULY = 0.0, dfLRX = 0.0, dfLRY = 0.0;
	double	z_min = 1e50, z_max = -1e50;
	GDALDatasetH	hDataset;
	GDALRasterBandH	hBand;
	GDALDriverH	hDriver;
	GInt16	*tmpI16, *outI16;
	GUInt16	*tmpUI16, *outUI16;
	GInt32	*tmpI32, *outI32;
	GUInt32	*tmpUI32, *outUI32;

	anSrcWin[0] = anSrcWin[1] = anSrcWin[2] = anSrcWin[3] = 0;

	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}
	argc++;			/* to account for the program's name to be inserted in argv[0] */

	/* get the length of the input string */
	argv = (char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "gdalread";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'B':	/* We have a selected bands request */
					for (n = 2, n_commas = 0; argv[i][n]; n++) if (argv[i][n] == ',') n_commas = n;
					for (n = 2, n_dash = 0; argv[i][n]; n++) if (argv[i][n] == '-') n_dash = n;
					nn = MAX(n_commas, n_dash);
					if (nn)
						nn = atoi(&argv[i][nn+1]);
					else
						nn = atoi(&argv[i][2]);
					whichBands = mxCalloc(nn, sizeof(int));
					nReqBands = decode_columns (&argv[i][2], whichBands, nn);
					break;
				case 'C':
					correct_bounds = TRUE;
					break;
				case 'c':
					p = &argv[i][2];
					if ((p = strchr(p, '/')) == NULL) {
						mexPrintf("Error using option -c. It must be on the form -c<key>/value\n");
					}
					else {
						*p = '\0';	*p++;
						CPLSetConfigOption( &argv[i][2], p ); 
						mexPrintf("key = %s\tval = %s\n",&argv[i][2], p);
					}
					break;
				case 'F':
					pixel_reg = TRUE;
					break;
				case 'I':
					insitu = TRUE;
					break;
				case 'L':
					fliplr = TRUE;
					break;
				case 'M':
					metadata_only = TRUE;
					break;
				case 'R':
					error += decode_R (argv[i], &dfULX, &dfLRX, &dfLRY, &dfULY);
					got_R = TRUE;
					break;
				case 'r':	/* Region is given in pixels */
					error += decode_R (argv[i], &dfULX, &dfLRX, &dfLRY, &dfULY);
					got_r = TRUE;
					break;
				case 'P':
					jump = atoi(&argv[i][2]);
					break;
				case 'S':
					scale_range = TRUE;
					break;
				case 's':
					forceSingle = TRUE;
					break;
				case 'U':
					flipud = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}
	
	if (error || nrhs < 1 || nlhs > 2) {
		mexPrintf ("usage(s): z = gdalread('filename',['-C'],['-F'],['-I'],['-Rw/e/s/n']);\n");
		mexPrintf ("                       ['-S'],['-U'], ['-c<key>/<value>'], ['-s']);\n");
		mexPrintf (" 	 [z,attrib] = gdalread('filename', ...);\n");
		mexPrintf (" 	 attrib     = gdalread('','-M');\n\n");
		mexPrintf ("\t   attrib is a structure with metadata.\n");
		mexPrintf ("\t-C correct the grid bounds reported by GDAL (it thinks that all grids are pixel registered)\n");
		mexPrintf ("\t-c Sets the named configuration keyword to the given value. Some common configuration\n");
		mexPrintf ("\t   keywords are GDAL_CACHEMAX (memory used internally for caching in megabytes) and\n");
		mexPrintf ("\t   GDAL_DATA (path of the GDAL 'data' directory)\n");
		mexPrintf ("\t-F force pixel registration in attrib.GMT_hdr\n");
		mexPrintf ("\t-I force importing via the 'insitu' transposition (about 10 times slower)\n");
		mexPrintf ("\t-L flip the grid LeftRight (For the time being, ignored if -U).\n");
		mexPrintf ("\t-M ouputs only the metadata structure\n");
		mexPrintf ("\t-S scale ouptut into the [0-255] range\n");
		mexPrintf ("\t-s Force the output 'z' array to be of float type (singles)\n");
		mexPrintf ("\t-R read only the sub-region enclosed by <west/east/south/north>\n");
		mexPrintf ("\t-U flip the grid UpDown (needed for all DEM grids in Mirone)\n");
		return;
	}

	/* Load the file name into a char string */

	if ( (gdal_filename = (char *)mxArrayToString(prhs[0])) == NULL) {
		mexPrintf ("%s\n", gdal_filename);
		mexErrMsgTxt ("gdalread: failure to decode gdal_filename string \n");
	}

	if (metadata_only && nlhs == 2)
		mexErrMsgTxt ("gdalread: -M option implies only one output\n");

	if (nlhs == 2)		/* Output also the metadata struct */
		gdal_dump = TRUE;

	/* Open gdal - If it has not been opened before (in a previous encarnation) */

	if (!runed_once)		/* Do next call only at first time this MEX is loaded */
		GDALAllRegister();

	if (metadata_only) {
		plhs[0] = populate_metadata_struct (gdal_filename, correct_bounds, pixel_reg, got_R, 
						nXSize, nYSize, dfULX, dfULY, dfLRX, dfLRY, z_min, z_max);
		runed_once = TRUE;	/* Signals that next call won't need to call GDALAllRegister() again */
		return;
	}

	hDataset = GDALOpen(gdal_filename, GA_ReadOnly);

	if (hDataset == NULL) {
		mexPrintf ("GDALOpen failed %s\n", CPLGetLastErrorMsg());
		plhs[0] = mxCreateNumericMatrix (0,0,mxDOUBLE_CLASS, mxREAL);
		mexErrMsgTxt ("\n");
	}

	/* Some formats (tipically DEMs) have their origin at Bottom Left corner.
	   For those we have to flip the data matrix to be in accord with matlab (Top Left) */

	hDriver = GDALGetDatasetDriver(hDataset);
	format = GDALGetDriverShortName(hDriver);

	if (!strcmp(format,"SDTS"))
		flipud = TRUE;
	else if (!strcmp(format,"USGSDEM"))
		flipud = TRUE;

	if (!strcmp(format,"ESAT"))	/* ENVISAT data are flipped left-right */
		fliplr = TRUE;

	if (GDAL_VERSION_NUM < 1700 && !strcmp(format,"netCDF"))
		flipud = FALSE;

	if (got_R || got_r) {
		/* -------------------------------------------------------------------- */
		/*      Compute the source window from the projected source window      */
		/*      if the projected coordinates were provided.  Note that the      */
		/*      projected coordinates are in ulx, uly, lrx, lry format,         */
		/*      while the anSrcWin is xoff, yoff, xsize, ysize with the         */
		/*      xoff,yoff being the ulx, uly in pixel/line.                     */
		/* -------------------------------------------------------------------- */
		double	adfGeoTransform[6];

		GDALGetGeoTransform( hDataset, adfGeoTransform );

		if( adfGeoTransform[2] != 0.0 || adfGeoTransform[4] != 0.0 ) {
			mexPrintf("The -projwin option was used, but the geotransform is\n"
					"rotated. This configuration is not supported.\n");
			GDALClose( hDataset );
			GDALDestroyDriverManager();
			return;
		}

		if (got_R) {	/* Region in map coordinates */
			anSrcWin[0] = (int) ((dfULX - adfGeoTransform[0]) / adfGeoTransform[1] + 0.001);
			anSrcWin[1] = (int) ((dfULY - adfGeoTransform[3]) / adfGeoTransform[5] + 0.001);
			anSrcWin[2] = (int) ((dfLRX - dfULX) / adfGeoTransform[1] + 0.5);
			anSrcWin[3] = (int) ((dfLRY - dfULY) / adfGeoTransform[5] + 0.5);
			if (GDAL_VERSION_NUM < 1700 && !strcmp(format,"netCDF")) {
				/* PATCH against the never ending GDAL bug of reading netCDF files */
				anSrcWin[1] = GDALGetRasterYSize(hDataset) - (anSrcWin[1] + anSrcWin[3]) - 1;
			}
		}
		else {		/* Region in pixel/line */
			anSrcWin[0] = (int) (dfULX);
			anSrcWin[1] = (int) (dfLRY);
			anSrcWin[2] = (int) (dfLRX - dfULX);
			anSrcWin[3] = (int) (dfULY - dfLRY);
		}

		if( anSrcWin[0] < 0 || anSrcWin[1] < 0
			|| anSrcWin[0] + anSrcWin[2] > GDALGetRasterXSize(hDataset)
			|| anSrcWin[1] + anSrcWin[3] > GDALGetRasterYSize(hDataset) ) {
			mexPrintf("Computed -srcwin falls outside raster size of %dx%d.\n",
			GDALGetRasterXSize(hDataset),
			GDALGetRasterYSize(hDataset) );
			return;
		}
		xOrigin = anSrcWin[0];
		yOrigin = anSrcWin[1];
		nXSize = anSrcWin[2];
		nYSize = anSrcWin[3];
		if (correct_bounds) {	/* Patch for the bug reading GMT grids */
			nXSize++;
			nYSize++;
		}
	}
	else {			/* Use entire dataset */
		xOrigin = yOrigin = 0;
		nXSize = GDALGetRasterXSize(hDataset);
		nYSize = GDALGetRasterYSize(hDataset);
	}

	/* The following assumes that all bands have the same PixelSize. Otherwise ... */
	hBand = GDALGetRasterBand(hDataset,1);
	nPixelSize = GDALGetDataTypeSize(GDALGetRasterDataType(hBand)) / 8;	/* /8 because return value is in BITS */

	if (jump) {
		nBufXSize = GDALGetRasterXSize(hDataset) / jump;
		nBufYSize = GDALGetRasterYSize(hDataset) / jump;
	}
	else {
		nBufXSize = nXSize;	nBufYSize = nYSize;
	}

	nBands = GDALGetRasterCount(hDataset);
	nXYSize = nBufXSize * nBufYSize;

	if (nReqBands) nBands = MIN(nBands,nReqBands);	/* If a band selection was made */

	if ( scale_range && (GDALGetRasterDataType(hBand) == GDT_Byte) )	/* Just sanitizing */
		scale_range = FALSE;

	/* Create a matrix for the return array */
	if (nBands == 1)
		ndims = 2;
	else if (nBands > 1) {
		ndims = 3;	dims[2] = nBands;
	}

	dims[0] = nBufYSize;	dims[1] = nBufXSize;

	if (scale_range) {
		plhs[0] = mxCreateNumericArray (ndims,dims,mxUINT8_CLASS, mxREAL);
		outByte = mxGetData(plhs[0]);
		if (!got_R)		/* Otherwise, Min/Max will be computed bellow */
			GDALComputeRasterMinMax(hBand, FALSE, adfMinMax);
	}
	else if (forceSingle) {
		plhs[0] = mxCreateNumericArray (ndims,dims,mxSINGLE_CLASS, mxREAL);
		outF32 = mxGetData(plhs[0]);
		nPixelSize = 4;
		insitu = FALSE;		/* Just to be sure */
	}
	else {
		switch( GDALGetRasterDataType(hBand) ) {
			case GDT_Byte:
				plhs[0] = mxCreateNumericArray (ndims,dims,mxUINT8_CLASS, mxREAL);
				outByte = mxGetData(plhs[0]);
				insitu = FALSE;		/* Just to be sure */
				break;
			case GDT_Int16:
				plhs[0] = mxCreateNumericArray (ndims,dims,mxINT16_CLASS, mxREAL);
				outI16 = mxGetData(plhs[0]);
				insitu = FALSE;		/* Just to be sure */
				break;
			case GDT_UInt16:
				plhs[0] = mxCreateNumericArray (ndims,dims,mxUINT16_CLASS, mxREAL);
				outUI16 = mxGetData(plhs[0]);
				insitu = FALSE;		/* Just to be sure */
				break;
			case GDT_Int32:
				plhs[0] = mxCreateNumericArray (ndims,dims,mxINT32_CLASS, mxREAL);
				outI32 = mxGetData(plhs[0]);
				break;
			case GDT_UInt32:
				plhs[0] = mxCreateNumericArray (ndims,dims,mxUINT32_CLASS, mxREAL);
				outUI32 = mxGetData(plhs[0]);
				break;
			case GDT_Float32:
				plhs[0] = mxCreateNumericArray (ndims,dims,mxSINGLE_CLASS, mxREAL);
				outF32 = mxGetData(plhs[0]);
				break;
		}
	}

	if (!insitu) {		/* EXPLICAR PORQUE */
		tmp = mxCalloc(nBufYSize * nBufXSize, nPixelSize);
		if (tmp == NULL)
			mexErrMsgTxt ("gdalread: failure to allocate enough memory\n");
	}
	else
		tmp = mxGetData(plhs[0]);	/* Apropriate type was already selected above */

	/* ------ compute two vectors indices that will be used inside loops below --------- */
	/* In the "Preview" mode those guys bellow are different and what we need is the BufSize */
	if (jump) {
		nX = nBufXSize;
		nY = nBufYSize;
	}
	else {
		nX = nXSize;
		nY = nYSize;
	}
	nVector = mxCalloc(nX, sizeof(int));
	mVector = mxCalloc(nY, sizeof(int));
	for (m = 0; m < nY; m++) mVector[m] = m*nX;
	if (flipud)
		for (n = 0; n < nX; n++) nVector[n] = n*nY-1 + nY;
	else {		/* For now I do the fliplr test only here because of the ENVISAT case */
		if (fliplr)
			for (n = 0, m = nX-1; n < nX; m--, n++) nVector[n] = m*nY;
		else
			for (n = 0; n < nX; n++) nVector[n] = n*nY;
	}
	/* --------------------------------------------------------------------------------- */

	for (i = 0; i < nBands; i++) {
		if (!nReqBands)		/* No band selection, read them sequentialy */
			hBand = GDALGetRasterBand( hDataset, i+1 );
		else			/* Band selection. Read only the requested ones */
			hBand = GDALGetRasterBand( hDataset, whichBands[i] );

		dataType = GDALGetRasterDataType(hBand);
		if (forceSingle)
			dataType = GDT_Float32;

		GDALRasterIO(hBand, GF_Read, xOrigin, yOrigin, nXSize, nYSize,
			tmp, nBufXSize, nBufYSize, dataType, 0, 0 );

        	dfNoData = GDALGetRasterNoDataValue(hBand, &bGotNodata);
		/* If we didn't computed it yet, its time to do it now */
		if (got_R) ComputeRasterMinMax(tmp, hBand, adfMinMax, nXSize, nYSize, z_min, z_max);
		if (nBands > 1 && scale_range) {	/* got_R && scale_range && nBands > 1 Should never be true */
			adfMinMax[0] = GDALGetRasterMinimum( hBand, &bGotMin );
			adfMinMax[1] = GDALGetRasterMaximum( hBand, &bGotMax );
			if(!(bGotMin && bGotMax))
				GDALComputeRasterMinMax(hBand, FALSE, adfMinMax);
		}

		/* In the "Preview" mode those guys bellow are different and what we need is the BufSize */
		if (jump) {
			nXSize = nBufXSize;
			nYSize = nBufYSize;
			i_x_nXYSize = i*nXSize*nYSize;		/* We don't need to recompute this everytime */
		}
		else
			i_x_nXYSize = i*nXYSize;

		switch (dataType) {
			case GDT_Byte:
				if (scale_range) {	 /* Scale data into the [0 255] range */
					range = adfMinMax[1] - adfMinMax[0];
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (bGotNodata) {
							if (tmp[mVector[m]+n] != (GByte)dfNoData) {
								aux = tmp[mVector[m]+n];
								aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
							}
							else	/* Put the NoData value to a cte */
								aux = 0;
						}
						else {
							aux = tmp[mVector[m]+n];
							aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
						}
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outByte[nn] = (char)aux;
					}
				}
				else {	/* No scaling */
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outByte[nn] = tmp[mVector[m]+n];
					}
				}
				break;
			case GDT_Int16:
				tmpI16 = (GInt16 *) tmp;
				if (scale_range) {	 /* Scale data into the [0 255] range */
					range = adfMinMax[1] - adfMinMax[0];
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (bGotNodata) {
							if (tmpI16[mVector[m]+n] != (GInt16)dfNoData) {
								aux = tmpI16[mVector[m]+n];
								aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
							}
							else	/* Put the NoData value to a cte */
								aux = 0;
						}
						else {
							aux = tmpI16[mVector[m]+n];
							aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
						}
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outByte[nn] = (char)aux;
					}
				}
				else {	/* No scaling */
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outI16[nn] = tmpI16[mVector[m]+n];
					}
				}
				break;
			case GDT_UInt16:
				tmpUI16 = (GUInt16 *) tmp;
				if (scale_range) {	 /* Scale data into the [0 255] range */
					range = adfMinMax[1] - adfMinMax[0];
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (bGotNodata) {
							if (tmpUI16[mVector[m]+n] != (GUInt16)dfNoData) {
								aux = tmpUI16[mVector[m]+n];
								aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
							}
							else	/* Put the NoData value to a cte */
								aux = 0;
						}
						else {
							aux = tmpUI16[mVector[m]+n];
							aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
						}
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outByte[nn] = (char)aux;
					}
				}
				else {	/* No scaling */
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outUI16[nn] = tmpUI16[mVector[m]+n];
					}
				}
				break;
			case GDT_Int32:
				tmpI32 = (GInt32 *) tmp;
				if (scale_range) {	 /* Scale data into the [0 255] range */
					range = adfMinMax[1] - adfMinMax[0];
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (bGotNodata) {
							if (tmpI32[mVector[m]+n] != (GInt32)dfNoData) {
								aux = (float)tmpI32[mVector[m]+n];
								aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
							}
							else	/* Put the NoData value to a cte */
								aux = 0;
						}
						else {
							aux = (float)tmpI32[mVector[m]+n];
							aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
						}
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outByte[nn] = (char)aux;
					}
				}
				else 	/* No scaling */
					to_col_majorI32(tmpI32, outI32, nXSize, nYSize, flipud, insitu, nVector, mVector, i_x_nXYSize);
				break;
			case GDT_UInt32:
				tmpUI32 = (GUInt32 *) tmp;
				if (scale_range) {	 /* Scale data into the [0 255] range */
					range = adfMinMax[1] - adfMinMax[0];
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (bGotNodata) {
							if (tmpUI32[mVector[m]+n] != (GUInt32)dfNoData) {
								aux = (float)tmpUI32[mVector[m]+n];
								aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
							}
							else	/* Put the NoData value to a cte */
								aux = 0;
						}
						else {
							aux = (float)tmpUI32[mVector[m]+n];
							aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
						}
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outByte[nn] = (char)aux;
					}
				}
				else 	/* No scaling */
					to_col_majorUI32(tmpUI32, outUI32, nXSize, nYSize, flipud, insitu, nVector, mVector, i_x_nXYSize);
				break;
			case GDT_Float32:
				tmpF32 = (float *) tmp;
				if (scale_range) {	 /* Scale data into the [0 255] range */
					range = adfMinMax[1] - adfMinMax[0];
					for (m = 0; m < nYSize; m++) for (n = 0; n < nXSize; n++) {
						if (bGotNodata) {
							if (tmpF32[mVector[m]+n] > (float)dfNoData) {
								aux = tmpF32[mVector[m]+n];
								aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
							}
							else	/* Put the NoData value to a cte */
								aux = 0;
						}
						else {
							aux = tmpF32[mVector[m]+n];
							aux = ((aux - adfMinMax[0]) / range ) * 254 + 1;
						}
						if (flipud)	nn = nVector[n]-m + i_x_nXYSize;
						else		nn = nVector[n]+m + i_x_nXYSize;
						outByte[nn] = (char)aux;
					}
				}
				else 	/* No scaling */
					to_col_majorF32(tmpF32, outF32, nXSize, nYSize, flipud, insitu, nVector, mVector, i_x_nXYSize);
				break;
			default:
				CPLAssert( FALSE );
		}
	}

	mxFree(argv);
	mxFree(nVector);
	mxFree(mVector);
	if (!insitu) 		/* EXPLICAR PORQUE */
		mxFree(tmp);

	GDALClose(hDataset);

	if (gdal_dump)
		plhs[1] = populate_metadata_struct (gdal_filename, correct_bounds, pixel_reg, got_R, 
						nXSize, nYSize, dfULX, dfULY, dfLRX, dfLRY, z_min, z_max);

	mxFree(gdal_filename);
	runed_once = TRUE;	/* Signals that next call won't need to call GDALAllRegister() again */
}

/* =============================================================================================== */
void to_col_majorI32(int *in, int *out, int n_col, int n_row, int flipud, int insitu, int *nVector, int *mVector, int offset) {
	int nn, m, n;
	if (!insitu) {
		for (m = 0; m < n_row; m++) for (n = 0; n < n_col; n++) {
			if (flipud)	nn = nVector[n]-m + offset;
			else		nn = nVector[n]+m + offset;
			out[nn] = in[mVector[m]+n];
		}
	}
	else {
		grd_FLIPUD_I32(out, n_col, n_row);
		troca_insituI32(out, n_col, n_row);
	}
}

void to_col_majorUI32(unsigned int *in, unsigned int *out, int n_col, int n_row, int flipud, int insitu, int *nVector, int *mVector, int offset) {
	int nn, m, n;
	if (!insitu) {
		for (m = 0; m < n_row; m++) for (n = 0; n < n_col; n++) {
			if (flipud)	nn = nVector[n]-m + offset;
			else		nn = nVector[n]+m + offset;
			out[nn] = in[mVector[m]+n];
		}
	}
	else {
		grd_FLIPUD_UI32(out, n_col, n_row);
		troca_insituUI32(out, n_col, n_row);
	}
}

void to_col_majorF32(float *in, float *out, int n_col, int n_row, int flipud, int insitu, int *nVector, int *mVector, int offset) {
	int nn, m, n;
	if (!insitu) {
		for (m = 0; m < n_row; m++) for (n = 0; n < n_col; n++) {
			if (flipud)	nn = nVector[n]-m + offset;
			else		nn = nVector[n]+m + offset;
			out[nn] = in[mVector[m]+n];
		}
	}
	else {
		grd_FLIPUD_F32(out, n_col, n_row);
		troca_insituF32(out, n_col, n_row);
	}
}

void troca_insituI32(int *a, int n, int m) {
	int i,row,column,current, size = m*n;
	for(i = 1, size -= 2; i < size; i++) {
		current = i;
		do      {
			column = current / m;
			if (column == 0) {current *= n;	continue;}
			row = current % m;
			current = n*row + column;
		} while(current < i);

		if (current > i)
			i32_swap(a[i],a[current]);
        }
}

void troca_insituUI32(unsigned int *a, int n, int m) {
	int i,row,column,current, size = m*n;
	for(i = 1, size -= 2; i < size; i++) {
		current = i;
		do      {
			column = current / m;
			if (column == 0) {current *= n;	continue;}
			row = current % m;
			current = n*row + column;
		} while(current < i);

		if (current > i)
			ui32_swap(a[i],a[current]);
        }
}

void troca_insituF32(float *a, int n, int m) {
	int i,row,column,current, size = m*n;
	for(i = 1, size -= 2; i < size; i++) {
		current = i;
		do      {
			column = current / m;
			if (column == 0) {current *= n;	continue;}
			row = current % m;
			current = n*row + column;
		} while(current < i);

		if (current > i)
			f_swap(a[i],a[current]);
        }
}

void grd_FLIPUD_I32(int data[], int nx, int ny) {
	int i, j, k, ny1, ny_half;
	/* Reverse order of all columns */
	ny_half = ny / 2;
	ny1 = ny - 1;
	for (i = 0; i < nx; i++) {
		for (j = 0, k = ny1; j < ny_half; j++, k--) {	/* Do this to all rows */
			i32_swap (data[j*nx+i], data[k*nx+i]);
		}
	}
}

void grd_FLIPUD_UI32(unsigned int data[], int nx, int ny) {
	int i, j, k, ny1, ny_half;
	/* Reverse order of all columns */
	ny_half = ny / 2;
	ny1 = ny - 1;
	for (i = 0; i < nx; i++) {
		for (j = 0, k = ny1; j < ny_half; j++, k--) {	/* Do this to all rows */
			ui32_swap (data[j*nx+i], data[k*nx+i]);
		}
	}
}

void grd_FLIPUD_F32(float data[], int nx, int ny) {
	int i, j, k, ny1, ny_half;
	/* Reverse order of all columns */
	ny_half = ny / 2;
	ny1 = ny - 1;
	for (i = 0; i < nx; i++) {
		for (j = 0, k = ny1; j < ny_half; j++, k--) {	/* Do this to all rows */
			f_swap (data[j*nx+i], data[k*nx+i]);
		}
	}
}

/*
 * POPULATE_METADATA_STRUCT
 *
 * This routine just queries the GDAL raster file for all the metadata
 * that can be squeezed out of it.
 *
 * The resulting matlab structure is by necessity nested.  Each raster
 * file can have several bands, e.g. PNG files usually have 3, a red, a
 * blue, and a green channel.  Each band can have several overviews (tiffs
 * come to mind here).
 *
 * Fields:
 *    ProjectionRef:  a string describing the projection.  Not parsed.
 *    GeoTransform:
 *        a 6-tuple.  Entries are as follows.
 *            [0] --> top left x
 *            [1] --> w-e pixel resolution
 *            [2] --> rotation, 0 if image is "north up"
 *            [3] --> top left y
 *            [4] --> rotation, 0 if image is "north up"
 *            [5] --> n-s pixel resolution
 *
 *    DriverShortName:  describes the driver used to query *this* raster file
 *    DriverLongName:  describes the driver used to query *this* raster file
 *    RasterXSize, RasterYSize:
 *        These are the primary dimensions of the raster.  See "Overview", though.
 *    RasterCount:
 *        Number of raster bands present in the file.
 *    Driver:
 *        This itself is a structure array.  Each element describes a driver
 *        that the locally compiled GDAL library has available.  So you recompile
 *        GDAL with new format support, this structure will change.
 *
 *        Fields:
 *            DriverShortName, DriverLongName:
 *                Same as fields in top level structure with same name.
 *
 *    Band:
 *        Also a structure array.  One element for each raster band present in
 *        the GDAL file.  See "RasterCount".
 *
 *        Fields:
 *            XSize, YSize:
 *                Dimensions of the current raster band.
 *            Overview:
 *                A structure array, one element for each overview present. If
 *                empty, then there are no overviews.
 *            NoDataValue:
 *                When passed back to MATLAB, one can set pixels with this value to NaN.
 *            MinMax:
 *                [min max] vector (one per band) with the band's MinMax. If not known [NaN NaN]
 *            ScaleOffset:
 *                [scale_factor offset] vector (one per band) with the band's Scale_factor Offset. If not known [1 0]
 *            ColorMap:
 *                A Mx3 double array with the colormap, or empty if it does not exists
 *
 *    ColorInterp:
 *
 *    Corners:
 *        Also a structure array with the Fields:
 *            LL, UL, UR, LR. Each of this contains the (x,y) coords (or pixel/line coords
 *            if image is not referenced) of the LL - LowerLeft - point, and so on.
 *
 *    GEOGCorners:
 *        A 4x4 cell array with the corner in geographical coordinates (if that is possible,
 *        or otherwise a single empty cell).
 *	  First two elements in a row contain a string version of Lon,Lat, the other two a numeric
 *	  version. The corner's order is the same as the 'Corners' field. That is LL, UL, UR, LR.
 *
 *    GMT_hdr:
 *        contains a row vector with (xmin, xmax, ymin, ymax, zmin, zmax, pixel_reg, xinc, yinc)
 *        for this data set. Pixel_reg is 1 for pixel registration and 0 for grid node
 *        registration, but it is actually set the -F option [default is 0]
 *
 *    GCPprojection:  a string describing the GCPs projection system. Not parsed.
 *
 *    GCPvalues:
 *        A Nx5 double array. First 2 columns are image (x,y) pixels. Remaining 3 contain
 *        the X,Y,Z coords in the projection system. Z is normaly zero.
 *
 *    Metadata:
 *        A cell array with the metadata associated to the dataset. In most cases it will be empty
 *
 *    Name:
 *        A string with the name of was is being open, which can be the main file or a subdataset
 *
 * */
mxArray *populate_metadata_struct (char *gdal_filename , int correct_bounds, int pixel_reg, int got_R, 
				int nXSize, int nYSize, double dfULX, double dfULY, double dfLRX, 
				double dfLRY, double z_min, double z_max) {
	static int driverCount = 0;	/* Number of available drivers for the version of GDAL we are using. */

	/* These are used to define the metadata structure about available GDAL drivers. */
	char *driver_fieldnames[100];
	const char	*format;
	int num_driver_fields;

	mxArray *driver_struct;
	mxArray *mxtmp;
	mxArray *mxProjectionRef;
	mxArray *mxGeoTransform;
	mxArray *mxGDALDriverShortName;
	mxArray *mxGDALDriverLongName;
	mxArray *mxGDALRasterCount;
	mxArray *mxGDALRasterXSize;
	mxArray *mxGDALRasterYSize;
	mxArray *mxCorners;
	mxArray *mxGMT_header;
	mxArray *mxCMap, *toMinMax, *toScaleOff;

	/* These will be matlab structures that hold the metadata.
	 * "metadata_struct" actually encompasses "band_struct",
	 * which encompasses "overview_struct" */
	mxArray *metadata_struct;
	mxArray *band_struct;
	mxArray *overview_struct;
	mxArray *colormap_struct;
	mxArray *corner_struct;

	int	i,j, overview, band_number;	/* Loop indices */
	int	n_colors;		/* Number of colors in the eventual Color Table */ 
	double	*dptr, *dptr2, *dptrMM,	/* short cut to the mxArray data */
		*dptrSO;

	GDALDriverH hDriver;		/* This is the driver chosen by the GDAL library to query the dataset. */
	GDALDatasetH hDataset;		/* pointer structure used to query the gdal file. */
	GDALRasterBandH hBand, overview_hBand;
	GDALColorTableH	hTable;
	GDALColorEntry	sEntry;
	const	GDAL_GCP *psGCP;	/* Pointer to a GCP structure */

	char	error_msg[500];		/* Construct error and warning messages using this buffer. */
	int	num_overview_fields;	/* Number of metadata items for each overview structure. */
	int	status;			/* success or failure */
	double 	adfGeoTransform[6];	/* bounds on the dataset */
	int	num_overviews;		/* Number of overviews in the current band. */
	int	num_struct_fields;	/* number of fields in the metadata structures. */
	int	num_band_fields;
	int	num_corner_fields;
	int	nCounter;		/* Generic counter */
	char	*fieldnames[20];	/* this array contains the names of the fields of the metadata structure. */
	char	*band_fieldnames[20];
	char	*corner_fieldnames[20];
	char	*overview_fieldnames[20];
	char	*colormap_fieldnames[1];
	char	**papszMetadata;
	int	xSize, ySize, raster_count; /* Dimensions of the dataset */
	int	gdal_type;		/* Datatype of the bands. */
	double	nan, tmpdble;		/* NaN & temporary value */
	double	xy_c[2], xy_geo[4][2];	/* Corner coordinates in the local coords system and geogs (if it exists) */
	int	dims[2], bSuccess;
	int	bGotMin, bGotMax;	/* To know if driver transmited Min/Max */
	double	adfMinMax[2];		/* Dataset Min Max */

	/* ------------------------------------------------------------------------- */
	/* Retrieve information on all drivers. */
	/* ------------------------------------------------------------------------- */
	if (driverCount == 0)		/* Do the next call only at the first time MEX is loaded */
		driverCount = GDALGetDriverCount();

	/* ------------------------------------------------------------------------- */
	/* Create the metadata structure. Just one element, with XXX fields. */
	/* ------------------------------------------------------------------------- */
	num_struct_fields = 19;
	fieldnames[0] = strdup ("ProjectionRef");
	fieldnames[1] = strdup ("GeoTransform");
	fieldnames[2] = strdup ("DriverShortName");
	fieldnames[3] = strdup ("DriverLongName");
	fieldnames[4] = strdup ("RasterXSize");
	fieldnames[5] = strdup ("RasterYSize");
	fieldnames[6] = strdup ("RasterCount");
	fieldnames[7] = strdup ("ColorInterp");
	fieldnames[8] = strdup ("Driver" );
	fieldnames[9] = strdup ("Band");
	fieldnames[10] = strdup ("Corners");
	fieldnames[11] = strdup("GEOGCorners");
	fieldnames[12] = strdup("GMT_hdr");
	fieldnames[13] = strdup("GCPprojection");
	fieldnames[14] = strdup("GCPvalues");
	fieldnames[15] = strdup("Metadata");
	fieldnames[16] = strdup("Subdatasets");
	fieldnames[17] = strdup("ImageStructure");
	fieldnames[18] = strdup("Name");
	metadata_struct = mxCreateStructMatrix ( 1, 1, num_struct_fields, (const char **)fieldnames );

	num_driver_fields = 2;
	driver_fieldnames[0] = strdup( "DriverLongName" );
	driver_fieldnames[1] = strdup( "DriverShortName" );
	driver_struct = mxCreateStructMatrix ( driverCount, 1, num_driver_fields, (const char **)driver_fieldnames );
	for ( j = 0; j < driverCount; ++j ) {
		hDriver = GDALGetDriver( j );
		mxtmp = mxCreateString ( GDALGetDriverLongName(hDriver) );
		mxSetField ( driver_struct, j, (const char *) "DriverLongName", mxtmp );

		mxtmp = mxCreateString ( GDALGetDriverShortName(hDriver) );
		mxSetField ( driver_struct, j, (const char *) "DriverShortName", mxtmp );

	}
	mxSetField( metadata_struct, 0, "Driver", driver_struct );

	if (strlen(gdal_filename) == 0)
		return(metadata_struct);

	/* ------------------------------------------------------------------------- */
	/* Open the file (if we can). */
	/* ------------------------------------------------------------------------- */
	hDataset = GDALOpen ( gdal_filename, GA_ReadOnly );
	if ( hDataset == NULL ) {
		sprintf ( error_msg, "Unable to open %s.\n", gdal_filename );
		mexErrMsgTxt ( error_msg );
	}

	/* ------------------------------------------------------------------------- */
	/* Record the ProjectionRef. */
	/* ------------------------------------------------------------------------- */
	if( GDALGetProjectionRef( hDataset ) != NULL ) {
		OGRSpatialReferenceH  hSRS;
		char		      *pszProjection;
		pszProjection = (char *) GDALGetProjectionRef( hDataset );

		hSRS = OSRNewSpatialReference(NULL);
		if( OSRImportFromWkt( hSRS, &pszProjection ) == CE_None ) {
			char	*pszPrettyWkt = NULL;
			OSRExportToPrettyWkt( hSRS, &pszPrettyWkt, FALSE );
			mxProjectionRef = mxCreateString ( pszPrettyWkt );
			CPLFree( pszPrettyWkt );
		}
		else
			mxProjectionRef = mxCreateString ( GDALGetProjectionRef( hDataset ) );

		OSRDestroySpatialReference( hSRS );
	}
	mxSetField ( metadata_struct, 0, "ProjectionRef", mxProjectionRef );

	/* ------------------------------------------------------------------------- */
	/* Record the geotransform. */
	/* ------------------------------------------------------------------------- */
	mxGeoTransform = mxCreateNumericMatrix ( 6, 1, mxDOUBLE_CLASS, mxREAL );
	dptr = mxGetPr ( mxGeoTransform );

	status = record_geotransform ( gdal_filename, hDataset, adfGeoTransform );
	if ( status == 0 ) {
		dptr[0] = adfGeoTransform[0];
		dptr[1] = adfGeoTransform[1];
		dptr[2] = adfGeoTransform[2];
		dptr[3] = adfGeoTransform[3];
		dptr[4] = adfGeoTransform[4];
		dptr[5] = adfGeoTransform[5];
		if (got_R) {			/* Need to change those */
			dptr[0] = dfULX;
			dptr[3] = dfULY;
		}
		mxSetField ( metadata_struct, 0, "GeoTransform", mxGeoTransform );
	} 
	else {
		sprintf ( error_msg, "No internal georeferencing exists for %s, and could not find a suitable world file either.\n", gdal_filename );	
		/*mexWarnMsgTxt ( error_msg );*/
	}

	/* ------------------------------------------------------------------------- */
	/* Get driver information */
	/* ------------------------------------------------------------------------- */
	hDriver = GDALGetDatasetDriver( hDataset );

	mxGDALDriverShortName = mxCreateString ( GDALGetDriverShortName( hDriver ) );
	mxSetField ( metadata_struct, 0, (const char *) "DriverShortName", mxGDALDriverShortName );

	mxGDALDriverLongName = mxCreateString ( GDALGetDriverLongName( hDriver ) );
	mxSetField ( metadata_struct, 0, (const char *) "DriverLongName", mxGDALDriverLongName );

	if (!got_R) {		/* That is, if we are using the entire dataset */
		xSize = GDALGetRasterXSize( hDataset );
		ySize = GDALGetRasterYSize( hDataset );
	}
	else if (got_R && nXSize == 0 && nYSize == 0) {		/* metadata_only */
		int	anSrcWin[4];
		anSrcWin[0] = anSrcWin[1] = anSrcWin[2] = anSrcWin[3] = 0;
		if( adfGeoTransform[2] != 0.0 || adfGeoTransform[4] != 0.0 ) {
			mexPrintf("The -projwin option was used, but the geotransform is\n"
					"rotated. This configuration is not supported.\n");
			GDALClose( hDataset );
			GDALDestroyDriverManager();
			mexErrMsgTxt ("Quiting with error\n");
		}

		anSrcWin[0] = (int) ((dfULX - adfGeoTransform[0]) / adfGeoTransform[1] + 0.001);
		anSrcWin[1] = (int) ((dfULY - adfGeoTransform[3]) / adfGeoTransform[5] + 0.001);
		anSrcWin[2] = (int) ((dfLRX - dfULX) / adfGeoTransform[1] + 0.5);
		anSrcWin[3] = (int) ((dfLRY - dfULY) / adfGeoTransform[5] + 0.5);

		if( anSrcWin[0] < 0 || anSrcWin[1] < 0
			|| anSrcWin[0] + anSrcWin[2] > GDALGetRasterXSize(hDataset)
			|| anSrcWin[1] + anSrcWin[3] > GDALGetRasterYSize(hDataset) ) {
			mexPrintf("Computed -srcwin falls outside raster size of %dx%d.\n",
			GDALGetRasterXSize(hDataset),
			GDALGetRasterYSize(hDataset) );
			mexErrMsgTxt ("Quiting with error\n");
		}
		nXSize = anSrcWin[2];	nYSize = anSrcWin[3];
		if (correct_bounds) {	/* Patch for the bug reading GMT grids */
			nXSize++;	nYSize++;
		}

		xSize = nXSize;
		ySize = nYSize;
	}
	else {
		xSize = nXSize;
		ySize = nYSize;
	}

	mxGDALRasterXSize = mxCreateDoubleScalar ( (double) xSize );
	mxSetField ( metadata_struct, 0, (const char *) "RasterXSize", mxGDALRasterXSize );

	mxGDALRasterYSize = mxCreateDoubleScalar ( (double) ySize );
	mxSetField ( metadata_struct, 0, (const char *) "RasterYSize", mxGDALRasterYSize );

	raster_count = GDALGetRasterCount( hDataset );
	mxGDALRasterCount = mxCreateDoubleScalar ( (double)raster_count );
	mxSetField ( metadata_struct, 0, (const char *) "RasterCount", mxGDALRasterCount );

	/* ------------------------------------------------------------------------- */
	/* Get the Color Interpretation Name */
	/* ------------------------------------------------------------------------- */
	hBand = GDALGetRasterBand( hDataset, 1 );
	if (raster_count > 0) {
		mxtmp = mxCreateString ( GDALGetColorInterpretationName(
			GDALGetRasterColorInterpretation(hBand)) );
		mxSetField ( metadata_struct, 0, (const char *) "ColorInterp", mxtmp );
	}
	else
		mxSetField ( metadata_struct, 0, (const char *) "ColorInterp", mxCreateString("nikles") );

	/* ------------------------------------------------------------------------- */
	/* Get the metadata for each band. */
	/* ------------------------------------------------------------------------- */
	num_band_fields = 8;
	band_fieldnames[0] = strdup ( "XSize" );
	band_fieldnames[1] = strdup ( "YSize" );
	band_fieldnames[2] = strdup ( "Overview" );
	band_fieldnames[3] = strdup ( "NoDataValue" );
	band_fieldnames[4] = strdup ( "MinMax" );
	band_fieldnames[5] = strdup ( "ScaleOffset" );
	band_fieldnames[6] = strdup ( "DataType" );
	band_fieldnames[7] = strdup ( "ColorMap" );
	band_struct = mxCreateStructMatrix ( raster_count, 1, num_band_fields, (const char **)band_fieldnames );

	num_overview_fields = 2;
	overview_fieldnames[0] = strdup ( "XSize" );
	overview_fieldnames[1] = strdup ( "YSize" );
	colormap_fieldnames[0] = strdup ( "CMap" );
	nan = mxGetNaN();

	/* ==================================================================== */
	/*      Loop over bands.                                                */
	/* ==================================================================== */
	for ( band_number = 0; band_number < raster_count; ++band_number ) {

		hBand = GDALGetRasterBand( hDataset, band_number+1 );

		if (!got_R) 		/* Not sure about what will realy happen in this case */
			mxtmp = mxCreateDoubleScalar ( (double) GDALGetRasterBandXSize( hBand ) );
		else
			mxtmp = mxCreateDoubleScalar ( (double)nXSize );
		mxSetField ( band_struct, band_number, "XSize", mxtmp );

		if (!got_R)
			mxtmp = mxCreateDoubleScalar ( (double) GDALGetRasterBandYSize( hBand ) );
		else
			mxtmp = mxCreateDoubleScalar ( (double)nYSize );
		mxSetField ( band_struct, band_number, "YSize", mxtmp );
	
		gdal_type = GDALGetRasterDataType ( hBand );
		mxtmp = mxCreateString ( GDALGetDataTypeName ( gdal_type ) );
		mxSetField ( band_struct, band_number, (const char *) "DataType", mxtmp );

		mxtmp = mxCreateDoubleScalar ( (double) (GDALGetRasterNoDataValue ( hBand, &status ) ) );
		mxSetField ( band_struct, band_number, "NoDataValue", mxtmp );
	
		/* Get band's Min/Max. If the band has no record of it we won't compute it */
		toMinMax = mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL);	/* To hold band's min/max */
		dptrMM = mxGetPr(toMinMax);
		adfMinMax[0] = GDALGetRasterMinimum( hBand, &bGotMin );
		adfMinMax[1] = GDALGetRasterMaximum( hBand, &bGotMax );
		if((bGotMin && bGotMax)) {
			dptrMM[0] = adfMinMax[0];
			dptrMM[1] = adfMinMax[1];
		}
		else {		/* No Min/Max, we won't computer it either */
			dptrMM[0] = dptrMM[1] = nan;
		}
		mxSetField ( band_struct, band_number, "MinMax", toMinMax );

		/* Get band's Scale/Offset. If the band does not have them use the neutral 1/0 values */
		toScaleOff = mxCreateNumericMatrix(1, 2, mxDOUBLE_CLASS, mxREAL);	/* To hold band's min/max */
		dptrSO = mxGetPr(toScaleOff);
		if (GDALGetRasterScale( hBand, &bSuccess ) != 1 || GDALGetRasterOffset( hBand, &bSuccess ) != 0) {
			dptrSO[0] = GDALGetRasterScale ( hBand, &bSuccess );
			dptrSO[1] = GDALGetRasterOffset( hBand, &bSuccess );
		}
		else {
			dptrSO[0] = 1.0;
			dptrSO[1] = 0.0;
		}
		mxSetField ( band_struct, band_number, "ScaleOffset", toScaleOff );

		num_overviews = GDALGetOverviewCount( hBand );	/* Can have multiple overviews per band. */
		if ( num_overviews > 0 ) {

			overview_struct = mxCreateStructMatrix ( num_overviews, 1, num_overview_fields,
						(const char **)overview_fieldnames );
	
			for ( overview = 0; overview < num_overviews; ++overview ) {
				overview_hBand = GDALGetOverview ( hBand, overview );
	
				xSize = GDALGetRasterBandXSize ( overview_hBand );
				mxtmp = mxCreateDoubleScalar ( xSize );
				mxSetField ( overview_struct, overview, "XSize", mxtmp );
	
				ySize = GDALGetRasterBandYSize ( overview_hBand );
				mxtmp = mxCreateDoubleScalar ( ySize );
				mxSetField ( overview_struct, overview, "YSize", mxtmp );
			}
			mxSetField ( band_struct, band_number, "Overview", overview_struct );
		}

		/* See if we have Color Tables(s) */
		colormap_struct = mxCreateStructMatrix(raster_count, 1, 1, (const char **)colormap_fieldnames);
		if( GDALGetRasterColorInterpretation(hBand) == GCI_PaletteIndex 
			&& (hTable = GDALGetRasterColorTable( hBand )) != NULL ) {

			n_colors = GDALGetColorEntryCount( hTable );
			dims[0] = n_colors;	 dims[1] = 4;
			mxCMap = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
			dptr = mxGetPr(mxCMap);

			for( i = 0; i < n_colors; i++ ) {
				GDALGetColorEntryAsRGB( hTable, i, &sEntry );
				dptr[i] = (double)sEntry.c1 / 255;
				dptr[i+n_colors] = (double)sEntry.c2 / 255;
				dptr[i+2*n_colors] = (double)sEntry.c3 / 255;
				dptr[i+3*n_colors] = (double)sEntry.c4 / 255;
			}
			mxSetField ( colormap_struct, band_number, "CMap", mxCMap );
			mxSetField ( band_struct, band_number, "ColorMap", colormap_struct );
		}
	}
	mxSetField ( metadata_struct, 0, "Band", band_struct );

	/* ------------------------------------------------------------------------- */
	/* Record the GMT header. This will be interleaved with "corners" because they share somes values */
	/* ------------------------------------------------------------------------- */
	mxGMT_header = mxCreateNumericMatrix(1, 9, mxDOUBLE_CLASS, mxREAL);
	dptr2 = mxGetPr(mxGMT_header);

	/* ------------------------------------------------------------------------- */
	/* Record corners. */
	/* ------------------------------------------------------------------------- */
	num_corner_fields = 4;
	corner_fieldnames[0] = strdup ("LL");
	corner_fieldnames[1] = strdup ("UL");
	corner_fieldnames[2] = strdup ("UR");
	corner_fieldnames[3] = strdup ("LR");

	corner_struct = mxCreateStructMatrix(1, 1, num_corner_fields, (const char **)corner_fieldnames);
	dims[0] = 1;	 dims[1] = 2;

	if (!got_R)					/* Lower Left */ 
		ReportCorner(hDataset, 0.0, GDALGetRasterYSize(hDataset), xy_c, xy_geo[0]);
	else
		{xy_c[0] = dfULX;	xy_c[1] = dfLRY;}
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	dptr2[0] = xy_c[0];		dptr2[2] = xy_c[1];	/* xmin, ymin */
	mxSetField(corner_struct, 0, "LL", mxCorners );

	if (!got_R)					/* Upper Left */
		ReportCorner(hDataset, 0.0, 0.0, xy_c, xy_geo[1]);				
	else
		{xy_c[0] = dfULX;	xy_c[1] = dfULY;}
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	mxSetField(corner_struct, 0, "UL", mxCorners );

	if (!got_R)					/* Upper Rigt */
		ReportCorner(hDataset, GDALGetRasterXSize(hDataset), 0.0, xy_c, xy_geo[2]);
	else
		{xy_c[0] = dfLRX;	xy_c[1] = dfULY;}
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	dptr2[1] = xy_c[0];		dptr2[3] = xy_c[1];	/* xmax, ymax */
	mxSetField(corner_struct, 0, "UR", mxCorners );

	if (!got_R)					/* LR */
		ReportCorner(hDataset, GDALGetRasterXSize(hDataset), GDALGetRasterYSize(hDataset), xy_c, xy_geo[3]);
	else
		{xy_c[0] = dfLRX;	xy_c[1] = dfLRY;}
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	mxSetField(corner_struct, 0, "LR", mxCorners );

	mxSetField (metadata_struct, 0, "Corners", corner_struct);

	/* --------------------------------------------------------------------------------------
	 * Record Geographical corners (if they exist) in a 4x4 cell array. First two
	 * elements in a row contain a string version of Lon,Lat, the other a numeric
	 * version. The corner's order is the same as the 'Corners' field. That is LL, UL, UR, LR.
	 * -------------------------------------------------------------------------------------- */
	if( !got_R ) {
		if( !mxIsNaN(xy_geo[0][0]) ) {
			mxtmp = mxCreateCellMatrix(4, 4);
			for( i = 0; i < 4; i++ ) {
				mxSetCell( mxtmp,i,mxDuplicateArray( mxCreateString(GDALDecToDMS( xy_geo[i][0], "Long", 2 )) ) );
				mxSetCell( mxtmp,i+4,mxDuplicateArray( mxCreateString(GDALDecToDMS( xy_geo[i][1], "Lat", 2 )) ) );
				mxSetCell( mxtmp,i+8,mxDuplicateArray( mxCreateDoubleScalar(xy_geo[i][0])) );
				mxSetCell( mxtmp,i+12,mxDuplicateArray( mxCreateDoubleScalar(xy_geo[i][1])) );
			}
		}
		else
			mxtmp = mxCreateCellMatrix(0, 1);
	}
	else
		mxtmp = mxCreateCellMatrix(0, 1);

	mxSetField (metadata_struct, 0, "GEOGCorners", mxtmp);

	/* ------------------------------------------------------------------------- */
	/* Fill in the rest of the GMT header values (If ...) */
	if (raster_count > 0) {
		if (z_min == 1e50) {		/* We don't know yet the dataset Min/Max */
			/*adfMinMax[0] = GDALGetRasterMinimum( hBand, &bGotMin );
			adfMinMax[1] = GDALGetRasterMaximum( hBand, &bGotMax );
			if(!(bGotMin && bGotMax))
				GDALComputeRasterMinMax( hBand, FALSE, adfMinMax );*/
			/* If file is a "VRT/Virtual Raster" do NOT try to compute min/max and trust on XML info */
			if (strcmp(GDALGetDriverShortName(hDriver), "VRT"))
				GDALComputeRasterMinMax( hBand, FALSE, adfMinMax );		/* NO VRT, scan file to compute min/max */
			else
				GDALComputeRasterMinMax( hBand, TRUE, adfMinMax );		/* VRT, believe in metadata info */

			dptr2[4] = adfMinMax[0];
			dptr2[5] = adfMinMax[1];
		}
		else {
			dptr2[4] = z_min;
			dptr2[5] = z_max;
		}	

		/* See if we want grid or pixel registration */
		if (correct_bounds == FALSE && pixel_reg == FALSE)	/* Neither -C or -F so use GDAL defaults */
			dptr2[6] = 1;
		else if (pixel_reg)
			dptr2[6] = 1;		/* pixel registration */
		else
			dptr2[6] = 0;		/* grid registration */

		dptr2[7] = adfGeoTransform[1];
		dptr2[8] = fabs(adfGeoTransform[5]);

		if (dptr2[2] > dptr2[3]) {	/* Sometimes GDAL does it: y_min > y_max. If so, revert it */
			tmpdble = dptr2[2];
			dptr2[2] = dptr2[3];
			dptr2[3] = tmpdble;
		}

		if (correct_bounds && !got_R) {
			dptr2[0] += dptr2[7] / 2;	dptr2[1] -= dptr2[7] / 2;
			dptr2[2] += dptr2[8] / 2;	dptr2[3] -= dptr2[8] / 2;
		}
		else if (got_R) {
			dptr2[0] = dfULX;	dptr2[1] = dfLRX;
			dptr2[2] = dfLRY;	dptr2[3] = dfULY;
		}
	}
	mxSetField (metadata_struct, 0, "GMT_hdr", mxGMT_header);

	/* ------------------------------------------------------------------------- */
	/* Record the GCPprojection. */
	/* ------------------------------------------------------------------------- */
	mxProjectionRef = mxCreateString ( GDALGetGCPProjection( hDataset ) );
	mxSetField ( metadata_struct, 0, "GCPprojection", mxProjectionRef );

	/* ------------------------------------------------------------------------- */
	/* Record the GCPvalues. */
	/* ------------------------------------------------------------------------- */
	nCounter = GDALGetGCPCount( hDataset );
	mxtmp = mxCreateNumericMatrix(nCounter, 5, mxDOUBLE_CLASS, mxREAL);	/* Empy matrix if no GCPs */
	if( nCounter > 0 ) {
		format = GDALGetDriverShortName(hDriver);
		dptr = mxGetPr(mxtmp);
		for( i = 0; i < nCounter; i++ ) {
			psGCP = GDALGetGCPs( hDataset ) + i;
			if (!strcmp(format,"ESAT"))	/* ENVISAT GCPs are left-right fliped */ 
				dptr[i] = xSize - psGCP->dfGCPPixel;
			else
				dptr[i] = psGCP->dfGCPPixel;
			dptr[i+nCounter] = psGCP->dfGCPLine;
			dptr[i+2*nCounter] = psGCP->dfGCPX;
			dptr[i+3*nCounter] = psGCP->dfGCPY;
			dptr[i+4*nCounter] = psGCP->dfGCPZ;
		}
	}
	mxSetField (metadata_struct, 0, "GCPvalues", mxtmp);

	/* ------------------------------------------------------------------------- */
	/* Record Metadata (if any).
	/* ------------------------------------------------------------------------- */
	papszMetadata = GDALGetMetadata( hDataset, NULL );
	nCounter = CSLCount(papszMetadata);
	mxtmp = mxCreateCellMatrix(nCounter, 1);
	if( nCounter > 0 ) {
		for( i = 0; i < nCounter; i++ )
			mxSetCell( mxtmp,i,mxDuplicateArray( mxCreateString(papszMetadata[i]) ) );
	}
	mxSetField (metadata_struct, 0, "Metadata", mxtmp);

	/* ------------------------------------------------------------------------- */
	/* Record Subdatasets (if any).
	/* ------------------------------------------------------------------------- */
	papszMetadata = GDALGetMetadata( hDataset, "SUBDATASETS" );
	nCounter = CSLCount(papszMetadata);
	mxtmp = mxCreateCellMatrix(nCounter, 1);
	if( nCounter > 0 ) {
		for( i = 0; i < nCounter; i++ )
			mxSetCell( mxtmp,i,mxDuplicateArray( mxCreateString(papszMetadata[i]) ) );
	}
	mxSetField (metadata_struct, 0, "Subdatasets", mxtmp);

	/* ------------------------------------------------------------------------- */
	/* Record ImageStructure (if any).
	/* ------------------------------------------------------------------------- */
	papszMetadata = GDALGetMetadata( hDataset, "IMAGE_STRUCTURE" );
	nCounter = CSLCount(papszMetadata);
	mxtmp = mxCreateCellMatrix(nCounter, 1);
	if( nCounter > 0 ) {
		for( i = 0; i < nCounter; i++ )
			mxSetCell( mxtmp,i,mxDuplicateArray( mxCreateString(papszMetadata[i]) ) );
	}
	mxSetField (metadata_struct, 0, "ImageStructure", mxtmp);

	/* ------------------------------------------------------------------------- */
	/* Record Name
	/* ------------------------------------------------------------------------- */
	mxSetField( metadata_struct, 0, "Name", mxCreateString( gdal_filename ) );

	GDALClose ( hDataset );

	return ( metadata_struct );
}

/************************************************************************/
/*                        ReportCorner()                                */
/************************************************************************/

int ReportCorner(GDALDatasetH hDataset, double x, double y, double *xy_c, double *xy_geo) {
	double	dfGeoX, dfGeoY;
	const char  *pszProjection;
	double	adfGeoTransform[6];
	OGRCoordinateTransformationH hTransform = NULL;
	OGRSpatialReferenceH hProj, hLatLong = NULL;
	char *str_lon, *str_lat; 

	xy_geo[0] = xy_geo[1] = mxGetNaN();		/* Default return values */
/* -------------------------------------------------------------------- */
/*      Transform the point into georeferenced coordinates.             */
/* -------------------------------------------------------------------- */
	if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None ) {
        	pszProjection = GDALGetProjectionRef(hDataset);
        	dfGeoX = adfGeoTransform[0] + adfGeoTransform[1] * x + adfGeoTransform[2] * y;
        	dfGeoY = adfGeoTransform[3] + adfGeoTransform[4] * x + adfGeoTransform[5] * y;
	}
	else {
		xy_c[0] = x;	xy_c[1] = y;
        	return FALSE;
	}

/* -------------------------------------------------------------------- */
/*      Report the georeferenced coordinates.                           */
/* -------------------------------------------------------------------- */
	xy_c[0] = dfGeoX;	xy_c[1] = dfGeoY;

/* -------------------------------------------------------------------- */
/*      Setup transformation to lat/long.                               */
/* -------------------------------------------------------------------- */
    if( pszProjection != NULL && strlen(pszProjection) > 0 ) {

        hProj = OSRNewSpatialReference( pszProjection );
        if( hProj != NULL )
            hLatLong = OSRCloneGeogCS( hProj );

        if( hLatLong != NULL ) {
            CPLPushErrorHandler( CPLQuietErrorHandler );
            hTransform = OCTNewCoordinateTransformation( hProj, hLatLong );
            CPLPopErrorHandler();
            OSRDestroySpatialReference( hLatLong );
        }

        if( hProj != NULL )
            OSRDestroySpatialReference( hProj );
    }
/*
/* -------------------------------------------------------------------- */
/*      Transform to latlong and report.                                */
/* -------------------------------------------------------------------- */
	if( hTransform != NULL && OCTTransform(hTransform,1,&dfGeoX,&dfGeoY,NULL) )
		xy_geo[0] = dfGeoX;	xy_geo[1] = dfGeoY;

	if( hTransform != NULL )
		OCTDestroyCoordinateTransformation( hTransform );

	return TRUE;
}

/* ---------------------------------------------------------------------------------
 * record_geotransform:
 *
 * If the gdal file is not internally georeferenced, try to get the world file.
 * Returns -1 in case no world file is found.
 */
int record_geotransform ( char *gdal_filename, GDALDatasetH hDataset, double *adfGeoTransform ) {
	int status = -1;
	char generic_buffer[5000];

	if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None )
		return (0);

	/* Try a world file.  First the generic extension. If the gdal_filename
	   is, say, "a.tif", then this will look for "a.wld". */
	if ( GDALReadWorldFile ( gdal_filename, "wld", adfGeoTransform ) )
		return (0);

	/* Try again, but try "a.tif.wld" instead. */
	sprintf ( generic_buffer, "%s.xxx", gdal_filename );
	status = GDALReadWorldFile ( generic_buffer, "wld", adfGeoTransform );
	if (status == 1)
		return (0);

	return (-1);
}

/* -------------------------------------------------------------------- */
void ComputeRasterMinMax(char *tmp, GDALRasterBandH hBand, double adfMinMax[2], 
			int nXSize, int nYSize, double z_min, double z_max) {
	/* Compute Min/Max of a sub-region. I'm forced to do this because the
	GDALComputeRasterMinMax works only on the entire dataset */
	int	i, bGotNoDataValue;
	GInt16	*tmpI16, *outI16;
	GUInt16	*tmpUI16, *outUI16;
	GInt32	*tmpI32, *outI32;
	GUInt32	*tmpUI32, *outUI32;
	float	*tmpF32;
	double	dfNoDataValue;

        dfNoDataValue = GDALGetRasterNoDataValue(hBand, &bGotNoDataValue);
	switch( GDALGetRasterDataType(hBand) ) {
		case GDT_Byte:
			for (i = 0; i < nXSize*nYSize; i++) {
				if( bGotNoDataValue && tmp[i] == dfNoDataValue ) continue;
				z_min = MIN(tmp[i], z_min);
				z_max = MAX(tmp[i], z_max);
			}
			break;
		case GDT_Int16:
			tmpI16 = (GInt16 *) tmp;
			for (i = 0; i < nXSize*nYSize; i++) {
				if( bGotNoDataValue && tmpI16[i] == dfNoDataValue ) continue;
				z_min = MIN(tmpI16[i], z_min);
				z_max = MAX(tmpI16[i], z_max);
			}
			break;
		case GDT_UInt16:
			tmpUI16 = (GUInt16 *) tmp;
			for (i = 0; i < nXSize*nYSize; i++) {
				if( bGotNoDataValue && tmpUI16[i] == dfNoDataValue ) continue;
				z_min = MIN(tmpUI16[i], z_min);
				z_max = MAX(tmpUI16[i], z_max);
			}
			break;
		case GDT_Int32:
			tmpI32 = (GInt32 *) tmp;
			for (i = 0; i < nXSize*nYSize; i++) {
				if( bGotNoDataValue && tmpI32[i] == dfNoDataValue ) continue;
				z_min = MIN(tmpI32[i], z_min);
				z_max = MAX(tmpI32[i], z_max);
			}
			break;
		case GDT_UInt32:
			tmpUI32 = (GUInt32 *) tmp;
			for (i = 0; i < nXSize*nYSize; i++) {
				if( bGotNoDataValue && tmpUI32[i] == dfNoDataValue ) continue;
				z_min = MIN(tmpUI32[i], z_min);
				z_max = MAX(tmpUI32[i], z_max);
			}
			break;
		case GDT_Float32:
			tmpF32 = (float *) tmp;
			for (i = 0; i < nXSize*nYSize; i++) {
				if( bGotNoDataValue && tmpF32[i] == dfNoDataValue ) continue;
				z_min = MIN(tmpF32[i], z_min);
				z_max = MAX(tmpF32[i], z_max);
			}
			break;
	}
	adfMinMax[0] = z_min;
	adfMinMax[1] = z_max;
}

/* -------------------------------------------------------------------- */
int decode_R (char *item, double *w, double *e, double *s, double *n) {
	char *text, string[BUFSIZ];
	
	/* Minimalist code to decode option -R extracted from GMT_get_common_args */
	
	int i, error = 0;
	double *p[4];
	
	p[0] = w;	p[1] = e;	p[2] = s;	p[3] = n;
			
	i = 0;
	strcpy (string, &item[2]);
	text = strtok (string, "/");
	while (text) {
		*p[i] = ddmmss_to_degree (text);
		i++;
		text = strtok (CNULL, "/");
	}
	if (item[strlen(item)-1] == 'r')	/* Rectangular box given, but valid here */
		error++;
	if (i != 4 || check_region (*p[0], *p[1], *p[2], *p[3]))
		error++;
	w = p[0];	e = p[1];
	s = p[2];	n = p[3];
	return (error);
}

/* -------------------------------------------------------------------- */
int check_region (double w, double e, double s, double n) {
	/* If region is given then we must have w < e and s < n */
	return ((w >= e || s >= n));
}

/* -------------------------------------------------------------------- */
double ddmmss_to_degree (char *text) {
	int i, colons = 0, suffix;
	double degree, minute, degfrac, second;

	for (i = 0; text[i]; i++) if (text[i] == ':') colons++;
	suffix = (int)text[i-1];	/* Last character in string */
	if (colons == 2) {	/* dd:mm:ss format */
		sscanf (text, "%lf:%lf:%lf", &degree, &minute, &second);
		degfrac = degree + Loc_copysign (minute / 60.0 + second / 3600.0, degree);
	}
	else if (colons == 1) {	/* dd:mm format */
		sscanf (text, "%lf:%lf", &degree, &minute);
		degfrac = degree + Loc_copysign (minute / 60.0, degree);
	}
	else
		degfrac = atof (text);
	if (suffix == 'W' || suffix == 'w' || suffix == 'S' || suffix == 's') degfrac = -degfrac;	/* Sign was given implicitly */
	return (degfrac);
}

/* -------------------------------------------------------------------- */
int decode_columns (char *txt, int *whichBands, int n_col) {
	int i, start, stop, pos, n = 0;
	char p[1024];

	pos = 0;
	while ((GMT_strtok (txt, ",", &pos, p))) {
		if (strchr (p, '-'))
			sscanf (p, "%d-%d", &start, &stop);
		else {
			sscanf (p, "%d", &start);
			stop = start;
		}
		stop = MIN (stop, n_col-0);
		for (i = start; i <= stop; i++) {
			whichBands[n] = i;
			n++;
		}
	}
	return (n);
}

/* -------------------------------------------------------------------- */
int GMT_strtok (const char *string, const char *sep, int *pos, char *token) {
        /* Reentrant replacement for strtok that uses no static variables.
         * Breaks string into tokens separated by one of more separator
         * characters (in sep).  Set *pos to 0 before first call.  Unlike
         * strtok, always pass the original string as first argument.
         * Returns 1 if it finds a token and 0 if no more tokens left.
         * pos is updated and token is returned.  char *token must point
         * to memory of length >= strlen (string).
         */

        int i, j, string_len;

        string_len = strlen (string);

        /* Wind up *pos to first non-separating character: */
        while (string[*pos] && strchr (sep, (int)string[*pos])) (*pos)++;

        if (*pos >= string_len || string_len == 0) return 0;    /* Got NULL string or no more string left to search */

        /* Search for next non-separating character */
        for (i = *pos; string[i] && !strchr (sep, (int)string[i]); i++);

        /* Copy token */
        j = i - *pos;
        strncpy (token, &string[*pos], j);
        token[j] = 0;

        /* Wind up *pos to next non-separating character */
        while (string[i] && strchr (sep, (int)string[i])) i++;
        *pos = i;

        return 1;
}
