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

/* Program:	gdawrite.c
 * Purpose:	matlab callable routine to write files supported by gdal
 *
 * Revision 5.0  07/12/2008 Accept GCPs in the "gcp" field of the hdr structure
 * Revision 4.0  10/11/2008 Added a "meta" field to the hdr structure
 * Revision 3.0  06/09/2007 Start to add a nodata option (not finished)
 * Revision 2.0  12/06/2007 Was not aware of receiving a WKT proj string
 * Revision 1.0  24/06/2006 Joaquim Luis
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
#define debug	0

#include "mex.h"
#include "gdal.h"
#include "ogr_srs_api.h"
#include "cpl_string.h"
#include "cpl_conv.h"

void DEBUGA(int n);

/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	ns, flipud = FALSE, i_x_nXYSize;
	char **papszOptions = NULL;
	char *pszFormat = "GTiff", *projWKT = NULL, *metaString = NULL; 
	double adfGeoTransform[6] = {0,1,0,0,0,1}; 
	double dfNoData;
	char *pszSRS_WKT = NULL;
	OGRSpatialReferenceH hSRS;
	GDALDatasetH hDstDS;
	GDALDriverH	hDriver;
	GDALRasterBandH hBand;
	GDALColorTableH	hColorTable = NULL;
	GDALColorEntry	sEntry;
	GDAL_GCP	*pasGCPs = NULL;
	int	nGCPCount = 0;

	const int *dim_array;
	int	nx, ny, i, j, m, n, nn, n_bands_in, registration = 1;
	int	n_dims, typeCLASS, mx_CLASS, nBytes, nColors;
	int	*nVector, *mVector, is_geog = 0;
	const char *fname;
	mxArray	*mx_ptr;
	void	*in_data;

	unsigned char *tmpByte, *outByte;
	unsigned short int *tmpUI16, *outUI16;
	short int *tmpI16, *outI16;
	int	*tmpI32, *outI32;
	unsigned int *tmpUI32, *outUI32;
	float	*tmpF32, *outF32;
	double	*tmpF64, *outF64, *ptr_d;
	char	**papszMetadataOptions = NULL;

	if (nrhs == 3) {
		if(!mxIsChar(prhs[1]))
			mexErrMsgTxt("GDALWRITE Second argument must be a valid string!");
		else 
			fname = (char *)mxArrayToString(prhs[1]);

		if(!mxIsChar(prhs[2]))
			mexErrMsgTxt("GDALWRITE Third argument must be a valid string!");
	}
	else if (nrhs == 2 && mxIsStruct(prhs[1])) {
		mx_ptr = mxGetField(prhs[1], 0, "driver");
		if (mx_ptr == NULL)
			mexErrMsgTxt("GDALWRITE Driver name not provided");
		pszFormat = (char *)mxArrayToString(mx_ptr);

		mx_ptr = mxGetField(prhs[1], 0, "name");
		if (mx_ptr == NULL)
			mexErrMsgTxt("GDALWRITE Output file name not provided");
		fname = (char *)mxArrayToString(mx_ptr);

		ptr_d = mxGetPr(mxGetField(prhs[1], 0, "ULx"));
		if (ptr_d == NULL)
			mexErrMsgTxt("GDALWRITE 'ULx' field not provided");
		adfGeoTransform[0] = *ptr_d;

		ptr_d = mxGetPr(mxGetField(prhs[1], 0, "Xinc"));
		if (ptr_d == NULL)
			mexErrMsgTxt("GDALWRITE 'Xinc' field not provided");
		adfGeoTransform[1] = *ptr_d;

		ptr_d = mxGetPr(mxGetField(prhs[1], 0, "ULy"));
		if (ptr_d == NULL)
			mexErrMsgTxt("GDALWRITE 'ULy' field not provided");
		adfGeoTransform[3] = *ptr_d;

		ptr_d = mxGetPr(mxGetField(prhs[1], 0, "Yinc"));
		if (ptr_d == NULL)
			mexErrMsgTxt("GDALWRITE 'Yinc' field not provided");
		adfGeoTransform[5] = -*ptr_d;

		/*ptr_d = mxGetPr(mxGetField(prhs[1], 0, "nodata"));
		if (ptr_d != NULL)
			dfNoData = *ptr_d;*/

		mx_ptr = mxGetField(prhs[1], 0, "Reg");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			registration = (int)ptr_d[0];
		}

		mx_ptr = mxGetField(prhs[1], 0, "Flip");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			flipud = (int)ptr_d[0];
		}

		mx_ptr = mxGetField(prhs[1], 0, "Geog");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			is_geog = (int)ptr_d[0];
		}

		mx_ptr = mxGetField(prhs[1], 0, "projWKT");
		if (mx_ptr != NULL)
			projWKT = (char *)mxArrayToString(mx_ptr);

		mx_ptr = mxGetField(prhs[1], 0, "meta");
		if (mx_ptr != NULL)
			metaString = (char *)mxArrayToString(mx_ptr);

		mx_ptr = mxGetField(prhs[1], 0, "Cmap");
		if (mx_ptr != NULL) {
			nColors = mxGetM(mx_ptr);
			ptr_d = mxGetPr(mx_ptr);
			hColorTable = GDALCreateColorTable (GPI_RGB); 
			for (i = 0; i < nColors; i++) {
				sEntry.c1 = (short)(ptr_d[i] * 255); 
				sEntry.c2 = (short)(ptr_d[i+nColors] * 255); 
				sEntry.c3 = (short)(ptr_d[i+2*nColors] * 255); 
				sEntry.c4 = (short)255; 
				GDALSetColorEntry( hColorTable, i, &sEntry ); 
			}
		}

		/* -------- Do we have GCPs? ----------------------- */
		mx_ptr = mxGetField(prhs[1], 0, "gcp");
		if (mx_ptr != NULL) {
			nGCPCount = mxGetM(mx_ptr);
			if (mxGetN(mx_ptr) != 4)
				mexErrMsgTxt("GDALWRITE: GCPs must be a Mx4 array");
			ptr_d = mxGetPr(mx_ptr);
			pasGCPs = (GDAL_GCP *) mxCalloc( nGCPCount, sizeof(GDAL_GCP) * nGCPCount );
			GDALInitGCPs( 1, pasGCPs + nGCPCount - 1 );
			for (i = 0; i < nGCPCount; i++) {
				pasGCPs[i].dfGCPPixel = ptr_d[i];
				pasGCPs[i].dfGCPLine = ptr_d[i+nGCPCount];
				pasGCPs[i].dfGCPX = ptr_d[i+2*nGCPCount];
				pasGCPs[i].dfGCPY = ptr_d[i+3*nGCPCount];
				pasGCPs[i].dfGCPZ = 0;
			}
		}
		/* -------------------------------------------------- */

		/* If grid limits were in grid registration, convert them to pixel reg */
		if (registration == 0) {
			adfGeoTransform[0] -= adfGeoTransform[1]/2.;
			adfGeoTransform[3] -= adfGeoTransform[5]/2.;
		}
	}
	else {
		mexPrintf("\tUsage: gdalwrite(data,hdr_struct)\n\n");
		GDALAllRegister();
        	mexPrintf( "The following format drivers are configured and support Create() method:\n" );
        	for( i = 0; i < GDALGetDriverCount(); i++ ) {
			hDriver = GDALGetDriver(i);
			if( GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATE, NULL ) != NULL)
				mexPrintf("%s: %s\n", GDALGetDriverShortName(hDriver), GDALGetDriverLongName(hDriver));
		}
		return;
	}

	n_dims = mxGetNumberOfDimensions(prhs[0]);
	dim_array=mxGetDimensions(prhs[0]);
	ny = dim_array[0];
	nx = dim_array[1];
	n_bands_in = dim_array[2];

	if (n_dims == 2)	/* Otherwise it would stay undefined */
		n_bands_in = 1;


	/* Find out in which data type was given the input array */
	if (mxIsUint8(prhs[0])) {
		typeCLASS = GDT_Byte;		nBytes = 1;
		outByte = (unsigned char *)mxCalloc (nx*ny, sizeof (unsigned char));
	}
	else if (mxIsUint16(prhs[0])) {
		typeCLASS = GDT_UInt16;		nBytes = 2;
		outUI16 = (unsigned short int *)mxCalloc (nx*ny, sizeof (short int));
	}
	else if (mxIsInt16(prhs[0])) {
		typeCLASS = GDT_Int16;		nBytes = 2;
		outI16 = (short int *)mxCalloc (nx*ny, sizeof (short int));
	}
	else if (mxIsInt32(prhs[0])) {
		typeCLASS = GDT_Int32;		nBytes = 4;
		outI32 = (int *)mxCalloc (nx*ny, sizeof (int));
	}
	else if (mxIsUint32(prhs[0])) {
		typeCLASS = GDT_UInt32;		nBytes = 4;
		outUI32 = (unsigned int *)mxCalloc (nx*ny, sizeof (int));
	}
	else if (mxIsSingle(prhs[0])) {
		typeCLASS = GDT_Float32;	nBytes = 4;
		outF32 = (float *)mxCalloc (nx*ny, sizeof (float));
	}
	else if (mxIsDouble(prhs[0])) {
		typeCLASS = GDT_Float64;	nBytes = 8;
		outF64 = (double *)mxCalloc (nx*ny, sizeof (double));
	}
	else
		mexErrMsgTxt("GDALWRITE Unknown input data class!");


	in_data = (void *)mxGetData(prhs[0]);

	GDALAllRegister();

	hDriver = GDALGetDriverByName( pszFormat ); 

	/* Use LZW compression with GeoTiff driver */
	if (!strcmp(pszFormat,"GTiff")) {
		papszOptions = CSLAddString( papszOptions, "COMPRESS=DEFLATE" ); 
		papszOptions = CSLAddString( papszOptions, "TILED=YES" ); 
		/*papszOptions = CSLAddString( papszOptions, "INTERLEAVE=BAND" ); */
	}

	if (metaString)
		papszOptions = CSLAddString( papszOptions, metaString ); 

	hDstDS = GDALCreate( hDriver, fname, nx, ny, n_bands_in, typeCLASS, papszOptions );

	if (hDstDS == NULL) {
		mexPrintf ("GDALOpen failed - %d\n%s\n",
                CPLGetLastErrorNo(), CPLGetLastErrorMsg());
		return;
	}
	GDALSetGeoTransform( hDstDS, adfGeoTransform ); 

	/* This was the only trick I found to set a "projection". The docs still have a long way to go */
	if (is_geog || projWKT) {
		hSRS = OSRNewSpatialReference( NULL );
		if (is_geog && !projWKT)	/* Only thing we know is that it is Geog */
			OSRSetFromUserInput( hSRS, "+proj=latlong +datum=WGS84" );
		else				/* Even if is_geog == TRUE, use the WKT string */ 
			OSRSetFromUserInput( hSRS, projWKT );
		OSRExportToWkt( hSRS, &pszSRS_WKT );
		OSRDestroySpatialReference( hSRS );
		GDALSetProjection( hDstDS, pszSRS_WKT );
		if ( nGCPCount == 0 )		/* Otherwise we still need this for setting the GCPs */
			CPLFree( pszSRS_WKT );
	}

	if ( nGCPCount != 0 ) {
DEBUGA(1);
		if (GDALSetGCPs( hDstDS, nGCPCount, pasGCPs, "" ) != CE_None)
			mexPrintf("WARNING: writing GCPs failed.\n");
		CPLFree( pszSRS_WKT );
	}

	nVector = mxCalloc(nx, sizeof(int));
	mVector = mxCalloc(ny, sizeof(int));
	for (n = 0; n < nx; n++) nVector[n] = n*ny-1 + ny;
	if (flipud)
		for (m = 0; m < ny; m++) mVector[m] = (ny-1)*nx - m*nx;
	else
		for (m = 0; m < ny; m++) mVector[m] = m*nx;

	for (i = 1; i <= n_bands_in; i++) {
		hBand = GDALGetRasterBand( hDstDS, i ); 
		if( i == 1 && hColorTable != NULL ) {
			if (GDALSetRasterColorTable( hBand, hColorTable ) == CE_Failure)
				mexPrintf("\tERROR creating Color Table");
			GDALDestroyColorTable( hColorTable );
		}
		/*GDALSetRasterNoDataValue(hBand,dfNoData);	test, which worked in a geotiff image */
		i_x_nXYSize = (i-1)*nx*ny;		/* We don't need to recompute this everytime */
		switch( typeCLASS ) {
			case GDT_Byte:
			 	tmpByte = (unsigned char *)in_data;	
				for (m = 0; m < ny; m++) for (n = 0; n < nx; n++) {
					nn = nVector[n]-m;
					outByte[mVector[m]+n] = tmpByte[nn + i_x_nXYSize];
				}
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outByte, nx, ny, typeCLASS, 0, 0 );
				break;
			case GDT_UInt16:
			 	tmpUI16 = (unsigned short int *)in_data;	
				for (m = 0; m < ny; m++) for (n = 0; n < nx; n++) {
					nn = nVector[n]-m;
					outUI16[mVector[m]+n] = tmpUI16[nn + i_x_nXYSize];
				}
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outUI16, nx, ny, typeCLASS, 0, 0 );
				break;
			case GDT_Int16:
			 	tmpI16 = (short int *)in_data;	
				for (m = 0; m < ny; m++) for (n = 0; n < nx; n++) {
					nn = nVector[n]-m;
					outI16[mVector[m]+n] = tmpI16[nn + i_x_nXYSize];
				}
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outI16, nx, ny, typeCLASS, 0, 0 );
				break;
			case GDT_UInt32:
			 	tmpUI32 = (unsigned int *)in_data;	
				for (m = 0; m < ny; m++) for (n = 0; n < nx; n++) {
					nn = nVector[n]-m;
					outUI32[mVector[m]+n] = tmpUI32[nn + i_x_nXYSize];
				}
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outUI32, nx, ny, typeCLASS, 0, 0 );
				break;
			case GDT_Int32:
			 	tmpI32 = (int *)in_data;	
				for (m = 0; m < ny; m++) for (n = 0; n < nx; n++) {
					nn = nVector[n]-m;
					outI32[mVector[m]+n] = tmpI32[nn + i_x_nXYSize];
				}
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outI32, nx, ny, typeCLASS, 0, 0 );
				break;
			case GDT_Float32:
			 	tmpF32 = (float *)in_data;	
				for (m = 0; m < ny; m++) for (n = 0; n < nx; n++) {
					nn = nVector[n]-m;
					outF32[mVector[m]+n] = tmpF32[nn + i_x_nXYSize];
				}
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outF32, nx, ny, typeCLASS, 0, 0 );
				break;
			case GDT_Float64:
			 	tmpF64 = (double *)in_data;	
				for (m = 0; m < ny; m++) for (n = 0; n < nx; n++) {
					nn = nVector[n]-m;
					outF64[mVector[m]+n] = tmpF64[nn + i_x_nXYSize];
				}
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outF64, nx, ny, typeCLASS, 0, 0 );
				break;
		}
	}

	switch( typeCLASS ) {
		case GDT_Byte:		mxFree((void *)outByte);	break;
		case GDT_UInt16:	mxFree((void *)outUI16);	break; 
		case GDT_Int16:		mxFree((void *)outI16);		break; 
		case GDT_UInt32:	mxFree((void *)outUI32);	break; 
		case GDT_Int32:		mxFree((void *)outI32);		break; 
		case GDT_Float32:	mxFree((void *)outF32);		break; 
		case GDT_Float64:	mxFree((void *)outF64);		break; 
	}

	GDALClose( hDstDS );
	mxFree(nVector);
	mxFree(mVector);
	if (nGCPCount) {
		GDALDeinitGCPs( nGCPCount, pasGCPs );
		mxFree((void *) pasGCPs );
	}
}

void DEBUGA(int n) {
#if debug
	mexPrintf("Merda %d\n",n);
#endif
}
