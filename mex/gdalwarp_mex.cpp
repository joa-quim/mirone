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

/* Program:	gdalwarp_mex.c
 * Purpose:	matlab callable routine to reproject files using gdal/proj4
 *
 * Revision 4.0  07/12/2008 Warp with GCPs Joaquim Luis
 * Revision 3.0  07/04/2008 Joaquim Luis
 * Revision 2.0  23/07/2007 Joaquim Luis
 * Revision 1.0  24/06/2006 Joaquim Luis
 *
 */

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#define debug	0

#include "mex.h"
#include "gdal.h"
#include "cpl_string.h"
#include "ogr_srs_api.h"
#include "cpl_conv.h"
#include "gdalwarper.h"
#include "ogr_spatialref.h"

mxArray *populate_metadata_struct (GDALDatasetH hDataset, int correct_bounds);
int ReportCorner(GDALDatasetH hDataset, double x, double y, double *xy_c);
int getNK(const mxArray *p, int which);
void DEBUGA(int n);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	nXYSize;
	double	adfGeoTransform[6] = {0,1,0,0,0,1}, adfDstGeoTransform[6];
	char	*pszSRS_WKT = NULL;
	char	**papszWarpOptions = NULL;
	GDALDatasetH	hSrcDS, hDstDS;
	GDALDriverH	hDriver;
	GDALRasterBandH hBand;
	GDALColorTableH	hColorTable = NULL;
	OGRSpatialReference oSrcSRS, oDstSRS; 
	GDALResampleAlg	interpMethod = GRA_NearestNeighbour;
	GDALTransformerFunc pfnTransformer = NULL;
	CPLErr		eErr;
	GDAL_GCP	*pasGCPs = NULL;
	static int runed_once = FALSE;	/* It will be set to true if reaches end of main */

	const int *dim_array;
	int	nx, ny, i, j, m, n, c, nBands, registration = 1;
	int	n_dims, typeCLASS, nBytes;
	char	*pszSrcSRS = NULL, *pszSrcWKT = NULL;
	char	*pszDstSRS = NULL, *pszDstWKT = NULL;
	void	*in_data;
	mxArray	*mx_ptr;

	unsigned char *tmpByte, *outByte;
	unsigned short int *tmpUI16, *outUI16;
	short int *tmpI16, *outI16;
	int	*tmpI32, *outI32;
	int	nPixels=0, nLines=0, nForceWidth=0, nForceHeight=0;
	int	nGCPCount = 0, nOrder = 0;
	unsigned int *tmpUI32, *outUI32;
	float	*tmpF32, *outF32;
	double	*tmpF64, *outF64, *ptr_d;
	double	dfMinX=0, dfMaxX=0, dfMinY=0, dfMaxY=0, dfResX=0, dfResY=0;
	double	adfExtent[4];
	double	dfXRes=0.0, dfYRes=0.0;
	double	dfWarpMemoryLimit = 0.0;
	double	*pdfDstNodata = NULL; 
	char	**papszMetadataOptions = NULL;
	char	*tmp, *txt;


	if (nrhs == 2 && mxIsStruct(prhs[1])) {
		mx_ptr = mxGetField(prhs[1], 0, "ULx");
		if (mx_ptr == NULL)
			mexErrMsgTxt("GDALWARP 'ULx' field not provided");
		ptr_d = mxGetPr(mx_ptr);
		adfGeoTransform[0] = *ptr_d;

		mx_ptr = mxGetField(prhs[1], 0, "Xinc");
		if (mx_ptr == NULL)
			mexErrMsgTxt("GDALWARP 'Xinc' field not provided");
		ptr_d = mxGetPr(mx_ptr);
		adfGeoTransform[1] = *ptr_d;

		mx_ptr = mxGetField(prhs[1], 0, "ULy");
		if (mx_ptr == NULL)
			mexErrMsgTxt("GDALWARP 'ULy' field not provided");
		ptr_d = mxGetPr(mx_ptr);
		adfGeoTransform[3] = *ptr_d;

		mx_ptr = mxGetField(prhs[1], 0, "Yinc");
		if (mx_ptr == NULL)
			mexErrMsgTxt("GDALWARP 'Yinc' field not provided");
		ptr_d = mxGetPr(mx_ptr);
		adfGeoTransform[5] = -*ptr_d;

		/* -------- See for resolution requests ------------ */
		mx_ptr = mxGetField(prhs[1], 0, "t_size");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			if (mxGetN(mx_ptr) == 2) {
				nForceWidth  = (int)ptr_d[0];
				nForceHeight = (int)ptr_d[1];
			}
			else if (mxGetN(mx_ptr) == 1) {	/* pick max(nrow,ncol) */
				if (mxGetM(prhs[0]) > getNK(prhs[0],1))
					nForceHeight = mxGetM(prhs[0]);
				else
					nForceWidth  = getNK(prhs[0], 1);
			}
			else {
				nForceHeight = mxGetM(prhs[0]);
				nForceWidth  = getNK(prhs[0], 1);
			}
		}

		mx_ptr = mxGetField(prhs[1], 0, "t_res");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			if (mxGetN(mx_ptr) == 2) {
				dfXRes = ptr_d[0];
				dfYRes = ptr_d[1];
			}
			else if (mxGetN(mx_ptr) == 1) {
				dfXRes = dfYRes = ptr_d[0];
			}
		}
		/* -------------------------------------------------- */

		/* -------- Change Warping cache size?  ------------ */
		mx_ptr = mxGetField(prhs[1], 0, "wm");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			dfWarpMemoryLimit = *ptr_d * 1024 * 1024;
		}
		/* -------------------------------------------------- */

		/* -------- Have a nodata value order? -------------- */
		mx_ptr = mxGetField(prhs[1], 0, "nodata");
		if (mx_ptr != NULL) {
			pdfDstNodata = mxGetPr(mx_ptr);
		}
		/* -------------------------------------------------- */

		/* -------- See for projection stuff ---------------- */
		mx_ptr = mxGetField(prhs[1], 0, "SrcProjSRS");
		if (mx_ptr != NULL)
			pszSrcSRS = (char *)mxArrayToString(mx_ptr);

		mx_ptr = mxGetField(prhs[1], 0, "SrcProjWKT");
		if (mx_ptr != NULL)
			pszSrcWKT = (char *)mxArrayToString(mx_ptr);

		mx_ptr = mxGetField(prhs[1], 0, "DstProjSRS");
		if (mx_ptr != NULL)
			pszDstSRS = (char *)mxArrayToString(mx_ptr);

		mx_ptr = mxGetField(prhs[1], 0, "DstProjWKT");
		if (mx_ptr != NULL)
			pszDstWKT = (char *)mxArrayToString(mx_ptr);
		/* -------------------------------------------------- */

		/* -------- Do we have GCPs? ----------------------- */
		mx_ptr = mxGetField(prhs[1], 0, "gcp");
		if (mx_ptr != NULL) {
			nGCPCount = mxGetM(mx_ptr);
			if (mxGetN(mx_ptr) != 4)
				mexErrMsgTxt("GDALWARP: GCPs must be a Mx4 array");
			ptr_d = mxGetPr(mx_ptr);
			pasGCPs = (GDAL_GCP *) mxCalloc( nGCPCount, sizeof(GDAL_GCP) );
			GDALInitGCPs( 1, pasGCPs + nGCPCount - 1 );
			for (i = 0; i < nGCPCount; i++) {
				pasGCPs[i].dfGCPPixel = ptr_d[i];
				pasGCPs[i].dfGCPLine = ptr_d[i+nGCPCount];
				pasGCPs[i].dfGCPX = ptr_d[i+2*nGCPCount];
				pasGCPs[i].dfGCPY = ptr_d[i+3*nGCPCount];
				pasGCPs[i].dfGCPZ = 0;
			}
		}
			/* ---- Have we an order request? --- */
		mx_ptr = mxGetField(prhs[1], 0, "order");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			nOrder = (int)*ptr_d;
			if (nOrder != -1 || nOrder != 0 || nOrder != 1 || nOrder != 2 || nOrder != 3)
				nOrder = 0;
		}
		/* -------------------------------------------------- */

		mx_ptr = mxGetField(prhs[1], 0, "ResampleAlg");
		if (mx_ptr != NULL) {
			txt = (char *)mxArrayToString(mx_ptr);
			if (!strcmp(txt,"nearest"))
				interpMethod = GRA_NearestNeighbour;
			else if (!strcmp(txt,"bilinear"))
				interpMethod = GRA_Bilinear;
			else if (!strcmp(txt,"cubic") || !strcmp(txt,"bicubic"))
				interpMethod = GRA_Cubic;
			else if (!strcmp(txt,"spline"))
				interpMethod = GRA_CubicSpline;
		}

		/* If grid limits were in grid registration, convert them to pixel reg */
		mx_ptr = mxGetField(prhs[1], 0, "Reg");
		if (mx_ptr != NULL) {
			ptr_d = mxGetPr(mx_ptr);
			registration = (int)ptr_d[0];
		}

		if (registration == 0) {
			adfGeoTransform[0] -= adfGeoTransform[1]/2.;
			adfGeoTransform[3] -= adfGeoTransform[5]/2.;
		}
	}
	else {
		mexPrintf("Usage: B = gdalwarp_mex(IMG,HDR_STRUCT)\n\n");
		mexPrintf("\tIMG -> is a Mx2 or Mx3 array with an grid/image data to reproject\n");
		mexPrintf("\tHDR_STRUCT -> is a structure with the following fields:\n");
		mexPrintf("\t\t'ULx' X coordinate of the uper left corner\n");
		mexPrintf("\t\t'ULy' Y coordinate of the uper left corner\n");
		mexPrintf("\t\t'Xinc' distance between columns in target grid/image coordinates\n");
		mexPrintf("\t\t'Yinc' distance between rows in target grid/image coordinates\n");
		mexPrintf("\t\t'SrcProjSRS', 'SrcProjWKT' -> Source projection string\n");
		mexPrintf("\t\t'DstProjSRS', 'DstProjWKT' -> Target projection string\n");
		mexPrintf("\t\t\tSRS stands for a string of the type used by proj4\n");
		mexPrintf("\t\t\tWKT stands for a string on the 'Well Known Text' format\n\n");
		mexPrintf("\t\t\tIf one of the Src or Dst fields is absent a GEOGRAPHIC WGS84 is assumed\n");
		mexPrintf("\nOPTIONS\n");
		mexPrintf("\t\t'gcp' a [Mx4] array with Ground Control Points\n");
		mexPrintf("\t\t't_size' a [width height] vector to set output file size in pixels\n");
		mexPrintf("\t\t't_res' a [xres yres] vector to set output file resolution (in target georeferenced units)\n");
		mexPrintf("\t\t'wm' amount of memory (in megabytes) that the warp API is allowed to use for caching\n");
		mexPrintf("\t\t'nodata' Set nodata values for output bands.\n");
		mexPrintf("\t\t'ResampleAlg' To set up the algorithm used during warp operation. Options are: \n");
		mexPrintf("\t\t\t'nearest' Use nearest neighbour resampling (default, fastest algorithm, worst interpolation quality).\n");
		mexPrintf("\t\t\t'bilinear' Use bilinear resampling.\n");
		mexPrintf("\t\t\t'cubic' Use cubic resampling.\n");
		mexPrintf("\t\t\t'spline' Use cubic spline resampling.\n\n");

		if (!runed_once)		/* Do next call only at first time this MEX is loaded */
			GDALAllRegister();

        	mexPrintf( "The following format drivers are configured and support Create() method:\n" );
        	for( i = 0; i < GDALGetDriverCount(); i++ ) {
			hDriver = GDALGetDriver(i);
			if( GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATE, NULL ) != NULL)
				mexPrintf("%s: %s\n", GDALGetDriverShortName(hDriver), 
							GDALGetDriverLongName(hDriver));
		}
		return;
	}

	n_dims = mxGetNumberOfDimensions(prhs[0]);
	dim_array=mxGetDimensions(prhs[0]);
	ny = dim_array[0];
	nx = dim_array[1];
	nBands = dim_array[2];

	if (n_dims == 2)	/* Otherwise it would stay undefined */
		nBands = 1;

	/* Find out in which data type was given the input array */
	if (mxIsUint8(prhs[0])) {
		typeCLASS = GDT_Byte;		nBytes = 1;
		outByte = (unsigned char *)mxMalloc (nx*ny * sizeof(unsigned char));
	}
	else if (mxIsUint16(prhs[0])) {
		typeCLASS = GDT_UInt16;		nBytes = 2;
		outUI16 = (unsigned short int *)mxMalloc (nx*ny * sizeof(short int));
	}
	else if (mxIsInt16(prhs[0])) {
		typeCLASS = GDT_Int16;		nBytes = 2;
		outI16 = (short int *)mxMalloc (nx*ny * sizeof(short int));
	}
	else if (mxIsInt32(prhs[0])) {
		typeCLASS = GDT_Int32;		nBytes = 4;
		outI32 = (int *)mxMalloc (nx*ny * sizeof(int));
	}
	else if (mxIsUint32(prhs[0])) {
		typeCLASS = GDT_UInt32;		nBytes = 4;
		outUI32 = (unsigned int *)mxMalloc (nx*ny * sizeof(int));
	}
	else if (mxIsSingle(prhs[0])) {
		typeCLASS = GDT_Float32;	nBytes = 4;
		outF32 = (float *)mxMalloc (nx*ny * sizeof(float));
	}
	else if (mxIsDouble(prhs[0])) {
		typeCLASS = GDT_Float64;	nBytes = 8;
		outF64 = (double *)mxMalloc (nx*ny * sizeof(double));
	}
	else
		mexErrMsgTxt("GDALWARP Unknown input data class!");


	in_data = (void *)mxGetData(prhs[0]);

	if (!runed_once)		/* Do next call only at first time this MEX is loaded */
		GDALAllRegister();

	hDriver = GDALGetDriverByName( "MEM" ); 

	hSrcDS = GDALCreate( hDriver, "mem", nx, ny, nBands, (GDALDataType)typeCLASS, NULL );
	if (hSrcDS == NULL) {
		mexPrintf ("GDALOpen failed - %d\n%s\n", CPLGetLastErrorNo(), CPLGetLastErrorMsg());
		return;
	}
	GDALSetGeoTransform( hSrcDS, adfGeoTransform ); 

	/* ---------- Set the Source projection ---------------------------- */
	/* If it was not provided assume it is Geog WGS84 */
	if (pszSrcSRS == NULL && pszSrcWKT == NULL)
		oSrcSRS.SetWellKnownGeogCS( "WGS84" ); 
	else if (pszSrcWKT != NULL)
		oSrcSRS.importFromWkt( &pszSrcWKT );

	else {
		if( oSrcSRS.SetFromUserInput( pszSrcSRS ) != OGRERR_NONE )
			mexErrMsgTxt("GDAL_WARP_MEX: Translating source SRS failed.");
	}
	if (pszSrcWKT == NULL)
		oSrcSRS.exportToWkt( &pszSrcWKT );

	GDALSetProjection( hSrcDS, pszSrcWKT );	
	//pszSrcWKT = (char *)GDALGetProjectionRef( hSrcDS );
	CPLAssert( pszSrcWKT != NULL && strlen(pszSrcWKT) > 0 );
	/* ------------------------------------------------------------------ */


	/* -------------- Copy input data into the hSrcDS dataset ----------- */
	for (i = 1; i <= nBands; i++) {
		hBand = GDALGetRasterBand( hSrcDS, i ); 
		nXYSize = (i-1)*nx*ny;
		switch( typeCLASS ) {
			case GDT_Byte:
			 	tmpByte = (unsigned char *)in_data;	
				for (m = ny-1, c = 0; m >= 0; m--) for (n = 0; n < nx; n++)
					outByte[c++] = tmpByte[m + n*ny + nXYSize];
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outByte, nx, ny, (GDALDataType)typeCLASS, 0, 0 );
				break;
			case GDT_UInt16:
			 	tmpUI16 = (unsigned short int *)in_data;	
				for (m = ny-1, c = 0; m >= 0; m--) for (n = 0; n < nx; n++)
					outUI16[c++] = tmpUI16[m + n*ny + nXYSize];
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outUI16, nx, ny, (GDALDataType)typeCLASS, 0, 0 );
				break;
			case GDT_Int16:
			 	tmpI16 = (short int *)in_data;	
				for (m = ny-1, c = 0; m >= 0; m--) for (n = 0; n < nx; n++)
					outI16[c++] = tmpI16[m + n*ny + nXYSize];
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outI16, nx, ny, (GDALDataType)typeCLASS, 0, 0 );
				break;
			case GDT_UInt32:
			 	tmpUI32 = (unsigned int *)in_data;	
				for (m = ny-1, c = 0; m >= 0; m--) for (n = 0; n < nx; n++)
					outUI32[c++] = tmpUI32[m + n*ny + nXYSize];
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outUI32, nx, ny, (GDALDataType)typeCLASS, 0, 0 );
				break;
			case GDT_Int32:
			 	tmpI32 = (int *)in_data;	
				for (m = ny-1, c = 0; m >= 0; m--) for (n = 0; n < nx; n++)
					outI32[c++] = tmpI32[m + n*ny + nXYSize];
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outI32, nx, ny, (GDALDataType)typeCLASS, 0, 0 );
				break;
			case GDT_Float32:
			 	tmpF32 = (float *)in_data;	
				for (m = ny-1, c = 0; m >= 0; m--) for (n = 0; n < nx; n++)
					outF32[c++] = tmpF32[m + n*ny + nXYSize];
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outF32, nx, ny, (GDALDataType)typeCLASS, 0, 0 );
				break;
			case GDT_Float64:
			 	tmpF64 = (double *)in_data;	
				for (m = ny-1, c = 0; m >= 0; m--) for (n = 0; n < nx; n++)
					outF64[c++] = tmpF64[m + n*ny + nXYSize];
				GDALRasterIO( hBand, GF_Write, 0, 0, nx, ny,outF64, nx, ny, (GDALDataType)typeCLASS, 0, 0 );
				break;
		}
	}

	/* ---------- Set up the Target coordinate system ------------------- */
	/* If it was not provided assume it is Geog WGS84 */
	CPLErrorReset();
	if (pszDstSRS == NULL && pszDstWKT == NULL)
		oDstSRS.SetWellKnownGeogCS( "WGS84" ); 
	else if (pszDstWKT != NULL)
		oDstSRS.importFromWkt( &pszDstWKT );
	else {
		if( oDstSRS.SetFromUserInput( pszDstSRS ) != OGRERR_NONE )
			mexErrMsgTxt("GDAL_WARP_MEX: Translating target SRS failed.");
	}
	if (pszDstWKT == NULL)
		oDstSRS.exportToWkt( &pszDstWKT );
	/* ------------------------------------------------------------------ */

	if ( nGCPCount != 0 ) {
		if (GDALSetGCPs(hSrcDS, nGCPCount, pasGCPs, "") != CE_None)
			mexPrintf("GDALWARP WARNING: writing GCPs failed.\n");
	}

	/* Create a transformer that maps from source pixel/line coordinates
	   to destination georeferenced coordinates (not destination pixel line) 
	   We do that by omitting the destination dataset handle (setting it to NULL). */

	void *hTransformArg;

	hTransformArg = GDALCreateGenImgProjTransformer(hSrcDS, pszSrcWKT, NULL, pszDstWKT, 
											nGCPCount == 0 ? FALSE : TRUE, 0, nOrder);
	if( hTransformArg == NULL )
		mexErrMsgTxt("GDALTRANSFORM: Generating transformer failed.");

	GDALTransformerInfo *psInfo = (GDALTransformerInfo*)hTransformArg;

	/* -------------------------------------------------------------------------- */
	/*      Get approximate output georeferenced bounds and resolution for file
	/* -------------------------------------------------------------------------- */
	if (GDALSuggestedWarpOutput2(hSrcDS, GDALGenImgProjTransform, hTransformArg, 
	                             adfDstGeoTransform, &nPixels, &nLines, adfExtent,
	                             0) != CE_None ) {
	    GDALClose(hSrcDS);
		mexErrMsgTxt("GDALWARP: GDALSuggestedWarpOutput2 failed.");
	}

	if (CPLGetConfigOption( "CHECK_WITH_INVERT_PROJ", NULL ) == NULL) {
		double MinX = adfExtent[0];
		double MaxX = adfExtent[2];
		double MaxY = adfExtent[3];
		double MinY = adfExtent[1];
		int bSuccess = TRUE;
            
		/* Check that the the edges of the target image are in the validity area */
		/* of the target projection */
#define N_STEPS 20
		for (i = 0; i <= N_STEPS && bSuccess; i++) {
			for (j = 0; j <= N_STEPS && bSuccess; j++) {
				double dfRatioI = i * 1.0 / N_STEPS;
				double dfRatioJ = j * 1.0 / N_STEPS;
				double expected_x = (1 - dfRatioI) * MinX + dfRatioI * MaxX;
				double expected_y = (1 - dfRatioJ) * MinY + dfRatioJ * MaxY;
				double x = expected_x;
				double y = expected_y;
				double z = 0;
				/* Target SRS coordinates to source image pixel coordinates */
				if (!psInfo->pfnTransform(hTransformArg, TRUE, 1, &x, &y, &z, &bSuccess) || !bSuccess)
					bSuccess = FALSE;
				/* Source image pixel coordinates to target SRS coordinates */
				if (!psInfo->pfnTransform(hTransformArg, FALSE, 1, &x, &y, &z, &bSuccess) || !bSuccess)
					bSuccess = FALSE;
				if (fabs(x - expected_x) > (MaxX - MinX) / nPixels ||
					fabs(y - expected_y) > (MaxY - MinY) / nLines)
					bSuccess = FALSE;
			}
		}
            
		/* If not, retry with CHECK_WITH_INVERT_PROJ=TRUE that forces ogrct.cpp */
		/* to check the consistency of each requested projection result with the */
		/* invert projection */
		if (!bSuccess) {
			CPLSetConfigOption( "CHECK_WITH_INVERT_PROJ", "TRUE" );
			CPLDebug("WARP", "Recompute out extent with CHECK_WITH_INVERT_PROJ=TRUE");

			if (GDALSuggestedWarpOutput2(hSrcDS, GDALGenImgProjTransform, hTransformArg, 
			                             adfDstGeoTransform, &nPixels, &nLines, adfExtent,
			                              0) != CE_None ) {
			    GDALClose(hSrcDS);
				mexErrMsgTxt("GDALWARO: GDALSuggestedWarpOutput2 failed.");
			}
		}
	}

	/* -------------------------------------------------------------------- */
	/*      Expand the working bounds to include this region, ensure the    */
	/*      working resolution is no more than this resolution.             */
	/* -------------------------------------------------------------------- */
	if( dfMaxX == 0.0 && dfMinX == 0.0 ) {
		dfMinX = adfExtent[0];
		dfMaxX = adfExtent[2];
		dfMaxY = adfExtent[3];
		dfMinY = adfExtent[1];
		dfResX = adfDstGeoTransform[1];
		dfResY = ABS(adfDstGeoTransform[5]);
	}
	else {
		dfMinX = MIN(dfMinX,adfExtent[0]);
		dfMaxX = MAX(dfMaxX,adfExtent[2]);
		dfMaxY = MAX(dfMaxY,adfExtent[3]);
		dfMinY = MIN(dfMinY,adfExtent[1]);
		dfResX = MIN(dfResX,adfDstGeoTransform[1]);
		dfResY = MIN(dfResY,ABS(adfDstGeoTransform[5]));
	}

	GDALDestroyGenImgProjTransformer( hTransformArg );

	/* -------------------------------------------------------------------- */
	/*      Turn the suggested region into a geotransform and suggested     */
	/*      number of pixels and lines.                                     */
	/* -------------------------------------------------------------------- */

	adfDstGeoTransform[0] = dfMinX;
	adfDstGeoTransform[1] = dfResX;
	adfDstGeoTransform[2] = 0.0;
	adfDstGeoTransform[3] = dfMaxY;
	adfDstGeoTransform[4] = 0.0;
	adfDstGeoTransform[5] = -1 * dfResY;

	nPixels = (int) ((dfMaxX - dfMinX) / dfResX + 0.5);
	nLines  = (int) ((dfMaxY - dfMinY) / dfResY + 0.5);

	/* -------------------------------------------------------------------- */
	/*      Did the user override some parameters?                          */
	/* -------------------------------------------------------------------- */
	if( dfXRes != 0.0 && dfYRes != 0.0 ) {
		dfMinX = adfDstGeoTransform[0];
		dfMaxX = adfDstGeoTransform[0] + adfDstGeoTransform[1] * nPixels;
		dfMaxY = adfDstGeoTransform[3];
		dfMinY = adfDstGeoTransform[3] + adfDstGeoTransform[5] * nLines;

		nPixels = (int) ((dfMaxX - dfMinX + (dfXRes/2.0)) / dfXRes);
		nLines = (int) ((dfMaxY - dfMinY + (dfYRes/2.0)) / dfYRes);
		adfDstGeoTransform[0] = dfMinX;
		adfDstGeoTransform[3] = dfMaxY;
		adfDstGeoTransform[1] = dfXRes;
		adfDstGeoTransform[5] = -dfYRes;
	}
	else if( nForceWidth != 0 && nForceHeight != 0 ) {
		dfXRes = (dfMaxX - dfMinX) / nForceWidth;
		dfYRes = (dfMaxY - dfMinY) / nForceHeight;

		adfDstGeoTransform[0] = dfMinX;
		adfDstGeoTransform[3] = dfMaxY;
		adfDstGeoTransform[1] = dfXRes;
		adfDstGeoTransform[5] = -dfYRes;

		nPixels = nForceWidth;
		nLines = nForceHeight;
	}
	else if( nForceWidth != 0) {
		dfXRes = (dfMaxX - dfMinX) / nForceWidth;
		dfYRes = dfXRes;

		adfDstGeoTransform[0] = dfMinX;
		adfDstGeoTransform[3] = dfMaxY;
		adfDstGeoTransform[1] = dfXRes;
		adfDstGeoTransform[5] = -dfYRes;

		nPixels = nForceWidth;
		nLines = (int) ((dfMaxY - dfMinY + (dfYRes/2.0)) / dfYRes);
	}
	else if( nForceHeight != 0) {
		dfYRes = (dfMaxY - dfMinY) / nForceHeight;
		dfXRes = dfYRes;

		adfDstGeoTransform[0] = dfMinX;
		adfDstGeoTransform[3] = dfMaxY;
		adfDstGeoTransform[1] = dfXRes;
		adfDstGeoTransform[5] = -dfYRes;

		nPixels = (int) ((dfMaxX - dfMinX + (dfXRes/2.0)) / dfXRes);
		nLines = nForceHeight;
	}

	/* --------------------- Create the output --------------------------- */
	hDstDS = GDALCreate( hDriver, "mem", nPixels, nLines, 
			GDALGetRasterCount(hSrcDS), (GDALDataType)typeCLASS, NULL );
    
	CPLAssert( hDstDS != NULL );

	/* -------------- Write out the projection definition ---------------- */
	GDALSetProjection( hDstDS, pszDstWKT );
	GDALSetGeoTransform( hDstDS, adfDstGeoTransform );

	/* --------------------- Setup warp options -------------------------- */
	GDALWarpOptions *psWO = GDALCreateWarpOptions();

	psWO->hSrcDS = hSrcDS;
	psWO->hDstDS = hDstDS;

	psWO->nBandCount = nBands;
	psWO->panSrcBands = (int *) CPLMalloc(psWO->nBandCount * sizeof(int) );
	psWO->panDstBands = (int *) CPLMalloc(psWO->nBandCount * sizeof(int) );
	for( i = 0; i < nBands; i++ ) {
		psWO->panSrcBands[i] = i+1;
		psWO->panDstBands[i] = i+1;
	}

	if( dfWarpMemoryLimit != 0.0 )
		psWO->dfWarpMemoryLimit = dfWarpMemoryLimit;

	/* --------------------- Setup the Resampling Algo ------------------- */
	psWO->eResampleAlg = interpMethod;


	/* --------------------- Setup NODATA options ------------------------ */
	papszWarpOptions = CSLSetNameValue(papszWarpOptions, "INIT_DEST", "NO_DATA" );

	if ( pdfDstNodata == NULL && (typeCLASS == GDT_Float32 || typeCLASS == GDT_Float64) ) {
		pdfDstNodata = (double *) mxCalloc((size_t)1, sizeof(double));
		*pdfDstNodata = mxGetNaN();
	}
	else if (pdfDstNodata != NULL) {
#define CLAMP(val,type,minval,maxval) \
    do { if (val < minval) { val = minval; } \
    else if (val > maxval) { val = maxval; } \
    else if (val != (type)val) { val = (type)(val + 0.5); } } \
    while(0)
		switch( typeCLASS ) {
			case GDT_Byte:
				CLAMP(pdfDstNodata[0], GByte, 0.0, 255.0);
				break;
			case GDT_UInt16:
				CLAMP(pdfDstNodata[0], GInt16, -32768.0, 32767.0);
				break;
			case GDT_Int16:
				CLAMP(pdfDstNodata[0], GUInt16, 0.0, 65535.0);
				break;
			case GDT_UInt32:
				CLAMP(pdfDstNodata[0], GInt32, -2147483648.0, 2147483647.0);
				break;
			case GDT_Int32:
				CLAMP(pdfDstNodata[0], GUInt32, 0.0, 4294967295.0);
				break;
			default:
				break;
		}
	}

	psWO->papszWarpOptions = CSLDuplicate(papszWarpOptions);

	if (pdfDstNodata != NULL) {
		psWO->padfDstNoDataReal = (double *) CPLMalloc(psWO->nBandCount*sizeof(double));
		psWO->padfDstNoDataImag = (double *) CPLMalloc(psWO->nBandCount*sizeof(double));
		for (i = 0; i < nBands; i++) {
                        psWO->padfDstNoDataReal[i] = pdfDstNodata[0];
                        psWO->padfDstNoDataImag[i] = 0.0;
			GDALSetRasterNoDataValue( GDALGetRasterBand(hDstDS, i+1), pdfDstNodata[0]);
		}
	}

	/* ------------ Establish reprojection transformer ------------------- */
	psWO->pTransformerArg = GDALCreateGenImgProjTransformer( hSrcDS, GDALGetProjectionRef(hSrcDS), 
							hDstDS, GDALGetProjectionRef(hDstDS), 
							nGCPCount == 0 ? FALSE : TRUE, 0.0, nOrder );
	psWO->pfnTransformer = GDALGenImgProjTransform;

	/* ----------- Initialize and execute the warp operation ------------- */
	GDALWarpOperation oOperation;

	oOperation.Initialize( psWO );
	eErr = oOperation.ChunkAndWarpImage( 0, 0, GDALGetRasterXSize( hDstDS ),
						GDALGetRasterYSize( hDstDS ) );
	CPLAssert( eErr == CE_None );

	GDALDestroyGenImgProjTransformer( psWO->pTransformerArg );
	GDALDestroyWarpOptions( psWO );
	GDALClose( hSrcDS );

	/* ------------ Free memory used to fill the hSrcDS dataset ---------- */
	switch( typeCLASS ) {
		case GDT_Byte:		mxFree((void *)outByte);	break;
		case GDT_UInt16:	mxFree((void *)outUI16);	break; 
		case GDT_Int16:		mxFree((void *)outI16);		break; 
		case GDT_UInt32:	mxFree((void *)outUI32);	break; 
		case GDT_Int32:		mxFree((void *)outI32);		break; 
		case GDT_Float32:	mxFree((void *)outF32);		break; 
		case GDT_Float64:	mxFree((void *)outF64);		break; 
	}

	int out_dims[3];
	out_dims[0] = nLines;
	out_dims[1] = nPixels;
	out_dims[2] = nBands;
	plhs[0] = mxCreateNumericArray (n_dims,out_dims,mxGetClassID(prhs[0]), mxREAL);
	tmp = (char *)mxCalloc(nPixels * nLines, nBytes);

	/* ------ Allocate memory to be used in filling the hDstDS dataset ---- */
	switch( typeCLASS ) {
		case GDT_Byte:
			outByte = (unsigned char *)mxGetData(plhs[0]);		break;
		case GDT_UInt16:
			outUI16 = (unsigned short int *)mxGetData(plhs[0]);	break;
		case GDT_Int16:
			outI16 = (short int *)mxGetData(plhs[0]);		break;
		case GDT_UInt32:
			outUI32 = (unsigned int *)mxGetData(plhs[0]);		break;
		case GDT_Int32:
			outI32 = (int *)mxGetData(plhs[0]);			break;
		case GDT_Float32:
			outF32 = (float *)mxGetData(plhs[0]);			break;
		case GDT_Float64:
			outF64 = (double *)mxGetData(plhs[0]);			break;
	}

	/* ----------- Copy the output hSrcDS dataset data into plhs  ---------- */
	for (i = 1; i <= nBands; i++) {
		hBand = GDALGetRasterBand( hDstDS, i ); 
		GDALRasterIO( hBand, GF_Read, 0, 0, nPixels, nLines, tmp, nPixels, nLines,
				(GDALDataType)typeCLASS, 0, 0 );
		nXYSize = (i-1) * nPixels * nLines;
		switch( typeCLASS ) {
			case GDT_Byte:
				for (m = nLines-1, c = 0; m >= 0; m--) for (n = 0; n < nPixels; n++)
					outByte[m + n*nLines + nXYSize] = tmp[c++];
				break;
			case GDT_UInt16:
				tmpUI16 = (GUInt16 *) tmp;
				for (m = nLines-1, c = 0; m >= 0; m--) for (n = 0; n < nPixels; n++)
					outUI16[m + n*nLines + nXYSize] = tmpUI16[c++];
				break;
			case GDT_Int16:
				tmpI16 = (GInt16 *) tmp;
				for (m = nLines-1, c = 0; m >= 0; m--) for (n = 0; n < nPixels; n++)
					outI16[m + n*nLines + nXYSize] = tmpI16[c++];
				break;
			case GDT_UInt32:
				tmpUI32 = (GUInt32 *) tmp;
				for (m = nLines-1, c = 0; m >= 0; m--) for (n = 0; n < nPixels; n++)
					outUI32[m + n*nLines + nXYSize] = tmpUI32[c++];
				break;
			case GDT_Int32:
				tmpI32 = (GInt32 *) tmp;
				for (m = nLines-1, c = 0; m >= 0; m--) for (n = 0; n < nPixels; n++)
					outI32[m + n*nLines + nXYSize] = tmpI32[c++];
				break;
			case GDT_Float32:
				tmpF32 = (float *) tmp;
				for (m = nLines-1, c = 0; m >= 0; m--) for (n = 0; n < nPixels; n++)
					outF32[m + n*nLines + nXYSize] = tmpF32[c++];
				break;
			case GDT_Float64:
				tmpF64 = (double *) tmp;
				for (m = nLines-1, c = 0; m >= 0; m--) for (n = 0; n < nPixels; n++)
					outF64[m + n*nLines + nXYSize] = tmpF64[c++];
				break;
		}
	}

	mxFree(tmp);
	if (nGCPCount) {
		GDALDeinitGCPs( nGCPCount, pasGCPs );	/* makes this mex crash in the next call - Is it still true??? */
		mxFree((void *) pasGCPs );
	}

	if (nlhs == 2)
		plhs[1] = populate_metadata_struct (hDstDS, 1);

	runed_once = TRUE;	/* Signals that next call won't need to call GDALAllRegister() again */

	/*GDALDestroyDriverManager();
	OGRFree(pszDstWKT);*/
	GDALClose( hDstDS );
	CSLDestroy( papszWarpOptions );
	if (pszDstWKT && strlen(pszDstWKT) > 1 ) OGRFree(pszDstWKT);	
	if (pszSrcWKT && strlen(pszSrcWKT) > 1 ) OGRFree(pszSrcWKT);
}
	
/* ---------------------------------------------------------------------- */
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
 *                A structure array, one element for each overview present.  If
 *                empty, then there are no overviews.
 *            NoDataValue:
 *                When passed back to MATLAB, one can set pixels with this value to NaN.
 *            ColorMap:
 *                A Mx3 double array with the colormap, or empty if it does not exists
 *
 * */
mxArray *populate_metadata_struct (GDALDatasetH hDataset, int correct_bounds) {

	/* These are used to define the metadata structure about available GDAL drivers. */
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

	/*
	 * These will be matlab structures that hold the metadata.
	 * "metadata_struct" actually encompasses "band_struct",
	 * which encompasses "overview_struct"
	 * */
	mxArray *metadata_struct;
	mxArray *band_struct;
	mxArray *overview_struct;
	mxArray *corner_struct;

	int	overview, band_number;	/* Loop indices */
	double	*dptr;			/* short cut to the mxArray data */
	double	*dptr2;			/*        ""           */

	GDALDriverH hDriver;		/* This is the driver chosen by the GDAL library to query the dataset. */
	GDALRasterBandH hBand, overview_hBand;

	int num_overview_fields;	/* Number of metadata items for each overview structure. */
	int status;			/* success or failure */
	double 	adfGeoTransform[6];	/* bounds on the dataset */
	int num_overviews;		/* Number of overviews in the current band. */
	int num_struct_fields;		/* number of fields in the metadata structures. */
	int num_band_fields;
	int num_corner_fields;
	char *fieldnames[100];		/* this array contains the names of the fields of the metadata structure. */
	char *band_fieldnames[100];
	char *corner_fieldnames[100];
	char *overview_fieldnames[100];
	int xSize, ySize, raster_count; /* Dimensions of the dataset */
	int gdal_type;			/* Datatype of the bands. */
	double tmpdble;			/* temporary value */
	double	xy_c[2];		/* Corner coordinates */
	int	dims[2];
	int	bGotMin, bGotMax;	/* To know if driver transmited Min/Max */
	double	adfMinMax[2];		/* Dataset Min Max */
	double	z_min = 1e50, z_max = -1e50;

	/* Create the metadata structure. Just one element, with XXX fields. */
	num_struct_fields = 10;
	fieldnames[0] = strdup ("ProjectionRef");
	fieldnames[1] = strdup ("GeoTransform");
	fieldnames[2] = strdup ("DriverShortName");
	fieldnames[3] = strdup ("DriverLongName");
	fieldnames[4] = strdup ("RasterXSize");
	fieldnames[5] = strdup ("RasterYSize");
	fieldnames[6] = strdup ("RasterCount");
	fieldnames[7] = strdup ("Band");
	fieldnames[8] = strdup ("Corners");
	fieldnames[9] = strdup("GMT_hdr");
	metadata_struct = mxCreateStructMatrix ( 1, 1, num_struct_fields, (const char **)fieldnames );

	/* Record the ProjectionRef. */
	mxProjectionRef = mxCreateString ( GDALGetProjectionRef( hDataset ) );
	mxSetField ( metadata_struct, 0, "ProjectionRef", mxProjectionRef );

	/* Record the geotransform. */
	mxGeoTransform = mxCreateNumericMatrix ( 6, 1, mxDOUBLE_CLASS, mxREAL );
	dptr = mxGetPr ( mxGeoTransform );

	if( GDALGetGeoTransform( hDataset, adfGeoTransform ) == CE_None ) {
		dptr[0] = adfGeoTransform[0];
		dptr[1] = adfGeoTransform[1];
		dptr[2] = adfGeoTransform[2];
		dptr[3] = adfGeoTransform[3];
		dptr[4] = adfGeoTransform[4];
		dptr[5] = adfGeoTransform[5];
		mxSetField ( metadata_struct, 0, "GeoTransform", mxGeoTransform );
	}

	/* Get driver information */
	hDriver = GDALGetDatasetDriver( hDataset );

	mxGDALDriverShortName = mxCreateString ( GDALGetDriverShortName( hDriver ) );
	mxSetField ( metadata_struct, 0, (const char *) "DriverShortName", mxGDALDriverShortName );

	mxGDALDriverLongName = mxCreateString ( GDALGetDriverLongName( hDriver ) );
	mxSetField ( metadata_struct, 0, (const char *) "DriverLongName", mxGDALDriverLongName );

	xSize = GDALGetRasterXSize( hDataset );
	ySize = GDALGetRasterYSize( hDataset );
	mxGDALRasterXSize = mxCreateDoubleScalar ( (double) xSize );
	mxSetField ( metadata_struct, 0, (const char *) "RasterXSize", mxGDALRasterXSize );

	mxGDALRasterYSize = mxCreateDoubleScalar ( (double) ySize );
	mxSetField ( metadata_struct, 0, (const char *) "RasterYSize", mxGDALRasterYSize );

	raster_count = GDALGetRasterCount( hDataset );
	mxGDALRasterCount = mxCreateDoubleScalar ( (double)raster_count );
	mxSetField ( metadata_struct, 0, (const char *) "RasterCount", mxGDALRasterCount );

	/* Get the metadata for each band. */
	num_band_fields = 5;
	band_fieldnames[0] = strdup ( "XSize" );
	band_fieldnames[1] = strdup ( "YSize" );
	band_fieldnames[2] = strdup ( "Overview" );
	band_fieldnames[3] = strdup ( "NoDataValue" );
	band_fieldnames[4] = strdup ( "DataType" );
	band_struct = mxCreateStructMatrix ( raster_count, 1, num_band_fields, (const char **)band_fieldnames );

	num_overview_fields = 2;
	overview_fieldnames[0] = strdup ( "XSize" );
	overview_fieldnames[1] = strdup ( "YSize" );
	for ( band_number = 1; band_number <= raster_count; ++band_number ) {	/* Loop over bands */

		hBand = GDALGetRasterBand( hDataset, band_number );

		mxtmp = mxCreateDoubleScalar ( (double) GDALGetRasterBandXSize( hBand ) );
		mxSetField ( band_struct, 0, "XSize", mxtmp );

		mxtmp = mxCreateDoubleScalar ( (double) GDALGetRasterBandYSize( hBand ) );
		mxSetField ( band_struct, 0, "YSize", mxtmp );
	
		gdal_type = GDALGetRasterDataType ( hBand );
		mxtmp = mxCreateString ( GDALGetDataTypeName ( (GDALDataType)gdal_type ) );
		mxSetField ( band_struct, 0, (const char *) "DataType", mxtmp );

		tmpdble = GDALGetRasterNoDataValue ( hBand, &status );
		mxtmp = mxCreateDoubleScalar ( (double) (GDALGetRasterNoDataValue ( hBand, &status ) ) );
		mxSetField ( band_struct, 0, "NoDataValue", mxtmp );
	
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
			mxSetField ( band_struct, 0, "Overview", overview_struct );
		}

	}
	mxSetField ( metadata_struct, 0, "Band", band_struct );

	/* Record the GMT header. This will be interleaved with "corners" because they share somes values */
	mxGMT_header = mxCreateNumericMatrix(1, 9, mxDOUBLE_CLASS, mxREAL);
	dptr2 = mxGetPr(mxGMT_header);

	/* Record corners. */
	num_corner_fields = 4;
	corner_fieldnames[0] = strdup ("LL");
	corner_fieldnames[1] = strdup ("UL");
	corner_fieldnames[2] = strdup ("UR");
	corner_fieldnames[3] = strdup ("LR");

	corner_struct = mxCreateStructMatrix(1, 1, num_corner_fields, (const char **)corner_fieldnames);
	dims[0] = 1;	 dims[1] = 2;

	ReportCorner(hDataset, 0.0, GDALGetRasterYSize(hDataset), xy_c);	/* Lower Left */
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	dptr2[0] = xy_c[0];		dptr2[2] = xy_c[1];	/* xmin, ymin */
	mxSetField(corner_struct, 0, "LL", mxCorners );

	ReportCorner(hDataset, 0.0, 0.0, xy_c);					/* Upper Left */
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	mxSetField(corner_struct, 0, "UL", mxCorners );

	ReportCorner(hDataset, GDALGetRasterXSize(hDataset), 0.0, xy_c);	/* Upper Rigt */
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	dptr2[1] = xy_c[0];		dptr2[3] = xy_c[1];	/* xmax, ymax */
	mxSetField(corner_struct, 0, "UR", mxCorners );

	ReportCorner(hDataset, GDALGetRasterXSize(hDataset), GDALGetRasterYSize(hDataset), xy_c);	/* Lower Rigt */
	mxCorners = mxCreateNumericArray (2,dims,mxDOUBLE_CLASS, mxREAL);
	dptr = mxGetPr(mxCorners);
	dptr[0] = xy_c[0];		dptr[1] = xy_c[1];
	mxSetField(corner_struct, 0, "LR", mxCorners );

	mxSetField (metadata_struct, 0, "Corners", corner_struct);

	/* Fill in the rest of the GMT header values */
	if (z_min == 1e50) {		/* We don't know yet the dataset Min/Max */
		adfMinMax[0] = GDALGetRasterMinimum( hBand, &bGotMin );
		adfMinMax[1] = GDALGetRasterMaximum( hBand, &bGotMax );
		if(!(bGotMin && bGotMax))
			GDALComputeRasterMinMax( hBand, TRUE, adfMinMax );
		dptr2[4] = adfMinMax[0];
		dptr2[5] = adfMinMax[1];
	}
	else {
		dptr2[4] = z_min;
		dptr2[5] = z_max;
	}	

	dptr2[6] = 0;
	dptr2[7] = adfGeoTransform[1];
	dptr2[8] = fabs(adfGeoTransform[5]);

	if (correct_bounds) {
		dptr2[0] += dptr2[7] / 2;	dptr2[1] -= dptr2[7] / 2;
		dptr2[2] += dptr2[8] / 2;	dptr2[3] -= dptr2[8] / 2;
	}

	mxSetField (metadata_struct, 0, "GMT_hdr", mxGMT_header);

	return ( metadata_struct );
}

/************************************************************************/
/*                        ReportCorner()                                */
/************************************************************************/

int ReportCorner(GDALDatasetH hDataset, double x, double y, double *xy_c) {
    double	dfGeoX, dfGeoY;
    const char  *pszProjection;
    double	adfGeoTransform[6];

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

    return TRUE;
}

void DEBUGA(int n) {
#if debug
	mexPrintf("Merda %d\n",n);
#endif
}

/* ---------------------------------------------------------------------- */
int getNK(const mxArray *p, int which) {
	/* Get number of columns or number of bands of a mxArray */
	int nx, nBands, nDims;
	const int *dim_array;

	nDims     = mxGetNumberOfDimensions(p);
	dim_array = mxGetDimensions(p);
	nx = dim_array[1];
	nBands = dim_array[2];
	if (nDims == 2) 	/* Otherwise it would stay undefined */
		nBands = 1;

	if (which == 1)
		return(nx);
	else if (which == 2)
		return(nBands);
	else
		mexErrMsgTxt("getNK: Bad dimension number!");
	return (-1);
}
