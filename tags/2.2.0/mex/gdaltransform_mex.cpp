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

/* Program:	gdaltranslate_mex.c
 * Purpose:	matlab callable routine to reproject a list of coordinates into 
 *		any supported projection,including GCP-based transformations.
 * NOTE:	This is a crude first version that works when reprojecting
 *		with GCPs. The other features were not even tested.
 *
 * Revision 1.0  08/12/2008 Joaquim Luis
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

void DEBUGA(int n);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	char	*pszSRS_WKT = NULL;
	GDALDatasetH	hSrcDS, hDstDS;
	GDALDriverH	hDriver;
	GDALRasterBandH hBand;
	OGRSpatialReference oSrcSRS, oDstSRS; 
	GDALResampleAlg	interpMethod = GRA_NearestNeighbour;
	GDALTransformerFunc pfnTransformer = NULL;
	CPLErr		eErr;
	GDAL_GCP	*pasGCPs = NULL;
	void		*hTransformArg;
	static int runed_once = FALSE;	/* It will be set to true if reaches end of main */

	int	nx, ny, i, j, m, n, c = 0, n_pts, n_fields;
	char	*pszSrcSRS = NULL, *pszSrcWKT = NULL;
	char	*pszDstSRS = NULL, *pszDstWKT = NULL;
	double	*in_data, *ptr_d, *x, *y, *z;
	mxArray	*mx_ptr;

	int	nGCPCount = 0, nOrder = 0;
	char	**papszMetadataOptions = NULL;
	char	*tmp, *txt;
	int	bInverse = FALSE;
	int	*bSuccess;
	char	**papszTO = NULL;


	if (nrhs == 2 && mxIsStruct(prhs[1])) {

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
			if (nOrder != -1 || nOrder != 0 || nOrder != 1 ||
				nOrder != 2 || nOrder != 3)
				nOrder = 0;
		}
		/* -------------------------------------------------- */
		mx_ptr = mxGetField(prhs[1], 0, "inverse");
		if (mx_ptr != NULL)
			bInverse = TRUE;

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

	}
	else {
		mexPrintf("Usage: out = gdaltransform_mex(IN,PAR_STRUCT)\n\n");
		mexPrintf("      IN is a Mx2 or Mx3 array of doubles with X, Y (,Z)\n");
		mexPrintf("      PAR_STRUCT is a structure with at most two of the next fields:\n");
		mexPrintf("\t\t'SrcProjSRS', 'SrcProjWKT' -> Source projection string\n");
		mexPrintf("\t\t'DstProjSRS', 'DstProjWKT' -> Target projection string\n");
		mexPrintf("\t\t\tSRS stands for a string of the type used by proj4\n");
		mexPrintf("\t\t\tWKT stands for a string on the 'Well Known Text' format\n\n");
		mexPrintf("\t\t\tIf one of the Src or Dst fields is absent a GEOGRAPHIC WGS84 is assumed\n");
		mexPrintf("\nOPTIONS\n");
		mexPrintf("\t\t'gcp' a [Mx4] array with Ground Control Points\n");
		mexPrintf("\t\t'inverse' reverse transformation. e.g from georef to rows-columns\n");

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


	if (!runed_once) 		/* Do next call only at first time this MEX is loaded */
		GDALAllRegister();


	/* ----- Check that first argument contains at least a mx2 table */
	n_pts = mxGetM (prhs[0]);
	n_fields = mxGetN(prhs[0]);
	if (!mxIsNumeric(prhs[0]) || (n_fields < 2)) {
		mexPrintf("GDALTRANSFORM ERROR: first argument must contain a mx2 (or mx3) table\n");
		mexErrMsgTxt("               with the x,y (,z) positions to convert.\n");
	}

DEBUGA(1);
	/* ---------- Set the Source projection ---------------------------- */
	/* If it was not provided assume it is Geog WGS84 */
	if (pszSrcSRS == NULL && pszSrcWKT == NULL)
		oSrcSRS.SetWellKnownGeogCS( "WGS84" ); 
	else if (pszSrcWKT != NULL)
		oSrcSRS.importFromWkt( &pszSrcWKT );

	else {
		if( oSrcSRS.SetFromUserInput( pszSrcSRS ) != OGRERR_NONE )
			mexErrMsgTxt("GDALTRANSFORM: Translating source SRS failed.");
	}
	if (pszSrcWKT == NULL)
		oSrcSRS.exportToWkt( &pszSrcWKT );

	/* ------------------------------------------------------------------ */


DEBUGA(2);
	/* ---------- Set up the Target coordinate system ------------------- */
	/* If it was not provided assume it is Geog WGS84 */
	CPLErrorReset();
	if (pszDstSRS == NULL && pszDstWKT == NULL)
		oDstSRS.SetWellKnownGeogCS( "WGS84" ); 
	else if (pszDstWKT != NULL)
		oDstSRS.importFromWkt( &pszDstWKT );
	else {
		if( oDstSRS.SetFromUserInput( pszDstSRS ) != OGRERR_NONE )
			mexErrMsgTxt("GDALTRANSFORM: Translating target SRS failed.");
	}
	if (pszDstWKT == NULL)
		oDstSRS.exportToWkt( &pszDstWKT );
	/* ------------------------------------------------------------------ */



	/* -------------------------------------------------------------------- */
	/*      Create a transformation object from the source to               */
	/*      destination coordinate system.                                  */
	/* -------------------------------------------------------------------- */
	if( nGCPCount != 0 && nOrder == -1 ) {
		pfnTransformer = GDALTPSTransform;
		hTransformArg = GDALCreateTPSTransformer( nGCPCount, pasGCPs, FALSE );
	}
	else if( nGCPCount != 0 ) {
DEBUGA(3);
		pfnTransformer = GDALGCPTransform;
		hTransformArg = GDALCreateGCPTransformer( nGCPCount, pasGCPs, nOrder, FALSE );
	}
	else {
		pfnTransformer = GDALGenImgProjTransform;
		//hTransformArg = GDALCreateGenImgProjTransformer2( hSrcDS, hDstDS, papszTO );
		hTransformArg = 
    			GDALCreateGenImgProjTransformer( NULL, pszSrcWKT, NULL, pszDstWKT, 
						 nGCPCount == 0 ? FALSE : TRUE, 0, nOrder );
	}

	//CSLDestroy( papszTO );
DEBUGA(4);

	if( hTransformArg == NULL )
		mexErrMsgTxt("GDALTRANSFORM: Generating transformer failed.");

	in_data = (double *)mxGetData(prhs[0]);

	if (n_pts == 1) {
		x = &in_data[0];
		y = &in_data[1];
		bSuccess = &c;
	}
	else {
		x = (double *)mxCalloc(n_pts, sizeof(double));
		y = (double *)mxCalloc(n_pts, sizeof(double));
		bSuccess = (int *)mxCalloc(n_pts, sizeof(int));

		for (j = 0; j < n_pts; j++) {
			x[j] = in_data[j];
			y[j] = in_data[j+n_pts];
		}
	}
DEBUGA(5);

	z = (double *)mxCalloc(n_pts, sizeof(double));
	if (n_fields == 3)
		for (j = 0; j < n_pts; j++)
			z[j] = in_data[j+2*n_pts];

	if ( !pfnTransformer( hTransformArg, bInverse, n_pts, x, y, z, bSuccess ) )
		mexPrintf( "Transformation failed.\n" );


DEBUGA(6);
	if( nGCPCount != 0 && nOrder == -1 )
		GDALDestroyTPSTransformer(hTransformArg);
	else if( nGCPCount != 0 )
		GDALDestroyGCPTransformer(hTransformArg);
	else
		GDALDestroyGenImgProjTransformer(hTransformArg);


	/* -------------- Copy the result into plhs  --------------------------------- */
	plhs[0] = mxCreateDoubleMatrix (n_pts,n_fields, mxREAL);
	ptr_d = mxGetPr(plhs[0]);

	for (j = 0; j < n_pts; j++) {
		ptr_d[j] = x[j];
		ptr_d[j+n_pts] = y[j];
	}

	if (n_pts > 1) {
		mxFree((void *)x);	mxFree((void *)y);	mxFree((void *)bSuccess);
	}

	if (n_fields == 3) {
		for (j = 0; j < n_pts; j++)
			ptr_d[j+2*n_pts] = z[j];
	}
	mxFree((void *)z);


	if (pszDstWKT && strlen(pszDstWKT) > 1 ) OGRFree(pszDstWKT);	
	if (pszSrcWKT && strlen(pszSrcWKT) > 1 ) OGRFree(pszSrcWKT);
	//OGRFree(pszSrcWKT);
	//OGRFree(pszDstWKT);
	if (nGCPCount > 0) {
		GDALDeinitGCPs( nGCPCount, pasGCPs );	// makes this mex crash in the next call
		mxFree((void *) pasGCPs );
	}

	runed_once = TRUE;	/* Signals that next call won't need to call GDALAllRegister() again */

}

void DEBUGA(int n) {
#if debug
	mexPrintf("Merda %d\n",n);
#endif
}
