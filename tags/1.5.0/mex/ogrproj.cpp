/*
 *      Coffeeright (c) 2002-2003 by J. Luis
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: w3.ualg.pt/~jluis/m_gmt
 *--------------------------------------------------------------------*/
/* Program:	ogrprof.c
 * Purpose:	matlab callable routine to do vector data projection via gdal
 *
 * Revision 1.0  24/6/2006 Joaquim Luis
 *
 */

#define CNULL	((char *)NULL)
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#include "mex.h"
#include "gdal.h"
#include "ogr_spatialref.h"

char *SanitizeSRS( const char *pszUserInput );
void Usage();

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	i, j, n_pts, n_fields;
	double	*in_data, *ptr_d, *x, *y, *z;
	char	*pszSrcSRS = NULL, *pszSrcWKT = NULL;
	char	*pszDstSRS = NULL, *pszDstWKT = NULL;
	mxArray	*mx_ptr;
	OGRSpatialReference oSrcSRS, oDstSRS; 
	OGRCoordinateTransformation *poCT; 

	if (nrhs == 0) { Usage(); return; }

	/* Case of just translate an SRS string into a WKT form */
	if (nrhs == 1 && mxIsChar(prhs[0])) {
		pszSrcSRS = (char *)mxArrayToString(prhs[0]);

		//pszSrcWKT = SanitizeSRS(pszSrcSRS);

		if( oSrcSRS.SetFromUserInput( pszSrcSRS ) != OGRERR_NONE )
			mexErrMsgTxt("OGRPROJ: Translating SRS string failed.");

		if (pszSrcSRS[0] == '+')	// from Proj4 to WKT
			oSrcSRS.exportToPrettyWkt( &pszSrcWKT, FALSE );
		else				// from others to Proj4
			oSrcSRS.exportToProj4( &pszSrcWKT );

		if (nlhs == 1)
			plhs[0] = mxCreateString(pszSrcWKT);
		else
			mexPrintf("%s",pszSrcWKT);

		OGRFree(pszSrcWKT);
		return;
	}

	if (nrhs == 2 && mxIsStruct(prhs[1])) {
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
	}
	else if (nrhs == 2 && mxIsChar(prhs[1]))
		pszSrcSRS = (char *)mxArrayToString(prhs[1]);
	else
		mexErrMsgTxt("OGRPROJ: Wrong number/type of arguments.");

	/* Check that first argument contains at least a mx2 table */
	n_pts = mxGetM (prhs[0]);
	n_fields = mxGetN(prhs[0]);
	if (!mxIsNumeric(prhs[0]) || (n_fields < 2)) {
		mexPrintf("OGRPROJ ERROR: first argument must contain a mx2 (or mx3) table\n");
		mexErrMsgTxt("               with the x,y (,z) positions to convert.\n");
	}

	/* ---------- Set the Source projection ---------------------------- */
	/* If it was not provided assume it is Geog WGS84 */
	if (pszSrcSRS == NULL && pszSrcWKT == NULL)
		oSrcSRS.SetWellKnownGeogCS( "WGS84" ); 
	else if (pszSrcWKT != NULL)
		oSrcSRS.importFromWkt( &pszSrcWKT );
	else {
		/*if ((pszSrcWKT = SanitizeSRS(pszSrcSRS)) == CNULL)
			mexErrMsgTxt("OGRPROJ: Translating source SRS failed.");
		oSrcSRS.importFromWkt( &pszSrcWKT );*/

		if( oSrcSRS.SetFromUserInput( pszSrcSRS ) != OGRERR_NONE )
			mexErrMsgTxt("OGRPROJ: Translating source SRS failed.");
	}
	/* ------------------------------------------------------------------ */

	/* ---------- Set the Target projection ---------------------------- */
	/* If it was not provided assume it is Geog WGS84 */
	CPLErrorReset();
	if (pszDstSRS == NULL && pszDstWKT == NULL)
		oDstSRS.SetWellKnownGeogCS( "WGS84" ); 
	else if (pszDstWKT != NULL)
		oDstSRS.importFromWkt( &pszDstWKT );
	else {
		//oDstSRS = new OGRSpatialReference();
		if( oDstSRS.SetFromUserInput( pszDstSRS ) != OGRERR_NONE )
			mexErrMsgTxt("OGRPROJ: Translating target SRS failed.");
	}
	/* ------------------------------------------------------------------ */

	poCT = OGRCreateCoordinateTransformation( &oSrcSRS, &oDstSRS );
	if( poCT == NULL ) {
		mexPrintf("Failed to create coordinate transformation between the\n"
			"following coordinate systems.  This may be because they\n"
			"are not transformable, or because projection services\n"
			"(PROJ.4 DLL/.so) could not be loaded.\n" );
		oSrcSRS.exportToPrettyWkt( &pszSrcWKT, FALSE );
		oDstSRS.exportToPrettyWkt( &pszDstWKT, FALSE );
		mexPrintf( "Source:\n%s\n%s\n", pszSrcWKT, pszDstWKT );
		OGRFree(pszSrcWKT);
		OGRFree(pszDstWKT);
		mexErrMsgTxt("");
	}

	in_data = (double *)mxGetData(prhs[0]);

	if (n_pts == 1) {
		x = &in_data[0];
		y = &in_data[1];
	}
	else {
		x = (double *)mxCalloc(n_pts, sizeof(double));
		y = (double *)mxCalloc(n_pts, sizeof(double));

		for (j = 0; j < n_pts; j++) {
			x[j] = in_data[j];
			y[j] = in_data[j+n_pts];
		}
	}

	if (n_fields == 3) {
		z = (double *)mxCalloc(n_pts, sizeof(double));
		for (j = 0; j < n_pts; j++)
			z[j] = in_data[j+2*n_pts];

		poCT->Transform( n_pts, x, y, z );
	}
	else {
		poCT->Transform( n_pts, x, y );
		//if( !poCT->Transform( n_pts, x, y ) )
			//mexPrintf( "Transformation failed.\n" );
	}

	delete (poCT);
	OGRFree(pszSrcWKT);
	OGRFree(pszDstWKT);

	/* -------------- Copy the result into plhs  --------------------------------- */
	plhs[0] = mxCreateDoubleMatrix (n_pts,n_fields, mxREAL);
	ptr_d = mxGetPr(plhs[0]);

	for (j = 0; j < n_pts; j++) {
		ptr_d[j] = x[j];
		ptr_d[j+n_pts] = y[j];
	}

	if (n_pts > 1) {
		mxFree(x);	mxFree(y);
	}

	if (n_fields == 3) {
		for (j = 0; j < n_pts; j++)
			ptr_d[j+2*n_pts] = z[j];

		mxFree(z);
	}
}

/************************************************************************/
/*                             SanitizeSRS                              */
/************************************************************************/

char *SanitizeSRS( const char *pszUserInput ) {
	OGRSpatialReferenceH hSRS = NULL;
	char *pszResult = NULL;

	CPLErrorReset();

	hSRS = OSRNewSpatialReference( NULL );

	if( OSRSetFromUserInput( hSRS, pszUserInput ) == OGRERR_NONE )
		OSRExportToWkt( hSRS, &pszResult );
	else
		return CNULL;
    
	OSRDestroySpatialReference( hSRS );

	return pszResult;
}

/* ------------------------------------------------------------------------- */
void Usage() {
	mexPrintf("out = ogrproj(IN,PAR_STRUCT)\n");
	mexPrintf("      IN is a Mx2 or Mx3 array of doubles with lon, lat (,z)\n");
	mexPrintf("      PAR_STRUCT is a structure with at most two of the next fields:\n");
	mexPrintf("      SrcProjSRS, SrcProjWKT -> Source projection string\n");
	mexPrintf("      DstProjSRS, DstProjWKT -> Target projection string\n");
	mexPrintf("      SRS stands for a string of the type used by proj4\n");
	mexPrintf("      WKT stands for a string on the 'Well Known Text' format\n\n");
	mexPrintf("      If one of the Src or Dst fields is absent a GEOGRAPHIC WGS84 is assumed\n");
	mexPrintf("\nout = ogrproj('SrcProjSRS')\n");
	mexPrintf("      converts the SRS Proj4 string into its WKT form,\n");
	mexPrintf("      or from others (see 'SetFromUserInput' method info) into a Proj4 string.\n");
}
