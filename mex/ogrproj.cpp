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

void Usage();

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	i, j, n_pts, n_fields, nc, nr;
	double	*in_data, *ptr_d, *x = NULL, *y = NULL, *z = NULL;
	char	*pszSrcSRS = NULL, *pszSrcWKT = NULL;
	char	*pszDstSRS = NULL, *pszDstWKT = NULL;
	char	*t;
	mxArray	*mx_ptr;
	OGRSpatialReference oSrcSRS, oDstSRS; 
	OGRCoordinateTransformation *poCT; 

	if (nrhs == 0 && nlhs == 0) { Usage(); return; }

	/* Case of just translate an SRS string into a WKT form */
	if (nrhs == 1 && mxIsChar(prhs[0])) {
		t = (char *)mxArrayToString(prhs[0]);
		pszSrcSRS = strdup(t);
		mxFree(t);

		if ( oSrcSRS.SetFromUserInput( pszSrcSRS ) != OGRERR_NONE )
			mexErrMsgTxt("OGRPROJ: Translating SRS string failed.");

		if (pszSrcSRS[0] == '+')	/* from Proj4 to WKT */
			oSrcSRS.exportToPrettyWkt( &pszSrcWKT, FALSE );
		else				/* from others to Proj4 */
			oSrcSRS.exportToProj4( &pszSrcWKT );

		if (nlhs == 1)
			plhs[0] = mxCreateString(pszSrcWKT);
		else
			mexPrintf("%s",pszSrcWKT);

		/*OGRSpatialReference::DestroySpatialReference ( &oSrcSRS );*/
		OGRFree(pszSrcWKT);
		free((void *)pszSrcSRS);
		return;
	}

	if (!mxIsDouble(prhs[0])) 
		mexErrMsgTxt("Data points must be of type double (SORRY).\n");

	if (mxIsNumeric(prhs[0]) && mxIsNumeric(prhs[1])) {
		if (!mxIsDouble(prhs[1])) 
			mexErrMsgTxt("All data points must be of type double (SORRY).\n");
		if (mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1])) 
			mexErrMsgTxt("First two args (x and y) have not the same size.\n");
	}

	if (nrhs >= 2 && mxIsStruct(prhs[nrhs-1])) {
		/* Here we must be careful to not allow mixture between Matlab and GDAL pointers 
		   accessing. mxArrayToString are hold by ML, but they would later be attempt to
		   be free by GDAL. That is why we go through the "neutral" strdup */
		mx_ptr = mxGetField(prhs[nrhs-1], 0, "SrcProjSRS");
		if (mx_ptr != NULL) {
			t = (char *)mxArrayToString(mx_ptr);
			pszSrcSRS = strdup(t);
			mxFree(t);
		}

		mx_ptr = mxGetField(prhs[nrhs-1], 0, "SrcProjWKT");
		if (mx_ptr != NULL) {
			t = (char *)mxArrayToString(mx_ptr);
			pszSrcWKT = strdup(t);
			mxFree(t);
		}

		mx_ptr = mxGetField(prhs[nrhs-1], 0, "DstProjSRS");
		if (mx_ptr != NULL) {
			t = (char *)mxArrayToString(mx_ptr);
			pszDstSRS = strdup(t);
			mxFree(t);
		}

		mx_ptr = mxGetField(prhs[nrhs-1], 0, "DstProjWKT");
		if (mx_ptr != NULL) {
			t = (char *)mxArrayToString(mx_ptr);
			pszDstWKT = strdup(t);
			mxFree(t);
		}
	}
	else if (nrhs >= 2 && mxIsChar(prhs[nrhs-1]))
		pszSrcSRS = (char *)mxArrayToString(prhs[nrhs-1]);
	else
		mexErrMsgTxt("OGRPROJ: Wrong number/type of arguments.");

	/* Try to guess if input is row or column arrays and if it's Mx2|3 or a x, y */
	nr = mxGetM(prhs[0]);
	if (mxIsNumeric(prhs[0]) && mxIsNumeric(prhs[1]))
		nc = 2;
	else
		nc = mxGetN(prhs[0]);

	if (nr > 3) {				/* Easy case. Column vector */
		n_pts = nr;	n_fields = nc;
	}
	else if (nc > 3) {			/* Easy case. Row vector */
		n_pts = nc;	n_fields = nr;
	}
	else if ((nr == 2 && nc == 2) || (nr == 3 && nc == 3)) {	/* AMBIGUOUS case. Assume column vectors */
		n_pts = n_fields = nr;
	}
	else if (nr == 1 && (nc == 2 || nc == 3)) {	/* Easy case. One point only */
		n_pts = nr;	n_fields = nc;
	}
	else {
		mexPrintf("Case not foreseen in input array size guessing. Boing??\n");
		n_pts = nr;	n_fields = nc;
	}

	if (!mxIsNumeric(prhs[0]) || (n_fields < 1) || (n_fields > 3) || (n_pts == 0)) {
		mexPrintf("OGRPROJ ERROR: first argument must contain a Mx1, Mx2 or Mx3 table\n");
		mexErrMsgTxt("               with the x,y (,z) positions to convert.\n");
	}

	/* ---------- Set the Source projection ---------------------------- */
	/* If it was not provided assume it is Geog WGS84 */
	if (pszSrcSRS == NULL && pszSrcWKT == NULL)
		oSrcSRS.SetWellKnownGeogCS( "WGS84" ); 
	else if (pszSrcWKT != NULL)
		oSrcSRS.importFromWkt( &pszSrcWKT );
	else {
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
		if (pszSrcWKT && strlen(pszSrcWKT) > 1 ) free((void *)pszSrcWKT);
		if (pszDstWKT && strlen(pszDstWKT) > 1 ) free((void *)pszDstWKT);
		mexErrMsgTxt("");
	}

	if (nlhs) {			/* If not in-place conversion */
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
		else
			poCT->Transform( n_pts, x, y );
	}
	else {
		x = mxGetPr(prhs[0]);	y = mxGetPr(prhs[1]);
		poCT->Transform( n_pts, x, y );
	}

	OGRCoordinateTransformation::DestroyCT(poCT);
	/*OGRSpatialReference::DestroySpatialReference ( &oSrcSRS );
	OGRSpatialReference::DestroySpatialReference ( &oDstSRS );*/
	if (pszSrcWKT && strlen(pszSrcWKT) > 1 ) OGRFree(pszSrcWKT);
	if (pszDstWKT && strlen(pszDstWKT) > 1 ) OGRFree(pszDstWKT);
	if (pszSrcWKT && strlen(pszSrcWKT) > 1 ) free((void *)pszSrcWKT);
	if (pszDstWKT && strlen(pszDstWKT) > 1 ) free((void *)pszDstWKT);


	if (nlhs) {	/* --------- Copy the result into plhs  ------------- */
		plhs[0] = mxCreateDoubleMatrix (n_pts,n_fields, mxREAL);
		ptr_d = mxGetPr(plhs[0]);

		for (j = 0; j < n_pts; j++) {
			ptr_d[j] = x[j];
			ptr_d[j+n_pts] = y[j];
		}

		if (n_pts > 1) {
			mxFree((void *)x);	mxFree((void *)y);
		}

		if (n_fields == 3) {
			for (j = 0; j < n_pts; j++)
				ptr_d[j+2*n_pts] = z[j];

			mxFree((void *)z);
		}
	}
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
	mexPrintf("      If one of the Src or Dst fields is absent a GEOGRAPHIC WGS84 is assumed\\nn");

	mexPrintf("out = ogrproj(X,Y,PAR_STRUCT)\n");
	mexPrintf("      Same as above but X and Y are sent in different variables.\n\n");

	mexPrintf("      ogrproj(X,Y,PAR_STRUCT)\n");
	mexPrintf("      Same as above but conversion is done in place (no Z).\n\n");

	mexPrintf("\nout = ogrproj('SrcProjSRS')\n");
	mexPrintf("      converts the SRS Proj4 string into its WKT form,\n");
	mexPrintf("      or from others (see 'SetFromUserInput' method info) into a Proj4 string.\n");
}
