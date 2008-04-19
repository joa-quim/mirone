/*--------------------------------------------------------------------
 *	$Id: grdtrack.c,v 1.19 2004/07/16 17:07:02 pwessel Exp $
 *
 *	Copyright (c) 1991-2004 by P. Wessel and W. H. F. Smith
 *	See COPYING file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * grdtrack reads a xyfile, opens the 2d binary gridded grdfile, 
 * and samples the dataset at the xy positions with a bilinear or bicubic
 * interpolant.  This new data is added to the input as an extra column
 * and printed to standard output.  In order to evaluate derivatives along
 * the edges of the grdfile region, we assume natural bicubic spline
 * boundary conditions (d2z/dn2 = 0, n being the normal to the edge;
 * d2z/dxdy = 0 in the corners).  Rectangles of size x_inc by y_inc are 
 * mapped to [0,1] x [0,1] by affine transformation, and the interpolation
 * done on the normalized rectangle.
 *
 * Author:	Walter H F Smith
 * Date:	23-SEP-1993
 * 
 * Based on the original grdtrack, which had this authorship/date/history:
 *
 * Author:	Paul Wessel
 * Date:	29-JUN-1988
 * Revised:	5-JAN-1990	PW: Updated to v.2.0
 *		4-AUG-1993	PW: Added -Q
 *		14-AUG-1998	PW: GMT 3.1
 *  Modified:	10 Jul 2000 3.3.5  by PW to allow plain -L to indicate geographic coordinates
 * Version:	4
 *
 *--------------------------------------------------------------------*/
/*
 * Mexified version of grdtrack
 * Author:	J. Luis
 * Date: 	17 Aug 2004
 *
 * Usage
 * out = grdtrack_m(Zin,head,xy,'options');
 *
 * where	Z is the array containing the input file
 *		head is the header descriptor in the format of grdread/grdwrite mexs
 *		xy is an array with at least two columns (first column x and second y)
 *		and options may be for example: '-Z'
 *
 * IMPORTANT NOTE. The data type of Zin can be any of: double, single, Int32, Int16,
 * Uint16 or Uint8 but "out" will allways be of type double
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		03/03/08 J Luis, Adapted for 4.2.1
 */

#include "gmt.h"
#include "mex.h"

BOOLEAN GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, nx, ny, n_read = 0, n_points = 0, one_or_zero;
	int n_output = 0, n_fields, n_pts, ii, jj;
	
	BOOLEAN error = FALSE, suppress = FALSE, node = FALSE, z_only = FALSE;
	BOOLEAN is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	BOOLEAN is_uint16 = FALSE, is_uint8 = FALSE;
	
	double value, west, east, south, north, threshold = 1.0, i_dx, i_dy, half, *in, *out;
	float *f;

	int	i2, argc = 0, nc_h, nr_h, mx, n_arg_no_char = 0, *i_4, interpolant = BCR_BICUBIC;
	short int *i_2;
	unsigned short int *ui_2;
	char	**argv;
	unsigned char *ui_1;
	float	*z_4;
	double	*pdata_d, *z_8, *head;
	
	struct GRD_HEADER grd;
	struct GMT_EDGEINFO edgeinfo;
	struct GMT_BCR bcr;

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
	argv[0] = "grdtrack";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	west = east = south = north = 0.0;
	
	if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);

	GMT_boundcond_init (&edgeinfo);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				/* Common parameters */
			
				case 'R':
				case 'V':
				case 'f':
				case '\0':
					error += GMT_get_common_args (argv[i], &west, &east, &south, &north);
					break;

				/* Supplemental parameters */				
				case 'L':
					if (argv[i][2]) {
						error += GMT_boundcond_parse (&edgeinfo, &argv[i][2]);
						if (edgeinfo.gn) {
							GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
							GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
						}
					}
					else {
						GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
						GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
						mexPrintf ("%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
					}
					break;
				case 'M':
					GMT_multisegment (&argv[i][2]);
					break;
				case 'N':
					node = TRUE;
					break;
				case 'Q':
					interpolant = BCR_BILINEAR;
					threshold = (argv[i][2]) ? atof (&argv[i][2]) : 1.0;
					break;
				case 'S':
					suppress = TRUE;
					break;
				case 'Z':
					z_only = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}
	
	if (argc == 1 || error) {
		mexPrintf ("grdtrack %s - Sampling of a 2-D gridded netCDF grdfile along 1-D trackline\n\n", GMT_VERSION);
		mexPrintf ("usage: out = grdtrack_m(grd,head,xyfile, ['-L<flag>'], ['-M'], ['-N']\n"); 
		mexPrintf ("\t['-Q[<value>]'], ['-R<west/east/south/north>[r]'] ['-S'] ['-Z'] ['-f[i|o]<colinfo>']\n");
				
		mexPrintf ("\t<xyfile> is an multicolumn file with (lon,lat) in the first two columns\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-L sets boundary conditions.  <flag> can be either\n");
		mexPrintf ("\t   g for geographic boundary conditions\n");
		mexPrintf ("\t   or one or both of\n");
		mexPrintf ("\t   x for periodic boundary conditions on x\n");
		mexPrintf ("\t   y for periodic boundary conditions on y\n");
		mexPrintf ("\t-M Input file contain multiple segments separated by a record with NaNs\n");
		mexPrintf ("\t-N Report value at nearest node instead of interpolating\n");
		mexPrintf ("\t-Q Quick mode, use bilinear rather than bicubic interpolation.\n");
		mexPrintf ("\t   Optionally, append <value> in the 0 < value <= 1 range.\n");
		mexPrintf ("\t   [Default = 1 requires all 4 nodes to be non-NaN.], <value> = 0.5\n");
		mexPrintf ("\t   will interpolate about 1/2 way from a non-NaN to a NaN node, while\n");
		mexPrintf ("\t   0.1 will go about 90%% of the way, etc.\n");
		mexPrintf ("\t-R specifies a subregion [Default is old region]\n");
		mexPrintf ("\t-S Suppress output when result equals NaN\n");
		mexPrintf ("\t-Z only output z-values [Default gives all columns]\n");
		mexPrintf ("\t-f Special formatting of input/output columns\n");
		return;
	}

	if (threshold <= 0.0 || threshold > 1.0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -Q:  threshold must be in <0,1] range\n", "grdtrack");
		error++;
	}
	if (error) return;

	if (nlhs == 0) {
		mexPrintf("ERROR: Must provide an output.\n");
		return;
	}

	/* Find out in which data type was given the input array */
	if (mxIsDouble(prhs[0])) {
		z_8 = mxGetPr(prhs[0]);
		is_double = TRUE;
	}
	else if (mxIsSingle(prhs[0])) {
		z_4 = mxGetData(prhs[0]);
		is_single = TRUE;
	}
	else if (mxIsInt32(prhs[0])) {
		i_4 = mxGetData(prhs[0]);
		is_int32 = TRUE;
	}
	else if (mxIsInt16(prhs[0])) {
		i_2 = mxGetData(prhs[0]);
		is_int16 = TRUE;
	}
	else if (mxIsUint16(prhs[0])) {
		ui_2 = mxGetData(prhs[0]);
		is_uint16 = TRUE;
	}
	else if (mxIsUint8(prhs[0])) {
		ui_1 = mxGetData(prhs[0]);
		is_uint8 = TRUE;
	}
	else {
		mexPrintf("GRDTRACK ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, In32, In16, UInt16 and Uint8.\n");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2)
		mexErrMsgTxt("First argument must contain a decent array\n");

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("Second argument must contain a valid header of the input array.\n");
	
	head  = mxGetPr(prhs[1]);		/* Get header info */

	/* Check that thirth argument contains at least a mx2 table */
	n_pts = mxGetM (prhs[2]);
	n_fields = mxGetN(prhs[2]);
	if (!mxIsNumeric(prhs[2]) || (n_fields < 2))
		mexErrMsgTxt("GRDTRACK ERROR: thirth argument must contain the x,y positions where to interpolate.\n");
	if (z_only)
		n_fields = 0;

	/* Read the interpolation points and convert them to double */
	if (mxIsDouble(prhs[2]))
		in = mxGetPr(prhs[2]);
	else if (mxIsSingle(prhs[2]))
		in = mxGetData(prhs[2]);

	grd.x_min = head[0];	grd.x_max = head[1];
	grd.y_min = head[2];	grd.y_max = head[3];
	grd.z_min = head[4];	grd.z_max = head[5];
	grd.x_inc = head[7];	grd.y_inc = head[8];
	grd.nx = nx;			grd.ny = ny;
	grd.node_offset = irint(head[6]);
	mx = nx + 4;

	if (west == east) {	/* No subset asked for */
		west = grd.x_min;
		east = grd.x_max;
		south = grd.y_min;
		north = grd.y_max;
	}
	one_or_zero = (grd.node_offset) ? 0 : 1;
	half = (grd.node_offset) ? 0.5 : 0.0;
	nx = irint ( (east - west) / grd.x_inc) + one_or_zero;
	ny = irint ( (north - south) / grd.y_inc) + one_or_zero;
	i_dx = 1.0 / grd.x_inc;
	i_dy = 1.0 / grd.y_inc;

	f = mxCalloc ((nx+4)*(ny+4), sizeof (float));
	GMT_end_for_mex (argc, argv);

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) f[i2*mx + j + 2*mx + 2] = (float)z_8[j*ny+i];
	}
	else if (is_single) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) f[i2*mx + j + 2*mx + 2] = z_4[j*ny+i];
	}
	else if (is_int32) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) f[i2*mx + j + 2*mx + 2] = (float)i_4[j*ny+i];
	}
	else if (is_int16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) f[i2*mx + j + 2*mx + 2] = (float)i_2[j*ny+i];
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) f[i2*mx + j + 2*mx + 2] = (float)ui_2[j*ny+i];
	}
	else if (is_uint8) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) f[i2*mx + j + 2*mx + 2] = (float)ui_1[j*ny+i];
	}

	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;

	GMT_boundcond_param_prep (&grd, &edgeinfo);
	
	project_info.w = west;	project_info.e = east;
	project_info.s = south;	project_info.n = north;
	
	/* Initialize bcr structure:  */

	/*GMT_bcr_init (&grd, GMT_pad, bilinear, threshold, &bcr);*/
	GMT_bcr_init (&grd, GMT_pad, interpolant, threshold, &bcr);

	/* Set boundary conditions  */
	
	GMT_boundcond_set (&grd, &edgeinfo, GMT_pad, f);
	
	if ((out = mxCalloc(n_pts * (n_fields+1), sizeof (double))) == 0)
		mexErrMsgTxt("GRDTRACK ERROR: Could not allocate memory\n");

	for (i = 0; i < n_pts; i++) {
		while ( (GMT_is_dnan(in[i]) || GMT_is_dnan(in[i+n_pts])) && !z_only) {
			for (j = 0; j < n_fields; j++) out[j*n_pts+i] = in[j*n_pts+i];
			out[j*n_pts+i] = mxGetNaN();
			i++;
		}

		/* If point is outside grd area, shift it using periodicity or skip if not periodic. */
		while ( (in[i+n_pts] < grd.y_min) && (edgeinfo.nyp > 0) ) in[i+n_pts] += (grd.y_inc * edgeinfo.nyp);
		if (in[i+n_pts] < grd.y_min) continue;

		while ( (in[i+n_pts] > grd.y_max) && (edgeinfo.nyp > 0) ) in[i+n_pts] -= (grd.y_inc * edgeinfo.nyp);
		if (in[i+n_pts] > grd.y_max) continue;

		while ( (in[i] < grd.x_min) && (edgeinfo.nxp > 0) ) in[i] += (grd.x_inc * edgeinfo.nxp);
		if (in[i] < grd.x_min) continue;

		while ( (in[i] > grd.x_max) && (edgeinfo.nxp > 0) ) in[i] -= (grd.x_inc * edgeinfo.nxp);
		if (in[i] > grd.x_max) continue;
		
		if (node) {
			ii = irint ((in[i] - grd.x_min) * i_dx - half) + one_or_zero;
			jj = irint ((grd.y_max - in[i+n_pts]) * i_dy - half) + one_or_zero;
			value = f[(jj+GMT_pad[3])*mx+ii+GMT_pad[0]];
		}
		else
			value = GMT_get_bcr_z(&grd, in[i], in[i+n_pts], f, &edgeinfo, &bcr);

		if (suppress && GMT_is_dnan (value)) continue;

		if (z_only) {	/* Simply print out value */
			out[i] = value;
		}
		else {	/* Simply copy other columns, append value, and output */
			for (j = 0; j < n_fields; j++) out[j*n_pts+i] = in[j*n_pts+i];
			out[j*n_pts+i] = value;
		}
	}

	plhs[0] = mxCreateDoubleMatrix (n_pts,n_fields+1, mxREAL);
	pdata_d = mxGetPr(plhs[0]);
	memcpy(pdata_d, out, n_pts*(n_fields+1)*8);
	mxFree(out);
}
