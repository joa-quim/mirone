/*--------------------------------------------------------------------
 *	$Id: grdsample.c,v 1.13 2004/07/16 17:07:02 pwessel Exp $
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
 * grdsample reads a grdfile and evaluates the grid at new grid positions
 * specified by new dx/dy values using a 2-D Taylor expansion of order 3.
 * In order to evaluate derivatives along the edges of the surface, I assume 
 * natural bi-cubic spline conditions, i.e. both the second and third normal 
 * derivatives are zero, and that the dxdy derivative in the corners are zero, too.
 *
 * Author:	Paul Wessel
 * Date:	19-JUL-1989
 * Revised:	6-JAN-1990	PW: Updated to v.2.0
 * Revised:	16-JUN-1998	PW: Updated to v.3.1
 * Version:	4
 *
 *--------------------------------------------------------------------*/
/*
 * Mexified version of grdsample
 * Author:	J. Luis
 * Date: 	17 Aug 2004
 *
 * Usage
 * Zout = grdsample(Zin,head,'options');
 *
 * where	Z is the array containing the input file
 *			head is the header descriptor in the format of grdread/grdwrite mexs
 *			and options may be for example: '-N100/100',-Q'
 *
 * IMPORTANT NOTE. The data type of Zin is preserved in Zout. That means you can send Zin
 * as a double, single, In32, In16, UInt16 or Uint8 and receive Zout in one of those types
 *
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */

#include "gmt.h"
#include "mex.h"

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, ij, one, i2, nx, ny, nc_h, nr_h, mx;
	int	argc = 0, n_arg_no_char = 0, *i_4, *o_i4, *pdata_i4;
	short int *i_2, *o_i2, *pdata_i2;
	unsigned short int *ui_2, *o_ui2, *pdata_ui2;
	char	**argv;
	unsigned char *ui_1, *o_ui1, *pdata_ui1;
	int error = FALSE, greenwich = FALSE, offset = FALSE, bilinear = FALSE;
	int area_set = FALSE, n_set = FALSE, inc_set = FALSE, toggle = FALSE;
	int is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int is_uint16 = FALSE, is_uint8 = FALSE;

	double *lon, lat, dx2, dy2, threshold = 1.0, *pdata_d, *z_8, *head, *o_d;
	float *a, *b, *z_4, *pdata_s, *o_s;
	struct GRD_HEADER grd_a, grd_b;
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
	argv[0] = "grdsample";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);*/
	argc = GMT_begin (argc, argv);

	GMT_grd_init (&grd_b, argc, argv, FALSE);
	GMT_boundcond_init (&edgeinfo);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
		
				/* Common parameters */
			
				case 'R':
				case 'f':
				case '\0':
					error += GMT_get_common_args (argv[i], &grd_b.x_min, &grd_b.x_max, &grd_b.y_min, &grd_b.y_max);
					break;
				
				/* Supplemental parameters */
			
				case 'F':
					offset = TRUE;
					break;
				case 'N':
					sscanf (&argv[i][2], "%d/%d", &grd_b.nx, &grd_b.ny);
					if (grd_b.ny == 0) grd_b.ny = grd_b.nx;
					n_set = TRUE;
					break;
				case 'I':
					GMT_getinc (&argv[i][2], &grd_b.x_inc, &grd_b.y_inc);
					inc_set = TRUE;
					break;
				case 'L':
					error += GMT_boundcond_parse (&edgeinfo, &argv[i][2]);
					break;
				case 'Q':
					bilinear = TRUE;
					threshold = (argv[i][2]) ? atof (&argv[i][2]) : 1.0;
					break;
				case 'T':	/* Convert from pixel file <-> gridfile */
					toggle = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
	}
	
	if (argc == 1 || error) {
		mexPrintf ("grdsample - Resample a gridded file onto a new grid\n\n");
		mexPrintf ("usage: new_grdfile = grdsample(old_grdfile,head, '[-F]', '[-I<dx>[m|c][/<dy>[m|c]]]', '[-L<flag>]',\n");
		mexPrintf ("\t'[-N<nx/ny>]', '[-Q[<value>]]', '[-R<west/east/south/north>]', '[-T]', '[-f[i|o]<colinfo>]')\n");
		
		mexPrintf ("\t<old_grdfile> is data set to be resampled\n");
		mexPrintf ("\t<new_grdfile> is the interpolated output data set\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-F force pixel node registration  [Default is centered]\n");
		mexPrintf ("\t-I sets the grid spacing (dx, dy) for the new grid\n");
		mexPrintf ("\t-L sets boundary conditions.  <flag> can be either\n");
		mexPrintf ("\t   g for geographic boundary conditions\n");
		mexPrintf ("\t   or one or both of\n");
		mexPrintf ("\t   x for periodic boundary conditions on x\n");
		mexPrintf ("\t   y for periodic boundary conditions on y\n");
		mexPrintf ("\t-N specifies number of columns (nx) and rows (ny) of new grid\n");
		mexPrintf ("\t-Q Quick mode, use bilinear rather than bicubic interpolation.\n");
		mexPrintf ("\t   Optionally, append <value> in the 0 < value <= 1 range.\n");
		mexPrintf ("\t   [Default = 1 requires all 4 nodes to be non-NaN.], <value> = 0.5\n");
		mexPrintf ("\t   will interpolate about 1/2 way from a non-NaN to a NaN node, while\n");
		mexPrintf ("\t   0.1 will go about 90%% of the way, etc.\n");
		mexPrintf ("\t-R specifies a subregion [Default is old region]\n");
		mexPrintf ("\t-T Toggles between grid registration and pixel registration\n");
		mexPrintf ("\t   In this case, -N is ignored.\n");
		mexPrintf ("\t   One only of -N -I must be specified\n");
		return;
	}
	if (bilinear && (threshold <= 0.0 || threshold > 1.0)) {
		mexPrintf ("%s: GMT SYNTAX ERROR -Q:  threshold must be in <0,1] range\n", "grdsample");
		error++;
	}
	if (!toggle) {
		if (grd_b.x_min != grd_b.x_max && grd_b.y_min != grd_b.y_max) area_set = TRUE;
		if (inc_set && n_set) {
			mexPrintf ("%s: GMT SYNTAX ERROR:  Only one of -I, -N may be specified\n", "grdsample");
			error++;
		}
		if (n_set && (grd_b.nx <= 0 || grd_b.ny <= 0)) {
			mexPrintf ("%s: GMT SYNTAX ERROR -N:  Must specify positive integers\n", "grdsample");
			error++;
		}
		if (inc_set && (grd_b.x_inc <= 0.0 || grd_b.y_inc <= 0.0)) {
			mexPrintf ("%s: GMT SYNTAX ERROR -I:  Must specify positive increments\n", "grdsample");
			error++;
		}
		if (!(inc_set || n_set)) {
			mexPrintf ("%s: GMT SYNTAX ERROR:  One of -I, -N must be specified\n", "grdsample");
			error++;
		}
	}

	if (error) return;

	if (nlhs == 0)
		mexErrMsgTxt("GRDSAMPLE ERROR: Must provide an output.\n");

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
		mexPrintf("GRDSAMPLE ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, In32, In16, UInt16 and Uint8.\n");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2)
		mexErrMsgTxt("First argument must contain a decent array\n");

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("Second argument must contain a valid header of the input array\n");

	head  = mxGetPr(prhs[1]);		/* Get header info */

	grd_a.x_min = head[0];	grd_a.x_max = head[1];
	grd_a.y_min = head[2];	grd_a.y_max = head[3];
	grd_a.z_min = head[4];	grd_a.z_max = head[5];
	grd_a.x_inc = head[7];	grd_a.y_inc = head[8];
	grd_a.nx = nx;			grd_a.ny = ny;
	grd_a.node_offset = irint(head[6]);
	mx = nx + 4;

	GMT_boundcond_param_prep (&grd_a, &edgeinfo);

	if (toggle) {
		offset = !grd_a.node_offset;	/* Change to the opposite of what it is */
		grd_b.nx = (offset) ? grd_a.nx - 1 : grd_a.nx + 1;
		grd_b.ny = (offset) ? grd_a.ny - 1 : grd_a.ny + 1;
		area_set = inc_set = FALSE;
	}
	
	if (area_set) {
		if ((grd_b.y_min + grd_b.y_inc) < grd_a.y_min || (grd_b.y_max - grd_b.y_inc) > grd_a.y_max)
			mexErrMsgTxt ("GRDSAMPLE:  Selected region exceeds the boundaries of the grdfile!\n");
		else if (!edgeinfo.nxp && ((grd_b.x_min + grd_b.x_inc) < grd_a.x_min || (grd_b.x_max - grd_b.x_inc) > grd_a.x_max))
			mexErrMsgTxt ("GRDSAMPLE:  Selected region exceeds the boundaries of the grdfile!\n");
	}
	
	a = mxCalloc ((nx+4)*(ny+4), sizeof (float));

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) a[i2*mx + j + 2*mx + 2] = (float)z_8[j*ny+i];
	}
	else if (is_single) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) a[i2*mx + j + 2*mx + 2] = z_4[j*ny+i];
	}
	else if (is_int32) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) a[i2*mx + j + 2*mx + 2] = (float)i_4[j*ny+i];
	}
	else if (is_int16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) a[i2*mx + j + 2*mx + 2] = (float)i_2[j*ny+i];
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) a[i2*mx + j + 2*mx + 2] = (float)ui_2[j*ny+i];
	}
	else if (is_uint8) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) a[i2*mx + j + 2*mx + 2] = (float)ui_1[j*ny+i];
	}

	if (!offset && !toggle) offset = grd_a.node_offset;
	one = !offset;
	grd_b.node_offset = offset;
	
	if (!area_set) {
		grd_b.x_min = grd_a.x_min;
		grd_b.x_max = grd_a.x_max;
		grd_b.y_min = grd_a.y_min;
		grd_b.y_max = grd_a.y_max;
	}
	
	if (edgeinfo.nxp && grd_b.x_min < 0.0 && grd_b.x_max > 0.0)
		greenwich = TRUE;
	else if (edgeinfo.nxp && grd_b.x_max > 360.0) {
		greenwich = TRUE;
		grd_b.x_min -= 360.0;
		grd_b.x_max -= 360.0;
	}
	if (inc_set) {
		grd_b.nx = irint ((grd_b.x_max - grd_b.x_min) / grd_b.x_inc) + one;
		grd_b.ny = irint ((grd_b.y_max - grd_b.y_min) / grd_b.y_inc) + one;
		grd_b.x_inc = (grd_b.x_max - grd_b.x_min) / (grd_b.nx - one);
		grd_b.y_inc = (grd_b.y_max - grd_b.y_min) / (grd_b.ny - one);
	}
	else {
		grd_b.x_inc = (grd_b.x_max - grd_b.x_min) / (grd_b.nx - one);
		grd_b.y_inc = (grd_b.y_max - grd_b.y_min) / (grd_b.ny - one);
	}
	
	/*GMT_grd_RI_verify (&grd_b, 1);*/	/* IF (IVAN == TRUE)  ==> Matlab = BOOM */

	b = mxCalloc (grd_b.nx * grd_b.ny, sizeof (float));
	
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;	/* Leave room for 2 empty boundary rows/cols */
	
	/* Initialize bcr structure:  */

	GMT_bcr_init (&grd_a, GMT_pad, bilinear, threshold, &bcr);

	/* Set boundary conditions  */
	
	GMT_boundcond_set (&grd_a, &edgeinfo, GMT_pad, a);

	/* Precalculate longitudes */
	
	dx2 = 0.5 * grd_b.x_inc;
	dy2 = 0.5 * grd_b.y_inc;
	lon = (double *) GMT_memory (VNULL, (size_t)grd_b.nx, sizeof (double), GMT_program);
	for (i = 0; i < grd_b.nx; i++) {
		lon[i] = grd_b.x_min + (i * grd_b.x_inc) + ((offset) ? dx2 : 0.0);
		if (edgeinfo.nxp && greenwich && lon[i] > 180.0) lon[i] -= 360.0;
	}

	for (j = ij = 0; j < grd_b.ny; j++) {
		lat = grd_b.y_max - (j * grd_b.y_inc);
		if (offset) lat -= dy2;
		for (i = 0; i < grd_b.nx; i++, ij++) b[ij] = (float)GMT_get_bcr_z (&grd_a, lon[i], lat, a, &edgeinfo, &bcr);
	}
	
	GMT_free ((void *)lon);
	mxFree(a);
	GMT_end (argc, argv);

	/* Transpose from gmt grd orientation to Matlab orientation */
	/* Because we need to do the transposition and also a type conversion, we need a extra array */
	nx = grd_b.nx;		ny = grd_b.ny;

	if (is_double) {
		o_d = mxCalloc (nx*ny, sizeof (double));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_d[j*ny+ny-i-1] = (double)b[i*nx+j];
		plhs[0] = mxCreateDoubleMatrix (ny,nx, mxREAL);
		pdata_d = mxGetPr(plhs[0]);
		memcpy(pdata_d, o_d, ny*nx * 8);
	}
	else if (is_single) {
		o_s = mxCalloc (nx*ny, sizeof (float));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_s[j*ny+ny-i-1] = b[i*nx+j];
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		pdata_s = mxGetData(plhs[0]);
		memcpy(pdata_s, o_s, ny*nx * 4);
	}
	else if (is_int32) {
		o_i4 = mxCalloc (nx*ny, sizeof (int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_i4[j*ny+ny-i-1] = irint(b[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT32_CLASS,mxREAL);
		pdata_i4 = mxGetData(plhs[0]);
		memcpy(pdata_i4, o_i4, ny*nx * 4);
	}
	else if (is_int16) {
		o_i2 = mxCalloc (nx*ny, sizeof (short int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_i2[j*ny+ny-i-1] = (short int)irint(b[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT16_CLASS,mxREAL);
		pdata_i2 = mxGetData(plhs[0]);
		memcpy(pdata_i2, o_i2, ny*nx * 2);
	}
	else if (is_uint16) {
		o_ui2 = mxCalloc (nx*ny, sizeof (short int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_ui2[j*ny+ny-i-1] = (unsigned short int)irint(b[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT16_CLASS,mxREAL);
		pdata_ui2 = mxGetData(plhs[0]);
		memcpy(pdata_ui2, o_ui2, ny*nx * 2);
	}
	else if (is_uint8) {
		o_ui1 = mxCalloc (nx*ny, sizeof (char));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_ui1[j*ny+ny-i-1] = (unsigned char)b[i*nx+j];
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT8_CLASS ,mxREAL);
		pdata_ui1 = mxGetData(plhs[0]);
		memcpy(pdata_ui1, o_ui1, ny*nx * 1);
	}

	mxFree(b);

}
