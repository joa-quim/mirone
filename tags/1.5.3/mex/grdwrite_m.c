/*
 *	$Id: grdwrite.c,v 1.3 2003/10/20 17:43:41 pwessel Exp $
 *
 *      Copyright (c) 1999-2001 by P. Wessel
 *      See COPYING file for copying and redistribution conditions.
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
 *      Contact info: www.soest.hawaii.edu/pwessel
 *--------------------------------------------------------------------*/
/* Program:	grdwrite.c
 * Purpose:	matlab callable routine to write a grd file
 * Author:	P Wessel, modified from D Sandwell's original version
 * Date:	07/01/93
 * Updates:	06/04/96: Now can take [x,y,z,file,title] instead, assuming node-format
 *		09/15/97 Phil Sharfstein: modified to Matlab 5 API
 *		10/06/98 P Wessel, upgrade to GMT 3.1 function calls
 *		11/12/98 P Wessel, ANSI-C and calls GMT_begin()
 *		10/20/03 P Wessel, longer path names [R Mueller]
 *		09/09/04 J Luis, Also accepts input Z array as single (float).
 *			 Due to that the grdwrite subroutine was integrated in the gateway routine
 *		04/06/06 J Luis, Replaced call to GMT_grdio_init to GMT_grd_init
 *				 as required by version 4.1.?
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		25/08/07 J Luis, Input Z array may be of any type (except int64)
*/
 
#include "gmt.h"
#include "mex.h"

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
 /* Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
/* z_8[]	:	double array for input */
/* x[], y[]	:	arrays for x/y */
/* info[]	:	array for xmin, xmax, ymin, ymax, zmin, zmax, node-offset dx dy */
/* *fileout	:	output filename */
/* *title	:	title */
/* nx		:	number of x points */
/* ny		:	number of y points */
/* pix		:	1 if pixel reg, 0 if gridline registered */

	GMT_LONG	i, j, i2, pad[4], error = 0, argc = 0;
	GMT_LONG	ns, ssz, nx, ny, k, pix = 0;	/* If no info we assume gridline reg */
	int	is_double = 0, is_single = 0, is_int32 = 0, is_uint32 = 0;
	int	is_int16 = 0, is_uint16 = 0, is_int8 = 0, is_uint8 = 0;
	float *z_4;           /* real array for output */

	double		*z_8, *info = (double *)NULL, *x, *y;
	float		*z_4i;
	int		*i_4;
	unsigned int	*ui_4;
	short int 	*i_2;
	unsigned short int *ui_2;
	char 		*i_1;
	unsigned char	*ui_1;
	char *fileout, title[80], *argv = "grdwrite-mex";
	struct GRD_HEADER grd; 

	pad[0] = pad[1] = pad[2] = pad[3] = 0;

	if (nrhs < 3 || nrhs > 6) {
		mexPrintf ("usage: grdwrite(Z, D, 'filename' [,'title']);\n");
		mexPrintf ("       grdwrite(X, Y, Z, 'filename', 'title', [1]);\n");
		return;
	}

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, &argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, &argv);*/
	argc = GMT_begin (argc, &argv);

	GMT_grd_init (&grd, 0, &argv, FALSE);	/* New in version 4.?.? */
 
	/* Load the file name into a char string */

	k = (nrhs >= 5) ? 3 : 2;
	ns = mxGetN (prhs[k]) + 1;
	ssz = ns * sizeof (mxChar);

	if (ssz > BUFSIZ)
		mexErrMsgTxt ("grdread: filename too long\n");

	fileout = mxMalloc (ssz);

	if (mxGetString (prhs[k], fileout, ns + 1) ) {
		mexPrintf ("%s\n", fileout);
		mexErrMsgTxt ("grdwrite: failure to decode string\n");
	}
	
	title[0] = 0;
	if (nrhs >= 4) { /* Load the title into a char string */

		k = (nrhs >= 5) ? 4 : 3;
		ns = mxGetN (prhs[k]) + 1;

		if (mxGetString (prhs[k], title, ns + 1) ) {
			mexPrintf ("%s\n", title);
			mexErrMsgTxt (" *** grdwrite  failure to decode string \n");
		}
	}

	/*  get the data and info pointers */
    
	if (nrhs >= 5) {
		x   = mxGetPr (prhs[0]);
		y   = mxGetPr (prhs[1]);
		/* Find out in which data type was given the input array */
		if (mxIsDouble(prhs[0])) {
			z_8 = mxGetPr (prhs[2]);
			is_double = 1;
		}
		else if (mxIsSingle(prhs[0])) {
			z_4i = mxGetData(prhs[2]);
			is_single = 1;
		}
		else if (mxIsInt32(prhs[0])) {
			i_4 = (int *)mxGetData(prhs[0]);
			is_int32 = 1;
		}
		else if (mxIsUint32(prhs[0])) {
			ui_4 = (unsigned int *)mxGetData(prhs[0]);
			is_uint32 = 1;
		}
		else if (mxIsInt16(prhs[0])) {
			i_2 = (short int *)mxGetData(prhs[0]);
			is_int16 = 1;
		}
		else if (mxIsUint16(prhs[0])) {
			ui_2 = (unsigned short int *)mxGetData(prhs[0]);
			is_uint16 = 1;
		}
		else if (mxIsInt8(prhs[0])) {
			i_1 = (char *)mxGetData(prhs[0]);
			is_int8 = 1;
		}
		else if (mxIsUint8(prhs[0])) {
			ui_1 = (unsigned char *)mxGetData(prhs[0]);
			is_uint8 = 1;
		}
		else {
			mexPrintf("grdwrite error: Unknown input data type.\n");
			mexErrMsgTxt("Valid types are:double, single, Int32, Int16, UInt16, Int8, UInt8.\n");
		}
		nx  = mxGetN (prhs[2]);
		ny  = mxGetM (prhs[2]);
	}
	else {
		/* Find out in which data type was given the input array */
		if (mxIsDouble(prhs[0])) {
			z_8 = mxGetPr (prhs[0]);
			is_double = 1;
		}
		else if (mxIsSingle(prhs[0])) {
			z_4i = mxGetData(prhs[0]);
			is_single = 1;
		}
		else if (mxIsInt32(prhs[0])) {
			i_4 = (int *)mxGetData(prhs[0]);
			is_int32 = 1;
		}
		else if (mxIsUint32(prhs[0])) {
			ui_4 = (unsigned int *)mxGetData(prhs[0]);
			is_uint32 = 1;
		}
		else if (mxIsInt16(prhs[0])) {
			i_2 = (short int *)mxGetData(prhs[0]);
			is_int16 = 1;
		}
		else if (mxIsUint16(prhs[0])) {
			ui_2 = (unsigned short int *)mxGetData(prhs[0]);
			is_uint16 = 1;
		}
		else if (mxIsInt8(prhs[0])) {
			i_1 = (char *)mxGetData(prhs[0]);
			is_int8 = 1;
		}
		else if (mxIsUint8(prhs[0])) {
			ui_1 = (unsigned char *)mxGetData(prhs[0]);
			is_uint8 = 1;
		}
		else {
			mexPrintf("grdwrite error: Unknown input data type.\n");
			mexErrMsgTxt("Valid types are:double, single, Int32, Int16, UInt16, Int8, UInt8.\n");
		}
		nx = mxGetN (prhs[0]);
		ny = mxGetM (prhs[0]);
		info = mxGetPr (prhs[1]);
	}

	GMT_grd_init (&grd, 0, (char **)NULL, 0);
	
	grd.nx = nx;
	grd.ny = ny;
	if (info) {
		grd.x_min = info[0];
		grd.x_max = info[1];
		grd.y_min = info[2];
		grd.y_max = info[3];
		grd.node_offset = (int) info[6];
		if (grd.node_offset) {
			grd.x_inc = (grd.x_max - grd.x_min) / grd.nx;
			grd.y_inc = (grd.y_max - grd.y_min) / grd.ny;
		}
		else {
			grd.x_inc = (grd.x_max - grd.x_min) / (grd.nx - 1);
			grd.y_inc = (grd.y_max - grd.y_min) / (grd.ny - 1);
		}
	}
	else {
		grd.x_inc = x[1] - x[0];
		grd.y_inc = y[1] - y[0];
		for (i = 2; !error && i < grd.nx; i++) if ((x[i] - x[i-1]) != grd.x_inc) error = 1;
		for (j = 2; !error && j < grd.ny; j++) if ((y[j] - y[j-1]) != grd.y_inc) error = 1;
		if (error)
			mexErrMsgTxt ("grdwrite: x and/or y not equidistant\n"); 
		grd.x_min = (pix) ? x[0] - 0.5 * grd.x_inc : x[0];
		grd.x_max = (pix) ? x[nx-1] + 0.5 * grd.x_inc : x[nx-1];
		grd.y_min = (pix) ? y[0] - 0.5 * grd.y_inc : y[0];
		grd.y_max = (pix) ? y[ny-1] + 0.5 * grd.y_inc : y[ny-1];
		grd.node_offset = pix;
	}
		
	strcpy (grd.remark, "File written from Matlab with grdwrite");
	
	/*  Allocate memory */

	if ((z_4 = (float *) malloc (sizeof (float) * nx * ny)) == (float *)NULL)
		mexErrMsgTxt ("grdwrite: failure to allocate memory\n");

	/* Transpose from Matlab orientation to grd orientation */

	if (is_double)
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = (float)z_8[j*ny+i];
	else if (is_single)
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = z_4i[j*ny+i];
	else if (is_int32)
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = (float)i_4[j*ny+i];
	else if (is_uint32)
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = (float)ui_4[j*ny+i];
	else if (is_int16)
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = (float)i_2[j*ny+i];
	else if (is_uint16)
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = (float)ui_2[j*ny+i];
	else if (is_int8)
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = (float)i_1[j*ny+i];
	else		/* mut be (is_uint8) */
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) for (j = 0; j < nx; j++) z_4[i2*nx+j] = (float)ui_1[j*ny+i];
     
	/* Update the header using values passed */

	strncpy (grd.title, title, 80);
 
	/*  Write the file */

	if (GMT_write_grd (fileout, &grd, z_4, 0.0, 0.0, 0.0, 0.0, pad, 0)) {
		free ((void *)z_4);
		mexErrMsgTxt ("grdwrite: failure to write file\n"); 
	}
	
	/*  Free memory */

	free ((void *)z_4);
	mxFree (fileout);
	GMT_end (argc, &argv);
	return;
}
