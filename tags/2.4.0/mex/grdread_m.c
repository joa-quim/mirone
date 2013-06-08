/*
 *	$Id: grdread.c,v 1.2 2003/10/20 17:43:41 pwessel Exp $
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
 *      Contact info: www.soest.hawaii.edu/wessel
 *--------------------------------------------------------------------*/
/* Program:	grdread.c
 * Purpose:	matlab callable routine to read a grd file
 * Author:	P Wessel, modified from D Sandwell's original version
 * Date:	07/01/93
 * Update:	06/04/96: P Wessel: Now can return [x,y,z] as option.
 *		09/15/97 Phil Sharfstein: modified to Matlab 5 API
 *		10/06/98 P Wessel, upgrade to GMT 3.1 function calls
 *		11/12/98 P Wessel, ANSI-C and calls GMT_begin()
 *		10/07/99 P Wessel, Did not set x,y if [x,y,z,d] was used
 *		10/20/03 P Wessel, longer path names [R Mueller]
 *		03/22/04 J Luis, allow for recognition of non standard GMT grid formats
 *		09/09/04 J Luis, Also accepts output z array as single (float).
 *			 Due to that the grdwrite subroutine was integrated in the gateway routine
 *		25/02/05 J Luis, Added the "insitu" option
 *		04/09/05 J Luis, Extended "info" to return also GMT_grd_i_format
 *		This is as temporary solution until a propor mechanism exists to tell us
 *		what is the data type we are reading. I consider that extremly important
 *		now that data size grows indefnitely and short ints are a "good salvation" 
 *
 *		04/06/06 J Luis, Replaced call to GMT_grdio_init to GMT_grd_init
 *				 as required by version 4.1.?
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */
 
#include "gmt.h"
#include "mex.h"

#define f_swap(x, y) {float tmp; tmp = x, x = y, y = tmp;}
#define ij(i,j) ((i) + (j)*n_cols)

void grd_FLIPLR(float data[], int nx, int ny);
int troca_insitu(float *a, int nx, int ny);

int n_cols;
/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	struct GRD_HEADER grd;
	float *z_4, *z_4o;          /* z_4o if real array for output is selected */
	double *z_8, *info = (double *)NULL, *x, *y, off;
	char *filein, *argv = "grdread-mex", *text, string[BUFSIZ];
	char first_opt[10], second_opt[10];
	GMT_LONG	ns, ssz, pz, i, j, pad[4];
	int is_double = 0;	/* Defaults to output in double */
	int insitu = 0;		/* Defaults to no in-situ transposition */
	int grd_min_size = 0, argc = 0;
 
	/* info contains xmin, xmax, ymin, ymax, zmin, zmax, node-offset, dx, dy */

	pad[0] = pad[1] = pad[2] = pad[3] = 0;

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, &argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, &argv);*/
	argc = GMT_begin (argc, &argv);

	GMT_grd_init (&grd, 0, &argv, FALSE);	/* New in version 4.?.? */

	if (nrhs == 0 || nrhs > 3 || nlhs < 1 || nlhs > 4) {
		mexPrintf ("usage: z = grdread('filename');\n");
		mexPrintf (" 	[z,info] = grdread('filename');\n");
		mexPrintf ("	[x,y,z] = grdread('filename');\n");
		mexPrintf ("	[x,y,z,info] = grdread('filename');\n");
		mexPrintf ("	[...] = grdread('filename','single','insitu');\n");
		return;
	}

	/* Load the file name into a char string */

	ns = mxGetN (prhs[0]) + 1;
	ssz = ns * sizeof (mxChar);

	if (ssz > BUFSIZ)
		mexErrMsgTxt ("grdread: filename too long\n");

	filein = mxMalloc (ssz);

	if (mxGetString (prhs[0], filein, ns + 1)) {
		mexPrintf ("%s\n", filein);
		mexErrMsgTxt ("grdread: failure to decode string \n");
	}

	if (nrhs == 2) is_double = 0;

	if (nrhs == 3 && mxIsChar(prhs[2])) {
		ns = mxGetN (prhs[2]) + 1;
		mxGetString (prhs[2], second_opt, ns + 1);
		if (strcmp(second_opt,"insitu") == 0) insitu = TRUE;
	}
	else if (nrhs == 3 && !mxIsChar(prhs[2]))
		grd_min_size = (int)(*mxGetPr(prhs[2]));

	/* Read the header */
	if (GMT_read_grd_info (filein, &grd))
		mexErrMsgTxt ("grdread: failure to read header\n");

	/* See if in-situ transposition was required on basis of the grid size */
	if (grd_min_size > 0)
		if (grd.nx * grd.ny > grd_min_size) insitu = TRUE; 

	/* Create a matrix for the return array */

	pz = (nlhs >= 3) ? 2 : 0;

	if (is_double) {
		plhs[pz] = mxCreateDoubleMatrix(grd.ny, grd.nx, mxREAL);
		z_8 = mxGetPr(plhs[pz]);
	}
	else {
		plhs[pz] = mxCreateNumericMatrix(grd.ny,grd.nx,mxSINGLE_CLASS,mxREAL);
		z_4o = mxGetData(plhs[pz]);
	}
    
	/* Create scalars for return arguments */

	if (nlhs == 2) {	/* Also return info array */
		plhs[1] = mxCreateDoubleMatrix(1, 10, mxREAL);
		info = mxGetPr(plhs[1]);
	}
	else if (nlhs >= 3) {	/* Return x,y arrays instead */
		plhs[0] = mxCreateDoubleMatrix(1, grd.nx, mxREAL);
		plhs[1] = mxCreateDoubleMatrix(1, grd.ny, mxREAL);
		x = mxGetPr(plhs[0]);
		y = mxGetPr(plhs[1]);
		if (nlhs == 4) {	/* Also return info array */
			plhs[3] = mxCreateDoubleMatrix(1, 10, mxREAL);
			info = mxGetPr(plhs[3]);
		}
	}
 
	if (info) {
		info[0] = grd.x_min;	info[1] = grd.x_max;
		info[2] = grd.y_min;	info[3] = grd.y_max;
		info[4] = grd.z_min;	info[5] = grd.z_max;
		info[6] = grd.node_offset;
		info[7] = grd.x_inc;	info[8] = grd.y_inc;
		/* This is temporary, until a better mechanisms exists to tell us the data type */
		info[9] = grd.type;
	}
 
 	/* Check for file access here since the exit returned by the read routine shuts down Matlab... */
	
	/* In order that non standard GMT grid formats can be recognized, we have to strip
	   the chars "=..." from the file name (in case that they do exist) */
	strcpy (string, filein);
	if ((text = strtok (string, "=")) == NULL) {
		if (access (filein, R_OK));
			mexErrMsgTxt("grdread0: failure to read file\n");
	}
	else {
		if (access (text, R_OK))
			mexErrMsgTxt("grdread: failure to read file\n");
	}

	/*  Allocate memory */
	if (!insitu) {
		if ((z_4 = (float *) mxMalloc (sizeof (float) * grd.nx * grd.ny)) == (float *)NULL)
			mexErrMsgTxt("grdread: failure to allocate memory\n");

		/*  Read the grid */
		if (GMT_read_grd(filein, &grd, z_4, 0.0, 0.0, 0.0, 0.0, pad, 0)) {
			mxFree ((void *)z_4); 
			mexErrMsgTxt("grdread: failure to read file\n");
		}
	}
	else {
		/*  Read the grid */
		if (GMT_read_grd(filein, &grd, z_4o, 0.0, 0.0, 0.0, 0.0, pad, 0)) {
			mxFree ((void *)z_4o); 
			mexErrMsgTxt("grdread: failure to read file\n");
		}
	}
 
	/*  Load the real grd array into a double (or single) matlab array
	    by transposing from grd format to matlab format */
    
	n_cols = grd.nx;
	if (!insitu) {
		if (is_double)
			for (i = 0; i < grd.ny; i++) for (j = 0; j < grd.nx; j++) z_8[j*grd.ny+grd.ny-i-1] = z_4[i*grd.nx+j];
		else
			for (i = 0; i < grd.ny; i++) for (j = 0; j < grd.nx; j++) z_4o[j*grd.ny+grd.ny-i-1] = z_4[i*grd.nx+j];
	}
	else {
		troca_insitu(z_4o, grd.nx, grd.ny);
		grd_FLIPLR(z_4o,grd.ny,grd.nx);
	}

	/*  Free memory */
	if (!insitu)
		mxFree ((void *)z_4);

	if (nlhs >= 3) {	/* Fill in the x and y arrayx */
		off = (grd.node_offset) ? 0.5 : 0.0;
		for (i = 0; i < grd.nx; i++) x[i] = grd.x_min + (i + off) * grd.x_inc;
		for (i = 0; i < grd.ny; i++) y[i] = grd.y_min + (i + off) * grd.y_inc;
	}

	mxFree (filein);
	GMT_end (argc, &argv);
	return;
}

int troca_insitu(float *a, int n, int m) {
	int i,row,column,current, size = m*n;
	float temp;
	for(i = 1, size -= 2; i < size; i++) {
		current = i;
		do      {
			column = current / m;
			row = current % m;
			current = n*row +  column;
		} while(current < i);

		if (current > i)
			f_swap(a[i],a[current]);
        }
	return (0);
}

void grd_FLIPLR (float data[], int nx, int ny) {
	int i, j, k, k0, nx1, nx_half;

	/* Reverse order of all rows */
	nx_half = nx / 2;
	nx1 = nx - 1;
	for (j = k0 = 0; j < ny; j++, k0 += nx) {	/* Do this to all rows */
		for (i = 0, k = nx1; i < nx_half; i++, k--) {
			f_swap (data[k0+i], data[k0+k]);
		}
	}
}
