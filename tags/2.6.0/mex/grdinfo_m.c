/*
 *	$Id: grdinfo2.c,v 1.3 2003/10/20 17:43:41 pwessel Exp $
 *
 *      Copyright (c) 1999-2010 by P. Wessel
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
/* Program:	grdinfo.c
 * Purpose:	matlab callable routine to read a grd file header
 * Author:	P Wessel, modified from D Sandwell's original version
 * Date:	07/01/93
 * Update:	09/15/97 Phil Sharfstein: modified to Matlab 5 API
 *		10/06/98 P Wessel, upgrade to GMT 3.1 function calls
 *		11/12/98 P Wessel, ANSI-C and calls GMT_begin()
 *		10/20/03 P Wessel, longer path names [R Mueller]
 *		03/22/04 J Luis, allow for recognition of non standard GMT grid formats
 *		02/23/05 J Luis, added a 'silent' option
 *		10/08/05 J Luis, added the 'hdr_struct' option
 *			 In order to implement this option I removed the grdinfo routine
 *		04/06/06 J Luis, Replaced call to GMT_grdio_init to GMT_grd_init
 *				 as required by version 4.1.?
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */
 
#include "mex.h"
#include "gmt.h"

/*int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

 /* Gateway routine */
   
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	struct GRD_HEADER grd;
	double *info = (double *)NULL, *pdata, tmp2[2], tmp3[3], tmp4[4];
	char *filein, *argv = "grdinfo-mex";
	char *text, string[BUFSIZ], *opt_char;
	int silent = 0, struct_hdr = 0;	/* Defaults to original verbose behaviour */
	int ns, ssz, argc = 0;
	char	*fieldnames[8];		/* this array contains the names of the fields of the hdr structure. */
	mxArray *hdr_struct;
	mxArray *mxTitle;
	mxArray *mxCommand;
	mxArray *mxRemark;
	mxArray *mxRegistration;
	mxArray *mxXinfo;
	mxArray *mxYinfo;
	mxArray *mxZinfo;
	mxArray *mxScale;

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, &argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, &argv);*/
	argc = GMT_begin (argc, &argv);

	/* Load the file name into a char string */

	ns = mxGetN (prhs[0]) + 1;
	ssz = ns * sizeof (mxChar);

	if (ssz > BUFSIZ)
		mexErrMsgTxt ("grdread: filename too long\n");

	filein = mxMalloc (ssz);
  
	if (nrhs == 0 || nrhs > 2) {
		mexPrintf (" usage: [info =] grdinfo('filename'[,'silent']); \n");
		mexPrintf (" usage: info = grdinfo('filename'[,'hdr_struct']); \n");
		mexPrintf (" If 'hdr_struct' is passed as the second argument; \n");
		mexPrintf (" info is a structure with the following fields:\n");
		mexPrintf (" Title (string)\n");
		mexPrintf (" Command (string)\n");
		mexPrintf (" Remark (string)\n");
		mexPrintf (" Registration (string)\n");
		mexPrintf (" X_info (a row vector with [x_min x_max x_inc nx])\n");
		mexPrintf (" Y_info (a row vector with [x_min y_max y_inc ny])\n");
		mexPrintf (" Z_info (a row vector with [z_min z_max])\n");
		mexPrintf (" Scale (a row vector with [scale_factor add_offset grdtype])\n");
		return;
	}

	if (mxGetString (prhs[0], filein, ns + 1)) {
		mexPrintf ("%s\n", filein);
		mexErrMsgTxt ("grdinfo: Failure to decode string \n");
	}

	/* Create scalars for return arguments */

	if (nlhs == 1 && nrhs == 1) {		/* If only one input and output */
		plhs[0] = mxCreateDoubleMatrix (1, 9, mxREAL);
		info = mxGetPr (plhs[0]);
 	}
	else if (nlhs == 1 && nrhs == 2) {
		silent = 1;
		opt_char = (char *)mxArrayToString(prhs[1]); 
		if (!strcmp(opt_char,"hdr_struct")) {	/* header output is a structure */
			struct_hdr = 1;
			fieldnames[0] = strdup ("Title");
			fieldnames[1] = strdup("Command");
			fieldnames[2] = strdup("Remark");
			fieldnames[3] = strdup("Registration");
			fieldnames[4] = strdup("X_info");
			fieldnames[5] = strdup("Y_info");
			fieldnames[6] = strdup("Z_info");
			fieldnames[7] = strdup("Scale");
			hdr_struct = mxCreateStructMatrix(1, 1, 8, (const char **)fieldnames);
		}
		else {			/* headr output is a numeric vector */
			plhs[0] = mxCreateDoubleMatrix (1, 9, mxREAL);
			info = mxGetPr (plhs[0]);
		}
	}
 	
	/* Make sure file is readable */
	
	/* In order that non standard GMT grid formats can be recognized, we have to strip
	   the chars "=..." from the file name (in case that they do exist) */

	strcpy (string, filein);
	if ((text = strtok (string, "=")) == NULL) {
		if (access (filein, R_OK))
			mexErrMsgTxt ("grdinfo: Cannot find or open file\n");
	}
	else {
		if (access (text, R_OK))
			mexErrMsgTxt ("grdinfo: Cannot find or open file\n");
	}

	/* Read the header */
 
	if (GMT_read_grd_info (filein, &grd))
		mexErrMsgTxt ("grdinfo: Failure to read header\n");


	if (!silent) {		/* This is used for backward compatibilitie */
		mexPrintf("%s: Title: %s\n", filein, grd.title);
		mexPrintf("%s: Command: %s\n", filein, grd.command);
		mexPrintf("%s: Remark: %s\n", filein, grd.remark);
		if (grd.node_offset)
			mexPrintf("%s: Pixel registration used\n", filein);
		else
			mexPrintf("%s: Normal node registration used\n", filein);
		mexPrintf("%s: x_min: %lg x_max: %lg x_inc: %lg units: %s nx: %d\n",
			filein, grd.x_min, grd.x_max, grd.x_inc, grd.x_units, grd.nx);
		mexPrintf("%s: y_min: %lg y_max: %lg y_inc: %lg units: %s ny: %d\n",
			filein, grd.y_min, grd.y_max, grd.y_inc, grd.y_units, grd.ny);
		mexPrintf("%s: z_min: %lg z_max: %lg units: %s\n",
			filein, grd.z_min, grd.z_max, grd.z_units);
		mexPrintf("%s: scale_factor: %lg add_offset: %lg\n",
			filein, grd.z_scale_factor, grd.z_add_offset);
	}

	if (!struct_hdr) {	/* Return info in a numeric vector */
		info[0] = grd.x_min;
		info[1] = grd.x_max;
		info[2] = grd.y_min;
		info[3] = grd.y_max;
		info[4] = grd.z_min;
		info[5] = grd.z_max;
		info[6] = grd.node_offset;
		info[7] = grd.x_inc;
		info[8] = grd.y_inc;
	}
	else {
		/* First the string fields */
		mxTitle = mxCreateString(grd.title);
		mxSetField(hdr_struct, 0, "Title", mxTitle);
		mxCommand = mxCreateString(grd.command);
		mxSetField(hdr_struct, 0, "Command", mxCommand);
		mxRemark = mxCreateString(grd.remark);
		mxSetField(hdr_struct, 0, "Remark", mxRemark);
		if (grd.node_offset)
			mxRegistration = mxCreateString("Pixel registration used");
		else
			mxRegistration = mxCreateString("Normal node registration used");
		mxSetField(hdr_struct, 0, "Registration", mxRegistration);

		/* Now the numeric fields */
		mxXinfo = mxCreateDoubleMatrix (1, 4, mxREAL);
		pdata = mxGetPr(mxXinfo);
		tmp4[0] = grd.x_min;	tmp4[1] = grd.x_max;	
		tmp4[2] = grd.x_inc;	tmp4[3] = grd.nx;
		memcpy(pdata, tmp4, 4 * sizeof(double));
		mxSetField(hdr_struct, 0, "X_info" , mxXinfo);

		mxYinfo = mxCreateDoubleMatrix (1, 4, mxREAL);
		pdata = mxGetPr(mxYinfo);
		tmp4[0] = grd.y_min;	tmp4[1] = grd.y_max;
		tmp4[2] = grd.y_inc;	tmp4[3] = grd.ny;
		memcpy(pdata, tmp4, 4 * sizeof(double));
		mxSetField(hdr_struct, 0, "Y_info" , mxYinfo);

		mxZinfo = mxCreateDoubleMatrix (1, 2, mxREAL);
		pdata = mxGetPr(mxZinfo);
		tmp2[0] = grd.z_min;	tmp2[1] = grd.z_max;
		memcpy(pdata, tmp2, 2 * sizeof(double));
		mxSetField(hdr_struct, 0, "Z_info" , mxZinfo);

		mxScale = mxCreateDoubleMatrix (1, 3, mxREAL);
		pdata = mxGetPr(mxScale);
		tmp3[0] = grd.z_scale_factor;	tmp3[1] = grd.z_add_offset;
		tmp3[2] = grd.type;
		memcpy(pdata, tmp3, 3 * sizeof(double));
		mxSetField(hdr_struct, 0, "Scale" , mxScale);

		plhs[0] = hdr_struct;
	}
	
	mxFree (filein);
	GMT_end (argc, &argv);
	return;
}
