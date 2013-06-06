/*--------------------------------------------------------------------
 *	$Id: grd2cpt.c,v 1.14 2004/07/19 02:36:13 pwessel Exp $
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
 * Creates a cumulative distribution function f(z) describing the data
 * in the grdfile.  f(z) is sampled at z values supplied by the user
 * [with -S option] or guessed from the sample mean and standard deviation.
 * f(z) is then found by looping over the grd array for each z and counting
 * data values <= z.  Once f(z) is found then a master cpt table is resampled
 * based on a normalized f(z).
 *
 * Author:	Walter H. F. Smith
 * Date:	12-JAN-1994
 * Revised:	PW: 12-MAY-1998, for GMT 3.1
 *		PW: 08-MAR-1998, for GMT 3.2 to allow use of master cptfiles
 *		PW: 08-JUL-2000, for GMT 3.3.5
 *		JL: 27-APR-2003, added a -R option
 *		SE: 17-SEP-2003, added a -E option
 * Version:	4
 * 
 */

#include "gmt.h"
#include "mex.h"

struct CDF_CPT {
	double	z;	/* Data value  */
	double	f;	/* Cumulative distribution function f(z)  */
} *cdf_cpt = NULL;

void GMT_end_for_mex (int argc, char **argv);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, nxy, nx, ny, one_or_zero, nfound, ngood, ncdf, log_mode = 0;
	int error = FALSE, set_limits = FALSE, set_z_vals = FALSE, ok = FALSE;
	int equal_inc = FALSE, scale = TRUE;
	int is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int is_uint16 = FALSE;
	int	argc = 0, n_arg_no_char = 0, nc_h, nr_h, i2, *i_4, *o_i4;
	short int *i_2;
	unsigned short int *ui_2, *o_ui2;
	char	**argv;
	float	*zdata, *z_4;
	double *z, min_limit, max_limit, z_start, z_stop, z_inc, mean, sd;
	double	*z_8, *head, a, b;
	struct GRD_HEADER grd;

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
	argv[0] = "grd2cdf";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	argc = GMT_begin (argc, argv);

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'E':
					if (sscanf(&argv[i][2], "%d", &ncdf) != 1) {
						mexPrintf("%s: GMT SYNTAX ERROR -E option:  Cannot decode value\n", GMT_program);
						error++;
					}
					if (!error) equal_inc = TRUE;
					break;
				case 'L':
					if (sscanf(&argv[i][2], "%lf/%lf", &min_limit, &max_limit) != 2) {
						mexPrintf("%s: GMT SYNTAX ERROR -L option:  Cannot decode limits\n", GMT_program);
						error++;
					}
					else {
						if (min_limit >= max_limit) {
							mexPrintf("%s: GMT SYNTAX ERROR -L option:  min_limit must be less than max_limit.\n", GMT_program);
							error++;
						}
					}
					if (!error) set_limits = TRUE;
					break;
				case 'N':
					scale = FALSE;
					break;
				case 'Q':
					if (argv[i][2] == 'o')	/* Input data is z, but take log10(z) before interpolation colors */
						log_mode = 2;
					else			/* Input is log10(z) */
						log_mode = 1;
					break;
				case 'S':
					if (sscanf(&argv[i][2], "%lf/%lf/%lf", &z_start, &z_stop, &z_inc) != 3) {
						mexPrintf("%s: GMT SYNTAX ERROR -S option:  Cannot decode values\n", GMT_program);
						error++;
					}
					else {
						if (z_stop <= z_start || z_inc <= 0.0) {
							mexPrintf("%s: GMT SYNTAX ERROR -S option:  Bad arguments\n", GMT_program);
							error++;
						}
					}
					if (!error) set_z_vals = TRUE;
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
	}
	
	if (n_arg_no_char == 0 || error) {
		mexPrintf ("grd2cdf - Make a linear or histogram-equalized color palette table from a grdfile\n\n");
		mexPrintf ("usage: grd2cdf(indfile, head, ['-E<nlevels>'], ['-L<min_limit>/<max_limit>'],\n");
		mexPrintf ("\t['-Q[i|o]'], ['-S<z_start>/<z_stop>/<z_inc>'])\n");
		
		mexPrintf ("\t<infile> is name of input array\n");
		mexPrintf ("\t<head> is array header descriptor of the form\n");
		mexPrintf ("\t [x_min x_max y_min y_max z_min zmax 0 x_inc y_inc]\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t   ---------------------------------\n");
		mexPrintf ("\t-E nlevels equidistant color levels\n");
		mexPrintf ("\t-L Limit the range of the data [Default uses actual min,max of data].\n");
		mexPrintf ("\t-N do not scale the cdf table into the [0-1] range (default do)\n");
		mexPrintf ("\t-Q assign a logarithmic colortable [Default is linear]\n");
		mexPrintf ("\t   -Qi: z-values are actually log10(z). Assign colors and write z. [Default]\n");
		mexPrintf ("\t   -Qo: z-values are z, but take log10(z), assign colors and write z.\n");
		mexPrintf ("\t-S Sample points should Step from z_start to z_stop by z_inc [Default guesses some values].\n");
		mexErrMsgTxt("\n");
	}
	
	if (nlhs == 0) {
		mexPrintf("GRD2CDF ERROR: Must provide an output.\n");
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
	else {
		mexPrintf("GRD2CDF ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, In32, In16 and UInt16.\n");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2)
		mexErrMsgTxt("GRD2CDF ERROR: First argument must contain a decent array\n");

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("GRD2CDF ERROR: Second argument must contain a valid header of the input array.\n");
	
	head  = mxGetPr(prhs[1]);		/* Get header info */

	grd.x_min = head[0];	grd.x_max = head[1];
	grd.y_min = head[2];	grd.y_max = head[3];
	grd.z_min = head[4];	grd.z_max = head[5];
	grd.x_inc = head[7];	grd.y_inc = head[8];
	grd.nx = nx;			grd.ny = ny;
	grd.node_offset = irint(head[6]);

	nxy = grd.nx * grd.ny;
	zdata = mxCalloc (nxy, sizeof (float));

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) zdata[i2*nx + j] = (float)z_8[j*ny+i];
	}
	else if (is_single) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) zdata[i2*nx + j] = z_4[j*ny+i];
	}
	else if (is_int32) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) zdata[i2*nx + j] = (float)i_4[j*ny+i];
	}
	else if (is_int16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) zdata[i2*nx + j] = (float)i_2[j*ny+i];
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) zdata[i2*nx + j] = (float)ui_2[j*ny+i];
	}

	/* Loop over the file and find NaNs.  If set limits, may create more NaNs  */
	nfound = 0;
	mean = 0.0;
	sd = 0.0;
	if (set_limits) {
		/* Loop over the grdfile, and set anything outside the limiting values to NaN.  */
		grd.z_min = min_limit;
		grd.z_max = max_limit;
		for (i = 0; i < nxy; i++) {
			if (GMT_is_fnan (zdata[i]))
				nfound++;
			else {
				if (zdata[i] < min_limit || zdata[i] > max_limit) {
					nfound++;
					zdata[i] = GMT_f_NaN;
				}
				else {
					mean += zdata[i];
					sd += zdata[i] * zdata[i];
				}
			}
		}
	}
	else {
		min_limit = grd.z_max;	/* This is just to double check grd.z_min, grd.z_max  */
		max_limit = grd.z_min;
		for (i = 0; i < nxy; i++) {
			if (GMT_is_fnan (zdata[i]))
				nfound++;
			else {
				if (zdata[i] < min_limit) min_limit = zdata[i];
				if (zdata[i] > max_limit) max_limit = zdata[i];
				mean += zdata[i];
				sd += zdata[i] * zdata[i];
			}
		}
		grd.z_min = min_limit;
		grd.z_max = max_limit;
	}
	ngood = nxy - nfound;	/* This is the number of non-NaN points for the cdf function  */
	mean /= ngood;
	sd /= ngood;
	sd = sqrt(sd - mean*mean);

	/* Now the zdata are ready.  Decide how to make steps in z.  */
	if (set_z_vals) {
		ncdf =  (grd.z_min < z_start) ? 1 : 0;
		ncdf += (int)floor((z_stop - z_start)/z_inc) + 1;
		if (grd.z_max > z_stop) ncdf++;
		cdf_cpt = (struct CDF_CPT *)GMT_memory (VNULL, (size_t)ncdf, sizeof(struct CDF_CPT), GMT_program);
		if (grd.z_min < z_start) {
			cdf_cpt[0].z = grd.z_min;
			cdf_cpt[1].z = z_start;
			i = 2;
		}
		else {
			cdf_cpt[0].z = z_start;
			i = 1;
		}
		j = (grd.z_max > z_stop) ? ncdf - 1 : ncdf;
		while (i < j) {
			cdf_cpt[i].z = cdf_cpt[i-1].z + z_inc;
			i++;
		}
		if (j == ncdf-1) cdf_cpt[j].z = grd.z_max;
	}
	else {
		/* Make a equaldistant color map from grd.z_min to grd.z_max */
		if(equal_inc) {
			z_inc=(grd.z_max-grd.z_min)/(double)(ncdf-1);
			cdf_cpt = (struct CDF_CPT *)GMT_memory (VNULL, (size_t)ncdf, sizeof(struct CDF_CPT), GMT_program);
			for(i = 0; i < ncdf; i++) {
				cdf_cpt[i].z = grd.z_min+i*z_inc;
			}
		}
		else {
			/* This is completely ad-hoc.  It chooses z based on steps of 0.1 for a Gaussian CDF:  */
			ncdf = 11;
			cdf_cpt = (struct CDF_CPT *)GMT_memory (VNULL, (size_t)ncdf, sizeof(struct CDF_CPT), GMT_program);
                	/* Stupid bug fix here:  
                	        If (mean-1.28155*sd <= grd.z_min || mean+1.28155*sd >= grd.z_max) then
                	          reset mean and sd so they fit inside available range:  */
                          
			if ((mean - 1.28155*sd) <= grd.z_min || (mean + 1.28155*sd) >= grd.z_max) {
				mean = 0.5 * (grd.z_min + grd.z_max);
				sd = (grd.z_max - mean) / 1.5;
				if (sd <= 0.0) {
					mexPrintf ("%s:  ERROR.  Min and Max data values are equal.\n", GMT_program);
					mexErrMsgTxt("\n");
				}
			}	/* End of stupid bug fix  */
		
			cdf_cpt[0].z = grd.z_min;
			cdf_cpt[1].z = mean - 1.28155 * sd;
			cdf_cpt[2].z = mean - 0.84162 * sd;
			cdf_cpt[3].z = mean - 0.52440 * sd;
			cdf_cpt[4].z = mean - 0.25335 * sd;
			cdf_cpt[5].z = mean;
			cdf_cpt[6].z = mean + 0.25335 * sd;
			cdf_cpt[7].z = mean + 0.52440 * sd;
			cdf_cpt[8].z = mean + 0.84162 * sd;
			cdf_cpt[9].z = mean + 1.28155 * sd;
			cdf_cpt[10].z = grd.z_max;
		}
	}
	
	/* Get here when we are ready to go.  cdf_cpt[].z contains the sample points.  */

	for (j = 0; j < ncdf; j++) {
		if (cdf_cpt[j].z <= grd.z_min)
			cdf_cpt[j].f = 0.0;
		else if (cdf_cpt[j].z >= grd.z_max)
			cdf_cpt[j].f = 1.0;
		else {
			nfound = 0;
			for (i = 0; i < nxy; i++) {
				if (!GMT_is_fnan (zdata[i]) && zdata[i] <= cdf_cpt[j].z) nfound++;
			}
			cdf_cpt[j].f = (double)(nfound-1)/(double)(ngood-1);
		}
	}

	/* Now the cdf function has been found  */

	plhs[0] = mxCreateDoubleMatrix (ncdf,1, mxREAL);
	z = mxGetPr(plhs[0]);
	/*z = (double *) GMT_memory (VNULL, (size_t)ncdf, sizeof (double), GMT_program);*/
	for (i = 0; i < ncdf; i++) z[i] = cdf_cpt[i].z;
	if (log_mode == 2) for (i = 0; i < ncdf; i++) z[i] = d_log10 (z[i]);	/* Make log10(z) values for interpolation step */

	if (scale) {			 /* Normalize z-range to 0-1 */
		b = 1.0 / (z[ncdf-1] - z[0]);
		a = -z[0] * b;
		for (i = 0; i < ncdf; i++)
			z[i] = a + b * z[i];

		z[0] = 0.0;			/* Prevent roundoff errors */
		z[ncdf-1] = 1.0;
	}

	GMT_free ((void *)cdf_cpt);
	GMT_free ((void *)zdata);
	/*GMT_free ((void *)z);*/
	
	GMT_end_for_mex (argc, argv);
}

void GMT_end_for_mex (int argc, char **argv) {
	/* GMT_end will clean up after us. */
	int i;

	for (i = 0; i < N_UNIQUE; i++) if (GMT_oldargv[i]) GMT_free ((void *)GMT_oldargv[i]);
	for (i = 0; i < N_UNIQUE; i++) GMT_processed_option[i] = FALSE;
	/*if (GMT_lut) GMT_free ((void *)GMT_lut); */
	/*GMT_free_plot_array ();	/* Unknown function (?) */

#ifdef __FreeBSD__
	fpresetsticky (FP_X_DZ | FP_X_INV);
	fpsetmask (FP_X_DZ | FP_X_INV);
#endif

	/*fflush (GMT_stdout);	/* Make sure output buffer is flushed (in Windows, Matlab BOOMS)*/
}
