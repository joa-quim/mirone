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

/*--------------------------------------------------------------------
 * Read an existing GMT cpt table and convert it to a Matlab colormap.
 * That is, a mx3 matrix with RGB values ranging from 0 to 1
 *
 * Author:	Joaquim Luis
 * Date:	4-Nov-2004
 * Modified:	18-Feb-2005
 *		May return the corresponding colors grid Z slices
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		19/05/08 J Luis, Patch to not free GMT_lut in GMT_end_for_mex (a crash inducer)
 */

#include "gmt.h"
#include "mex.h"

void sample_cpt (double z[], int nz, int continuous, int reverse, int log_mode, double *pal);

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, nz, log_mode = 0, argc = 0, n_arg_no_char = 0, n_colors;
	int error = FALSE, ok = FALSE, continuous = FALSE, reverse = FALSE;
	char	*table = CNULL, CPT_file[BUFSIZ], **argv;
	double	*z, z_start = 0, z_stop = 256, z_inc = 1;
	double	*pal, *z_ints;

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
	argv[0] = "cpt2cmap";
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

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'C':
					table = &argv[i][2];
					break;
				case 'I':
					reverse = TRUE;
					break;
				case 'Q':
					if (argv[i][2] == 'o')	/* Input data is z, but take log10(z) before interpolation colors */
						log_mode = 2;
					else			/* Input is log10(z) */
						log_mode = 1;
					break;
				case 'T':
					if (sscanf(&argv[i][2], "%lf/%lf/%lf", &z_start, &z_stop, &z_inc) != 3) {
						mexPrintf("%s: GMT SYNTAX ERROR -T option:  Cannot decode values\n", GMT_program);
						error++;
					}
					else {
						if (z_stop <= z_start || z_inc <= 0.0) {
							mexPrintf("%s: GMT SYNTAX ERROR -T option:  Bad arguments\n", GMT_program);
							error++;
						}
					}
					break;
				case 'Z':
					continuous = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}
	
	if (argc == 1 || error) {
		mexPrintf ("cpt2cmap - Read and convert a GMT color palette to a Matlab colormap\n\n");
		mexPrintf ("usage: cmap = cpt2cmap('-C<table>', ['-Q[i|o]]', ['-I'], ['-T<z_start>/<z_stop>/<z_inc>'], ['-Z'])\n");
		mexPrintf ("\t-C Specify a colortable file name\n");
		
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t   ---------------------------------\n");
		mexPrintf ("\t-I reverses the sense of the color table\n");
		mexPrintf ("\t-Q assign a logarithmic colortable [Default is linear]\n");
		mexPrintf ("\t   -Qi: z-values are actually log10(z). Assign colors and write z. [Default]\n");
		mexPrintf ("\t   -Qo: z-values are z, but take log10(z), assign colors and write z.\n");
		mexPrintf ("\t-T Sample points should Step from z_start to z_stop by z_inc [Default 0/256/1].\n");
		mexPrintf ("\t-Z will create a continuous color palette.\n");
		mexErrMsgTxt("\n");
	}
	
	if (nlhs == 0) {
		mexPrintf("CPT2CMAP ERROR: Must provide an output.\n");
		mexErrMsgTxt("\n");
	}


	if (!table)
		mexErrMsgTxt("CPT2CMAP ERROR: Must provide a GMT color palette.\n");
	else {
		if (strstr (table, ".cpt"))
			strcpy (CPT_file, table);
		else
			sprintf (CPT_file, "%s.cpt", table);

		ok = !access (CPT_file, R_OK);
	}

	if (!ok) {
		mexPrintf ("%s: ERROR: Cannot find colortable %s\n", GMT_program, CPT_file);
		mexErrMsgTxt("\n"); 
	}

	if ((GMT_read_cpt (CPT_file)) == EXIT_FAILURE)
		mexErrMsgTxt("GMT Error while reading palette file");

	nz = irint ((z_stop - z_start) / z_inc) + 1;
	z = (double *) mxMalloc ((size_t)nz * sizeof(double));
	for (i = 0; i < nz; i++) z[i] = z_start + i * z_inc;	/* Desired z values */


	plhs[0] = mxCreateDoubleMatrix (nz-1,3, mxREAL);
	pal = mxGetPr(plhs[0]);
	sample_cpt (z, nz, continuous, reverse, log_mode, pal);

	if (nlhs == 2) {	/* Also returns the grid slice values */
		plhs[1] = mxCreateDoubleMatrix (GMT_n_colors,2, mxREAL);
		z_ints = mxGetPr(plhs[1]);
		for (i = 0; i < GMT_n_colors; i++) {
			z_ints[i] = GMT_lut[i].z_low;
			z_ints[i+GMT_n_colors] = GMT_lut[i].z_high;
		}
	}
	/*if (nlhs == 3) {	/* Also returns the BFN colors */
		/*plhs[2] = mxCreateDoubleMatrix (3,3, mxREAL);
		bfn = mxGetPr(plhs[2]);
		for (i = 0; i < 3; i++) {
			bfn[i] = GMT_bfn[i].rgb;
		}
	}*/

	mxFree ((void *)z);
	
	/* If we let GMT_end_for_mex free GMT_lut than a posterior call to grdgradient_m will crash
	   The trick is to temporarely set GMT_n_colors to 0 and than the if () case is not executed
	   This certainly results in one more memory leak */
	/*n_colors = GMT_n_colors;
	GMT_n_colors = FALSE;
	GMT_end_for_mex (argc, argv);
	GMT_n_colors = n_colors;*/
	GMT_end (argc, argv);
	GMT_n_colors = FALSE;		/* HAVE TO, because of the comment above */
}

void sample_cpt (double z[], int nz, int continuous, int reverse, int log_mode, double *pal) {
	/* Resamples the current cpt table based on new z-array.
	 * Old cpt is normalized to 0-1 range and scaled to fit new z range.
	 * New cpt may be continuous and/or reversed.
	 * This a cripled version of GMT_sample_cpt */

	int i, j, k, nx, upper, lower, rgb_low[3], rgb_high[3];
	double *x, *z_out, a, b, f;
	struct GMT_LUT *lut;

	if (!GMT_continuous && continuous) mexPrintf ("%s: Warning: Making a continous cpt from a discrete cpt may give unexpected results!\n", GMT_program);

	lut = (struct GMT_LUT *) GMT_memory (VNULL, (size_t)GMT_n_colors, sizeof (struct GMT_LUT), GMT_program);
	
	/* First normalize old cpt file so z-range is 0-1 */

	b = 1.0 / (GMT_lut[GMT_n_colors-1].z_high - GMT_lut[0].z_low);
	a = -GMT_lut[0].z_low * b;

	for (i = 0; i < GMT_n_colors; i++) {	/* Copy/normalize cpt file and reverse if needed */
		lut[i].z_low = a + b * GMT_lut[i].z_low;
		lut[i].z_high = a + b * GMT_lut[i].z_high;
		if (reverse) {
			j = GMT_n_colors - i - 1;
			memcpy ((void *)lut[i].rgb_high, (void *)GMT_lut[j].rgb_low,  (size_t)(3 * sizeof (int)));
			memcpy ((void *)lut[i].rgb_low,  (void *)GMT_lut[j].rgb_high, (size_t)(3 * sizeof (int)));
		}
		else {
			j = i;
			memcpy ((void *)lut[i].rgb_high, (void *)GMT_lut[j].rgb_high, (size_t)(3 * sizeof (int)));
			memcpy ((void *)lut[i].rgb_low,  (void *)GMT_lut[j].rgb_low,  (size_t)(3 * sizeof (int)));
		}
	}
	lut[0].z_low = 0.0;			/* Prevent roundoff errors */
	lut[GMT_n_colors-1].z_high = 1.0;
	
	/* Then set up normalized output locations x */

	nx = (continuous) ? nz : nz - 1;
	x = (double *) GMT_memory (VNULL, (size_t)nz, sizeof(double), GMT_program);
	if (log_mode) {	/* Our z values are actually log10(z), need array with z for output */
		z_out = (double *) GMT_memory (VNULL, (size_t)nz, sizeof(double), GMT_program);
		for (i = 0; i < nz; i++) z_out[i] = pow (10.0, z[i]);
	}
	else
		z_out = z;	/* Just point to the incoming z values */

	if (nx == 1) {	/* Want a single color point, assume 1/2 way */
		x[0] = 0.5;
	}
	else {	/* As with LUT, translate users z-range to 0-1 range */
		b = 1.0 / (z[nz-1] - z[0]);
		a = -z[0] * b;
		for (i = 0; i < nz; i++) x[i] = a + b * z[i];	/* Normalized z values 0-1 */
		x[nz-1] = 1.0;	/* To prevent bad roundoff */
	}

	for (i = 0; i < nz-1; i++) {

		lower = i;
		upper = i + 1;
		
		/* Interpolate lower color */

		for (j = 0; j < GMT_n_colors && x[lower] >= lut[j].z_high; j++);
		if (j == GMT_n_colors) j--;

		f = 1.0 / (lut[j].z_high - lut[j].z_low);
		
		for (k = 0; k < 3; k++) rgb_low[k] = lut[j].rgb_low[k] + irint ((lut[j].rgb_high[k] - lut[j].rgb_low[k]) * f * (x[lower] - lut[j].z_low));

		if (continuous) {	/* Interpolate upper color */

			while (j < GMT_n_colors && x[upper] > lut[j].z_high) j++;

			f = 1.0 / (lut[j].z_high - lut[j].z_low);
			
			for (k = 0; k < 3; k++) rgb_high[k] = lut[j].rgb_low[k] + irint ((lut[j].rgb_high[k] - lut[j].rgb_low[k]) * f * (x[upper] - lut[j].z_low));
		}
		else	/* Just copy lower */
			for (k = 0; k < 3; k++) rgb_high[k] = rgb_low[k];

		pal[i] = rgb_low[0] / 255.0;
		pal[i+nz-1] = rgb_low[1] / 255.0;
		pal[i+2*(nz-1)] = rgb_low[2] / 255.0;
	}

	GMT_free ((void *)x);
	GMT_free ((void *)lut);
	if (log_mode) GMT_free ((void *)z_out);
}
