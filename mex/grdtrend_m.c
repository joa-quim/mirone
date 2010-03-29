/*--------------------------------------------------------------------
 *	$Id: grdtrend.c,v 1.10 2004/01/02 22:45:13 pwessel Exp $
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
/* grdtrend <input.grd> -N<n_model>[r] [-T<trend.grd>] [-V]
	[-W<weight.grd] [-D<differences.grd]

Reads a grdfile and fits a trend surface.  Trend surface
is defined by:

m1 +m2*x + m3*y + m4*xy + m5*x*x + m6*y*y + m7*x*x*x
	+ m8*x*x*y + m9*x*y*y + m10*y*y*y.

n_model is set by the user to be an integer in [1,10]
which sets the number of model coefficients to fit.
Thus:
n_model = 1 gives the mean value of the surface,
n_model = 3 fits a plane,
n_model = 4 fits a bilinear surface,
n_model = 6 fits a biquadratic,
n_model = 10 fits a bicubic surface.

The user may write out grdfiles of the fitted surface
[-T<trend.grd>] and / or of the residuals (input data
minus fitted trend) [-D<differences.grd] and / or of
the weights used in iterative fitting [-W<weight.grd].
This last option applies only when the surface is fit
iteratively [-N<n>[r]].

A robust fit may be achieved by iterative fitting of
a weighted least squares problem, where the weights
are set according to a scale length based on the 
Median absolute deviation (MAD: Huber, 1982).  The
-N<n>r option achieves this.

Author:		W. H. F. Smith
Date:		21 May, 1991.
Version:	4
Calls:		uses the QR solution of the Normal
		equations furnished by Wm. Menke's
		C routine "gauss".  We gratefully
		acknowledge this contribution.
Revised:	12-JUN-1998 PW, for GMT 3.1

Remarks:

We adopt a translation and scaling of the x,y coordinates.
We choose x,y such that they are in [-1,1] over the range
of the grdfile.  If the problem is unweighted, all input
values are filled (no "holes" or NaNs in the input grdfile),
and n_model <= 4 (bilinear or simpler), then the normal
equations matrix (G'G in Menke notation) is diagonal under
this change of coordinates, and the solution is trivial.
In this case, it would be dangerous to try to accumulate
the sums which are the elements of the normal equations;
while they analytically cancel to zero, the addition errors
would likely prevent this.  Therefore we have written a
routine, grd_trivial_model(), to handle this case.

If the problem is more complex than the above trivial case,
(missing values, weighted problem, or n_model > 4), then
G'G is not trivial and we just naively accumulate sums in
the G'G matrix.  We hope that the changed coordinates in
[-1,1] will help the accuracy of the problem.  We also use
Legendre polynomials in this case so that the matrix elements
are conveniently sized near 0 or 1.

/*--------------------------------------------------------------------------------------
 * Mexified version of grdtrend. This version differs from the GMT's in the way that it
 * doesn't allow output of more than one file. That is, ouput of trend and residual is
 * not allowed. We have compute one at a time.
 *
 * Author:	J. Luis
 * Date: 	04 Sep 2004
 *
 * Usage
 * Zout = grdtrend_(Zin,head,'options');
 *
 * where	Z is the array containing the input file
 *			head is the header descriptor in the format of grdread/grdwrite mexs
 *			and options may be for example: '-N4',-T'
 *
 * IMPORTANT NOTE. The data type of Zin is preserved in Zout. That means you can send Zin
 * as a double, single, In32, In16, or UInt16 and receive Zout in one of those types
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */

#include "gmt.h"
#include "mex.h"

#define MAX_TABLE_COLS 10	/* Used by Menke routine gauss  */

char format[BUFSIZ];
/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

	void	gauss(double *a, double *vec, int n, int nstore, double test, int *ierror, int itriag);		/* QR solution of the Normal equations  */
	void	set_up_vals(double *val, int nval, double vmin, double vmax, double dv, int pixel_reg);		/* Store x[i], y[j] once for all to save time  */
	void	load_pstuff(double *pstuff, int n_model, double x, double y, int newx, int newy);		/* Compute Legendre polynomials of x[i],y[j] as needed  */
	void	grd_trivial_model(float *data, int nx, int ny, double *xval, double *yval, double *gtd, int n_model);	/* Fit trivial models.  See Remarks above.  */
	void	compute_trend(float *trend, int nx, int ny, double *xval, double *yval, double *gtd, int n_model, double *pstuff);	/* Find trend from a model  */
	void	compute_resid(float *data, float *trend, float *resid, int nxy);	/* Find residuals from a trend  */
	void	compute_chisq(float *resid, float *weight, int nxy, double *chisq, double scale);	/* Find Chi-Squared from weighted residuals  */
	void	compute_robust_weight(float *resid, float *weight, int nxy, double *scale);	/* Find weights from residuals  */
	void	write_model_parameters(double *gtd, int n_model);	/* Do reports if gmtdefs.verbose == TRUE  */
	void	load_gtg_and_gtd(float *data, int nx, int ny, double *xval, double *yval, double *pstuff, double *gtg, double *gtd, int n_model, float *weight, int weighted);		/* Fill normal equations matrices  */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	i, j, k, nx, ny, ierror = 0, iterations, nxy, n_model = 0;
	int	error = FALSE, robust = FALSE, trivial, weighted;
	int is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int is_uint16 = FALSE, is_uint8 = FALSE;
	int	d_filename = FALSE, t_filename = FALSE, w_filename = FALSE;
	double	chisq, old_chisq, zero_test = 1.0e-08, scale = 1.0;
	int		argc = 0, n_arg_no_char = 0, nc_h, nr_h, i2, *i_4, *o_i4, *pdata_i4;
	short int *i_2, *o_i2, *pdata_i2;
	unsigned short int *ui_2, *o_ui2, *pdata_ui2;
	char	**argv;
	float	*z_4, *pdata_s, *o_s, *out;
	double	*pdata_d, *z_8, *head, *o_d;

	float	*data;		/* Pointer for array from input grdfile  */
	float	*trend;		/* Pointer for array containing fitted surface  */
	float	*resid;		/* Pointer for array containing residual surface  */
	float	*weight;	/* Pointer for array containing data weights  */

	double	*xval;		/* Pointer for array of change of variable:  x[i]  */
	double	*yval;		/* Pointer for array of change of variable:  y[j]  */
	double	*gtg;		/* Pointer for array for matrix G'G normal equations  */
	double	*gtd;		/* Pointer for array for vector G'd normal equations  */
	double	*old;		/* Pointer for array for old model, used for robust sol'n  */
	double	*pstuff;	/* Pointer for array for Legendre polynomials of x[i],y[j]  */

	struct GRD_HEADER head_d;

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
	argv[0] = "grdtrend";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

/* Execution begins here with loop over arguments:  */

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);*/
	argc = GMT_begin (argc, argv);
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				/* Supplemental parameters */
			
				case 'D':
					d_filename = TRUE;
					break;
				case 'N':
					j = 2;
					if (argv[i][j] && (argv[i][j] == 'r' || argv[i][j] == 'r') ){
						robust = TRUE;
						j++;
					}
					if (argv[i][j]) n_model = atoi(&argv[i][j]);
					break;
				case 'T':
					t_filename = TRUE;
					break;
				case 'W':
					w_filename = TRUE;
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (argc == 1 || error) {
		mexPrintf("grdtrend %s - Fit trend surface to gridded data\n\n", "v4.0b");
		mexPrintf("usage:  grdtrend_m(infile,head, -N[r]<n_model>, '-D', '-T', '-W'\n\n");
		
		mexPrintf ("\t<infile> is name of input array\n");
		mexPrintf ("\t<head> is array header descriptor of the form\n");
		mexPrintf ("\t [x_min x_max y_min y_max z_min zmax 0 x_inc y_inc]\n");
		mexPrintf("\t-N # model parameters to fit; integer in [1,10].  Insert r for robust fit.\n");
		mexPrintf("\n\tOPTIONS:\n");
		mexPrintf("\t-D Output the array with differences (input - trend).\n");
		mexPrintf("\t-T Output the array with trend.\n");
		mexPrintf("\t-W If robust = TRUE, output the array with weights.\n");
		return;
	}

	if (n_model <= 0 || n_model > 10) {
		mexPrintf("%s: GMT SYNTAX ERROR -N option:  Specify 1-10 model parameters\n", "grdtrend_m");
		error++;
	}
	if ((d_filename + t_filename + w_filename) == 0) {
		mexPrintf("%s: GMT SYNTAX ERROR: Must select one of -D, -T, -W\n", "grdtrend_m");
		error++;
	}
	if ((d_filename + t_filename) == 2) {
		mexPrintf("%s: GMT WARNING: Residual and trend selected, only residual will be considered\n", "grdtrend_m");
		t_filename = FALSE;
	}
	if ((d_filename + w_filename) == 2) {
		mexPrintf("%s: GMT WARNING: Residual and weights selected, only residual will be considered\n", "grdtrend_m");
		w_filename = FALSE;
	}
	if ((t_filename + w_filename) == 2) {
		mexPrintf("%s: GMT WARNING: Trend and weights selected, only trend will be considered\n", "grdtrend_m");
		w_filename = FALSE;
	}

	if (error) return;

/* End of argument parsing.  */
	
	if (nlhs == 0)
		mexErrMsgTxt("ERROR: Must provide an output.\n");

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
		mexPrintf("GRDSAMPLE ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are: double, single, In32, In16 and UInt16.\n");
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

	head_d.x_min = head[0];		head_d.x_max = head[1];
	head_d.y_min = head[2];		head_d.y_max = head[3];
	head_d.z_min = head[4];		head_d.z_max = head[5];
	head_d.x_inc = head[7];		head_d.y_inc = head[8];
	head_d.nx = nx;			head_d.ny = ny;
	head_d.node_offset = irint(head[6]);
	nxy = nx * ny;

	trivial = ( (n_model < 5) && (!(robust)) && (!w_filename) );

	data = mxCalloc (nxy, sizeof (float));

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*nx + j] = (float)z_8[j*ny+i];
	}
	else if (is_single) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*nx + j] = z_4[j*ny+i];
	}
	else if (is_int32) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*nx + j] = (float)i_4[j*ny+i];
	}
	else if (is_int16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*nx + j] = (float)i_2[j*ny+i];
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*nx + j] = (float)ui_2[j*ny+i];
	}


	/* Check for NaNs:  */
	i = 0;
	while (trivial && i < nxy) {
		if (GMT_is_fnan (data[i])) trivial = FALSE;
		i++;
	}

/* End input read section.  */

/* Allocate other required arrays:  */

	trend = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), "grdtrend_m");
	resid = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), "grdtrend_m");
	xval = (double *) GMT_memory (VNULL, (size_t)head_d.nx, sizeof (double), "grdtrend_m");
	yval = (double *) GMT_memory (VNULL, (size_t)head_d.ny, sizeof (double), "grdtrend_m");
	gtg = (double *) GMT_memory (VNULL, (size_t)(n_model*n_model), sizeof (double), "grdtrend_m");
	gtd = (double *) GMT_memory (VNULL, (size_t)n_model, sizeof (double), "grdtrend_m");
	old = (double *) GMT_memory (VNULL, (size_t)n_model, sizeof (double), "grdtrend_m");
	pstuff = (double *) GMT_memory (VNULL, (size_t)n_model, sizeof (double), "grdtrend_m");
	

/* If a weight array is needed, get one:  */

	weighted = (robust || w_filename);
	if (weighted) {
		weight = (float *) GMT_memory (VNULL, (size_t)nxy, sizeof (float), "grdtrend_m");
		for (i = 0; i < nxy; i++) weight[i] = 1.0;
	}

/* Set up xval and yval lookup tables:  */

	set_up_vals(xval, head_d.nx, head_d.x_min, head_d.x_max, head_d.x_inc, head_d.node_offset);
	set_up_vals(yval, head_d.ny, head_d.y_min, head_d.y_max, head_d.y_inc, head_d.node_offset);

/* End of set up of lookup values.  */

/* Do the problem:  */

	if (trivial) {
		grd_trivial_model(data, head_d.nx, head_d.ny, xval, yval, gtd, n_model);
		compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, n_model, pstuff);
		compute_resid(data, trend, resid, nxy);
	}
	else {	/* Problem is not trivial  !!  */

		load_gtg_and_gtd(data, head_d.nx, head_d.ny, xval, yval, pstuff, gtg, gtd, n_model, weight, weighted);
		gauss(gtg, gtd, n_model, n_model, zero_test, &ierror, 1);
		if (ierror) {
			mexPrintf("%s:  Gauss returns error code %d\n", "grdtrend_m", ierror);
			return;
		}
		compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, n_model, pstuff);
		compute_resid(data, trend, resid, nxy);

		if (robust) {
			compute_chisq(resid, weight, nxy, &chisq, scale);
			iterations = 1;
			sprintf(format, "%%s Robust iteration %%d:  Old Chi Squared:  %s  New Chi Squared %s\n", gmtdefs.d_format, gmtdefs.d_format);
			do {
				old_chisq = chisq;
				for(k = 0; k < n_model; k++) old[k] = gtd[k];
				compute_robust_weight(resid, weight, nxy, &scale);
				load_gtg_and_gtd(data, head_d.nx, head_d.ny, xval, yval, pstuff, gtg, gtd, n_model, weight, weighted);
				gauss(gtg, gtd, n_model, n_model, zero_test, &ierror, 1);
				if (ierror) {
					mexPrintf("%s:  Gauss returns error code %d\n", "grdtrend_m", ierror);
					return;
				}
				compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, n_model, pstuff);
				compute_resid(data, trend, resid, nxy);
				compute_chisq(resid, weight, nxy, &chisq, scale);
				iterations++;
			} while ( (old_chisq / chisq) > 1.0001);

			/* Get here when new model not significantly better; use old one:  */

			for(k = 0; k < n_model; k++) gtd[k] = old[k];
			compute_trend(trend, head_d.nx, head_d.ny, xval, yval, gtd, n_model, pstuff);
			compute_resid(data, trend, resid, nxy);
		}
	}
	mxFree(data);

/* End of do the problem section.  */

/* Get here when ready to do output:  */
	
	/* Transpose from gmt grd orientation to Matlab orientation */
	/* Because we need to do the transposition and likely a type conversion, we need a extra array */
	if (t_filename)
		out = trend;
	else if (d_filename)
		out = resid;
	else if (w_filename && robust)
		out = weight;

	if (is_double) {
		o_d = mxCalloc (nx*ny, sizeof (double));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_d[j*ny+ny-i-1] = (double)out[i*nx+j];
		plhs[0] = mxCreateDoubleMatrix (ny,nx, mxREAL);
		pdata_d = mxGetPr(plhs[0]);
		memcpy(pdata_d, o_d, ny*nx * 8);	mxFree(o_d);
	}
	else if (is_single) {
		o_s = mxCalloc (nx*ny, sizeof (float));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_s[j*ny+ny-i-1] = out[i*nx+j];
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		pdata_s = mxGetData(plhs[0]);
		memcpy(pdata_s, o_s, ny*nx * 4);	mxFree(o_s);
	}
	else if (is_int32) {
		o_i4 = mxCalloc (nx*ny, sizeof (int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_i4[j*ny+ny-i-1] = irint(out[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT32_CLASS,mxREAL);
		pdata_i4 = mxGetData(plhs[0]);
		memcpy(pdata_i4, o_i4, ny*nx * 4);	mxFree(o_i4);
	}
	else if (is_int16) {
		o_i2 = mxCalloc (nx*ny, sizeof (short int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_i2[j*ny+ny-i-1] = (short int)irint(out[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT16_CLASS,mxREAL);
		pdata_i2 = mxGetData(plhs[0]);
		memcpy(pdata_i2, o_i2, ny*nx * 2);	mxFree(o_i2);
	}
	else if (is_uint16) {
		o_ui2 = mxCalloc (nx*ny, sizeof (short int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_ui2[j*ny+ny-i-1] = (unsigned short int)irint(out[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT16_CLASS,mxREAL);
		pdata_ui2 = mxGetData(plhs[0]);
		memcpy(pdata_ui2, o_ui2, ny*nx * 2);	mxFree(o_ui2);
	}

/* That's all, folks!  */

	if (weighted) GMT_free ((void *)weight);
	GMT_free ((void *)pstuff);
	GMT_free ((void *)gtd);
	GMT_free ((void *)gtg);
	GMT_free ((void *)yval);
	GMT_free ((void *)xval);
	GMT_free ((void *)resid);
	GMT_free ((void *)trend);

	GMT_end (argc, argv);
}		

void  set_up_vals(double *val, int nval, double vmin, double vmax, double dv, int pixel_reg)
{
	int	i;
	double  v, middle, drange, true_min, true_max;

	true_min = (pixel_reg) ? vmin + 0.5 * dv : vmin;
	true_max = (pixel_reg) ? vmax - 0.5 * dv : vmax;

	middle = 0.5 * (true_min + true_max);
	drange = 2.0 / (true_max - true_min);
	for (i = 0; i < nval; i++) {
		v = true_min + i * dv;
		val[i] = (v - middle) * drange;
	}
	/* Just to be sure no rounding outside:  */
	val[0] = -1.0;
	val[nval - 1] = 1.0;
	return;
}

void	load_pstuff(double *pstuff, int n_model, double x, double y, int newx, int newy)
{
	/* If either x or y has changed, compute new Legendre polynomials as needed  */

	if (newx) {
		if (n_model >= 2) pstuff[1] = x;
		if (n_model >= 5) pstuff[4] = 0.5*(3.0*pstuff[1]*pstuff[1] - 1.0);
		if (n_model >= 7) pstuff[6] = (5.0*pstuff[1]*pstuff[4] - 2.0*pstuff[1])/3.0;
	}
	if (newy) {
		if (n_model >= 3) pstuff[2] = y;
		if (n_model >= 6) pstuff[5] = 0.5*(3.0*pstuff[2]*pstuff[2] - 1.0);
		if (n_model >= 10) pstuff[9] = (5.0*pstuff[2]*pstuff[5] - 2.0*pstuff[2])/3.0;
	}
	/* In either case, refresh cross terms:  */

	if (n_model >= 4) pstuff[3] = pstuff[1]*pstuff[2];
	if (n_model >= 8) pstuff[7] = pstuff[4]*pstuff[2];
	if (n_model >= 9) pstuff[8] = pstuff[1]*pstuff[5];

	return;
}

void	compute_trend(float *trend, int nx, int ny, double *xval, double *yval, double *gtd, int n_model, double *pstuff)
{
	int	i, j, k, ij;

	for (ij = 0, j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++, ij++) {
			load_pstuff(pstuff, n_model, xval[i], yval[j], 1, (!(i)));
			trend[ij] = (float)gtd[0];
			for (k = 1; k < n_model; k++) {
				trend[ij] += (float)(pstuff[k]*gtd[k]);
			}
		}
	}
}

void	compute_resid(float *data, float *trend, float *resid, int nxy)
{
	int	i;

	for (i = 0; i < nxy; i++) {
		if (GMT_is_fnan (data[i])) {
			resid[i] = data[i];
		}
		else {
			resid[i] = data[i] - trend[i];
		}
	}
	return;
}

void	grd_trivial_model(float *data, int nx, int ny, double *xval, double *yval, double *gtd, int n_model)
{
	/* Routine to fit up elementary polynomial model of grd data, 
	model = gtd[0] + gtd[1]*x + gtd[2]*y + gtd[3] * x * y,
	where x,y are normalized to range [-1,1] and there are no
	NaNs in grdfile, and problem is unweighted least squares.  */

	int	i, j, ij;
	double	x2, y2, sumx2 = 0.0, sumy2 = 0.0, sumx2y2 = 0.0;

	/* First zero the model parameters to use for sums:  */

	for (i = 0; i < n_model; i++) gtd[i] = 0.0;

	/* Now accumulate sums:  */

	for (ij = 0, j = 0; j < ny; j++) {
		y2 = yval[j] * yval[j];
		for (i = 0; i < nx; i++, ij++) {
			x2 = xval[i] * xval[i];
			sumx2 += x2;
			sumy2 += y2;
			sumx2y2 += (x2 * y2);
			gtd[0] += data[ij];
			if (n_model >= 2) gtd[1] += data[ij] * xval[i];
			if (n_model >= 3) gtd[2] += data[ij] * yval[j];
			if (n_model == 4) gtd[3] += data[ij] * xval[i] * yval[j];
		}
	}

	/* See how trivial it is?  */

	gtd[0] /= (nx * ny);
	if (n_model >= 2) gtd[1] /= sumx2;
	if (n_model >= 3) gtd[2] /= sumy2;
	if (n_model == 4) gtd[3] /= sumx2y2;

	return;
}

void	compute_chisq(float *resid, float *weight, int nxy, double *chisq, double scale)
{
	int	i;
	double	tmp;

	*chisq = 0.0;
	for (i = 0; i < nxy; i++) {
		if (GMT_is_fnan (resid[i]))continue;
		tmp = resid[i];
		if (scale != 1.0) tmp /= scale;
		tmp *= tmp;

		if (weight[i] == 1.0) {
			*chisq += tmp;
		}
		else {
			/* Weight has already been squared  */
			*chisq += (tmp * weight[i]);
		}
	}
	return;
}

void	compute_robust_weight(float *resid, float *weight, int nxy, double *scale)
{
	int	i, j, j2;
	double	r, mad;

	for (i = j = 0; i < nxy; i++) {
		if (GMT_is_fnan (resid[i]))continue;
		weight[j] = (float)fabs((double)resid[i]);
		j++;
	}

	qsort ((void *)weight, (size_t)j, sizeof(float), GMT_comp_float_asc);

	j2 = j / 2;
	if (j%2) {
		mad = weight[j2];
	}
	else {
		mad = 0.5 *(weight[j2] + weight[j2 - 1]);
	}

	/* Adjust mad to equal Gaussian sigma:  */

	*scale = 1.4826 * mad;

	/* Use weight according to Huber (1981), but squared:  */

	for (i = 0; i < nxy; i++) {
		if (GMT_is_fnan (resid[i])) {
			weight[i] = resid[i];
			continue;
		}
		r = fabs(resid[i]) / (*scale);
		
		weight[i] = (float)((r <= 1.5) ? 1.0 : (3.0 - 2.25/r) / r);
	}
	return;
}

void	write_model_parameters(double *gtd, int n_model)
{
	int	i;
	char	pbasis[10][16];

	sprintf(pbasis[0], "Mean");
	sprintf(pbasis[1], "X");
	sprintf(pbasis[2], "Y");
	sprintf(pbasis[3], "X*Y");
	sprintf(pbasis[4], "P2(x)");
	sprintf(pbasis[5], "P2(y)");
	sprintf(pbasis[6], "P3(x)");
	sprintf(pbasis[7], "P2(x)*P1(y)");
	sprintf(pbasis[8], "P1(x)*P2(y)");
	sprintf(pbasis[9], "P3(y)");

	sprintf(format, "Coefficient fit to %%s:  %s\n", gmtdefs.d_format);
	for(i = 0; i < n_model; i++) mexPrintf( format, pbasis[i], gtd[i]);

	return;
}

void	load_gtg_and_gtd(float *data, int nx, int ny, double *xval, double *yval, double *pstuff, double *gtg, double *gtd, int n_model, float *weight, int weighted)
{
	/* Routine to load the matrix G'G (gtg) and vector G'd (gtd)
	for the normal equations.  Routine uses indices i,j to refer
	to the grdfile of data, and k,l to refer to the k_row, l_col
	of the normal equations matrix.  We need sums of [weighted]
	data and model functions in gtg and gtd.  We save time by
	loading only lower triangular part of gtg and then filling
	by symmetry after i,j loop.  */

	int	i, j, ij, k, l, n_used;

/*	First zero things out to start:  */

	n_used = 0;
	for (k = 0; k < n_model; k++) {
		gtd[k] = 0.0;
		for (l = 0; l < n_model; l++) {
			gtg[k*n_model+l] = 0.0;
		}
	}

/*  Now get going.  Have to load_pstuff separately in i and j,
	because it is possible that we skip data when i = 0.
	Loop over all data:  */

	for (ij = 0, j = 0; j < ny; j++ ) {
		load_pstuff(pstuff, n_model, xval[0], yval[j], 0, 1);
		for (i = 0; i < nx; i++, ij++) {

			if (GMT_is_fnan (data[ij]))continue;

			n_used++;
			load_pstuff(pstuff, n_model, xval[i], yval[j], 1, 0);

/* If weighted  */	if (weighted) {
				/* Loop over all gtg and gtd elements:  */
				gtd[0] += (data[ij] * weight[ij]);
				gtg[0] += (weight[ij]);
				for (k = 1; k < n_model; k++) {
					gtd[k] += (data[ij] * weight[ij] * pstuff[k]);
					gtg[k] += (weight[ij] * pstuff[k]);
					for (l = k; l < n_model; l++) {
						gtg[k + l*n_model] += (pstuff[k]*pstuff[l]*weight[ij]);
					}
				}
			}
/* If !weighted  */	else {
				/* Loop over all gtg and gtd elements:  */
				gtd[0] += data[ij];
				for (k = 1; k < n_model; k++) {
					gtd[k] += (data[ij] * pstuff[k]);
					gtg[k] += pstuff[k];
					for (l = k; l < n_model; l++) {
						gtg[k + l*n_model] += (pstuff[k]*pstuff[l]);
					}
				}
/* End if  */		}
		}
	}
/* End of loop over data i,j  */

/* Now if !weighted, use more accurate sum for gtg[0], and set symmetry:  */

	if (!(weighted)) gtg[0] = n_used;

	for (k = 0; k < n_model; k++) {
		for (l = 0; l < k; l++) {
			gtg[l + k*n_model] = gtg[k + l*n_model];
		}
	}
/* That is all there is to it!  */

	return;
}
				
void gauss (double *a, double *vec, int n, int nstore, double test, int *ierror, int itriag)
{
 
/* subroutine gauss, by william menke */
/* july 1978 (modified feb 1983, nov 85) */
 
/* a subroutine to solve a system of n linear equations in n unknowns*/
/* where n doesn't exceed MAX_TABLE_COLS */
/* gaussian reduction with partial pivoting is used */
/*      a               (sent, destroyed)       n by n matrix           */
/*      vec             (sent, overwritten)     n vector, replaced w/ solution*/
/*      nstore          (sent)                  dimension of a  */
/*      test            (sent)                  div by zero check number*/
/*      ierror          (returned)              zero on no error*/
/*      itriag          (sent)                  matrix triangularized only*/
/*                                               on TRUE useful when solving*/
/*                                               multiple systems with same a */
        static int isub[MAX_TABLE_COLS], l1;
        int line[MAX_TABLE_COLS], iet, ieb, i, j, k, l, j2;
        double big, testa, b, sum;
        

        iet=0;  /* initial error flags, one for triagularization*/
        ieb=0;  /* one for backsolving */

/* triangularize the matrix a*/
/* replacing the zero elements of the triangularized matrix */
/* with the coefficients needed to transform the vector vec */

        if (itriag) {   /* triangularize matrix */
 
                for( j=0; j<n; j++ ) {      /*line is an array of flags*/
                        line[j]=0; 
                        /* elements of a are not moved during pivoting*/
                        /* line=0 flags unused lines */
                        }    /*end for j*/
                        
                for( j=0; j<n-1; j++ ) {
                        /*  triangularize matrix by partial pivoting */
                       big = 0.0; /* find biggest element in j-th column*/
                                  /* of unused portion of matrix*/
                       for( l1=0; l1<n; l1++ ) {
                               if( line[l1]==0 ) {
                                       testa=(double) fabs(
                                                (double) (*(a+l1*nstore+j)) );
                                       if (testa>big) {
                                                i=l1;
                                                big=testa;
                                                } /*end if*/
                                        } /*end if*/
                                } /*end for l1*/
                       if( big<=test) {   /* test for div by 0 */
                               iet=1;
                               } /*end if*/
 
                       line[i]=1;  /* selected unused line becomes used line */
                       isub[j]=i;  /* isub points to j-th row of tri. matrix */
 
                       sum=1.0/(*(a+i*nstore+j)); 
                                /*reduce matrix towards triangle */
                       for( k=0; k<n; k++ ) {
                                if( line[k]==0 ) {
                                        b=(*(a+k*nstore+j))*sum;
                                        for( l=j+1; l<n; l++ ) {
                                               *(a+k*nstore+l)=
                                                        (*(a+k*nstore+l))
                                                        -b*(*(a+i*nstore+l));
                                               } /*end for l*/
                                       *(a+k*nstore+j)=b;
                                        } /*end if*/
                                } /*end for k*/
                        } /*end for j*/
 
               for( j=0; j<n; j++ ) {
                        /*find last unused row and set its pointer*/
                        /*  this row contains the apex of the triangle*/
                        if( line[j]==0) {
                                l1=j;   /*apex of triangle*/
                                isub[n-1]=j;
                                break;
                                } /*end if*/
                        } /*end for j*/
 
                } /*end if itriag true*/
                
        /*start backsolving*/
        
        for( i=0; i<n; i++ ) {  /* invert pointers. line(i) now gives*/
                                /* row no in triang matrix of i-th row*/
                                /* of actual matrix */
                line[isub[i]] = i;
                } /*end for i*/
 
        for( j=0; j<n-1; j++) { /*transform the vector to match triang. matrix*/
               b=vec[isub[j]];
               for( k=0; k<n; k++ ) {
                      if (line[k]>j) {  /* skip elements outside of triangle*/
                                vec[k]=vec[k]-(*(a+k*nstore+j))*b;
                                } /*end if*/
                        } /*end for k*/
                } /*end for j*/
 
      b = *(a+l1*nstore+(n-1));   /*apex of triangle*/
      if( ((double)fabs( (double) b))<=test) {
                /*check for div by zero in backsolving*/
                ieb=2;
                } /*end if*/
      vec[isub[n-1]]=vec[isub[n-1]]/b;
 
      for( j=n-2; j>=0; j-- ) { /* backsolve rest of triangle*/
                sum=vec[isub[j]];
                for( j2=j+1; j2<n; j2++ ) {
                        sum = sum - vec[isub[j2]] * (*(a+isub[j]*nstore+j2));
                        } /*end for j2*/
                        b = *(a+isub[j]*nstore+j);
               if( ((double)fabs((double)b))<=test) {
                        /* test for div by 0 in backsolving */
                        ieb=2;
                        } /*end if*/
                vec[isub[j]]=sum/b;   /*solution returned in vec*/
                } /*end for j*/

/*put the solution vector into the proper order*/

      for( i=0; i<n; i++ ) {    /* reorder solution */
                for( k=i; k<n; k++ ) {  /* search for i-th solution element */
                        if( line[k]==i ) {
                                j=k;
                                break;
                                } /*end if*/
                        } /*end for k*/
               b = vec[j];       /* swap solution and pointer elements*/
               vec[j] = vec[i];
               vec[i] = b;
               line[j] = line[i];
                } /*end for i*/
 
      *ierror = iet + ieb;   /* set final error flag*/
}
