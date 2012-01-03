/*--------------------------------------------------------------------
 *	$Id: trend1d.c,v 1.39 2008/03/24 08:58:32 guru Exp $
 *
 *	Copyright (c) 1991-2008 by P. Wessel and W. H. F. Smith
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
 * trend1d [<xy[w]file>] -F<output_flags> -N[f]<n_m_parameters>[r] 
 *	[-C<condition_#>] [-I[<confid>]] [-V] [-W]
 *
 * where:
 *	[<xy[w]file>] is an ascii file with x y in first 2 columns [or
 *		x y w in first 3 columns].  Default reads from GMT_stdin.
 *	-F<output_flags> is a string of at least one, up to five, in
 *		and order, from the set {x y m r w}.  x,y = input,
 *		m = model, r = residual = y-m, and w= weight used.
 *	-N[f]<n_m_parameters>[r]
 *		If iterative Robust fitting desired, use append r.
 *		To fit a Fourier model, use -Nf.
 *		Number of terms in the model is <n_m_parameters>.
 *		Example:  Robust quadratic polynomial:  -N2r.
 *	[-C<condition_#>] Cut off eigenvalue spectrum; use only eigen-
 *		values such that (lambda_max / lambda[i]) < condition_#.
 *	[-I[<confid>]] Iteratively Increment the number of model parameters,
 *		searching for the significant model size, up to a maximum
 *		set by <n_m_parameters>.  We start with a 1 parameter
 *		model and then iteratively increase the number of
 *		model parameters, m, while m <= <n_m_parameters> &&
 *		reduction in variance from i to i+1 is significant
 *		at the <confid> level according to F test.  If user sets
 *		-I without giving <confid> then <confid> = 0.95.
 *	[-V]	Verbose operation.
 *	[-W]	Weighted data are input.  Read 3 cols and use 3rd as weight.
 *
 *
 * Read GMT_stdin or file of x y pairs, or weighted pairs as x,y w data.  Fit 
 * a regression model y = f(x) + e, where e are error misfits and f(x) has
 * some user-prescribed functional form.  Presently available models are
 * polynomials and Fourier series.  The user may choose the number of terms
 * in the model to fit, whether to seek iterative refinement robust w.r.t.
 * outliers, and whether to seek automatic discovery of the significant
 * number of model parameters.
 *
 *
 * In trend1d I chose to construct the polynomial model using Chebyshev 
 * Polynomials so that the user may easily compare the sizes of the
 * coefficients (and compare with a Fourier series as well).  Tn(x)
 * is an n-degree polynomial with n zero-crossings in [-1,1] and n+1
 * extrema, at which the value of Tn(x) is +/- 1.  It is this property
 * which makes it easy to compare the size of the coefficients.
 *
 * During model fitting the data x coordinate is normalized into the domain
 * [-1, 1] for Chebyshev Polynomial fitting, or into the domain [-pi, pi]
 * for Fourier series fitting.  Before writing out the data the coordinate
 * is rescaled to match the original input values.
 *
 * An n degree polynomial can be written with terms of the form a0 + a1*x
 * + a2*x*x + ...  But it can also be written using other polynomial 
 * basis functions, such as a0*P0 + a1*P1 + a2*P2..., the Legendre
 * polynomials, and a0*T0 + a1*T1 + a2*T2..., the Chebyshev polynomials.
 * (The domain of the x values has to be in [-1, 1] in order to use P or T.)
 * It is well known that the ordinary polynomial basis 1, x, x*x, ... gives 
 * terribly ill- conditioned matrices.  The Ps and Ts do much better.
 * This is because the ordinary basis is far from orthogonal.  The Ps
 * are orthogonal on [-1,1] and the Ts are orthogonal on [-1,1] under a
 * simple weight function.
 * Because the Ps have ordinary orthogonality on [-1,1], I expected them
 * to be the best basis for a regression model; best meaning that they 
 * would lead to the most balanced G'G (matrix of normal equations) with
 * the smallest condition number and the most nearly diagonal model
 * parameter covariance matrix ((G'G)inverse).  It turns out, however, that
 * the G'G obtained from the Ts is very similar and usually has a smaller 
 * condition number than the Ps G'G.  Both of these are vastly superior to
 * the usual polynomials 1, x, x*x.  In a test with 1000 equally spaced
 * data and 8 model parameters, the Chebyshev system had a condition # = 10.6,
 * Legendre = 14.8, and traditional = 54722.7.  For 1000 randomly spaced data
 * and 8 model parameters, the results were C = 13.1, L = 15.6, and P = 54916.6.
 * As the number of model parameters approaches the number of data, the 
 * situation still holds, although all matrices get ill-conditioned; for 8 
 * random data and 8 model parameters, C = 1.8e+05, L = 2.6e+05, P = 1.1e+08.
 * I expected the Legendre polynomials to have a covariance matrix more nearly 
 * diagonal than that of the Chebyshev polynomials, but on this criterion also
 * the Chebyshev turned out to do better.  Only as ndata -> n_model_parameters
 * does the Legendre covariance matrix do better than the Chebyshev.   So for
 * all these reasons I use Chebyshev polynomials.
 *
 * Author:	W. H. F. Smith
 * Date:	25 February 1991-2000.
 * Revised:	11 June, 1991-2000 for v2.0 of GMT-SYSTEM.
 *		13-JUN-1998, for GMT 3.1 (PW)
 *		13-JUL-2000, for GMT 3.3.5 (PW)
 *		10-MAY-2001, PW: Use Numerical Recipes scheme to also output polynomial coefficients
 * Version:	3.4 18-APR-2001
 * Version:	4.1.x
 *
 *
 *	MEXIFIED by Joaquim Luis 09-Jul-2008
 *
 *			Added -L option to output linear trend parameters (slope & intercept)
 *	12-12-2008	Extended -L option  (not yet documented)
 *			Added -X option to output the Chi
 *	31-12-2008	Added -R option to output the R-squared (only linear)
 *	07-01-2009	Added -P option to output the statistical significance of R
 *
 *		NOTE: this MEX does not neeed to link against the GMT library
 */

#include "mex.h"
#include <string.h>
#include <math.h>

#define TREND1D_N_OUTPUT_CHOICES 5
#define	FALSE	0
#define	TRUE	1
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define d_log(x) ((x) <= 0.0 ? mxGetNaN() : log (x))

struct TREND1D_CTRL {
	struct C {	/* -C<condition_#> */
		int active;
		double value;
	} C;
	struct F {	/* -F<xymrw> */
		int active;
		char col[TREND1D_N_OUTPUT_CHOICES];	/* Character codes for desired output in the right order */
	} F;
	struct I {	/* -I[<confidence>] */
		int active;
		double value;
	} I;
	struct L {	/* -L When -N2[r] output slope & intercept*/
		int active;
		int chi;
		int mad;
		int Rsq;
		double value;
	} L;
	struct N {	/* -N[f]<n_model>[r] */
		int active;
		int robust;
		int mode;
		int value;
	} N;
	struct P {	/* -P */
		int active;
		double value;
	} P;
	struct R {	/* -R */
		int active;
	} R;
	struct X {	/* -X */
		int active;
	} X;
	struct W {	/* -W */
		int active;
	} W;
};

#define TREND1D_POLYNOMIAL 0
#define TREND1D_FOURIER 1
#define TREND1D_ANNUAL 2

struct	TREND1D_DATA {
	double	x;
	double	y;
	double	m;
	double	r;
	double	w;
};

int	GMT_ln_gamma_r(double x, double *lngam);
int	GMT_comp_double_asc (const void *p_1, const void *p_2);
int	GMT_jacobi (double *a, int *n, int *m, double *d, double *v, double *b, double *z, int *nrots);
int	GMT_inc_beta (double a, double b, double x, double *ibeta);
int     GMT_f_q (double chisq1, int nu1, double chisq2, int nu2, double *prob);
int	GMT_f_test_new (double chisq1, int nu1, double chisq2, int nu2, double *prob, int iside);
double	GMT_ln_gamma (double xx);
double	GMT_cf_beta (double a, double b, double x);
double	tpvalue (double x, double v);
void	GMT_cheb_to_pol (double c[], int n, double a, double b);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int i, j, k,n_outputs, n_model, significant, rank;
	
	int	n_data, np, n_fields;
	int	error = FALSE, weighted_output = FALSE;

	double	*gtg, *v, *gtd, *lambda, *workb, *workz, *c_model, *o_model, *w_model, *work, *in_m, *pdata;
	double	xmin, xmax, c_chisq, o_chisq, w_chisq, scale = 1.0, prob, Rsq, p = 0;
	double	get_chisq(struct TREND1D_DATA *data, int n_data, int n_model);

	int	argc = 0, n_arg_no_char = 0;
	char	**argv;

	struct	TREND1D_DATA *data;
	struct	TREND1D_CTRL *Ctrl;

	void read_data (struct TREND1D_DATA **data, double *in_m, int n_pts, double *xmin, double *xmax, int weighted_input, double **work);
	void write_output (struct TREND1D_DATA *data, int n_data, char *output_choice, int n_outputs, double *pdata);
	void transform_x (struct TREND1D_DATA *data, int n_data, int model_type, double xmin, double xmax);
	void untransform_x(struct TREND1D_DATA *data, int n_data, int model_type, double xmin, double xmax);
	void recompute_weights(struct TREND1D_DATA *data, int n_data, double *work, double *scale);
	void allocate_array_space(int np, double **gtg, double **v, double **gtd, double **lambda, double **workb, double **workz, double **c_model, double **o_model, double **w_model);
	void free_the_memory(double *gtg, double *v, double *gtd, double *lambda, double *workb, double *workz, double *c_model, double *o_model, double *w_model, struct TREND1D_DATA *data, double *work);
	void calc_m_and_r(struct TREND1D_DATA *data, int n_data, double *model, int n_model, int m_type, double *grow);
	void move_model_a_to_b(double *model_a, double *model_b, int n_model, double *chisq_a, double *chisq_b);
	void load_gtg_and_gtd(struct TREND1D_DATA *data, int n_data, double *gtg, double *gtd, double *grow, int n_model, int mp, int m_type);
	void solve_system(double *gtg, double *gtd, double *model, int n_model, int mp, double *lambda, double *v, double *b, double *z, double c_no, int *ir);
	int	GMT_sig_f (double chi1, int n1, double chi2, int n2, double level, double *prob);
	void	*New_Trend1d_Ctrl (), Free_Trend1d_Ctrl (struct TREND1D_CTRL *C);

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
	argv[0] = "trend1d";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	
	Ctrl = (struct TREND1D_CTRL *)New_Trend1d_Ctrl ();	/* Allocate and initialize a new control structure */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.value = atof(&argv[i][2]);
					break;
				case 'F':
					Ctrl->F.active = TRUE;
					for (j = 2, k = 0; argv[i][j]; j++, k++) {
						if (k < TREND1D_N_OUTPUT_CHOICES)
							Ctrl->F.col[k] = argv[i][j];
						else {
							error++;
							mexPrintf ("TREND1D: SYNTAX ERROR -F option: Too many output columns selected\n");
							mexPrintf ("TREND1D: SYNTAX ERROR -F option: Choose from -Fxymrw\n");
						}
					}
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					Ctrl->I.value = (argv[i][2]) ? atof(&argv[i][2]) : 0.51;
					break;
				case 'L':
					Ctrl->L.active = TRUE;
					if (argv[i][2] == 'x')		Ctrl->L.chi = TRUE;
					else if (argv[i][2] == 'm')	Ctrl->L.mad = TRUE;
					else if (argv[i][2] == 'r')	Ctrl->L.Rsq = TRUE;
					Ctrl->L.value = 1.5;		/* Default value, in case it is needed but not set */ 
					if (argv[i][3])
						Ctrl->L.value = atof(&argv[i][3]);
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					if (strchr (argv[i], 'r')) Ctrl->N.robust = TRUE;
					j = (argv[i][2] == 'r') ? 3 : 2;
					if (argv[i][j] == 'F' || argv[i][j] == 'f') {
						Ctrl->N.mode = TREND1D_FOURIER;
						j++;
					}
					else if (argv[i][j] == 'P' || argv[i][j] == 'p') {
						Ctrl->N.mode = TREND1D_POLYNOMIAL;
						j++;
					}
					else if (argv[i][j] == 'A' || argv[i][j] == 'a') {	/* Non doc. It was an experiment */
						Ctrl->N.mode = TREND1D_ANNUAL;
						j++;
					}
					if (argv[i][j])
						Ctrl->N.value = atoi(&argv[i][j]);
					else {
						error = TRUE;
 						mexPrintf ("%TREND1D: GMT SYNTAX ERROR -N option.  No model specified\n");
					}
					if (Ctrl->N.mode == TREND1D_ANNUAL)
						Ctrl->N.value = Ctrl->N.value * 2 + 2;
					break;
				case 'P':
					Ctrl->P.active = TRUE;
					Ctrl->P.value = 1.0;		/* Default value, in case it is needed but not set */ 
					if (argv[i][3])
						Ctrl->P.value = atof(&argv[i][3]);
					break;
				case 'R':
					Ctrl->R.active = TRUE;
					break;
				case 'X':
					Ctrl->X.active = TRUE;
					break;
				case 'W':
					Ctrl->W.active = TRUE;
					break;
				default:
					error = TRUE;
					mexPrintf("trend1d SYNTAX ERROR. %s\n", argv[i][1]);
					break;
			}
		}
	}

	if (argc == 1 || error) {
		mexPrintf("trend1d - Fit a [weighted] [robust] polynomial [or Fourier] model for y = f(x) to xy[w]\n\n");
		mexPrintf("usage:  [out, lin_params] = trend1d_m(in_array, '-F<xymrw>', '-N[f]<n_model>[r]', '[-C<condition_#>],'\n");
		mexPrintf("\t'[-I[<confidence>]]', '[-L]', '[-P]', '[-R]', '[-X]', '[-W])'\n\n");

		mexPrintf("\t-F Choose at least 1, up to 5, any order, of xymrw for ascii output to stdout.\n");
		mexPrintf("\t   x=x, y=y, m=model, r=residual=y-m, w=weight.  w determined iteratively if robust fit used.\n");
		mexPrintf("\t-N fit a Polynomial [Default] or Fourier (-Nf) model with <n_model> terms.\n");
		mexPrintf("\t   Append r for robust model.  E.g., robust quadratic = -N3r.\n");
		mexPrintf("\n\tOPTIONS:\n");
		mexPrintf("\tin_array <xy[w]>, first 2 cols = x y [3 cols = x y w].\n");
		mexPrintf("\t-C Truncate eigenvalue spectrum so matrix has <condition_#>.  [Default = 1.0e06].\n");
		mexPrintf("\t-I Iteratively Increase # model parameters, to a max of <n_model> so long as the\n");
		mexPrintf("\t   reduction in variance is significant at the <confidence> level.\n");
		mexPrintf("\t   Give -I without a number to default to 0.51 confidence level.\n");
		mexPrintf("\t-L When the fit is linear (-N2[r]) outputs slope & intercept.\n");
		mexPrintf("\t   If only one output is requested out = [slope, intercept].\n");
		mexPrintf("\t   If two outputs, first is determined by -F and the second is [slope, intercept].\n");
		mexPrintf("\t   Note, however, that by asking two outputs it is not necessary to use -L option.\n");
		mexPrintf("\t-P For linear (-N2[r]) Forces -L and lin_prams is [m b R^2 p] vector.\n");
		mexPrintf("\t-R For linear (-N2[r]) fits outputs the R-square. It forces -L and lin_prams is a 1x3 vector.\n");
		mexPrintf("\t-X For linear (-N2[r]) fits outputs SQRT(Chi-square). It forces -L and lin_prams is a 1x3 vector.\n");
		mexPrintf("\t-W Weighted input given, weights in 3rd column.  [Default is unweighted].\n");
		mexPrintf("\t   Default is 2 (or 3 if -W is set) input columns.\n");
		return;
	}

	if (Ctrl->C.value <= 1.0) {
		mexPrintf ("TREND1D: GMT SYNTAX ERROR -C option.  Condition number must be larger than unity\n");
		error++;
	}
	if (Ctrl->I.value < 0.0 || Ctrl->I.value > 1.0) {
		mexPrintf ("TREND1D: GMT SYNTAX ERROR -C option.  Give 0 < confidence level < 1.0\n");
		error++;
	}
	if (Ctrl->N.value <= 0.0) {
		mexPrintf ("TREND1D: GMT SYNTAX ERROR -N option.  A positive number of terms must be specified\n");
		error++;
	}
	for (k = n_outputs = 0; k < TREND1D_N_OUTPUT_CHOICES && Ctrl->F.col[k]; k++) {
		if (!strchr ("xymrw", Ctrl->F.col[k])) {
			mexPrintf ("TREND1D: SYNTAX ERROR -F option.  Unrecognized output choice %c\n", Ctrl->F.col[k]);
			error++;
		}
		else if (Ctrl->F.col[k] == 'w')
			weighted_output = TRUE;

		n_outputs++;
	}
        if (Ctrl->X.active)
		Ctrl->L.active = TRUE;

        if (n_outputs == 0 && !(Ctrl->L.active && (nlhs == 1))) {
                mexPrintf ("TREND1D: SYNTAX ERROR -F option.  Must specify at least one output columns \n");
                error++;
        }

	if (error) return;

	if (nlhs == 0)
		mexErrMsgTxt("TREND1D ERROR: Must provide an output.\n");


	/* Check that first argument contains at least a mx2 table */
	n_data = mxGetM (prhs[0]);
	n_fields = mxGetN(prhs[0]);
	if (!mxIsNumeric(prhs[0]) || (n_fields < 2))
		mexErrMsgTxt("TREND1D ERROR: first argument must contain at least a mx2 table with x,y[w].\n");

	if (Ctrl->W.active && n_fields == 2)
		mexPrintf("TREND1D WARNING: weighted fit asked but data has only 2 columns. Ignoring option -W\n");

	if (Ctrl->L.active && Ctrl->N.value != 2 && Ctrl->N.mode != TREND1D_ANNUAL) {
		mexPrintf("TREND1D WARNING: ouput model parameters is only possible with linear trend. Ignoring option -L\n");
		Ctrl->L.active = FALSE;
		Ctrl->X.active = FALSE;
	}

	if (nlhs == 2 && Ctrl->N.value == 2 && !Ctrl->L.active)		/* Two outputs and -N2 == -L */
		Ctrl->L.active = TRUE;

	/* Read the input points and convert them to double */
	in_m = (double *)mxGetData(prhs[0]);

	np = Ctrl->N.value;	/* Row dimension for matrices gtg and v  */
	allocate_array_space(np, &gtg, &v, &gtd, &lambda, &workb, &workz, &c_model, &o_model, &w_model);

	read_data(&data, in_m, n_data, &xmin, &xmax, Ctrl->W.active, &work);

	if (xmin == xmax)
		mexErrMsgTxt("TREND1D: Fatal error in input data.  X min = X max.\n");

	if (n_data < Ctrl->N.value) mexPrintf("TREND1D: Warning. Ill-posed problem.  n_data < n_model_max.\n");

	transform_x(data, n_data, Ctrl->N.mode, xmin, xmax);	/* Set domain to [-1, 1] or [-pi, pi]  */

	if (Ctrl->I.active) {
		n_model = 1;

		/* Fit first model  */
		load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
		solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
		calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
		c_chisq = get_chisq(data, n_data, n_model);
		if (Ctrl->N.robust) {
			do {
				recompute_weights(data, n_data, work, &scale);
				move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
				load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
				solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				c_chisq = get_chisq(data, n_data, n_model);
				significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
			} while (significant);
			/* Go back to previous model only if w_chisq < c_chisq  */
			if (w_chisq < c_chisq) {
				move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
			}
		}
		/* First [robust] model has been found  */

		significant = TRUE;
		while(n_model < Ctrl->N.value && significant) {
			move_model_a_to_b(c_model, o_model, n_model, &c_chisq, &o_chisq);
			n_model++;

			/* Fit next model  */
			load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
			solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
			calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
			c_chisq = get_chisq(data, n_data, n_model);
			if (Ctrl->N.robust) {
				do {
					recompute_weights(data, n_data, work, &scale);
					move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
					load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
					solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
					calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
					c_chisq = get_chisq(data, n_data, n_model);
					significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
				} while (significant);
				/* Go back to previous model only if w_chisq < c_chisq  */
				if (w_chisq < c_chisq) {
					move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
					calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
					if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
				}
			}
			/* Next [robust] model has been found  */
			significant = GMT_sig_f(c_chisq, n_data-n_model, o_chisq, n_data-n_model-1, Ctrl->I.value, &prob);
		}

		if (!(significant) ) {	/* Go back to previous [robust] model, stored in o_model  */
			n_model--;
			rank--;
			move_model_a_to_b(o_model, c_model, n_model, &o_chisq, &c_chisq);
			calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
			if (Ctrl->N.robust && weighted_output) recompute_weights(data, n_data, work, &scale);
		}
	}
	else {
		n_model = Ctrl->N.value;
		load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
		solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
		calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
		c_chisq = get_chisq(data, n_data, n_model);
		if (Ctrl->N.robust) {
			do {
				recompute_weights(data, n_data, work, &scale);
				move_model_a_to_b(c_model, w_model, n_model, &c_chisq, &w_chisq);
				load_gtg_and_gtd(data, n_data, gtg, gtd, workb, n_model, np, Ctrl->N.mode);
				solve_system(gtg, gtd, c_model, n_model, np, lambda, v, workb, workz, Ctrl->C.value, &rank);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				c_chisq = get_chisq(data, n_data, n_model);
				significant = GMT_sig_f(c_chisq, n_data-n_model, w_chisq, n_data-n_model, Ctrl->I.value, &prob);
			} while (significant);
			/* Go back to previous model only if w_chisq < c_chisq  */
			if (w_chisq < c_chisq) {
				move_model_a_to_b(w_model, c_model, n_model, &w_chisq, &c_chisq);
				calc_m_and_r(data, n_data, c_model, n_model, Ctrl->N.mode, workb);
				if (weighted_output && n_model == Ctrl->N.value) recompute_weights(data, n_data, work, &scale);
			}
		}
	}

	untransform_x(data, n_data, Ctrl->N.mode, xmin, xmax);

	/* -P, -R and -X options may be quite inconsistent still */
	if (Ctrl->R.active || Ctrl->P.active) {
		double SSR = 0, SST = 0, ymean = 0;
		for (i = 0; i < n_data; i++)
			ymean += data[i].y;

		ymean /= n_data;
		for (i = 0; i < n_data; i++){
			SST += ((data[i].y - ymean) * (data[i].y - ymean));
			SSR += ((data[i].m - ymean) * (data[i].m - ymean));
		}

		Rsq = SSR / SST;
		if (Ctrl->P.active) {
			p = 2 * tpvalue(-fabs(sqrt(Rsq * (n_data-2) / (1 - Rsq))),n_data-2);
		}
	}

	if (Ctrl->L.active || Ctrl->X.active || Ctrl->R.active || Ctrl->P.active) {	/* Outputs also slope and intercept */
		int	nLX;
		double	m, b, nan;

		nan = mxGetNaN();
		m = (data[n_data - 1].m - data[0].m) / (data[n_data - 1].x - data[0].x);
		b = -m * data[0].x + data[0].m;

		if (Ctrl->L.chi && (sqrt(c_chisq) > Ctrl->L.value)) {
			m = nan;	b = nan;
		}
		else if (Ctrl->L.mad && (scale > Ctrl->L.value)) {
			m = nan;	b = nan;
		}
		else if (Ctrl->L.Rsq && (scale > Ctrl->L.value)) {
			m = nan;	b = nan;
		}
		nLX = ((Ctrl->X.active || Ctrl->R.active) ? 3 : 2);
		if (Ctrl->P.active) nLX++;
		if (nlhs == 1) {
			plhs[0] = mxCreateDoubleMatrix (1, nLX, mxREAL);
			pdata = mxGetPr(plhs[0]);
			pdata[0] = m;	pdata[1] = b;
			if (Ctrl->X.active) pdata[2] = sqrt(c_chisq);
			if (Ctrl->R.active) pdata[2] = Rsq;
			if (Ctrl->P.active) pdata[3] = p;
		}
		else {		/* Output whatever -F plus slope and intercept */
			double *pdata_l;
			plhs[0] = mxCreateDoubleMatrix (n_data, n_outputs, mxREAL);
			pdata = mxGetPr(plhs[0]);
			write_output(data, n_data, Ctrl->F.col, n_outputs, pdata);

			plhs[1] = mxCreateDoubleMatrix (1, nLX, mxREAL);
			pdata_l = mxGetPr(plhs[1]);
			pdata_l[0] = m;	pdata_l[1] = b;
			if (Ctrl->X.active) pdata_l[2] = sqrt(c_chisq);
			if (Ctrl->R.active) pdata_l[2] = Rsq;
			if (Ctrl->P.active) pdata_l[3] = p;
		}
	}
	else {
		if (Ctrl->F.active) {
			plhs[0] = mxCreateDoubleMatrix (n_data, n_outputs, mxREAL);
			pdata = mxGetPr(plhs[0]);
			write_output(data, n_data, Ctrl->F.col, n_outputs, pdata);
			if (nlhs == 2) {
				double *pdata_c;
				GMT_cheb_to_pol (c_model, n_model, xmin, xmax);
				plhs[1] = mxCreateDoubleMatrix (1, n_model, mxREAL);
				pdata_c = mxGetPr(plhs[1]);
				for (i = 0, j = n_model-1; i < n_model; i++, j--)
					pdata_c[i] = c_model[j];	/* Output from higher to lower order */
			}
		}
		else {
			double *pdata_c;
			GMT_cheb_to_pol (c_model, n_model, xmin, xmax);
			plhs[0] = mxCreateDoubleMatrix (1, n_model, mxREAL);
			pdata_c = mxGetPr(plhs[1]);
			for (i = 0, j = n_model-1; i < n_model; i++, j--)
				pdata_c[i] = c_model[j];
		}
	}

	free_the_memory(gtg, v, gtd, lambda, workb, workz, c_model, o_model, w_model, data, work);

	Free_Trend1d_Ctrl (Ctrl);	/* Deallocate control structure */
}

void read_data (struct TREND1D_DATA **data, double *in_m, int n_pts, double *xmin, double *xmax, int weighted_input, double **work) {
	/* Copy the vectors transmited as inputs into the data structure. For large inputs,
	this implies a significant memory wasting.*/
	int n;

	(*data) = (struct TREND1D_DATA *) mxCalloc ((size_t)n_pts, sizeof(struct TREND1D_DATA));

	*xmin = in_m[0];
	*xmax = in_m[0];

	for (n = 0; n < n_pts; n++) {

		(*data)[n].x = in_m[n];
		(*data)[n].y = in_m[n + n_pts];
		(*data)[n].w = (weighted_input) ? in_m[n + 2*n_pts] : 1.0;

		if (*xmin > (*data)[n].x) *xmin = (*data)[n].x;
		if (*xmax < (*data)[n].x) *xmax = (*data)[n].x;
	}

	*work = (double *) mxCalloc ((size_t)n_pts, sizeof(double));
}

void write_output (struct TREND1D_DATA *data, int n_data, char *output_choice, int n_outputs, double *pdata) {
	int i, j;

	for (i = 0; i < n_data; i++) {
		for (j = 0; j < n_outputs; j++) {
			switch (output_choice[j]) {
				case 'x':
					pdata[i + j*n_data] = data[i].x;
					break;
				case 'y':
					pdata[i + j*n_data] = data[i].y;
					break;
				case 'm':
					pdata[i + j*n_data] = data[i].m;
					break;
				case 'r':
					pdata[i + j*n_data] = data[i].r;
					break;
				case 'w':
					pdata[i + j*n_data] = data[i].w;
					break;
			}
		}
	}

	/*m = (data[n_data - 1].m - data[0].m) / (data[n_data - 1].x - data[0].x);
	b = -m * data[0].x + data[0].y;

	y - y0 = m(x-x0) = mx - m.x0
	y = mx - m.x0 + y0
	y = -m.x0 + y0*/
}

void allocate_array_space (int np, double **gtg, double **v, double **gtd, double **lambda, double **workb, double **workz, double **c_model, double **o_model, double **w_model) {

	*gtg = (double *) mxCalloc ((size_t)(np*np), sizeof(double));
	*v = (double *) mxCalloc ((size_t)(np*np), sizeof(double));
	*gtd = (double *) mxCalloc ((size_t)np, sizeof(double));
	*lambda = (double *) mxCalloc ((size_t)np, sizeof(double));
	*workb = (double *) mxCalloc ((size_t)np, sizeof(double));
	*workz = (double *) mxCalloc ((size_t)np, sizeof(double));
	*c_model = (double *) mxCalloc ((size_t)np, sizeof(double));
	*o_model = (double *) mxCalloc ((size_t)np, sizeof(double));
	*w_model = (double *) mxCalloc ((size_t)np, sizeof(double));
}

void free_the_memory (double *gtg, double *v, double *gtd, double *lambda, double *workb, double *workz, double *c_model, double *o_model, double *w_model, struct TREND1D_DATA *data, double *work) {

	mxFree ((void *)w_model);
	mxFree ((void *)o_model);
	mxFree ((void *)c_model);
	mxFree ((void *)workz);
	mxFree ((void *)workb);
	mxFree ((void *)lambda);
	mxFree ((void *)gtd);
	mxFree ((void *)v);
	mxFree ((void *)gtg);
	mxFree ((void *)work);
	mxFree ((void *)data);
}

void transform_x (struct TREND1D_DATA *data, int n_data, int model_type, double xmin, double xmax) {
	int	i;
	double	offset, scale;

	offset = 0.5 * (xmin + xmax);	/* Mid Range  */
	scale = 2.0 / (xmax - xmin);	/* 1 / (1/2 Range)  */

	if (model_type == TREND1D_FOURIER || model_type == TREND1D_ANNUAL) 	/* Set Range to 1 period  */
		scale *= M_PI;

	for (i = 0; i < n_data; i++)
		data[i].x = (data[i].x - offset) * scale;
}

void untransform_x (struct TREND1D_DATA *data, int n_data, int model_type, double xmin, double xmax) {
	int	i;
	double	offset, scale;

	offset = 0.5 * (xmin + xmax);	/* Mid Range  */
	scale = 0.5 * (xmax - xmin);	/* 1/2 Range  */

	if (model_type == TREND1D_FOURIER || model_type == TREND1D_ANNUAL)
		scale /= M_PI;

	for (i = 0; i < n_data; i++)
		data[i].x = (data[i].x * scale) + offset;
}

double get_chisq (struct TREND1D_DATA *data, int n_data, int n_model) {
	int	i, nu;
	double	chi = 0.0;


	for (i = 0; i < n_data; i++) {	/* Weight is already squared  */
		if (data[i].w == 1.0)
			chi += (data[i].r * data[i].r);

		else
			chi += (data[i].r * data[i].r * data[i].w);
	}
	nu = n_data - n_model;
	if (nu > 1) return(chi/nu);
	return(chi);
}

void recompute_weights (struct TREND1D_DATA *data, int n_data, double *work, double *scale) {
	int	i;
	double	k, ksq, rr;

	/* First find median { fabs(data[].r) },
		estimate scale from this,
		and compute chisq based on this.  */ 

	for (i = 0; i < n_data; i++) {
		work[i] = fabs(data[i].r);
	}
	qsort((void *)work, (size_t)n_data, sizeof(double), GMT_comp_double_asc);

	if (n_data%2) {
		*scale = 1.4826 * work[n_data/2];
	}
	else {
		*scale = 0.7413 * (work[n_data/2 - 1] + work[n_data/2]);
	}

	k = 1.5 * (*scale);	/*  Huber[1964] weight; 95% efficient for Normal data  */
	ksq = k * k;

	for (i = 0; i < n_data; i++) {
		rr = fabs(data[i].r);
		if (rr <= k)
			data[i].w = 1.0;

		else
			data[i].w = (2*k/rr) - (ksq/(rr*rr) );	/* This is really w-squared  */
	}
}

void load_g_row (double x, int n, double *gr, int m) {
      	  	/* Current data position, appropriately normalized.  */
   	  	/* Number of model parameters, and elements of gr[]  */
      	     	/* Elements of row of G matrix.  */
   	  	/* Parameter indicating model type  */

	/* Routine computes the elements gr[j] in the ith row of the
		G matrix (Menke notation), where x is the ith datum's
		abscissa.  */

	int	j, k;
	double	xk;

	if (n) {

		gr[0] = 1.0;

		switch (m) {

			case TREND1D_POLYNOMIAL:
				/* Create Chebyshev polynomials  */
				if (n > 1) gr[1] = x;
				for (j = 2; j < n; j++) {
					gr[j] = 2 * x * gr[j-1] - gr[j-2];
				}
				break;

			case TREND1D_ANNUAL:
				gr[1] = x;
				for (j = 2, k = 1; j < n-1; j+=2, k++) {
					xk = x*k;
					gr[j] = cos(xk);
					gr[j+1] = sin(xk);
				}
				break;

			case TREND1D_FOURIER:
				for (j = 1; j < n; j++) {
					k = (j + 1)/2;
					if (k > 1) {
						if (j%2)	/* Odd */
							gr[j] = cos(k*x);
						else
							gr[j] = sin(k*x);
					}
					else {
						if (j%2)
							gr[j] = cos(x);
						else
							gr[j] = sin(x);
					}
				}
				break;
		}
	}
}

void calc_m_and_r (struct TREND1D_DATA *data, int n_data, double *model, int n_model, int m_type, double *grow) {
	/*	model[n_model] holds solved coefficients of m_type model.
		grow[n_model] is a vector for a row of G matrix.  */

	int	i, j;

	for (i = 0; i < n_data; i++) {
		load_g_row(data[i].x, n_model, grow, m_type);
		data[i].m = 0.0;
		for (j = 0; j < n_model; j++) {
			data[i].m += model[j]*grow[j];
		}
		data[i].r = data[i].y - data[i].m;
	}
}

void move_model_a_to_b (double *model_a, double *model_b, int n_model, double *chisq_a, double *chisq_b) {
	int	i;
	for(i = 0; i<  n_model; i++) {
		model_b[i] = model_a[i];
	}
	*chisq_b = *chisq_a;
}

void load_gtg_and_gtd (struct TREND1D_DATA *data, int n_data, double *gtg, double *gtd, double *grow, int n_model, int mp, int m_type) {
   	/* mp is row dimension of gtg  */

	int	i, j, k;
	double	wy;

	/* First zero the contents for summing:  */

	for (j = 0; j < n_model; j++) {
		for (k = 0; k < n_model; k++)
			gtg[j + k*mp] = 0.0;

		gtd[j] = 0.0;
	}

	/* Sum over all data  */
	for (i = 0; i < n_data; i++) {
		load_g_row(data[i].x, n_model, grow, m_type);
		if (data[i].w != 1.0) {
			wy = data[i].w * data[i].y;
			for (j = 0; j < n_model; j++) {
				for (k = 0; k < n_model; k++)
					gtg[j + k*mp] += (data[i].w * grow[j] * grow[k]);

				gtd[j] += (wy * grow[j]);
			}
		}
		else {
			for (j = 0; j < n_model; j++) {
				for (k = 0; k < n_model; k++)
					gtg[j + k*mp] += (grow[j] * grow[k]);

				gtd[j] += (data[i].y * grow[j]);
			}
		}
	}
}

void solve_system (double *gtg, double *gtd, double *model, int n_model, int mp, double *lambda, double *v, double *b, double *z, double c_no, int *ir) {

	int	i, j, k, n, m, rank = 0, nrots;
	double	c_test, temp_inverse_ij;

	if (n_model == 1) {
		model[0] = gtd[0] / gtg[0];
		*ir = 1;
	}
	else {
		n = (int)n_model;
		m = (int)mp;
		if(GMT_jacobi(gtg, &n, &m, lambda, v, b, z, &nrots))
			mexPrintf("TREND1D:  Warning:  Matrix Solver Convergence Failure.\n");

		c_test = fabs(lambda[0])/c_no;
		while(rank < n_model && lambda[rank] > 0.0 && lambda[rank] > c_test) rank++;
		for (i = 0; i < n_model; i++) {
			model[i] = 0.0;
			for (j = 0; j < n_model; j++) {
				temp_inverse_ij = 0.0;
				for (k = 0; k <  rank; k++) {
					temp_inverse_ij += (v[i + k*mp] * v[j + k*mp] / lambda[k]);
				}
				model[i] += (temp_inverse_ij * gtd[j]);
			}
		}
		*ir = rank;
	}
}

void GMT_cheb_to_pol (double c[], int n, double a, double b) {
	/* Convert from Chebyshev coefficients used on a t =  [-1,+1] interval
	 * to polynomial coefficients on the original x = [a b] interval.
	 * Modified from Numerical Miracles, ...eh Recipes */
	 
	 int j, k;
	 double sv, cnst, fac, *d, *dd;
	 
	 d  = mxCalloc ((size_t)n, sizeof (double));
	 dd = mxCalloc ((size_t)n, sizeof (double));
	 
	 /* First we generate coefficients for a polynomial in t */
	 
	 d[0] = c[n-1];
	 for (j = n - 2; j >= 1; j--) {
	 	for (k = n - j; k >= 1; k--) {
			sv = d[k];
			d[k] = 2.0 * d[k-1] - dd[k];
			dd[k] = sv;
		}
		sv = d[0];
		d[0] = -dd[0] + c[j];
		dd[0] = sv;
	}
	for (j = n - 1; j >= 1; j--) d[j] = d[j-1] - dd[j];
	/* d[0] = -dd[0] + 0.5 * c[0]; */	/* This is what Num. Rec. says, but we do not do the approx with 0.5 * c[0] */
	d[0] = -dd[0] + c[0];

	/* Next step is to undo the scaling so we can use coefficients with x */

	cnst = fac = 2.0 / (b - a);
	for (j = 1; j < n; j++) {
		d[j] *= fac;
		fac *= cnst;
	}
	cnst = 0.5 * (a + b);
	for (j = 0; j <= n - 2; j++) for (k = n - 2; k >= j; k--) d[k] -= cnst * d[k+1];

	/* Return the new coefficients via c */

	memcpy ((void *)c, (void *)d, (size_t)(n * sizeof (double)));

	mxFree ((void *) d);
	mxFree ((void *) dd);
}

void *New_Trend1d_Ctrl () {	/* Allocate and initialize a new control structure */
	struct TREND1D_CTRL *C;
	
	C = (struct TREND1D_CTRL *) mxCalloc ((size_t)1, sizeof (struct TREND1D_CTRL));
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->C.value = 1.0e06;		/* Condition number for matrix solution  */	
	C->I.value = 0.51;		/* Confidence interval for significance test  */
	C->N.mode = TREND1D_POLYNOMIAL;
	return ((void *)C);
}

void Free_Trend1d_Ctrl (struct TREND1D_CTRL *C) {	/* Deallocate control structure */
	mxFree ((void *)C);	
}

#define MAX_SWEEPS 50

int	GMT_jacobi (double *a, int *n, int *m, double *d, double *v, double *b, double *z, int *nrots) {
/*
 *
 * Find eigenvalues & eigenvectors of a square symmetric matrix by Jacobi's
 * method.  Given A, find V and D such that A = V * D * V-transpose, with
 * V an orthogonal matrix and D a diagonal matrix.  The eigenvalues of A
 * are on diag(D), and the j-th column of V is the eigenvector corresponding
 * to the j-th diagonal element of D.  Returns 0 if OK, -1 if it fails to
 * converge in MAX_SWEEPS.
 *
 * a is sent as a square symmetric matrix, of size n, and row dimension m.
 * Only the diagonal and super-diagonal elements of a will be used, so the
 * sub-diagonal elements could be used to preserve a, or could have been
 * destroyed by an earlier attempt to form the Cholesky decomposition of a.
 * On return, the super-diagonal elements are destroyed.  The diagonal and
 * sub-diagonal elements are unchanged.
 * d is returned as an n-vector containing the eigenvalues of a, sorted 
 * so that d[i] >= d[j] when i < j.  d = diag(D).
 * v is returned as an n by n matrix, V, with row dimension m, and the
 * columns of v are the eigenvectors corresponding to the values in d.
 * b is an n-vector of workspace, used to keep a copy of the diagonal
 * elements which is updated only after a full sweep.  
 * z is an n-vector of workspace, used to accumulate the updates to
 * the diagonal values of a during each sweep.  This reduces round-
 * off problems.
 * nrots is the number of rotations performed.  Bounds on round-off error
 * can be estimated from this if desired.
 *
 * Numerical Details:
 * The basic algorithms is in many textbooks.  The idea is to make an
 * infinite series (which turns out to be at quadratically convergent)
 * of steps, in each of which A_new = P-transpose * A_old * P, where P is
 * a plane-rotation matrix in the p,q plane, through an angle chosen to
 * zero A_new(p,q) and A_new(q,p).  The sum of the diagonal elements
 * of A is unchanged by these operations, but the sum of squares of 
 * diagonal elements of a is increased by 2 * |A_old(p,q)| at each step.
 * Although later steps make non-zero again the previously zeroed entries,
 * the sum of squares of diagonal elements increases with each rotation,
 * while the sum of squares of off-diagonals keeps decreasing, so that
 * eventually A_new is diagonal to machine precision.  This should 
 * happen in a few (3 to 7) sweeps.
 *
 * If only the eigenvalues are wanted then there are faster methods, but
 * if all eigenvalues and eigenvectors are needed, then this method is
 * only somewhat slower than the fastest method (Householder tri-
 * diagonalization followed by symmetric QR iterations), and this method
 * is numerically extremely stable.
 *
 * C G J Jacobi ("Ueber ein leichtes Vefahren, die in der Theorie der 
 * Saekularstoerungen vorkommenden Gelichungen numerisch aufzuloesen",
 * Crelle's Journal, v. 30, pp. 51--94, 1846) originally searched the
 * entire (half) matrix for the largest |A(p,q)| to select each step.
 * When the method was developed for machine computation (R T Gregory,
 * "Computing eigenvalues and eigenvectors of a symmetric matrix on
 * the ILLIAC", Math. Tab. and other Aids to Comp., v. 7, pp. 215--220,
 * 1953) it was done with a series of "sweeps" through the upper triangle,
 * visiting all p,q in turn.  Later, D A Pope and C Tompkins ("Maximizing
 * functions of rotations - experiments concerning speed of diagonalization
 * of symmetric matrices using Jacobi's method", J Assoc. Comput. Mach.
 * v. 4, pp. 459--466, 1957) introduced a variant that skips small
 * elements on the first few sweeps.  The algorithm here was given by 
 * Heinz Rutishauser (1918--1970) and published in Numer. Math. v. 9, 
 * pp 1--10, 1966, and in Linear Algebra (the Handbook for Automatic 
 * Computation, v. II), by James Hardy Wilkinson and C. Reinsch (Springer-
 * Verlag, 1971).  It also appears in Numerical Recipes.
 *
 * This algorithm takes care to avoid round-off error in several ways.
 * First, although there are four values of theta in (-pi, pi] that
 * would zero A(p,q), there is only one with magnitude <= pi/4.
 * This one is used.  This is most stable, and also has the effect
 * that, if A_old(p,p) >= A_old(q,q) then A_new(p,p) > A_new(q,q).
 * Two copies of the diagonal elements are maintained in d[] and b[].
 * d[] is updated immediately in each rotation, and each new rotation
 * is computed based on d[], so that each rotation gets the benefit
 * of the previous ones.  However, z[] is also used to accumulate
 * the sum of all the changes in the diagonal elements during one sweep, 
 * and z[] is used to update b[] after each sweep.  Then b is copied
 * to d.  In this way, at the end of each sweep, d is reset to avoid
 * accumulating round-off.
 *
 * This routine determines whether y is small compared to x by testing
 * if (fabs(y) + fabs(x) == fabs(x) ).  It is assumed that the
 * underflow which may occur here is nevertheless going to allow this
 * expression to be evaluated as TRUE or FALSE and execution to 
 * continue.  If the run environment doesn't allow this, the routine
 * won't work properly.
 *
 * programmer:	W. H. F. Smith, 7 June, 1991.
 * Revised:	PW: 12-MAR-1998 for GMT 3.1
 * Revision by WHF Smith, March 03, 2000, to speed up loop indexes.
 */
	int	p, q, pp, pq, mp1, pm, qm, nsweeps, j, jm, i, k;
	double	sum, threshold, g, h, t, theta, c, s, tau;


	/* Begin by initializing v, b, d, and z.  v = identity matrix,
		b = d = diag(a), and z = 0:  */
	
	memset ((void *)v, 0, (size_t)((*m)*(*n)*sizeof(double)) );
	memset ((void *)z, 0, (size_t)((*n)*sizeof(double)) );
	
	mp1 = (*m) + 1;
	
	for (p = 0, pp = 0; p < (*n); p++, pp+=mp1) {
		v[pp] = 1.0;
		b[p] = a[pp];
		d[p] = b[p];
	}

	/* End of initializations.  Set counters and begin:  */

	(*nrots) = 0;
	nsweeps = 0;

	while (nsweeps < MAX_SWEEPS) {
	
		/* Sum off-diagonal elements of upper triangle.  */
		sum = 0.0;
		for (q = 1, qm = (*m); q < (*n); q++, qm += (*m) ) {
			for (p = 0, pq = qm; p < q; p++, pq++) {
				sum += fabs(a[pq]);
			}
		}
		
		/* Exit this loop (converged) when sum == 0.0  */
		if (sum == 0.0) break;


		/* If (nsweeps < 3) do only bigger elements;  else all  */
		threshold =  (nsweeps < 3) ? 0.2 * sum / ( (*n) * (*n) ) : 0.0;

		/* Now sweep whole upper triangle doing Givens rotations:  */
		
		for (q = 1, qm = (*m); q < (*n); q++, qm += (*m) ) {
			for (p = 0, pm = 0, pq = qm; p < q; p++, pm += (*m), pq++) {
				/* In 3/2000 I swapped order of these loops,
					to allow simple incrementing of pq  */
			
				if (a[pq] == 0.0) continue;	/* New 3/2000  */
			
				g = 100.0 * fabs(a[pq]);
				
				/* After four sweeps, if g is small relative
					to a(p,p) and a(q,q), skip the 
					rotation and set a(p,q) to zero.  */

				if ( (nsweeps > 3) && ( (fabs(d[p])+g) == fabs(d[p]) ) && ( (fabs(d[q])+g) == fabs(d[q]) ) ) {
					a[pq] = 0.0;
				}
				else if (fabs(a[pq]) > threshold) {

					h = d[q] - d[p];
					
					if (h == 0.0) {
						t = 1.0;	/* This if block is new 3/2000  */
					}
					else if ( (fabs(h)+g) ==  fabs(h) ) {
						t = a[pq] / h;
					}
					else {
						theta = 0.5 * h / a[pq];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta) );
						if (theta < 0.0) t = -t;
					}

					c = 1.0 / sqrt(1.0 + t*t);
					s = t * c;
					tau = s / (1.0 + c);
					
					h = t * a[pq];
					z[p] -= h;
					z[q] += h;
					d[p] -= h;
					d[q] += h;
					a[pq] = 0.0;

					for (j = 0; j < p; j++) {
						g = a[j + pm];
						h = a[j + qm];
						a[j + pm] = g - s * (h + g * tau);
						a[j + qm] = h + s * (g - h * tau);
					}
					for (j = p+1, jm = (*m)*(p+1); j < q; j++, jm += (*m) ) {
						g = a[p + jm];
						h = a[j + qm];
						a[p + jm] = g - s * (h + g * tau);
						a[j + qm] = h + s * (g - h * tau);
					}
					for (j = q+1, jm = (*m)*(q+1); j < (*n); j++, jm += (*m) ) {
						g = a[p + jm];
						h = a[q + jm];
						a[p + jm] = g - s * (h + g * tau);
						a[q + jm] = h + s * (g - h * tau);
					}

					for (j = 0; j < (*n); j++) {
						g = v[j + pm];
						h = v[j + qm];
						v[j + pm] = g - s * (h + g * tau);
						v[j + qm] = h + s * (g - h * tau);
					}

					(*nrots)++;
				}
			}
		}
		
		/* End of one sweep of the upper triangle.  */
		
		nsweeps++;

		for (p = 0; p < (*n); p++) {
			b[p] += z[p];	/* Update the b copy of diagonal  */
			d[p] = b[p];	/* Replace d with b to reduce round-off error  */
			z[p] = 0.0;	/* Clear z.  */
		}
	}

	/* Get here via break when converged, or when nsweeps == MAX_SWEEPS.
		Sort eigenvalues by insertion:  */

	for (i = 0; i < (*n)-1; i++) {
		k = i;
		g = d[i];
		for (j = i+1; j < (*n); j++) {  /* Find max location  */
			if (d[j] >= g) {
				k = j;
				g = d[j];
			}
		}
		if (k != i) {  /*  Need to swap value and vector  */
			d[k] = d[i];
			d[i] = g;
			p = i * (*m);
			q = k * (*m);
			for (j = 0; j < (*n); j++) {
				g = v[j + p];
				v[j + p] = v[j + q];
				v[j + q] = g;
			}
		}
	}

	/* Return 0 if converged; else print warning and return -1:  */

	if (nsweeps == MAX_SWEEPS) {
		mexPrintf ("GMT_jacobi:  Failed to converge in %ld sweeps\n", nsweeps);
		return(-1);
	}
	return(0);
}


int	GMT_sig_f (double chi1, int n1, double chi2, int n2, double level, double *prob) {
	/* Returns TRUE if chi1/n1 significantly less than chi2/n2
		at the level level.  Returns FALSE if:
			error occurs in GMT_f_test_new();
			chi1/n1 not significantly < chi2/n2 at level.

			Changed 12 August 1999 to use GMT_f_test_new()  */

	int	trouble;

	trouble = GMT_f_test_new (chi1, n1, chi2, n2, prob, -1);
	if (trouble) return(0);
	return((*prob) >= level);
}


int	GMT_f_test_new (double chisq1, int nu1, double chisq2, int nu2, double *prob, int iside) {
	/* Given chisq1 and chisq2, random variables distributed as chi-square
		with nu1 and nu2 degrees of freedom, respectively, except that
		chisq1 is scaled by var1, and chisq2 is scaled by var2, let
		the null hypothesis, H0, be that var1 = var2.  This routine
		assigns prob, the probability that we can reject H0 in favor
		of a new hypothesis, H1, according to iside:
			iside=+1 means H1 is that var1 > var2
			iside=-1 means H1 is that var1 < var2
			iside=0  means H1 is that var1 != var2.
		This routine differs from the old GMT_f_test() by adding the
		argument iside and allowing one to choose the test.  The old
		routine in effect always set iside=0.
		This routine also differs from GMT_f_test() in that the former
		used the incomplete beta function and this one uses GMT_f_q().

		Returns 0 on success, -1 on failure.

		WHF Smith, 12 August 1999.
	*/

	double	q;	/* The probability from GMT_f_q(), which is the prob
				that H0 should be retained even though
				chisq1/nu1 > chisq2/nu2.  */

	if (chisq1 <= 0.0 || chisq2 <= 0.0 || nu1 < 1 || nu2 < 1) {
		*prob = mxGetNaN();
		mexPrintf ("GMT_f_test_new:  ERROR:  Bad argument(s).\n");
		return (-1);
	}

	GMT_f_q (chisq1, nu1, chisq2, nu2, &q);

	if (iside > 0) {
		*prob = 1.0 - q;
	}
	else if (iside < 0) {
		*prob = q;
	}
	else if ( (chisq1/nu1) <= (chisq2/nu2) ) {
		*prob = 2.0*q;
	}
	else {
		*prob = 2.0*(1.0 - q);
	}

	return (0);
}


int     GMT_f_q (double chisq1, int nu1, double chisq2, int nu2, double *prob) {
	/* Routine to compute Q(F, nu1, nu2) = 1 - P(F, nu1, nu2), where nu1
		and nu2 are positive integers, chisq1 and chisq2 are random
		variables having chi-square distributions with nu1 and nu2
		degrees of freedom, respectively (chisq1 and chisq2 >= 0.0),
		F = (chisq1/nu1)/(chisq2/nu2) has the F-distribution, and
		P(F, nu1, nu2) is the cumulative F-distribution, that is,
		the integral from 0 to (chisq1/nu1)/(chisq2/nu2) of the F-
		distribution.  Q = 1 - P is small when (chisq1/nu1)/(chisq2/nu2)
		is large with respect to 1.  That is, the value returned by
		this routine is the likelihood that an F >= (chisq1/nu1)/
		(chisq2/nu2) would occur by chance.

		Follows Abramowitz and Stegun.
		This is different from the method in Numerical Recipes, which
		uses the incomplete beta function but makes no use of the fact
		that nu1 and nu2 are known to be integers, and thus there is
		a finite limit on the sum for their expression.

		W H F Smith, August, 1999.

		REVISED by W H F Smith, October 27, 2000 after GMT 3.3.6 release.
		I found that the A&S methods overflowed for large nu1 and nu2, so
		I decided to go back to the GMT_inc_beta way of doing things.

	*/

	/* Check range of arguments:  */

	if (nu1 <= 0 || nu2 <= 0 || chisq1 < 0.0 || chisq2 < 0.0) {
		mexPrintf ("GMT_f_q:  Bad argument(s).\n");
		return (-1);
	}


	/* Extreme cases evaluate immediately:  */

	if (chisq1 == 0.0) {
		*prob = 1.0;
		return (0);
	}
	if (chisq2 == 0.0) {
		*prob = 0.0;
		return (0);
	}

	/* REVISION of Oct 27, 2000:  This inc beta call here returns
		the value.  All subsequent code is not used.  */

	if (GMT_inc_beta(0.5*nu2, 0.5*nu1, chisq2/(chisq2+chisq1), prob) ) {
		mexPrintf("GMT_q_p:  Trouble in GMT_inc_beta call.\n");
		return(-1);
	}
	return (0);
}


int	GMT_inc_beta (double a, double b, double x, double *ibeta) {
	double	bt, gama, gamb, gamab;

	if (a <= 0.0) {
		mexPrintf("GMT_inc_beta:  Bad a (a <= 0).\n");
		return(-1);
	}
	if (b <= 0.0) {
		mexPrintf("GMT_inc_beta:  Bad b (b <= 0).\n");
		return(-1);
	}
	if (x > 0.0 && x < 1.0) {
		GMT_ln_gamma_r(a, &gama);
		GMT_ln_gamma_r(b, &gamb);
		GMT_ln_gamma_r( (a+b), &gamab);
		bt = exp(gamab - gama - gamb
			+ a * d_log(x) + b * d_log(1.0 - x) );

		/* Here there is disagreement on the range of x which
			converges efficiently.  Abramowitz and Stegun
			say to use x < (a - 1) / (a + b - 2).  Editions
			of Numerical Recipes thru mid 1987 say
			x < ( (a + 1) / (a + b + 1), but the code has
			x < ( (a + 1) / (a + b + 2).  Editions printed
			late 1987 and after say x < ( (a + 1) / (a + b + 2)
			in text as well as code.  What to do ? */

		if (x < ( (a + 1) / (a + b + 2) ) ) {
			*ibeta = bt * GMT_cf_beta(a, b, x) / a;
		}
		else {
			*ibeta = 1.0 - bt * GMT_cf_beta(b, a, (1.0 - x) ) / b;
		}
		return(0);
	}
	else if (x == 0.0) {
		*ibeta = 0.0;
		return(0);
	}
	else if (x == 1.0) {
		*ibeta = 1.0;
		return(0);
	}
	else if (x < 0.0) {
		mexPrintf("GMT_inc_beta:  Bad x (x < 0).\n");
		*ibeta = 0.0;
	}
	else if (x > 1.0) {
		mexPrintf("GMT_inc_beta:  Bad x (x > 1).\n");
		*ibeta = 1.0;
	}
	return(-1);
}

int	GMT_ln_gamma_r(double x, double *lngam) {
	/* Get natural logrithm of Gamma(x), x > 0.
		To maintain full accuracy, this
		routine uses Gamma(1 + x) / x when
		x < 1.  This routine in turn calls
		GMT_ln_gamma(x), which computes the
		actual function value.  GMT_ln_gamma
		assumes it is being called in a
		smart way, and does not check the
		range of x.  */

	if (x > 1.0) {
		*lngam = GMT_ln_gamma(x);
		return(0);
	}
	if (x > 0.0 && x < 1.0) {
		*lngam = GMT_ln_gamma(1.0 + x) - d_log(x);
		return(0);
	}
	if (x == 1.0) {
		*lngam = 0.0;
		return(0);
	}
	mexPrintf("Ln Gamma:  Bad x (x <= 0).\n");
	return(-1);
}

double	GMT_ln_gamma (double xx) {
	/* Routine to compute natural log of Gamma(x)
		by Lanczos approximation.  Most accurate
		for x > 1; fails for x <= 0.  No error
		checking is done here; it is assumed
		that this is called by GMT_ln_gamma_r()  */

	static double	cof[6] = {
		76.18009173,
		-86.50532033,
		24.01409822,
		-1.231739516,
		0.120858003e-2,
		-0.536382e-5
	};

	static double	stp = 2.50662827465, half = 0.5, one = 1.0, fpf = 5.5;
	double	x, tmp, ser;

	int	i;

	x = xx - one;
	tmp = x + fpf;
	tmp = (x + half) * d_log(tmp) - tmp;
	ser = one;
	for (i = 0; i < 6; i++) {
		x += one;
		ser += (cof[i]/x);
	}
	return(tmp + d_log(stp*ser) );
}

double	GMT_cf_beta (double a, double b, double x) {
	/* Continued fraction method called by GMT_inc_beta.  */

	static int	itmax = 100;
	static double	eps = 3.0e-7;

	double	am = 1.0, bm = 1.0, az = 1.0;
	double	qab, qap, qam, bz, em, tem, d;
	double	ap, bp, app, bpp, aold;

	int	m = 0;

	qab = a + b;
	qap = a + 1.0;
	qam = a - 1.0;
	bz = 1.0 - qab * x / qap;

	do {
		m++;
		em = m;
		tem = em + em;
		d = em*(b-m)*x/((qam+tem)*(a+tem));
		ap = az+d*am;
		bp = bz+d*bm;
		d = -(a+m)*(qab+em)*x/((a+tem)*(qap+tem));
		app = ap+d*az;
		bpp = bp+d*bz;
		aold = az;
		am = ap/bpp;
		bm = bp/bpp;
		az = app/bpp;
		bz = 1.0;
	} while ( ( (fabs(az-aold) ) >= (eps * fabs(az) ) ) && (m < itmax) );

	if (m == itmax)
		mexPrintf("GMT_cf_beta:  A or B too big, or ITMAX too small.\n");

	return(az);
}

int GMT_comp_double_asc (const void *p_1, const void *p_2) {
	/* Returns -1 if point_1 is < that point_2,
	   +1 if point_2 > point_1, and 0 if they are equal
	*/
	int bad_1, bad_2;
	double *point_1, *point_2;

	point_1 = (double *)p_1;
	point_2 = (double *)p_2;
	bad_1 = mxIsNaN ((*point_1));
	bad_2 = mxIsNaN ((*point_2));

	if (bad_1 && bad_2) return (0);
	if (bad_1) return (1);
	if (bad_2) return (-1);

	if ( (*point_1) < (*point_2) )
		return (-1);
	else if ( (*point_1) > (*point_2) )
		return (1);
	else
		return (0);
}

double	tpvalue (double x, double v) {
	/* Compute p-value for t statistic 
	(from a function of the same name inside the MATLAB corrcoef) */

	double	p = 0, xx = 0;

	if (x == 0 && v != 1 && v > 0) return(0.5);

	/* Return NaN for invalid inputs. */
	if (v <= 0 || mxIsNaN(x) || mxIsNaN(v)) return(mxGetNaN());

	if (v == 1)
		p = 0.5 + atan(x) / M_PI;

	/* See Abramowitz and Stegun, formulas 26.5.27 and 26.7.1 */
	if (x != 0 && v != 1 && v > 0) {	/* first compute F(-|x|) */
		xx = v / (v + x * x);
		GMT_inc_beta(v/2, 0.5, xx, &p);
		p /= 2;
	}

	/* Adjust for x>0.  Right now p<0.5, so this is numerically safe. */
	if (x > 0 && v != 1 && v > 0)	p = 1 - p;

	return (p);
}
