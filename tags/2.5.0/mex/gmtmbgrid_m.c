/*--------------------------------------------------------------------
 *
 *	$Id: surface.c,v 1.18 2004/04/25 09:10:46 pwessel Exp $
 *	Copyright (c) 1991-2004 by P. Wessel and W. H. F. Smith
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
 *--------------------------------------------------------------------
 *
 * GMTMBGRID is an adaptation of the MB-SYSTEM MBGRID program to work with
 * the GMT package. This program uses a gaussian weighted mean algorithm to
 * grid regions covered by data and than fills gaps using a minimum curvature
 * algorithm. Two algorithms are available. The minimum curvature algorithm
 * from the GMT program surface and the IGPP/SIO zgrid routine for thin plate
 * spline interpolation.
 * Some of the MBGRID original options were removed (for example, it now works
 * only with xyz data), but it can now use both of the interpolation algorithms
 * on running time. It is also much more memory efficient than the original.
 *
 * With high density of data points the zgrid algorithm seams to be much
 * faster than surface, but extrapolation of large areas can result in
 * bizzare results.
 *
 *--------------------------------------------------------------------
 * MBGRID Author:	D. W. Caress
 * Date:	February 22, 1993
 * Rewrite:	May 2, 1994
 * Rerewrite:	April 25, 1995
 * Rererewrite:	January 2, 1996
 */

/*--------------------------------------------------------------------
 *    Coffeeright (c) 2004-2012 by J. Luis
 *
 * GMTMBGRID Author:	Joaquim Luis
 * Created:		04-JUN-2006 
 *
 * MEXIFIED:		11-JUL-2008 
 *
 * THIS IS A STANADLONE CODE. That is, there are no dependencies on either
 * MB-SYSTEM or the GMT libraries (several routines picked from there)
 *
 *	Contact info: jluis@ualg.pt
 *-------------------------------------------------------------------- */

#include "mex.h"
#include <math.h>
#include <string.h>

#define	FALSE	0
#define	TRUE	1
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define D2R (M_PI / 180.0)
#define cosd(x) cos ((x) * D2R)
#define CNULL	((char *)NULL)
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef rint
#define rint(x) (floor((x)+0.5))
#endif
#ifndef irint
#define irint(x) ((int)rint(x))
#endif

/* These macros calculate the number of nodes in x or y  or the increment dx, dy*/

#define GMT_get_n(min,max,inc,off) ((int)irint (((max) - (min)) / (inc)) + 1 - (off))
#define GMT_get_inc(min,max,n,off) (((max) - (min)) / ((n) + (off) - 1))

/* flag for no data in grid */
#define	NO_DATA_FLAG	9999999999999999

/* interpolation mode */
#define MBGRID_INTERP_NONE	0
#define MBGRID_INTERP_GAP	1
#define MBGRID_INTERP_NEAR	2
#define MBGRID_INTERP_ALL	3

/* maximum comment length in characters */
#define MB_COMMENT_MAXLINE 1944

#define INTERP_SURFACE		1
#define INTERP_ZGRID		0

#define OUTSIDE 2000000000	/* Index number indicating data is outside useable area */

/* Variables for surface */
struct SURFACE_DATA {
	float x;
	float y;
	float z;
	int index;
};

struct SURFACE_BRIGGS {
	double b[6];
};

struct GRD_HEADER {
/* Do not change the first three items. They are copied verbatim to the native grid header */
	int nx;				/* Number of columns */
	int ny;				/* Number of rows */
	int node_offset;		/* 0 for node grids, 1 for pixel grids */
/* This section is flexible. It is not copied to any grid header */
	int type;			/* Grid format */
	char name[256];			/* Actual name of the file after any ?<varname> and =<stuff> has been removed */
	char varname[80];		/* NetCDF: variable name */
	int y_order;			/* NetCDF: 1 if S->N, -1 if N->S */
	int z_id;			/* NetCDF: id of z field */
	int ncid;			/* NetCDF: file ID */
	int t_index[3];			/* NetCDF: index of higher coordinates */
	double nan_value;		/* Missing value as stored in grid file */
	double xy_off;			/* 0.0 (node_offset == 0) or 0.5 ( == 1) */
/* The following elements should not be changed. They are copied verbatim to the native grid header */
	double x_min;			/* Minimum x coordinate */
	double x_max;			/* Maximum x coordinate */
	double y_min;			/* Minimum y coordinate */
	double y_max;			/* Maximum y coordinate */
	double z_min;			/* Minimum z value */
	double z_max;			/* Maximum z value */
	double x_inc;			/* x increment */
	double y_inc;			/* y increment */
	double z_scale_factor;		/* grd values must be multiplied by this */
	double z_add_offset;		/* After scaling, add this */
	char x_units[80];		/* units in x-direction */
	char y_units[80];		/* units in y-direction */
	char z_units[80];		/* grid value units */
	char title[80];			/* name of data set */
	char command[320];		/* name of generating command */
	char remark[160];		/* comments re this data set */
}; 

int verbose = FALSE;
static int npoints=0;			/* Number of data points */
static int nx=0;			/* Number of nodes in x-dir. */
static int ny=0;			/* Number of nodes in y-dir. (Final grid) */
static int mx = 0;
static int my = 0;
static int ij_sw_corner, ij_se_corner, ij_nw_corner, ij_ne_corner;
static int block_nx;			/* Number of nodes in x-dir for a given grid factor */
static int block_ny;			/* Number of nodes in y-dir for a given grid factor */
static int max_iterations=250;		/* Max iter per call to iterate */
static int total_iterations = 0;
static int grid, old_grid;		/* Node spacings  */
static int grid_east;
static int n_fact = 0;			/* Number of factors in common (ny-1, nx-1) */
static int factors[32];		/* Array of common factors */
static int n_empty;			/* No of unconstrained nodes at initialization  */
static int set_low = 0;		/* 0 unconstrained,1 = by min data value, 2 = by user value */
static int set_high = 0;		/* 0 unconstrained,1 = by max data value, 2 = by user value */
static int constrained = FALSE;	/* TRUE if set_low or set_high is TRUE */
static double low_limit, high_limit;	/* Constrains on range of solution */
static double xmin, xmax, ymin, ymax;	/* minmax coordinates */
static float *lower, *upper;		/* arrays for minmax values, if set */
static double xinc, yinc;		/* Size of each grid cell (final size) */
static double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
static double r_xinc, r_yinc, r_grid_xinc, r_grid_yinc;	/* Reciprocals  */
static double converge_limit = 0.0;	/* Convergence limit */
static double radius = 0.0;			/* Search radius for initializing grid  */
static double	tension = 0.0;
static double	boundary_tension = 0.0;
static double	interior_tension = 0.0;
static double	a0_const_1, a0_const_2;	/* Constants for off grid point equation  */
static double	e_2, e_m2, one_plus_e2;
static double	eps_p2, eps_m2, two_plus_ep2, two_plus_em2;
static double	x_edge_const, y_edge_const;
static double	epsilon = 1.0;
static double	z_mean;
static double	z_scale = 1.0;		/* Root mean square range of z after removing planar trend  */
static double	r_z_scale = 1.0;	/* reciprocal of z_scale  */
static double	plane_c0, plane_c1, plane_c2;	/* Coefficients of best fitting plane to data  */
static double small;			/* Let data point coincide with node if distance < small */
static float *u;			/* Pointer to grid array */
static char *iu;			/* Pointer to grid info array */
static char mode_type[2] = {'I','D'};	/* D means include data points when iterating
				 	* I means just interpolate from larger grid */

static int	offset[25][12];		/* Indices of 12 nearby points in 25 cases of edge conditions  */
static double	coeff[2][12];		/* Coefficients for 12 nearby points, constrained and unconstrained  */

static double relax_old, relax_new = 1.4;	/* Coefficients for relaxation factor to speed up convergence */


int mb_zgrid(float *z, int *nx, int *ny, float *x1, float *y1, float *dx, float *dy, float *xyz, 
		int *n, float *zpij, int *knxt, int *imnew, float *cay, int *nrng);
int mb_surface(int verbose, int ndat, float *xdat, float *ydat, float *zdat,
		double xxmin, double xxmax, double yymin, double yymax, double xxinc, double yyinc,
		double ttension, float *sgrid);
int read_data(int ndat, float *xdat, float *ydat, float *zdat);
int gcd_euclid(int a, int b);	/* Finds the greatest common divisor  */
int get_prime_factors(int n, int *f), iterate(int mode);
int compare_points (const void *point_1v, const void *point_2v);
void set_grid_parameters(void), throw_away_unusables(void), remove_planar_trend(void), rescale_z_values(void);
void load_constraints(char *low, char *high), smart_divide(void), set_offset(void), set_index(void), initialize_grid(void), set_coefficients(void);
void find_nearest_point(void), fill_in_forecast(void), check_errors(void), replace_planar_trend(void);
void new_initialize_grid(void);
void get_output(float *sgrid);

int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (double w, double e, double s, double n);
double ddmmss_to_degree (char *text);
int GMT_getinc (char *line, double *dx, double *dy);
int GMT_getincn (char *line, double inc[], int n);
int GMT_strtok (const char *string, const char *sep, int *pos, char *token);
void GMT_RI_prepare (struct GRD_HEADER *h);

struct GRD_HEADER h;
static struct SURFACE_DATA *data;	/* Data point and index to node it currently constrains  */
static struct SURFACE_BRIGGS *briggs;	/* Coefficients in Taylor series for Laplacian(z) a la I. C. Briggs (1974)  */

int GMT_inc_code[2] = {0, 0};

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	error = 0;
	int	argc = 0, n_arg_no_char = 0;
	char	**argv;
	float	*pdata;

	/* MBIO read control parameters */
	double	bounds[4];
	char	*backgroundfile = NULL, *grdfile = NULL;

	/* mbgrid control variables */
	int	clip = 0;
	int	clipmode = MBGRID_INTERP_ALL;
	int	use_NaN = TRUE;
	double	clipvalue = NO_DATA_FLAG;
	float	outclipvalue = NO_DATA_FLAG;
	double	scale = 1.0;
	double	border = 0.0;
	double	extend = 0.0;
	
	int	grdrasterid = 0;

	/* grid variables */
	double	wbnd[4], xx, yy, xx2, factor, weight, zmin, zmax, zclip;
	int	gxdim, gydim, offx, offy, xtradim;
	double	*grid = NULL, *norm = NULL;
	float	*bdata = NULL, *sdata = NULL, *work1 = NULL;
	int	*work2 = NULL, *work3 = NULL;
	float	*bxdata = NULL, *bydata = NULL, *bzdata = NULL;
	float	*sxdata = NULL, *sydata = NULL, *szdata = NULL;
	float	*output = NULL, *sgrid = NULL;
	int	*num = NULL, *cnt = NULL;
	float	xmin, ymin, ddx, ddy, zflag = 5.0e34, cay;
	int	ndata, nbackground;
	int	nbinset, nbinzero, nbinspline, nbinbackground = 0;

	/* other variables */
	int	i, j, k, n, ii, jj, iii, jjj, kkk, ir, i1, i2, j1, j2;
	int	kgrid, kint, ix, iy, ix1, ix2, iy1, iy2, len;
	int	dmask[9], n_files = 0;
	int	interp_method = INTERP_SURFACE;
	int	n_pts, n_pts2, alloc_out = FALSE, got_R = FALSE;
	double	r, *in0, *in1, *in2, *info;
	char	fill[256];
	
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
	argv[0] = "gmtmbgrid";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
              
				/* Common parameters */
                      
				case 'V':
					verbose = TRUE;
					break;
				case 'R':
					error += decode_R (argv[i], &h.x_min, &h.x_max, &h.y_min, &h.y_max);
					got_R = TRUE;
					break;
                              
				case 'B':
					border = atof (&argv[i][2]);
					break;
				case 'C':
					len = strlen(argv[i]);
					strcpy (fill, &argv[i][2]);
					if (argv[i][len-1] == 'O' || argv[i][len-1] == 'o')
						clipmode = MBGRID_INTERP_NONE;
					else if (argv[i][len-1] == 'N' || argv[i][len-1] == 'n') {
						fill[len-1] = '\0';
						clip = atoi(fill);
						clipmode = MBGRID_INTERP_NEAR;
					}
					else if (argv[i][len-1] == 'G' || argv[i][len-1] == 'g') {
						fill[len-1] = '\0';
						clip = atoi(fill);
						clipmode = MBGRID_INTERP_GAP;
					}
					else {
						clip = atoi(&argv[i][2]);
						clipmode = MBGRID_INTERP_GAP;
					}
					break;
				case 'E':
					extend = atof (&argv[i][2]);
					break;
				case 'F':
					backgroundfile = &argv[i][2];
					if ((grdrasterid = atoi(backgroundfile)) <= 0)
						grdrasterid = -1;
					break;
				case 'I':
					if (GMT_getinc (&argv[i][2], &h.x_inc, &h.y_inc)) {
						mexPrintf ("GMTMBGRID: SYNTAX ERROR -I option.  Correct syntax:\n");
						mexPrintf ("\t-I<xinc>[m|c|e|k|i|n|+][=][/<yinc>[m|c|e|k|i|n|+][=]]\n");
						mexPrintf ("\t   Give increment and append unit (m)inute, se(c)ond, m(e)ter, (k)ilometer, m(i)les, (n)autical miles.\n");
						mexPrintf ("\t   (Note: m,c,e,k,i,n only apply to geographic regions specified in degrees)\n");
						mexPrintf ("\t   Append = to adjust the domain to fit the increment [Default adjusts increment to fit domain].\n");
						mexPrintf ("\t   Alternatively, specify number of nodes by appending +. Then, the increments are calculated\n");
						mexPrintf ("\t   from the given domain and node-registration settings.\n");
						error = TRUE;
					}
					break;
				case 'M':
					if (argv[i][2] == 'Z' || argv[i][2] == 'z')
						interp_method = INTERP_ZGRID;
					else if (argv[i][2] == 'S' || argv[i][2] == 's')
						interp_method = INTERP_SURFACE;
					else
						mexPrintf("GMTMBGRID: GMT SYNTAX ERROR -M option: Trash given - Ignored\n");
					break;
				case 'T':
					tension = atof (&argv[i][2]);
					break;
				case 'W':
					scale = atof (&argv[i][2]);
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	GMT_RI_prepare (&h);	/* Ensure -R -I consistency and set nx, ny */

	if (argc == 1 || error) {	/* Display usage */
		mexPrintf ("GMTMBGRID - Adjustable tension continuous curvature surface gridding\n\n");
		mexPrintf ("usage: [Z,head] = gmtmbgrid_m(x,y,z, '-I<xinc>[m|c][/<yinc>[m|c]]',\n");
		mexPrintf ("\t'-R<west>/<east>/<south>/<north>', '[-B<border>]', '[-C<clip>[g|n|o]]',\n");
		mexPrintf ("\t'[-E<extend>]', '[-F<background_grid>]', '[-W<scale>]'\n\n");

		mexPrintf ("\tGMTMBGRID will read from standard input or <xyz-file[s]>.\n");
		mexPrintf ("\tNotice that ONLY the first file is used in the main interpolation>.\n");
		mexPrintf ("\tA second file, if provided, will be used as 'background'.\n\n");
		mexPrintf ("\tRequired arguments to GMTMBGRID:\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-B<border> Sets the border of a smoothly interpolated grid to the value\n");
		mexPrintf ("\t-tborder wherever no data exist, provided border > 0.0.\n");
		mexPrintf ("\t-C<clip> clip value sets the distance from data (in grid cells) that\n");
		mexPrintf ("\t\tthe spline interpolation may be applied.\n");
		mexPrintf ("\t\tAppend n to fill all undefined cells within a distance of clip cells from data.\n");
		mexPrintf ("\t\tAppend g to fill data gaps up to two times clip grid cells in size.\n");
		mexPrintf ("\t\tUse -Co to not apply spline interpolation.\n");
		mexPrintf ("\t\tIf no letter is appended it defaults to n.\n");
		mexPrintf ("\t\tDefault (no -C) is to interpolate all grid cells.\n");
		mexPrintf ("\t-E<extend> uses data extending extend*nx|ny in the interpolation.\n");

		mexPrintf ("\t-F<background_grid> Enables filling in all undefined grid cells with\n");
		mexPrintf ("\t\tdata extracted from the <background_grid>. Alternatively, give a\n");
		mexPrintf ("\t\tsecond xyz-file in the command line. This will than be interpolated before use.\n");

		mexPrintf ("\t-M<s|z> Select interpolation algorithm. -Mz -> zgrid. -Ms -> surface [Default]\n");
		mexPrintf ("\t-T adds Tension to the gridding equation; use a value between 0 and 1.\n");
		mexPrintf ("\t\tOr a value between  0 and Inf if the -Mz was used.\n");
		mexPrintf ("\t\tdefault = 0 gives minimum curvature (smoothest; bicubic) solution.\n");
		mexPrintf ("\t-W<scale> Sets the width of the gaussian weighting function in terms of\n");
		mexPrintf ("\t\tthe grid spacing [Default: scale = 1.0].\n");
		mexPrintf ("\t\tDefault is 3 input columns.\n\n");
		return;
	}

	if (!got_R) {
		mexPrintf ("GMTMBGRID SYNTAX ERROR:  Must specify -R option\n");
		error++;
	}
	if (h.x_inc <= 0.0 || h.y_inc <= 0.0) {
		mexPrintf ("GMTMBGRID SYNTAX ERROR -I option.  Must specify positive increment(s)\n");
		error++;
	}
	if (n_files > 1 && backgroundfile) {
		mexPrintf ("GMTMBGRID ERROR. Make up your mind. Provide ONLY a background input file OR a background grid.\n");
		error++;
	}
	if (tension > 1 && interp_method == INTERP_SURFACE) {
		mexPrintf ("GMTMBGRID WARNING: option -T. Tension in the Surface algorithm must be between 0 and 1. Reset to 0.\n");
		tension = 0.;
	}

	if (error) return;



	/* --------------------  Test input numeric data consistency  ----------------------- */
	if (!mxIsNumeric(prhs[0]) || !mxIsNumeric(prhs[1]) || !mxIsNumeric(prhs[2]))
		mexErrMsgTxt("GMTMBGRID ERROR: first 3 input args must contain the x,y,z triplets to interpolate.\n");

	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("GMTMBGRID ERROR: (x,y,z) must be all doubles (sorry).\n");

	n_pts = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));
	k = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));
	n = MAX(mxGetM(prhs[2]), mxGetN(prhs[2]));

	if ( n_pts != k || n_pts != n)
		mexErrMsgTxt("GMTMBGRID ERROR: (x,y,z) triplets must all have number of points.\n");

	n_files = 1;

	if (mxIsNumeric(prhs[3])) {		/* We have a background dataset */
		if (!mxIsNumeric(prhs[4]) || !mxIsNumeric(prhs[5]))
			mexErrMsgTxt("GMTMBGRID ERROR: background data must come on the form of a secon x, y,z triplet.\n");

		n_pts2 = MAX(mxGetM(prhs[3]), mxGetN(prhs[3]));
		k = MAX(mxGetM(prhs[4]), mxGetN(prhs[4]));
		n = MAX(mxGetM(prhs[5]), mxGetN(prhs[5]));
		
		if ( n_pts2 != k || n_pts2 != n)
			mexErrMsgTxt("GMTMBGRID ERROR: background data (x,y,z) triplets must all have number of points.\n");

		n_files = 2;
	}
	/* --------------------------------------------------------------------------------------*/

	/* Since the mbgrid version used as base for gmtmbgrid was bugged regarding the
	   option of interpolating all nodes in the grid, I use this as a patch */
	if (clipmode == MBGRID_INTERP_ALL)
		clip = MAX(h.nx,h.ny) + 1;

	/* User selected a background file/grid but didn't specify -C option */
	if (clipmode == MBGRID_INTERP_ALL && (backgroundfile != NULL))
		clipmode = MBGRID_INTERP_NEAR;

	/* define NaN in case it's needed */
	if (use_NaN == TRUE)
		outclipvalue = (float)mxGetNaN();


	/* calculate other grid properties */
	factor = 4.0/(scale*scale*h.x_inc*h.y_inc);
	offx = offy = 0;
	if (extend > 0.0) {
		offx = (int) (extend*h.nx);
		offy = (int) (extend*h.ny);
	}
	xtradim = (int)scale + 2;
	gxdim = h.nx + 2*offx;
	gydim = h.ny + 2*offy;
	wbnd[0] = h.x_min - offx*h.x_inc;
	wbnd[1] = h.x_max + offx*h.x_inc;
	wbnd[2] = h.y_min - offy*h.y_inc;
	wbnd[3] = h.y_max + offy*h.y_inc;

	/* get data input bounds in lon lat */
	bounds[0] = wbnd[0];
	bounds[1] = wbnd[1];
	bounds[2] = wbnd[2];
	bounds[3] = wbnd[3];
		
	/* extend the bounds slightly to be sure no data gets missed */
	xx = MIN(0.05*(bounds[1] - bounds[0]), 0.1);
	yy = MIN(0.05*(bounds[3] - bounds[2]), 0.1);
	bounds[0] -= xx;
	bounds[1] += xx;
	bounds[2] -= yy;
	bounds[3] += yy;
	
	/* check interpolation parameters */
	if ((clipmode == MBGRID_INTERP_GAP || clipmode == MBGRID_INTERP_NEAR) && clip > h.nx && clip > h.ny)
		clipmode = MBGRID_INTERP_ALL;

	/***** do weighted mean gridding *****/
	/* allocate memory for arrays */
	if (offx == offy == 0) {
		plhs[0] = mxCreateNumericMatrix (h.ny,h.nx,mxDOUBLE_CLASS,mxREAL);
		grid = (double *)mxGetData(plhs[0]);
		alloc_out = TRUE;
	}
	else
		grid = (double *) mxCalloc ((size_t)(gxdim * gydim), sizeof(double));

	norm = (double *) mxCalloc ((size_t)(gxdim * gydim), sizeof(double));
	cnt = (int *) mxCalloc ((size_t)(gxdim * gydim), sizeof(int));
	num = (int *) mxCalloc ((size_t)(gxdim * gydim), sizeof(int));


	in0 = mxGetPr(prhs[0]);
	in1 = mxGetPr(prhs[1]);
	in2 = mxGetPr(prhs[2]);

	/* Read in data */
	ndata = 0;

	for (n = 0; n < n_pts; n++) {

		if (mxIsNaN (in2[n])) continue;

		/* get position in grid */
		ix = (int)((in0[n] - wbnd[0] + 0.5*h.x_inc)/h.x_inc);
		iy = (int)((in1[n] - wbnd[2] + 0.5*h.y_inc)/h.y_inc);

		/* process the data */
		if (ix >= -xtradim && ix < gxdim + xtradim && iy >= -xtradim && iy < gydim + xtradim) {
			ix1 = MAX(ix - xtradim, 0);
			ix2 = MIN(ix + xtradim, gxdim - 1);
			iy1 = MAX(iy - xtradim, 0);
			iy2 = MIN(iy + xtradim, gydim - 1);
			for (ii = ix1; ii <= ix2; ii++) {
				xx = wbnd[0] + ii*h.x_inc - in0[n];
				xx2 = xx * xx;
				for (jj = iy1; jj <= iy2; jj++) {
					kgrid = ii*gydim + jj;
					yy = wbnd[2] + jj*h.y_inc - in1[n];
					weight = exp(-(xx2 + yy*yy)*factor);
					norm[kgrid] += weight;
					grid[kgrid] += weight*in2[n];
					num[kgrid]++;
					if (ii == ix && jj == iy)
						cnt[kgrid]++;
				}
			}
			ndata++;
		}
		else if (ix >= 0 && ix < gxdim && iy >= 0 && iy < gydim) {
			kgrid = ix*gydim + iy;
			if (num[kgrid] <= 0) {
				norm[kgrid] = 1.0;
				grid[kgrid] = in2[n];
				num[kgrid] = 1;
				cnt[kgrid] = 1;
			}
			ndata++;
		}
	}	/* end of for loop */

	/* now loop over all points in the output grid */
	nbinset = 0;
	nbinspline = 0;
	for (i = 0; i < gxdim; i++) {
		for (j = 0; j < gydim; j++) {
			kgrid = i * gydim + j;
			if (cnt[kgrid] > 0) {
				grid[kgrid] = grid[kgrid]/norm[kgrid];
				nbinset++;
			}
			else
				grid[kgrid] = clipvalue;
		}
	}
	mxFree ((void *)norm);	/* Not needed anymore - JL */
	mxFree ((void *)cnt);
				
	/* if clip set do smooth interpolation */
	if (clipmode != MBGRID_INTERP_NONE && clip > 0 && nbinset > 0) {
		/* set up data vector */
		ndata = 0;
		if (border > 0.0)
			ndata = 2*gxdim + 2*gydim - 2;
		for (i = 0; i < gxdim; i++) {
			for (j = 0; j < gydim; j++) {
				kgrid = i * gydim + j;
				if (grid[kgrid] < clipvalue) ndata++;
			}
		}
		sgrid = (float *) mxCalloc ((size_t)(gxdim*gydim), sizeof(float));

		/* allocate and initialize sgrid */
		if (interp_method == INTERP_SURFACE) {
			sxdata = (float *) mxCalloc ((size_t)(ndata), sizeof(float));
			sydata = (float *) mxCalloc ((size_t)(ndata), sizeof(float));
			szdata = (float *) mxCalloc ((size_t)(ndata), sizeof(float));

			/* get points from grid */
			ndata = 0;
			for (i=0;i<gxdim;i++)
				for (j=0;j<gydim;j++) {
					kgrid = i * gydim + j;
					if (grid[kgrid] < clipvalue) {
						sxdata[ndata] = (float)(h.x_min + h.x_inc*i);
						sydata[ndata] = (float)(h.y_min + h.y_inc*j);
						szdata[ndata] = (float)(grid[kgrid]);
						ndata++;
					}
				}
			/* if desired set border */
			if (border > 0.0) {
				for (i=0;i<gxdim;i++) {
					j = 0;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sxdata[ndata] = (float)(h.x_min + h.x_inc*i);
						sydata[ndata] = (float)(h.y_min + h.y_inc*j);
						szdata[ndata] = (float)(border);
						ndata++;
					}
					j = gydim - 1;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sxdata[ndata] = (float)(h.x_min + h.x_inc*i);
						sydata[ndata] = (float)(h.y_min + h.y_inc*j);
						szdata[ndata] = (float)(border);
						ndata++;
					}
				}
				for (j=1;j<gydim-1;j++) {
					i = 0;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sxdata[ndata] = (float)(h.x_min + h.x_inc*i);
						sydata[ndata] = (float)(h.y_min + h.y_inc*j);
						szdata[ndata] = (float)(border);
						ndata++;
					}
					i = gxdim - 1;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sxdata[ndata] = (float)(h.x_min + h.x_inc*i);
						sydata[ndata] = (float)(h.y_min + h.y_inc*j);
						szdata[ndata] = (float)(border);
						ndata++;
					}
				}
			}

			/* do the interpolation */
			if (verbose)
				mexPrintf("\nDoing Surface spline interpolation with %d data points...\n",ndata);
			mb_surface(verbose,ndata,sxdata,sydata,szdata,
				h.x_min,h.x_max,h.y_min,h.y_max,h.x_inc,h.y_inc, tension,sgrid);
		}
		else {		/* ZGRID */
			sdata = (float *) mxCalloc ((size_t)(3*ndata), sizeof(float));
			work1 = (float *) mxCalloc ((size_t)(ndata), sizeof(float));
			work2 = (int *) mxCalloc ((size_t)(ndata), sizeof(int));
			work3 = (int *) mxCalloc ((size_t)(gxdim+gydim), sizeof(int));

			/* get points from grid */
			ndata = 0;
			for (i=0;i<gxdim;i++)
				for (j=0;j<gydim;j++) {
					kgrid = i * gydim + j;
					if (grid[kgrid] < clipvalue) {
						sdata[ndata++] = (float)(h.x_min + h.x_inc*i);
						sdata[ndata++] = (float)(h.y_min + h.y_inc*j);
						sdata[ndata++] = (float)(grid[kgrid]);
					}
				}
			/* if desired set border */
			if (border > 0.0) {
				for (i = 0; i < gxdim; i++) {
					j = 0;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sdata[ndata++] = (float)(h.x_min + h.x_inc*i);
						sdata[ndata++] = (float)(h.y_min + h.y_inc*j);
						sdata[ndata++] = (float)(border);
					}
					j = gydim - 1;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sdata[ndata++] = (float)(h.x_min + h.x_inc*i);
						sdata[ndata++] = (float)(h.y_min + h.y_inc*j);
						sdata[ndata++] = (float)(border);
					}
				}
				for (j=1;j<gydim-1;j++) {
					i = 0;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sdata[ndata++] = (float)(h.x_min + h.x_inc*i);
						sdata[ndata++] = (float)(h.y_min + h.y_inc*j);
						sdata[ndata++] = (float)(border);
					}
					i = gxdim - 1;
					kgrid = i * gydim + j;
					if (grid[kgrid] >= clipvalue) {
						sdata[ndata++] = (float)(h.x_min + h.x_inc*i);
						sdata[ndata++] = (float)(h.y_min + h.y_inc*j);
						sdata[ndata++] = (float)(border);
					}
				}
			}
			ndata = ndata/3;

			/* do the interpolation */
			if (verbose)
				mexPrintf("\nDoing Zgrid spline interpolation with %d data points...\n",ndata);
			cay = (float)tension;
			xmin = (float)h.x_min;
			ymin = (float)h.y_min;
			ddx = (float)h.x_inc;
			ddy = (float)h.y_inc;
			mb_zgrid(sgrid,&gxdim,&gydim,&xmin,&ymin, &ddx,&ddy,sdata,&ndata, work1,work2,work3,&cay,&clip);
		}

		if (verbose) {
			if (clipmode == MBGRID_INTERP_GAP)
				mexPrintf("Applying spline interpolation to fill gaps of %d cells or less...\n",clip);
			else if (clipmode == MBGRID_INTERP_NEAR)
				mexPrintf("Applying spline interpolation to fill %d cells from data...\n",clip);
			else if (clipmode == MBGRID_INTERP_ALL)
				mexPrintf("Applying spline interpolation to fill all undefined cells in the grid...\n");
		}

		/* translate the interpolation into the grid array filling only data gaps */
		if (clipmode == MBGRID_INTERP_GAP) {
			for (i = 0; i < gxdim; i++) {
				for (j = 0; j < gydim; j++) {
					kgrid = i * gydim + j;
					if (interp_method == INTERP_SURFACE)
						kint = i + (gydim -j - 1) * gxdim;
					else
						kint = i + j*gxdim;
					num[kgrid] = FALSE;
					if (grid[kgrid] >= clipvalue && sgrid[kint] < zflag) {
						/* initialize direction mask of search */
						for (ii = 0; ii < 9; ii++)
							dmask[ii] = FALSE;

						/* loop over rings around point, starting close */
						for (ir = 0; ir <= clip && num[kgrid] == FALSE; ir++) {
							/* set bounds of search */
							i1 = MAX(0, i - ir);
							i2 = MIN(gxdim - 1, i + ir);
							j1 = MAX(0, j - ir);
							j2 = MIN(gydim - 1, j + ir);
							jj = j1;
							for (ii = i1; ii <= i2 && num[kgrid] == FALSE;ii++) {
								if (grid[ii*gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = irint((ii - i)/r) + 1;
									jjj = irint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = TRUE;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										num[kgrid] = TRUE;
								}
							}
				      
							jj = j2;
							for (ii = i1; ii <= i2 && num[kgrid] == FALSE;ii++) {
								if (grid[ii*gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = irint((ii - i)/r) + 1;
									jjj = irint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = TRUE;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										num[kgrid] = TRUE;
								}
							}
				      
							ii = i1;
							for (jj=j1;jj<=j2 && num[kgrid] == FALSE;jj++) {
								if (grid[ii*gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = irint((ii - i)/r) + 1;
									jjj = irint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = TRUE;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										num[kgrid] = TRUE;
								}
							}
				      
							ii = i2;
							for (jj=j1;jj<=j2 && num[kgrid] == FALSE;jj++) {
								if (grid[ii*gydim+jj] < clipvalue) {
									r = sqrt((double)((ii-i)*(ii-i) + (jj-j)*(jj-j)));
									iii = irint((ii - i)/r) + 1;
									jjj = irint((jj - j)/r) + 1;
									kkk = iii * 3 + jjj;
									dmask[kkk] = TRUE;
									if ((dmask[0] && dmask[8]) || (dmask[3] && dmask[5])
										|| (dmask[6] && dmask[2]) || (dmask[1] && dmask[7]))
										num[kgrid] = TRUE;
								}
							}
						}
					}	/* end if */
				}
			}
			for (i=0;i<gxdim;i++) {
				for (j=0;j<gydim;j++) {
					kgrid = i * gydim + j;
					if (interp_method == INTERP_SURFACE)
						kint = i + (gydim -j - 1) * gxdim;
					else
						kint = i + j*gxdim;
					if (num[kgrid] == TRUE) {
						grid[kgrid] = sgrid[kint];
						nbinspline++;
					}
				}
			}
		}

		/* translate the interpolation into the grid array filling by proximity */
		else if (clipmode == MBGRID_INTERP_NEAR) {
			for (i = 0; i < gxdim; i++) {
				for (j = 0; j < gydim; j++) {
					kgrid = i * gydim + j;
					if (interp_method == INTERP_SURFACE)
						kint = i + (gydim -j - 1) * gxdim;
					else
						kint = i + j*gxdim;
					num[kgrid] = FALSE;
					if (grid[kgrid] >= clipvalue && sgrid[kint] < zflag) {
						/* loop over rings around point, starting close */
						for (ir = 0; ir <= clip && num[kgrid] == FALSE; ir++) {
							/* set bounds of search */
							i1 = MAX(0, i - ir);
							i2 = MIN(gxdim - 1, i + ir);
							j1 = MAX(0, j - ir);
							j2 = MIN(gydim - 1, j + ir);
				      
							jj = j1;
							for (ii=i1;ii<=i2 && num[kgrid] == FALSE;ii++) {
								if (grid[ii*gydim+jj] < clipvalue)
									num[kgrid] = TRUE;
							}
				      
							jj = j2;
							for (ii=i1;ii<=i2 && num[kgrid] == FALSE;ii++) {
								if (grid[ii*gydim+jj] < clipvalue)
									num[kgrid] = TRUE;
							}
				      
							ii = i1;
							for (jj=j1;jj<=j2 && num[kgrid] == FALSE;jj++) {
								if (grid[ii*gydim+jj] < clipvalue)
									num[kgrid] = TRUE;
							}
				      
							ii = i2;
							for (jj=j1;jj<=j2 && num[kgrid] == FALSE;jj++) {
								if (grid[ii*gydim+jj] < clipvalue)
									num[kgrid] = TRUE;
							}
						}
					}
				}
			}
			for (i=0;i<gxdim;i++) {
				for (j=0;j<gydim;j++) {
					kgrid = i * gydim + j;
					if (interp_method == INTERP_SURFACE)
						kint = i + (gydim -j - 1) * gxdim;
					else
						kint = i + j*gxdim;
					if (num[kgrid] == TRUE) {
						grid[kgrid] = sgrid[kint];
						nbinspline++;
					}
				}
			}
		}

		/* translate the interpolation into the grid array filling all empty bins */
		else {
			for (i=0;i<gxdim;i++) {
				for (j=0;j<gydim;j++) {
					kgrid = i * gydim + j;
					if (interp_method == INTERP_SURFACE)
						kint = i + (gydim -j - 1) * gxdim;
					else
						kint = i + j*gxdim;
					if (grid[kgrid] >= clipvalue && sgrid[kint] < zflag) {
						grid[kgrid] = sgrid[kint];
						nbinspline++;
					}
				}
			}
		}
			
		/* deallocate the interpolation arrays */
		if (interp_method == INTERP_SURFACE) {
			mxFree ((void *)sxdata);	mxFree ((void *)sydata);
			mxFree ((void *)szdata);
		}
		else {
			mxFree ((void *)sdata);		mxFree ((void *)work1);
			mxFree ((void *)work2);		mxFree ((void *)work3);
		}
		mxFree ((void *)sgrid);
	}
	mxFree ((void *)num);		/* Not needed anymore - JL */

	if (n_files == 2) {	/* Use background data from second input numeric package  */

		if (verbose) mexPrintf("Extracting background from file ...\n");

		/* allocate and initialize background data array */
		if (interp_method == INTERP_SURFACE) {
			bxdata = (float *) mxCalloc ((size_t)(n_pts2), sizeof(float));
			bydata = (float *) mxCalloc ((size_t)(n_pts2), sizeof(float));
			bzdata = (float *) mxCalloc ((size_t)(n_pts2), sizeof(float));
		}
		else
			bdata = (float *) mxCalloc ((size_t)(3*n_pts2), sizeof(float));

		in0 = mxGetPr(prhs[3]);
		in1 = mxGetPr(prhs[4]);
		in2 = mxGetPr(prhs[5]);

		nbackground = 0;
		for (n = 0; n < n_pts2; n++) {

			if (mxIsNaN (in2[n])) continue;

			if (in0[n] < h.x_min || in0[n] > h.x_max || in1[n] < h.y_min || in1[n] > h.y_max) continue;
			if (interp_method == INTERP_SURFACE) {
				bxdata[nbackground] = (float) in0[n];
				bydata[nbackground] = (float) in1[n];
				bzdata[nbackground] = (float) in2[n];
				nbackground++;
			}
			else {
				bdata[nbackground*3]   = (float) in0[n];
				bdata[nbackground*3+1] = (float) in1[n];
				bdata[nbackground*3+2] = (float) in2[n];
				nbackground++;
			}
		}
	}

	/* if grdrasterid set read background data extracted using grdraster then interpolate it onto internal grid */
	if (n_files == 2 && nbackground > 0) {		

		/* allocate grid and work arrays */
		sgrid = (float *) mxCalloc ((size_t)(gxdim*gydim), sizeof(float));
		if (interp_method == INTERP_ZGRID) {
			work1 = (float *) mxCalloc ((size_t)(nbackground), sizeof(float));
			work2 = (int *) mxCalloc ((size_t)(nbackground), sizeof(int));
			work3 = (int *) mxCalloc ((size_t)(gxdim+gydim), sizeof(int));
		}

		/* do the interpolation */
		if (verbose) mexPrintf("\nDoing spline interpolation with %d data points from background...\n",nbackground);
		if (interp_method == INTERP_SURFACE)
			mb_surface(verbose,nbackground,bxdata,bydata,bzdata,
				h.x_min,h.x_max,h.y_min,h.y_max,h.x_inc,h.y_inc, tension,sgrid);
		else {
			cay = (float)tension;
			xmin = (float)h.x_min;
			ymin = (float)h.y_min;
			ddx = (float)h.x_inc;
			ddy = (float)h.y_inc;
			clip = MAX(gxdim,gydim);
			mb_zgrid(sgrid,&gxdim,&gydim,&xmin,&ymin, &ddx,&ddy,bdata,
				&nbackground,work1,work2,work3,&cay,&clip);
		}

		/* translate the interpolation into the grid array 
		    - interpolate only to fill a data gap */
		for (i=0;i<gxdim;i++) {
			for (j=0;j<gydim;j++) {
				kgrid = i * gydim + j;
				if (interp_method == INTERP_SURFACE)
					kint = i + (gydim -j - 1) * gxdim;
				else
					kint = i + j*gxdim;
				if (grid[kgrid] >= clipvalue && sgrid[kint] < zflag) {
					grid[kgrid] = sgrid[kint];
					nbinbackground++;
				}
			}
		}
		if (interp_method == INTERP_SURFACE) {
			mxFree ((void *)bxdata);	mxFree ((void *)bydata);
			mxFree ((void *)bzdata);
		}
		else {
			mxFree ((void *)work1);	mxFree ((void *)work2);	
			mxFree ((void *)work3);	mxFree ((void *)bdata);
		}
		mxFree ((void *)sgrid);
	}

	/* get min max of data */
	zclip = clipvalue;
	zmin = zclip;
	zmax = zclip;
	for (i=0;i<gxdim;i++)
		for (j=0;j<gydim;j++) {
			kgrid = i * gydim + j;
			if (zmin == zclip && grid[kgrid] < zclip)
				zmin = grid[kgrid];
			if (zmax == zclip && grid[kgrid] < zclip)
				zmax = grid[kgrid];
			if (grid[kgrid] < zmin && grid[kgrid] < zclip)
				zmin = grid[kgrid];
			if (grid[kgrid] > zmax && grid[kgrid] < zclip)
				zmax = grid[kgrid];
		}
	if (zmin == zclip) zmin = 0.0;
	if (zmax == zclip) zmax = 0.0;

	if (verbose) {
		nbinzero = gxdim*gydim - nbinset - nbinspline - nbinbackground;
		mexPrintf("\nTotal number of bins:            %d\n",gxdim*gydim);
		mexPrintf("Bins set using data:             %d\n",nbinset);
		mexPrintf("Bins set using interpolation:    %d\n",nbinspline);
		mexPrintf("Bins set using background:       %d\n",nbinbackground);
		mexPrintf("Bins not set:                    %d\n",nbinzero);
		mexPrintf("Minimum value: %10.2f   Maximum value: %10.2f\n", zmin,zmax);
	}

	/* OUT - It is already in the Matlab orientation */

	if (alloc_out) {	/* plhs[0] is already allocated and grid points to it */
		for (i = 0; i < h.nx * h.ny; i++) {
			if (grid[i] == clipvalue)
				grid[i] = outclipvalue;
		}
	}
	else {			/* In this case output is SINGLE (a bit inconsistent but memory saver) */
		plhs[0] = mxCreateNumericMatrix (h.ny,h.nx,mxSINGLE_CLASS,mxREAL);
		pdata = (float *)mxGetData(plhs[0]);
		for (i = 0; i < h.nx; i++) {
			k = (i + offx)*gydim + offy;
			for (j = 0; j < h.ny; j++) {
				kgrid = k + j;
				pdata[kgrid] = (float) grid[kgrid];
				if (grid[kgrid] == clipvalue)
					pdata[kgrid] = outclipvalue;
			}
		}
		mxFree ((void *)grid);
	}

	if (nlhs == 2) {	/* User also wants the header */
		plhs[1] = mxCreateDoubleMatrix (1, 9, mxREAL);
		info = mxGetPr (plhs[1]);
		info[0] = h.x_min;
		info[1] = h.x_max;
		info[2] = h.y_min;
		info[3] = h.y_max;
		info[4] = h.z_min;
		info[5] = h.z_max;
		info[6] = h.node_offset;
		info[7] = h.x_inc;
		info[8] = h.y_inc;
	}
}


/*--------------------------------------------------------------------
 *    The MB-system:	mb_zgrid.c	    4/25/95
 *    $Id: mb_zgrid.c,v 5.0 2000/12/01 22:53:59 caress Exp $
 *
 *    Copyright (c) 1993, 1994, 1995, 2000 by
 *    David W. Caress (caress@mbari.org)
 *      Monterey Bay Aquarium Research Institute
 *      Moss Landing, CA 95039
 *    and Dale N. Chayes (dale@ldeo.columbia.edu)
 *      Lamont-Doherty Earth Observatory
 *      Palisades, NY 10964
 *
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/* 
 * This is a function to generate thin plate spline interpolation
 * of a data field. This code originated as fortran in the
 * 1960's and was used routinely at the Institute of 
 * Geophysics and Planetary Physics at the Scripps Institution
 * of Oceanography through the 1970's and 1980's. The Fortran
 * code was obtained from Professory Robert L. Parker at
 * IGPP in 1989.
 * It was converted to C in 1995 for use with the MB-System
 * software package.
 * 
 * The nature of the interpolation is controlled by the
 * parameters cay and nrng: cay sets the tension of the
 * interpolation such that cay=0.0 yields a pure Laplace
 * (minimum curvature) solution and cay=infinity yields
 * a pure thin plate spline solution. A cay=1e10 value
 * has commonly been used to yield spline solutions.
 * The nrng value sets the number of grid spaces from
 * data that will be interpolated; if nrng exceeds the
 * maximum dimension of the grid then the entire grid
 * will be interpolated.
 * 
 * The input parameters are:
 *     nx,ny = max subscripts of z in x and y directions . 
 *     x1,y1 = coordinates of z(1,1) 
 *     dx,dy = x and y increments . 
 *     xyz(3,*) = array giving x-y position and hgt of each data point. 
 *     n = length of xyz series. 
 *     zpij[n] = float work array
 *     knxt[n] = int work array
 *     imnew[MAX(nx, ny)+1] = int work array
 *     cay = k = amount of spline eqn (between 0 and inf.) 
 *     nrng...grid points more than nrng grid spaces from the nearest 
 *            data point are set to undefined. 
 *
 * Author:	Unknown, but "jdt", "ian crain",  and "dr t murty"
 *              obviously contributed.
 * Hacker:	D. W. Caress
 * Date:	April 25, 1995
 *
 *     The following are the original comments from the Fortran code:
 *
 *     sets up square grid for contouring , given arbitrarily placed 
 *     data points. laplace interpolation is used. 
 *     the method used here was lifted directly from notes left by 
 *     mr ian crain formerly with the comp.science div. 
 *     info on relaxation soln of laplace eqn supplied by dr t murty. 
 *     fortran ii   oceanography/emr   dec/68   jdt 
 *
 *     z = 2-d array of hgts to be set up. points outside region to be 
 *     contoured should be initialized to 10**35 . the rest should be 0.0 
 *
 *     modification feb/69   to get smoother results a portion of the 
 *     beam eqn  was added to the laplace eqn giving 
 *     delta2x(z)+delta2y(z) - k(delta4x(z)+delta4y(z)) = 0 . 
 *     k=0 gives pure laplace solution.  k=infinity gives pure spline 
 *     solution. 
 *
 *     nx,ny = max subscripts of z in x and y directions . 
 *     x1,y1 = coordinates of z(1,1) 
 *     dx,dy = x and y increments . 
 *     xyz(3,*) = array giving x-y position and hgt of each data point. 
 *     n = length of xyz series. 
 *     cay = k = amount of spline eqn (between 0 and inf.) 
 *     nrng...grid points more than nrng grid spaces from the nearest 
 *            data point are set to undefined. 
 *
 *     modification dec23/69   data pts no longer moved to grid pts. 
 *
 *     modification feb/85  common blocks work1 and work2 replaced by 
 *     dimension statement and parameters nwork, mwork introduced. 
 *
 *     modification feb/90  nwork and mwork replaced by maxdat and maxdim 
 *     for compatibility with command driven interface program 
 *     David W. Caress
 *
 *     zgrid.c -- translated by f2c (version 19950314) from zgrid.f.
 *     Work arrays zpij[n], knxt[n], imnew[MAX(nx, ny)+1] are now
 *     passed into the function.
 *     David W. Caress
 *     April 25,  1995
 *--------------------------------------------------------------------*/

/*----------------------------------------------------------------------- */
int mb_zgrid(float *z, int *nx, int *ny, float *x1, float *y1, float *dx, float *dy, float *xyz, 
		int *n, float *zpij, int *knxt, int *imnew, float *cay, int *nrng) {
    /* System generated locals */
    int z_dim1, z_offset, i__2;

    /* Local variables */
    int i, j, k, iter, nnew, itmax, jmnew;
    int kk, im, jm, npg, npt;
    float delz, zijn, zmin, zimm, zmax, zjmm, zipp, zjpp;
    float root, zsum, zpxy, a, b, c, d;
    float x, y, zbase, relax, delzm, derzm, dzmax, dzrms;
    float dzrms8, z00, dz, ze, hrange, zn, zs, zw, zrange, dzmaxf, 
	    relaxn, rootgs, dzrmsp, big, abz;
    float eps, zim, zjm, wgt, zip, zjp, tpy, zxy;

    /* Parameter adjustments */
    z_dim1 = *nx;
    z_offset = z_dim1 + 1;
    z -= z_offset;
    xyz -= 4;
    

    /* Function Body */
    itmax = 100;
    eps = (float).002;
    big = (float)9e29;

/*     get zbase which will make all zp values positive by 20*(zmax-zmin) 
 * ********************************************************************** 
*/

    zmin = xyz[6];
    zmax = xyz[6];
    for (k = 2; k <= *n; ++k) {
	if (xyz[k * 3 + 3] > zmax)
		zmax = xyz[k * 3 + 3];

	if (xyz[k * 3 + 3] < zmin)
		zmin = xyz[k * 3 + 3];
    }
    zrange = zmax - zmin;
    zbase = (float)(zrange * 20. - zmin);
    hrange = MIN(*dx * (*nx - 1), *dy * (*ny - 1));
    derzm = (float)(zrange * 2. / hrange);

/*     set pointer array knxt */
/* ********************************************************************** 
*/

    for (kk = 1; kk <= *n; ++kk) {
	k = *n + 1 - kk;
	knxt[k - 1] = 0;
	i = (int)((xyz[k * 3 + 1] - *x1) / *dx + 1.5);
	if (i * (*nx + 1 - i) <= 0)
		continue;

	j = (int)((xyz[k * 3 + 2] - *y1) / *dy + 1.5);
	if (j * (*ny + 1 - j) <= 0)
		continue;

	if (z[i + j * z_dim1] >= big)
		continue;

	knxt[k - 1] = *n + 1;

	if (z[i + j * z_dim1] > 0.)
		knxt[k - 1] = (int)(z[i + j * z_dim1] + .5);

	z[i + j * z_dim1] = (float) k;
    }

/*     affix each data point zp to its nearby grid point.  take avg zp if 
 *     more than one zp nearby the grid point. add zbase and complement. 
 * ********************************************************************** */

    for (k = 1; k <= *n; ++k) {
	if (knxt[k - 1] <= 0)
		continue;

	npt = 0;
	zsum = (float)0.;
	i = (int)((xyz[k * 3 + 1] - *x1) / *dx + 1.5);
	j = (int)((xyz[k * 3 + 2] - *y1) / *dy + 1.5);
	kk = k;
L70:
	++npt;
	zsum += xyz[kk * 3 + 3];
	knxt[kk - 1] = -knxt[kk - 1];
	kk = -knxt[kk - 1];
	if (kk <= *n)
		goto L70;
	z[i + j * z_dim1] = (float)(-zsum / npt - zbase);
    }

/*     initially set each unset grid point to value of nearest known pt. 
 * ********************************************************************** 
*/

    for (i = 1; i <= *nx; ++i) {
	for (j = 1; j <= *ny; ++j) {
	    if (z[i + j * z_dim1] != 0.)
		continue;
	    z[i + j * z_dim1] = (float)-1e35;
	}
    }
    for (iter = 1; iter <= *nrng; ++iter) {
	nnew = 0;
	for (i = 1; i <= *nx; ++i) {
	    for (j = 1; j <= *ny; ++j) {
		if (z[i + j * z_dim1] + big >= 0.)
		    goto L192;

		if (j - 1 <= 0)
		    goto L162;

		if (jmnew <= 0) goto L154;
		else
		    goto L162;

L154:
		zijn = (float)fabs(z[i + (j - 1) * z_dim1]);
		if (zijn >= big) goto L162;
		else
		    goto L195;

L162:
		if (i <= 1)
		    goto L172;

		if (imnew[j - 1] <= 0) goto L164;
		else
		    goto L172;

L164:
		zijn = (float)fabs(z[i - 1 + j * z_dim1]);
		if (zijn >= big) goto L172;
		else
		    goto L195;

L172:
		if (j >= *ny)
		    goto L182;

		zijn = (float)fabs(z[i + (j + 1) * z_dim1]);
		if (zijn >= big) goto L182;
		else
		    goto L195;

L182:
		if (i >= *nx)
		    goto L192;

		zijn = (float)fabs(z[i + 1 + j * z_dim1]);
		if (zijn >= big) goto L192;
		else
		    goto L195;

L192:
		imnew[j - 1] = 0;
		jmnew = 0;
		goto L197;
L195:
		imnew[j - 1] = 1;
		jmnew = 1;
		z[i + j * z_dim1] = zijn;
		++nnew;
L197:
		;
	    }
	}
	if (nnew <= 0)
	    goto L200;
    }
L200:
    for (i = 1; i <= *nx; ++i) {
	for (j = 1; j <= *ny; ++j) {
	    abz = (float)fabs(z[i + j * z_dim1]);
	    if (abz >= big)
	    	z[i + j * z_dim1] = abz;

	}
    }

/*     improve the non-data points by applying point over-relaxation */
/*     using the laplace-spline equation  (carres method is used) */
/* ********************************************************************** 
*/

    dzrmsp = zrange;
    relax = (float)1.;
    for (iter = 1; iter <= itmax; ++iter) {
	dzrms = (float)0.;
	dzmax = (float)0.;
	npg = 0;
	for (i = 1; i <= *nx; ++i) {
	    for (j = 1; j <= *ny; ++j) {
		z00 = z[i + j * z_dim1];
		if (z00 >= big)
		    continue;

		if (z00 < 0)
		    continue;

		wgt = (float)0.;
		zsum = (float)0.;

		im = 0;
		if (i <= 1)
		    goto L570;

		zim = (float)fabs(z[i - 1 + j * z_dim1]);
		if (zim >= big)
		    goto L570;

		im = 1;
		wgt += 1.;
		zsum += zim;
		if (i <= 2)
		    goto L570;

		zimm = (float)fabs(z[i - 2 + j * z_dim1]);
		if (zimm >= big)
		    goto L570;

		wgt += *cay;
		zsum -= *cay * (zimm - zim * 2);
L570:
		if (*nx <= i)
		    goto L700;

		zip = (float)fabs(z[i + 1 + j * z_dim1]);
		if (zip >= big)
		    goto L700;

		wgt += 1.;
		zsum += zip;
		if (im <= 0)
		    goto L620;

		wgt += *cay * 4;
		zsum += *cay * 2 * (zim + zip);
L620:
		if (*nx - 1 - i <= 0)
		    goto L700;

		zipp = (float)fabs(z[i + 2 + j * z_dim1]);
		if (zipp >= big)
		    goto L700;

		wgt += *cay;
		zsum -= *cay * (zipp - zip * 2);
L700:

		jm = 0;
		if (j <= 1)
		    goto L1570;

		zjm = (float)fabs(z[i + (j - 1) * z_dim1]);
		if (zjm >= big)
		    goto L1570;

		jm = 1;
		wgt += (float)1.;
		zsum += zjm;
		if (j <= 2)
		    goto L1570;

		zjmm = (float)fabs(z[i + (j - 2) * z_dim1]);
		if (zjmm >= big)
		    goto L1570;

		wgt += *cay;
		zsum -= (*cay * (zjmm - zjm * 2));
L1570:
		if (*ny <= j)
		    goto L1700;

		zjp = (float)fabs(z[i + (j + 1) * z_dim1]);
		if (zjp >= big)
		    goto L1700;

		wgt += 1.;
		zsum += zjp;
		if (jm <= 0)
		    goto L1620;

		wgt += (*cay * 4);
		zsum += (*cay * 2 * (zjm + zjp));
L1620:
		if (*ny - 1 - j <= 0)
		    goto L1700;

		zjpp = (float)fabs(z[i + (j + 2) * z_dim1]);
		if (zjpp >= big)
		    goto L1700;

		wgt += *cay;
		zsum -= (*cay * (zjpp - zjp * 2));
L1700:

		dz = zsum / wgt - z00;
		++npg;
		dzrms += dz * dz;
		dzmax = MAX((float)fabs(dz), dzmax);
		z[i + j * z_dim1] = z00 + dz * relax;
	    }
	}


/*     shift data points zp progressively back to their proper places as */
/*     the shape of surface z becomes evident. */
/* ********************************************************************** */

	if (iter - iter / 10 * 10 != 0)
	    goto L3600;

	for (k = 1; k <= *n; ++k) {
	    knxt[k - 1] = (i__2 = knxt[k - 1], abs(i__2));
	    if (knxt[k - 1] <= 0)
		continue;

	    x = (xyz[k * 3 + 1] - *x1) / *dx;
	    i = (int)(x + 1.5);
	    x = (float)(x + 1. - i);
	    y = (xyz[k * 3 + 2] - *y1) / *dy;
	    j = (int)(y + 1.5);
	    y = (float)(y + 1. - j);
	    zpxy = xyz[k * 3 + 3] + zbase;
	    z00 = (float)fabs(z[i + j * z_dim1]);

	    zw = (float)1e35;
	    if (i <= 1)
		goto L3120;

	    zw = (float)fabs(z[i - 1 + j * z_dim1]);
L3120:
	    ze = (float)1e35;
	    if (i >= *nx)
		goto L3140;

	    ze = (float)fabs(z[i + 1 + j * z_dim1]);
L3140:
	    if (ze >= big) {
		goto L3150;
	    } else {
		goto L3160;
	    }
L3150:
	    if (zw >= big) {
		goto L3170;
	    } else {
		goto L3180;
	    }
L3160:
	    if (zw >= big) {
		goto L3190;
	    } else {
		goto L3200;
	    }
L3170:
	    ze = z00;
	    zw = z00;
	    goto L3200;
L3180:
	    ze = (float)(z00 * 2. - zw);
	    goto L3200;
L3190:
	    zw = (float)(z00 * 2. - ze);

L3200:
	    zs = (float)1e35;
	    if (j <= 1)
		goto L3220;

	    zs = (float)fabs(z[i + (j - 1) * z_dim1]);
L3220:
	    zn = (float)1e35;
	    if (j >= *ny)
		goto L3240;

	    zn = (float)fabs(z[i + (j + 1) * z_dim1]);
L3240:
	    if (zn >= big) {
		goto L3250;
	    } else {
		goto L3260;
	    }
L3250:
	    if (zs >= big) {
		goto L3270;
	    } else {
		goto L3280;
	    }
L3260:
	    if (zs >= big) {
		goto L3290;
	    } else {
		goto L3300;
	    }
L3270:
	    zn = z00;
	    zs = z00;
	    goto L3300;
L3280:
	    zn = (float)(z00 * 2. - zs);
	    goto L3300;
L3290:
	    zs = (float)(z00 * 2. - zn);

L3300:
	    a = (float)((ze - zw) * .5);
	    b = (float)((zn - zs) * .5);
	    c = (float)((ze + zw) * .5 - z00);
	    d = (float)((zn + zs) * .5 - z00);
	    zxy = z00 + a * x + b * y + c * x * x + d * y * y;
	    delz = z00 - zxy;
	    delzm = (float)(derzm * ((float)fabs(x) * *dx + (float)fabs(y) * *dy) * .8);
	    if (delz > delzm)
	        delz = delzm;
	    if (delz + delzm >= 0.)
		goto L3365;

	    delz = -delzm;
L3365:
	    zpij[k - 1] = zpxy + delz;
	}

	for (k = 1; k <= *n; ++k) {
	    if (knxt[k - 1] > 0) {
		npt = 0;
		zsum = (float)0.;
		i = (int)((xyz[k * 3 + 1] - *x1) / *dx + 1.5);
		j = (int)((xyz[k * 3 + 2] - *y1) / *dy + 1.5);
		kk = k;
L3420:
		++npt;
		zsum += zpij[kk - 1];
		knxt[kk - 1] = -knxt[kk - 1];
		kk = -knxt[kk - 1];
		if (kk <= *n)
			goto L3420;
		z[i + j * z_dim1] = -zsum / npt;
	    }
	}
L3600:

/*     test for convergence */
/* ****************************************************************** */
/* all grid points assigned */
	if (npg <= 1)
	    goto L4010;

	dzrms = (float)sqrt(dzrms / npg);
	root = dzrms / dzrmsp;
	dzrmsp = dzrms;
	dzmaxf = dzmax / zrange;
	if (iter - iter / 10 * 10 - 2 != 0)
	    goto L3715;

	dzrms8 = dzrms;
L3715:
	if (iter - iter / 10 * 10 != 0)
	    goto L4000;

	root = (float)sqrt(sqrt(sqrt(dzrms / dzrms8)));
	if (root >= .9999)
	    goto L4000;

	if (dzmaxf / (1. - root) <= eps)
	    goto L4010;

/*     improve the relaxation factor. */
/* ******************************************************************** */

	if ((iter - 20) * (iter - 40) * (iter - 60) != 0)
	    goto L4000;

	if (relax - 1. - root >= 0.)
	    goto L4000;

	tpy = (float)((root + relax - 1.) / relax);
	rootgs = tpy * tpy / root;
	relaxn = (float)(2. / (sqrt(1. - rootgs) + 1.));
	if (iter != 60) {
	    goto L3780;
	} else {
	    goto L3785;
	}
L3780:
	relaxn -= (float)((2. - relaxn) * .25);
L3785:
	relax = (float)MAX(relax,relaxn);
L4000:
	;
	if (verbose)
		if (iter % 10 == 0)
			mexPrintf("Iteration %d of a maximum of %d\r", iter, itmax); 
    }
L4010:
	if (verbose) mexPrintf("\n");

/*     remove zbase from array z and return. */
/* ********************************************************************** 
*/

    for (i = 1; i <= *nx; ++i) {
	for (j = 1; j <= *ny; ++j) {
	    if (z[i + j * z_dim1] < big)
	    	z[i + j * z_dim1] = (float)fabs(z[i + j * z_dim1]) - zbase;
	}
    }
    return 0;
}


/*--------------------------------------------------------------------
 *    The MB-system:	mb_surface.c	5/2/94
 *    $Id: mb_surface.c,v 5.2 2005/03/25 04:09:53 caress Exp $
 *
 *    Inclusion in MB-System:
 *    Copyright (c) 1994, 2003 by
 *    David W. Caress (caress@mbari.org)
 *      Monterey Bay Aquarium Research Institute
 *      Moss Landing, CA 95039
 *    and Dale N. Chayes (dale@ldeo.columbia.edu)
 *      Lamont-Doherty Earth Observatory
 *      Palisades, NY 10964
 *
 *    Algorithm and original code:
 *    Copyright (c) 1991 by P. Wessel and W. H. F. Smith
 *
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * SURFUNC is a function for gridding data using a minimum curvature
 * algorithm developed by W.H.F. Smith and P. Wessel.  The source
 * code below is almost entirely taken from the source code
 * "surface.c" distributed as part of the GMT-system by Wessel
 * and Smith. The MB-System Copyright notice above applies only
 * to the inclusion of this code in MB-System and changes made to
 * the code as part of that inclusion.  The algorithm and the
 * bulk of the code remains Copyright (c) by P. Wessel and W. H. F. Smith.
 *
 *--------------------------------------------------------------------*
 *
 * The original Smith and Wessel comments follow:
 *
 * surface.c:  a gridding program.
 * reads xyz triples and fits a surface to the data.
 * surface satisfies (1 - T) D4 z - T D2 z = 0,
 * where D4 is the biharmonic operator,
 * D2 is the Laplacian,
 * and T is a "tension factor" between 0 and 1.
 * End member T = 0 is the classical minimum curvature
 * surface.  T = 1 gives a harmonic surface.  Use T = 0.25
 * or so for potential data; something more for topography.
 *
 * Program includes overrelaxation for fast convergence and
 * automatic optimal grid factorization.
 *
 * See Smith & Wessel (Geophysics, 3, 293-305, 1990) for details.
 *
 * Authors: Walter H. F. Smith and Paul Wessel
 * Date: April, 1988.
 *
 *--------------------------------------------------------------------*
 *
 * Particulars regarding turning the program "surface" version 4.3
 * (revision of 26 February, 1992) into a function "surfunc" follow:
 *
 * Author:	D. W. Caress
 * Date:	May 2, 1994
 *
 * $Log: mb_surface.c,v $
 * Revision 5.2  2005/03/25 04:09:53  caress
 * Problems with global variables in mb_surface.c stomping on similarly named global variables in some programs has been fixed.
 *
 * Revision 5.1  2003/04/29 20:27:48  caress
 * Fixed multiple definitions of "error".
 *
 * Revision 5.0  2003/03/22 03:10:36  caress
 * Reinserting this code into MB-System for first time in years.
 *
 * Revision 4.2  1994/10/21  13:02:31  caress
 * Release V4.0
 *
 * Revision 4.1  1994/06/04  02:02:01  caress
 * Fixed several bugs and made some stylistic changes to
 * the output.  Changed the data input bounds to be much
 * larger than the working grid bounds.
 *
 * Revision 4.0  1994/05/05  20:30:06  caress
 * First cut. This code derived from GMT program surface by
 * Walter Smith and Paul Wessel.
 *
 *
 */

int mb_surface(int verbose, int ndat, float *xdat, float *ydat, float *zdat,
		double xxmin, double xxmax, double yymin, double yymax, double xxinc, double yyinc,
		double ttension, float *sgrid) {

	/* int	size_query = FALSE; */
	int	serror = FALSE;
	char	low[100], high[100];

	/* copy parameters */
	xmin = xxmin;
	xmax = xxmax;
	ymin = yymin;
	ymax = yymax;
	xinc = xxinc;
	yinc = yyinc;
	tension = ttension;
	total_iterations = 0;

	/* New in v4.3:  Default to unconstrained:  */
	set_low = set_high = 0; 
	
	if (xmin >= xmax || ymin >= ymax) serror = TRUE;
	if (xinc <= 0.0 || yinc <= 0.0) serror = TRUE;

	if (tension != 0.0) {
		boundary_tension = tension;
		interior_tension = tension;
	}
	relax_old = 1.0 - relax_new;

	nx = irint((xmax - xmin)/xinc) + 1;
	ny = irint((ymax - ymin)/yinc) + 1;
	mx = nx + 4;
	my = ny + 4;
	r_xinc = 1.0 / xinc;
	r_yinc = 1.0 / yinc;

	/* New stuff here for v4.3:  Check out the grid dimensions:  */
	grid = gcd_euclid(nx-1, ny-1);

	/*
	if (gmtdefs.verbose || size_query || grid == 1) fprintf (stderr, "W: %.3lf E: %.3lf S: %.3lf N: %.3lf nx: %d ny: %d\n",
		xmin, xmax, ymin, ymax, nx, ny);
	if (grid == 1) fprintf(stderr,"surface:  WARNING:  Your grid dimensions are mutually prime.\n");
	if (grid == 1 || size_query) suggest_sizes_for_surface(nx-1, ny-1);
	if (size_query) exit(0);
	*/

	/* New idea: set grid = 1, read data, setting index.  Then throw
		away data that can't be used in end game, constraining
		size of briggs->b[6] structure.  */
	
	grid = 1;
	set_grid_parameters();
	read_data(ndat,xdat,ydat,zdat);
	throw_away_unusables();
	remove_planar_trend();
	rescale_z_values();
	load_constraints(low, high);
	
	/* Set up factors and reset grid to first value  */
	
	grid = gcd_euclid(nx-1, ny-1);
	n_fact = get_prime_factors(grid, factors);
	set_grid_parameters();
	while ( block_nx < 4 || block_ny < 4 ) {
		smart_divide();
		set_grid_parameters();
	}
	set_offset();
	set_index();
	/* Now the data are ready to go for the first iteration.  */

	/* Allocate more space  */
	
	briggs = (struct SURFACE_BRIGGS *) mxCalloc ((size_t)npoints, sizeof(struct SURFACE_BRIGGS));
	iu = (char *) mxCalloc ((size_t)(mx * my), sizeof(char));
	u = (float *) mxCalloc ((size_t)(mx * my), sizeof(float));

	if (radius > 0) initialize_grid(); /* Fill in nodes with a weighted avg in a search radius  */

	if (verbose) mexPrintf("Grid\tMode\tIteration\tMax Change\tConv Limit\tTotal Iterations\n");
	
	set_coefficients();
	
	old_grid = grid;
	find_nearest_point ();
	iterate (1);
	 
	while (grid > 1) {
		smart_divide ();
		set_grid_parameters();
		set_offset();
		set_index ();
		fill_in_forecast ();
		iterate(0);
		old_grid = grid;
		find_nearest_point ();
		iterate (1);
	}
	
	if (verbose) check_errors ();

	replace_planar_trend();

	get_output(sgrid);

	mxFree ((void *) data);
	mxFree ((void *) briggs);
	mxFree ((void *) iu);
	mxFree ((void *) u);
	if (set_low) mxFree ((void *)lower);
	if (set_high) mxFree ((void *)upper);

	return(0);
}

void	set_coefficients() {

	double	e_4, loose, a0;
	
	loose = 1.0 - interior_tension;
	e_2 = epsilon * epsilon;
	e_4 = e_2 * e_2;
	eps_p2 = e_2;
	eps_m2 = 1.0/e_2;
	one_plus_e2 = 1.0 + e_2;
	two_plus_ep2 = 2.0 + 2.0*eps_p2;
	two_plus_em2 = 2.0 + 2.0*eps_m2;
	
	x_edge_const = 4 * one_plus_e2 - 2 * (interior_tension / loose);
	e_m2 = 1.0 / e_2;
	y_edge_const = 4 * (1.0 + e_m2) - 2 * (interior_tension * e_m2 / loose);

	
	a0 = 1.0 / ( (6 * e_4 * loose + 10 * e_2 * loose + 8 * loose - 2 * one_plus_e2) + 4*interior_tension*one_plus_e2);
	a0_const_1 = 2 * loose * (1.0 + e_4);
	a0_const_2 = 2.0 - interior_tension + 2 * loose * e_2;
	
	coeff[1][4] = coeff[1][7] = -loose;
	coeff[1][0] = coeff[1][11] = -loose * e_4;
	coeff[0][4] = coeff[0][7] = -loose * a0;
	coeff[0][0] = coeff[0][11] = -loose * e_4 * a0;
	coeff[1][5] = coeff[1][6] = 2 * loose * one_plus_e2;
	coeff[0][5] = coeff[0][6] = (2 * coeff[1][5] + interior_tension) * a0;
	coeff[1][2] = coeff[1][9] = coeff[1][5] * e_2;
	coeff[0][2] = coeff[0][9] = coeff[0][5] * e_2;
	coeff[1][1] = coeff[1][3] = coeff[1][8] = coeff[1][10] = -2 * loose * e_2;
	coeff[0][1] = coeff[0][3] = coeff[0][8] = coeff[0][10] = coeff[1][1] * a0;
	
	e_2 *= 2;		/* We will need these in boundary conditions  */
	e_m2 *= 2;
	
	ij_sw_corner = 2 * my + 2;			/*  Corners of array of actual data  */
	ij_se_corner = ij_sw_corner + (nx - 1) * my;
	ij_nw_corner = ij_sw_corner + (ny - 1);
	ij_ne_corner = ij_se_corner + (ny - 1);

}

void	set_offset() {

	int	add_w[5], add_e[5], add_s[5], add_n[5], add_w2[5], add_e2[5], add_s2[5], add_n2[5];
	int	i, j, kase;
	
	add_w[0] = -my; add_w[1] = add_w[2] = add_w[3] = add_w[4] = -grid_east;
	add_w2[0] = -2 * my;  add_w2[1] = -my - grid_east;  add_w2[2] = add_w2[3] = add_w2[4] = -2 * grid_east;
	add_e[4] = my; add_e[0] = add_e[1] = add_e[2] = add_e[3] = grid_east;
	add_e2[4] = 2 * my;  add_e2[3] = my + grid_east;  add_e2[2] = add_e2[1] = add_e2[0] = 2 * grid_east;

	add_n[4] = 1; add_n[3] = add_n[2] = add_n[1] = add_n[0] = grid;
	add_n2[4] = 2;  add_n2[3] = grid + 1;  add_n2[2] = add_n2[1] = add_n2[0] = 2 * grid;
	add_s[0] = -1; add_s[1] = add_s[2] = add_s[3] = add_s[4] = -grid;
	add_s2[0] = -2;  add_s2[1] = -grid - 1;  add_s2[2] = add_s2[3] = add_s2[4] = -2 * grid;

	for (i = 0, kase = 0; i < 5; i++) {
		for (j = 0; j < 5; j++, kase++) {
			offset[kase][0] = add_n2[j];
			offset[kase][1] = add_n[j] + add_w[i];
			offset[kase][2] = add_n[j];
			offset[kase][3] = add_n[j] + add_e[i];
			offset[kase][4] = add_w2[i];
			offset[kase][5] = add_w[i];
			offset[kase][6] = add_e[i];
			offset[kase][7] = add_e2[i];
			offset[kase][8] = add_s[j] + add_w[i];
			offset[kase][9] = add_s[j];
			offset[kase][10] = add_s[j] + add_e[i];
			offset[kase][11] = add_s2[j];
		}
	}
}

void fill_in_forecast () {

	/* Fills in bilinear estimates into new node locations
	   after grid is divided.   
	 */

	int i, j, ii, jj, index_0, index_1, index_2, index_3;
	int index_new;
	double delta_x, delta_y, a0, a1, a2, a3;
	double old_size;
	
		
	old_size = 1.0 / (double)old_grid;

	/* first do from southwest corner */
	
	for (i = 0; i < nx-1; i += old_grid) {
		
		for (j = 0; j < ny-1; j += old_grid) {
			
			/* get indices of bilinear square */
			index_0 = ij_sw_corner + i * my + j;
			index_1 = index_0 + old_grid * my;
			index_2 = index_1 + old_grid;
			index_3 = index_0 + old_grid;
			
			/* get coefficients */
			a0 = u[index_0];
			a1 = u[index_1] - a0;
			a2 = u[index_3] - a0;
			a3 = u[index_2] - a0 - a1 - a2;
			
			/* find all possible new fill ins */
			
			for (ii = i;  ii < i + old_grid; ii += grid) {
				delta_x = (ii - i) * old_size;
				for (jj = j;  jj < j + old_grid; jj += grid) {
					index_new = ij_sw_corner + ii * my + jj;
					if (index_new == index_0) continue;
					delta_y = (jj - j) * old_size;
					u[index_new] = (float)(a0 + a1 * delta_x + delta_y * ( a2 + a3 * delta_x));
					iu[index_new] = 0;
				}
			}
			iu[index_0] = 5;
		}
	}
	
	/* now do linear guess along east edge */
	
	for (j = 0; j < (ny-1); j += old_grid) {
		index_0 = ij_se_corner + j;
		index_3 = index_0 + old_grid;
		for (jj = j;  jj < j + old_grid; jj += grid) {
			index_new = ij_se_corner + jj;
			delta_y = (jj - j) * old_size;
			u[index_new] = u[index_0] + (float)(delta_y * (u[index_3] - u[index_0]));
			iu[index_new] = 0;
		}
		iu[index_0] = 5;
	}
	/* now do linear guess along north edge */
	for (i = 0; i < (nx-1); i += old_grid) {
		index_0 = ij_nw_corner + i * my;
		index_1 = index_0 + old_grid * my;
		for (ii = i;  ii < i + old_grid; ii += grid) {
			index_new = ij_nw_corner + ii * my;
			delta_x = (ii - i) * old_size;
			u[index_new] = u[index_0] + (float)(delta_x * (u[index_1] - u[index_0]));
			iu[index_new] = 0;
		}
		iu[index_0] = 5;
	}
	/* now set northeast corner to fixed and we're done */
	iu[ij_ne_corner] = 5;
}

int compare_points (const void *point_1v, const void *point_2v) {
	/*  Routine for qsort to sort data structure for fast access to data by node location.
	    Sorts on index first, then on radius to node corresponding to index, so that index
	    goes from low to high, and so does radius.
	*/
	int block_i, block_j, index_1, index_2;
	double x0, y0, dist_1, dist_2;
	struct SURFACE_DATA *point_1, *point_2;
	
	point_1 = (struct SURFACE_DATA *)point_1v;
	point_2 = (struct SURFACE_DATA *)point_2v;
	index_1 = point_1->index;
	index_2 = point_2->index;
	if (index_1 < index_2)
		return (-1);
	else if (index_1 > index_2)
		return (1);
	else if (index_1 == OUTSIDE)
		return (0);
	else {	/* Points are in same grid cell, find the one who is nearest to grid point */
		block_i = point_1->index/block_ny;
		block_j = point_1->index%block_ny;
		x0 = xmin + block_i * grid_xinc;
		y0 = ymin + block_j * grid_yinc;
		dist_1 = (point_1->x - x0) * (point_1->x - x0) + (point_1->y - y0) * (point_1->y - y0);
		dist_2 = (point_2->x - x0) * (point_2->x - x0) + (point_2->y - y0) * (point_2->y - y0);
		if (dist_1 < dist_2)
			return (-1);
		if (dist_1 > dist_2)
			return (1);
		else
			return (0);
	}
}

void smart_divide () {
		/* Divide grid by its largest prime factor */
	grid /= factors[n_fact - 1];
	n_fact--;
}

void set_index () {
		/* recomputes data[k].index for new value of grid,
		   sorts data on index and radii, and throws away
		   data which are now outside the useable limits. */
	int i, j, k, k_skipped = 0;

	for (k = 0; k < npoints; k++) {
		i = (int)floor(((data[k].x-xmin)*r_grid_xinc) + 0.5);
		j = (int)floor(((data[k].y-ymin)*r_grid_yinc) + 0.5);
		if (i < 0 || i >= block_nx || j < 0 || j >= block_ny) {
			data[k].index = OUTSIDE;
			k_skipped++;
		}
		else
			data[k].index = i * block_ny + j;
	}
	
	qsort ((char *)data, npoints, sizeof (struct SURFACE_DATA), compare_points);
	
	npoints -= k_skipped;
	
}

void find_nearest_point() {
	int i, j, k, last_index, block_i, block_j, iu_index, briggs_index;
	double x0, y0, dx, dy, xys, xy1, btemp;
	double b0, b1, b2, b3, b4, b5;
	
	last_index = -1;
	small = 0.05 * ((grid_xinc < grid_yinc) ? grid_xinc : grid_yinc);

	for (i = 0; i < nx; i += grid)	/* Reset grid info */
		for (j = 0; j < ny; j += grid)
			iu[ij_sw_corner + i*my + j] = 0;
	
	briggs_index = 0;
	for (k = 0; k < npoints; k++) {	/* Find constraining value  */
		if (data[k].index != last_index) {
			block_i = data[k].index/block_ny;
			block_j = data[k].index%block_ny;
			last_index = data[k].index;
	 		iu_index = ij_sw_corner + (block_i * my + block_j) * grid;
	 		x0 = xmin + block_i*grid_xinc;
	 		y0 = ymin + block_j*grid_yinc;
	 		dx = (data[k].x - x0)*r_grid_xinc;
	 		dy = (data[k].y - y0)*r_grid_yinc;
	 		if (fabs(dx) < small && fabs(dy) < small) {
	 			iu[iu_index] = 5;
	 			u[iu_index] = data[k].z;
	 		}
	 		else {
	 			if (dx >= 0.0) {
	 				if (dy >= 0.0)
	 					iu[iu_index] = 1;
	 				else
	 					iu[iu_index] = 4;
	 			}
	 			else {
	 				if (dy >= 0.0)
	 					iu[iu_index] = 2;
	 				else
	 					iu[iu_index] = 3;
	 			}
	 			dx = fabs(dx);
	 			dy = fabs(dy);
	 			btemp = 2 * one_plus_e2 / ( (dx + dy) * (1.0 + dx + dy) );
	 			b0 = 1.0 - 0.5 * (dx + (dx * dx)) * btemp;
	 			b3 = 0.5 * (e_2 - (dy + (dy * dy)) * btemp);
	 			xys = 1.0 + dx + dy;
	 			xy1 = 1.0 / xys;
	 			b1 = (e_2 * xys - 4 * dy) * xy1;
	 			b2 = 2 * (dy - dx + 1.0) * xy1;
	 			b4 = b0 + b1 + b2 + b3 + btemp;
	 			b5 = btemp * data[k].z;
	 			briggs[briggs_index].b[0] = b0;
	 			briggs[briggs_index].b[1] = b1;
	 			briggs[briggs_index].b[2] = b2;
	 			briggs[briggs_index].b[3] = b3;
	 			briggs[briggs_index].b[4] = b4;
	 			briggs[briggs_index].b[5] = b5;
	 			briggs_index++;
	 		}
	 	}
	 }
}

						
void set_grid_parameters() {			
	block_ny = (ny - 1) / grid + 1;
	block_nx = (nx - 1) / grid + 1;
	grid_xinc = grid * xinc;
	grid_yinc = grid * yinc;
	grid_east = grid * my;
	r_grid_xinc = 1.0 / grid_xinc;
	r_grid_yinc = 1.0 / grid_yinc;
}

void initialize_grid() {
	/*
	 * For the initial gridsize, compute weighted averages of data inside the search radius
	 * and assign the values to u[i,j] where i,j are multiples of gridsize.
	 */
	 int	irad, jrad, i, j, imin, imax, jmin, jmax, index_1, index_2, k, ki, kj, k_index;
	 double	r, rfact, sum_w, sum_zw, weight, x0, y0;

	 irad = (int)ceil(radius/grid_xinc);
	 jrad = (int)ceil(radius/grid_yinc);
	 rfact = -4.5/(radius*radius);
	 
	 for (i = 0; i < block_nx; i ++ ) {
	 	x0 = xmin + i*grid_xinc;
	 	for (j = 0; j < block_ny; j ++ ) {
	 		y0 = ymin + j*grid_yinc;
	 		imin = i - irad;
	 		if (imin < 0) imin = 0;
	 		imax = i + irad;
	 		if (imax >= block_nx) imax = block_nx - 1;
	 		jmin = j - jrad;
	 		if (jmin < 0) jmin = 0;
	 		jmax = j + jrad;
	 		if (jmax >= block_ny) jmax = block_ny - 1;
	 		index_1 = imin*block_ny + jmin;
	 		index_2 = imax*block_ny + jmax + 1;
	 		sum_w = sum_zw = 0.0;
	 		k = 0;
	 		while (k < npoints && data[k].index < index_1) k++;
	 		for (ki = imin; k < npoints && ki <= imax && data[k].index < index_2; ki++) {
	 			for (kj = jmin; k < npoints && kj <= jmax && data[k].index < index_2; kj++) {
	 				k_index = ki*block_ny + kj;
	 				while (k < npoints && data[k].index < k_index) k++;
	 				while (k < npoints && data[k].index == k_index) {
	 					r = (data[k].x-x0)*(data[k].x-x0) + (data[k].y-y0)*(data[k].y-y0);
	 					weight = exp (rfact*r);
	 					sum_w += weight;
	 					sum_zw += weight*data[k].z;
	 					k++;
	 				}
	 			}
	 		}
	 		if (sum_w == 0.0) {
	 			mexPrintf ("surface: Warning: no data inside search radius at: %.8g %.8g\n", x0, y0);
	 			u[ij_sw_corner + (i * my + j) * grid] = (float)z_mean;
	 		}
	 		else {
	 			u[ij_sw_corner + (i*my+j)*grid] = (float)(sum_zw/sum_w);
	 		}
		}
	}
}


void new_initialize_grid(void) {
	/*
	 * For the initial gridsize, load constrained nodes with weighted avg of their data;
	 * and then do something with the unconstrained ones.
	 */
	 int	k, k_index, u_index, block_i, block_j;
	 double	sum_w, sum_zw, weight, x0, y0, dx, dy, dx_scale, dy_scale;

	dx_scale = 4.0 / grid_xinc;
	dy_scale = 4.0 / grid_yinc;
	n_empty = block_ny * block_nx;
	k = 0;
	while (k < npoints) {
		block_i = data[k].index / block_ny;
		block_j = data[k].index % block_ny;
		x0 = xmin + block_i*grid_xinc;
		y0 = ymin + block_j*grid_yinc;
		u_index = ij_sw_corner + (block_i*my + block_j) * grid;
		k_index = data[k].index;
		
		dy = (data[k].y - y0) * dy_scale;
		dx = (data[k].x - x0) * dx_scale;
		sum_w = 1.0 / (1.0 + dx*dx + dy*dy);
		sum_zw = data[k].z * sum_w;
		k++;

		while (k < npoints && data[k].index == k_index) {
			
			dy = (data[k].y - y0) * dy_scale;
			dx = (data[k].x - x0) * dx_scale;
			weight = 1.0 / (1.0 + dx*dx + dy*dy);
			sum_zw += data[k].z * weight;
			sum_w += weight;
			sum_zw += weight*data[k].z;
			k++;
	 	}
	 	u[u_index] = (float)(sum_zw/sum_w);
	 	iu[u_index] = 5;
	 	n_empty--;
	 }
}

/* This function rewritten by D.W. Caress 5/3/94 */
int read_data(int ndat, float *xdat, float *ydat, float *zdat) {

	int	i, j, k, kmax, kmin, idat;
	double	zmin = 1.0e38, zmax = -1.0e38;

	data = (struct SURFACE_DATA *) mxCalloc ((size_t)ndat, sizeof(struct SURFACE_DATA));
	
	/* Read in xyz data and computes index no and store it in a structure */
	k = 0;
	z_mean = 0;
	for (idat = 0; idat < ndat; idat++) {
		i = (int)floor(((xdat[idat]-xmin)*r_grid_xinc) + 0.5);
		j = (int)floor(((ydat[idat]-ymin)*r_grid_yinc) + 0.5);
		if (i >= 0 && i < block_nx && j >= 0 && j < block_ny) {
			data[k].index = i * block_ny + j;
			data[k].x = xdat[idat];
			data[k].y = ydat[idat];
			data[k].z = zdat[idat];
			if (zmin > zdat[idat]) {
				zmin = zdat[idat];
				kmin = k;
			}
			if (zmax < zdat[idat]) {
				zmax = zdat[idat];
				kmax = k;
			}
			k++;
			z_mean += zdat[idat];
		}
	}

	npoints = k;
	z_mean /= k;
	if( converge_limit == 0.0 ) {
		converge_limit = 0.001 * z_scale; /* c_l = 1 ppt of L2 scale */
	}
	if (verbose) {
		mexPrintf("surface: Minimum value of your dataset x,y,z at: %g %g %g\n",
			data[kmin].x, data[kmin].y, data[kmin].z);
		mexPrintf("surface: Maximum value of your dataset x,y,z at: %g %g %g\n",
			data[kmax].x, data[kmax].y, data[kmax].z);
	}
	
	if (set_low == 1)
		low_limit = data[kmin].z;
	else if (set_low == 2 && low_limit > data[kmin].z) {
	/*	low_limit = data[kmin].z;	*/
		/*
		fprintf (stderr, "surface: Warning:  Your lower value is > than min data value.\n");
		*/
	}
	if (set_high == 1)
		high_limit = data[kmax].z;
	else if (set_high == 2 && high_limit < data[kmax].z) {
	/*	high_limit = data[kmax].z;	*/
		/*
		fprintf (stderr, "surface: Warning:  Your upper value is < than max data value.\n");
		*/
	}
	return(0);
}

/* this function rewritten from write_output() by D.W. Caress 5/3/94 */
void get_output(float *sgrid) {
	
	int	index, i, j;

	index = ij_sw_corner;
	for(i = 0; i < nx; i++, index += my) 
		for (j = 0; j < ny; j++) 
			sgrid[j*nx+i] = u[index + ny - j - 1];
}
	
int	iterate(int mode) {

	int	i, j, k, ij, kase, briggs_index, ij_v2;
	int	x_case, y_case, x_w_case, x_e_case, y_s_case, y_n_case;
	int	iteration_count = 0;
	
	double	current_limit = converge_limit / grid;
	double	change, max_change = 0.0, busum, sum_ij;
	double	b0, b1, b2, b3, b4, b5;
	
	double	x_0_const = 4.0 * (1.0 - boundary_tension) / (2.0 - boundary_tension);
	double	x_1_const = (3 * boundary_tension - 2.0) / (2.0 - boundary_tension);
	double	y_denom = 2 * epsilon * (1.0 - boundary_tension) + boundary_tension;
	double	y_0_const = 4 * epsilon * (1.0 - boundary_tension) / y_denom;
	double	y_1_const = (boundary_tension - 2 * epsilon * (1.0 - boundary_tension) ) / y_denom;

	do {
		briggs_index = 0;	/* Reset the constraint table stack pointer  */
		
		max_change = -1.0;
		
		/* Fill in auxiliary boundary values (in new way) */
		
		/* First set d2[]/dn2 = 0 along edges:  */
		/* New experiment : (1-T)d2[]/dn2 + Td[]/dn = 0  */
		
		
		
		for (i = 0; i < nx; i += grid) {
			/* set d2[]/dy2 = 0 on south side:  */
			ij = ij_sw_corner + i * my;
			/* u[ij - 1] = 2 * u[ij] - u[ij + grid];  */
			u[ij - 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij + grid]);
			/* set d2[]/dy2 = 0 on north side:  */
			ij = ij_nw_corner + i * my;
			/* u[ij + 1] = 2 * u[ij] - u[ij - grid];  */
			u[ij + 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij - grid]);
			
		}
		
		for (j = 0; j < ny; j += grid) {
			/* set d2[]/dx2 = 0 on west side:  */
			ij = ij_sw_corner + j;
			/* u[ij - my] = 2 * u[ij] - u[ij + grid_east];  */
			u[ij - my] = (float)(x_1_const * u[ij + grid_east] + x_0_const * u[ij]);
			/* set d2[]/dx2 = 0 on east side:  */
			ij = ij_se_corner + j;
			/* u[ij + my] = 2 * u[ij] - u[ij - grid_east];  */
			u[ij + my] = (float)(x_1_const * u[ij - grid_east] + x_0_const * u[ij]);
		}
			
		/* Now set d2[]/dxdy = 0 at each corner:  */
		
		ij = ij_sw_corner;
		u[ij - my - 1] = u[ij + grid_east - 1] + u[ij - my + grid] - u[ij + grid_east + grid];
				
		ij = ij_nw_corner;
		u[ij - my + 1] = u[ij + grid_east + 1] + u[ij - my - grid] - u[ij + grid_east - grid];
				
		ij = ij_se_corner;
		u[ij + my - 1] = u[ij - grid_east - 1] + u[ij + my + grid] - u[ij - grid_east + grid];
				
		ij = ij_ne_corner;
		u[ij + my + 1] = u[ij - grid_east + 1] + u[ij + my - grid] - u[ij - grid_east - grid];
		
		/* Now set (1-T)dC/dn + Tdu/dn = 0 at each edge :  */
		/* New experiment:  only dC/dn = 0  */
		
		x_w_case = 0;
		x_e_case = block_nx - 1;
		for (i = 0; i < nx; i += grid, x_w_case++, x_e_case--) {
		
			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;
				
			/* South side :  */
			kase = x_case * 5;
			ij = ij_sw_corner + i * my;
			u[ij + offset[kase][11]] = 
				(float)(u[ij + offset[kase][0]] + eps_m2*(u[ij + offset[kase][1]] + u[ij + offset[kase][3]]
					- u[ij + offset[kase][8]] - u[ij + offset[kase][10]])
					+ two_plus_em2 * (u[ij + offset[kase][9]] - u[ij + offset[kase][2]]) );
				/*  + tense * eps_m2 * (u[ij + offset[kase][2]] - u[ij + offset[kase][9]]) / (1.0 - tense);  */
			/* North side :  */
			kase = x_case * 5 + 4;
			ij = ij_nw_corner + i * my;
			u[ij + offset[kase][0]] = 
				-(float)(-u[ij + offset[kase][11]] + eps_m2 * (u[ij + offset[kase][1]] + u[ij + offset[kase][3]]
					- u[ij + offset[kase][8]] - u[ij + offset[kase][10]])
					+ two_plus_em2 * (u[ij + offset[kase][9]] - u[ij + offset[kase][2]]) );
				/*  - tense * eps_m2 * (u[ij + offset[kase][2]] - u[ij + offset[kase][9]]) / (1.0 - tense);  */
		}
		
		y_s_case = 0;
		y_n_case = block_ny - 1;
		for (j = 0; j < ny; j += grid, y_s_case++, y_n_case--) {
				
			if(y_s_case < 2)
				y_case = y_s_case;
			else if(y_n_case < 2)
				y_case = 4 - y_n_case;
			else
				y_case = 2;
			
			/* West side :  */
			kase = y_case;
			ij = ij_sw_corner + j;
			u[ij+offset[kase][4]] = 
				u[ij + offset[kase][7]] + (float)(eps_p2 * (u[ij + offset[kase][3]] + u[ij + offset[kase][10]]
				-u[ij + offset[kase][1]] - u[ij + offset[kase][8]])
				+ two_plus_ep2 * (u[ij + offset[kase][5]] - u[ij + offset[kase][6]]));
				/*  + tense * (u[ij + offset[kase][6]] - u[ij + offset[kase][5]]) / (1.0 - tense);  */
			/* East side :  */
			kase = 20 + y_case;
			ij = ij_se_corner + j;
			u[ij + offset[kase][7]] = 
				- (float)(-u[ij + offset[kase][4]] + eps_p2 * (u[ij + offset[kase][3]] + u[ij + offset[kase][10]]
				- u[ij + offset[kase][1]] - u[ij + offset[kase][8]])
				+ two_plus_ep2 * (u[ij + offset[kase][5]] - u[ij + offset[kase][6]]) );
				/*  - tense * (u[ij + offset[kase][6]] - u[ij + offset[kase][5]]) / (1.0 - tense);  */
		}

			
			
		/* That's it for the boundary points.  Now loop over all data  */
		
		x_w_case = 0;
		x_e_case = block_nx - 1;
		for (i = 0; i < nx; i += grid, x_w_case++, x_e_case--) {
		
			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;
			
			y_s_case = 0;
			y_n_case = block_ny - 1;
			
			ij = ij_sw_corner + i * my;
			
			for (j = 0; j < ny; j += grid, ij += grid, y_s_case++, y_n_case--) {
	
				if (iu[ij] == 5) continue;	/* Point is fixed  */
				
				if(y_s_case < 2)
					y_case = y_s_case;
				else if(y_n_case < 2)
					y_case = 4 - y_n_case;
				else
					y_case = 2;
				
				kase = x_case * 5 + y_case;
				sum_ij = 0.0;

				if (iu[ij] == 0) {		/* Point is unconstrained  */
					for (k = 0; k < 12; k++) {
						sum_ij += (u[ij + offset[kase][k]] * coeff[0][k]);
					}
				}
				else {				/* Point is constrained  */
				
					b0 = briggs[briggs_index].b[0];
					b1 = briggs[briggs_index].b[1];
					b2 = briggs[briggs_index].b[2];
					b3 = briggs[briggs_index].b[3];
					b4 = briggs[briggs_index].b[4];
					b5 = briggs[briggs_index].b[5];
					briggs_index++;
					if (iu[ij] < 3) {
						if (iu[ij] == 1) {	/* Point is in quadrant 1  */
							busum = b0 * u[ij + offset[kase][10]]
								+ b1 * u[ij + offset[kase][9]]
								+ b2 * u[ij + offset[kase][5]]
								+ b3 * u[ij + offset[kase][1]];
						}
						else {			/* Point is in quadrant 2  */
							busum = b0 * u[ij + offset[kase][8]]
								+ b1 * u[ij + offset[kase][9]]
								+ b2 * u[ij + offset[kase][6]]
								+ b3 * u[ij + offset[kase][3]];
						}
					}
					else {
						if (iu[ij] == 3) {	/* Point is in quadrant 3  */
							busum = b0 * u[ij + offset[kase][1]]
								+ b1 * u[ij + offset[kase][2]]
								+ b2 * u[ij + offset[kase][6]]
								+ b3 * u[ij + offset[kase][10]];
						}
						else {		/* Point is in quadrant 4  */
							busum = b0 * u[ij + offset[kase][3]]
								+ b1 * u[ij + offset[kase][2]]
								+ b2 * u[ij + offset[kase][5]]
								+ b3 * u[ij + offset[kase][8]];
						}
					}
					for (k = 0; k < 12; k++) {
						sum_ij += (u[ij + offset[kase][k]] * coeff[1][k]);
					}
					sum_ij = (sum_ij + a0_const_2 * (busum + b5))
						/ (a0_const_1 + a0_const_2 * b4);
				}
				
				/* New relaxation here  */
				sum_ij = u[ij] * relax_old + sum_ij * relax_new;
				
				if (constrained) {	/* Must check limits.  Note lower/upper is v2 format and need ij_v2! */
					ij_v2 = (ny - j - 1) * nx + i;
					if (set_low /*&& !GMT_is_fnan((double)lower[ij_v2])*/ && sum_ij < lower[ij_v2])
						sum_ij = lower[ij_v2];
					else if (set_high /*&& !GMT_is_fnan((double)upper[ij_v2])*/ && sum_ij > upper[ij_v2])
						sum_ij = upper[ij_v2];
				}
					
				change = fabs(sum_ij - u[ij]);
				u[ij] = (float)sum_ij;
				if (change > max_change) max_change = change;
			}
		}
		iteration_count++;
		total_iterations++;
		max_change *= z_scale;	/* Put max_change into z units  */
		if (verbose > 1) 
			mexPrintf("%4d\t%c\t%8d\t%10g\t%10g\t%10d\n", grid, mode_type[mode], 
				iteration_count, max_change, current_limit, total_iterations);

	} while (max_change > current_limit && iteration_count < max_iterations);
	
	if (verbose) mexPrintf("%4d\t%c\t%8d\t%10g\t%10g\t%10d\n",
		grid, mode_type[mode], iteration_count, max_change, current_limit, total_iterations);

	return(iteration_count);
}

void check_errors () {

	int	i, j, k, ij, n_nodes, move_over[12];	/* move_over = offset[kase][12], but grid = 1 so move_over is easy  */
	
	double	x0, y0, dx, dy, mean_error, mean_squared_error, z_est, z_err, curvature, c;
	double	du_dx, du_dy, d2u_dx2, d2u_dxdy, d2u_dy2, d3u_dx3, d3u_dx2dy, d3u_dxdy2, d3u_dy3;
	
	double	x_0_const = 4.0 * (1.0 - boundary_tension) / (2.0 - boundary_tension);
	double	x_1_const = (3 * boundary_tension - 2.0) / (2.0 - boundary_tension);
	double	y_denom = 2 * epsilon * (1.0 - boundary_tension) + boundary_tension;
	double	y_0_const = 4 * epsilon * (1.0 - boundary_tension) / y_denom;
	double	y_1_const = (boundary_tension - 2 * epsilon * (1.0 - boundary_tension) ) / y_denom;
	
	
	move_over[0] = 2;
	move_over[1] = 1 - my;
	move_over[2] = 1;
	move_over[3] = 1 + my;
	move_over[4] = -2 * my;
	move_over[5] = -my;
	move_over[6] = my;
	move_over[7] = 2 * my;
	move_over[8] = -1 - my;
	move_over[9] = -1;
	move_over[10] = -1 + my;
	move_over[11] = -2;

	mean_error = 0;
	mean_squared_error = 0;
	
	/* First update the boundary values  */

	for (i = 0; i < nx; i ++) {
		ij = ij_sw_corner + i * my;
		u[ij - 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij + 1]);
		ij = ij_nw_corner + i * my;
		u[ij + 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij - 1]);
	}

	for (j = 0; j < ny; j ++) {
		ij = ij_sw_corner + j;
		u[ij - my] = (float)(x_1_const * u[ij + my] + x_0_const * u[ij]);
		ij = ij_se_corner + j;
		u[ij + my] = (float)(x_1_const * u[ij - my] + x_0_const * u[ij]);
	}

	ij = ij_sw_corner;
	u[ij - my - 1] = u[ij + my - 1] + u[ij - my + 1] - u[ij + my + 1];
	ij = ij_nw_corner;
	u[ij - my + 1] = u[ij + my + 1] + u[ij - my - 1] - u[ij + my - 1];
	ij = ij_se_corner;
	u[ij + my - 1] = u[ij - my - 1] + u[ij + my + 1] - u[ij - my + 1];
	ij = ij_ne_corner;
	u[ij + my + 1] = u[ij - my + 1] + u[ij + my - 1] - u[ij - my - 1];

	for (i = 0; i < nx; i ++) {
				
		ij = ij_sw_corner + i * my;
		u[ij + move_over[11]] = 
			(float)(u[ij + move_over[0]] + eps_m2*(u[ij + move_over[1]] + u[ij + move_over[3]]
				- u[ij + move_over[8]] - u[ij + move_over[10]])
				+ two_plus_em2 * (u[ij + move_over[9]] - u[ij + move_over[2]]) );
					
		ij = ij_nw_corner + i * my;
		u[ij + move_over[0]] = 
			-(float)(-u[ij + move_over[11]] + eps_m2 * (u[ij + move_over[1]] + u[ij + move_over[3]]
				- u[ij + move_over[8]] - u[ij + move_over[10]])
				+ two_plus_em2 * (u[ij + move_over[9]] - u[ij + move_over[2]]) );
	}
		
	for (j = 0; j < ny; j ++) {
			
		ij = ij_sw_corner + j;
		u[ij+move_over[4]] = 
			u[ij + move_over[7]] + (float)(eps_p2 * (u[ij + move_over[3]] + u[ij + move_over[10]]
			-u[ij + move_over[1]] - u[ij + move_over[8]])
			+ two_plus_ep2 * (u[ij + move_over[5]] - u[ij + move_over[6]]));
				
		ij = ij_se_corner + j;
		u[ij + move_over[7]] = 
			- (float)(-u[ij + move_over[4]] + eps_p2 * (u[ij + move_over[3]] + u[ij + move_over[10]]
			- u[ij + move_over[1]] - u[ij + move_over[8]])
			+ two_plus_ep2 * (u[ij + move_over[5]] - u[ij + move_over[6]]) );
	}

	/* That resets the boundary values.  Now we can test all data.  
		Note that this loop checks all values, even though only nearest were used.  */
	
	for (k = 0; k < npoints; k++) {
		i = data[k].index/ny;
		j = data[k].index%ny;
	 	ij = ij_sw_corner + i * my + j;
	 	if ( iu[ij] == 5 ) continue;
	 	x0 = xmin + i*xinc;
	 	y0 = ymin + j*yinc;
	 	dx = (data[k].x - x0)*r_xinc;
	 	dy = (data[k].y - y0)*r_yinc;
 
	 	du_dx = 0.5 * (u[ij + move_over[6]] - u[ij + move_over[5]]);
	 	du_dy = 0.5 * (u[ij + move_over[2]] - u[ij + move_over[9]]);
	 	d2u_dx2 = u[ij + move_over[6]] + u[ij + move_over[5]] - 2 * u[ij];
	 	d2u_dy2 = u[ij + move_over[2]] + u[ij + move_over[9]] - 2 * u[ij];
	 	d2u_dxdy = 0.25 * (u[ij + move_over[3]] - u[ij + move_over[1]]
	 			- u[ij + move_over[10]] + u[ij + move_over[8]]);
	 	d3u_dx3 = 0.5 * ( u[ij + move_over[7]] - 2 * u[ij + move_over[6]]
	 				+ 2 * u[ij + move_over[5]] - u[ij + move_over[4]]);
	 	d3u_dy3 = 0.5 * ( u[ij + move_over[0]] - 2 * u[ij + move_over[2]]
	 				+ 2 * u[ij + move_over[9]] - u[ij + move_over[11]]);
	 	d3u_dx2dy = 0.5 * ( ( u[ij + move_over[3]] + u[ij + move_over[1]] - 2 * u[ij + move_over[2]] )
	 				- ( u[ij + move_over[10]] + u[ij + move_over[8]] - 2 * u[ij + move_over[9]] ) );
	 	d3u_dxdy2 = 0.5 * ( ( u[ij + move_over[3]] + u[ij + move_over[10]] - 2 * u[ij + move_over[6]] )
	 				- ( u[ij + move_over[1]] + u[ij + move_over[8]] - 2 * u[ij + move_over[5]] ) );

	 	/* 3rd order Taylor approx:  */
	 		
	 	z_est = u[ij] + dx * (du_dx +  dx * ( (0.5 * d2u_dx2) + dx * (d3u_dx3 / 6.0) ) )
				+ dy * (du_dy +  dy * ( (0.5 * d2u_dy2) + dy * (d3u_dy3 / 6.0) ) )
	 			+ dx * dy * (d2u_dxdy) + (0.5 * dx * d3u_dx2dy) + (0.5 * dy * d3u_dxdy2);
	 		
	 	z_err = z_est - data[k].z;
	 	mean_error += z_err;
	 	mean_squared_error += (z_err * z_err);
	 }
	 mean_error /= npoints;
	 mean_squared_error = sqrt( mean_squared_error / npoints);
	 
	 curvature = 0.0;
	 n_nodes = nx * ny;
	 
	 for (i = 0; i < nx; i++) {
	 	for (j = 0; j < ny; j++) {
	 		ij = ij_sw_corner + i * my + j;
	 		c = u[ij + move_over[6]] + u[ij + move_over[5]]
	 			+ u[ij + move_over[2]] + u[ij + move_over[9]] - 4.0 * u[ij + move_over[6]];
			curvature += (c * c);
		}
	}

	if (verbose) {
		mexPrintf("\nSpline interpolation fit information:\n");
		mexPrintf("Data points   nodes    mean error     rms error     curvature\n");
		mexPrintf("%9d %9d   %10g   %10g  %10g\n",
			npoints, n_nodes, mean_error, mean_squared_error, curvature);
	}
 }

void	remove_planar_trend() {

	int	i;
	double	a, b, c, d, xx, yy, zz;
	double	sx, sy, sz, sxx, sxy, sxz, syy, syz;
	
	sx = sy = sz = sxx = sxy = sxz = syy = syz = 0.0;
	
	for (i = 0; i < npoints; i++) {

		xx = (data[i].x - xmin) * r_xinc;
		yy = (data[i].y - ymin) * r_yinc;
		zz = data[i].z;
		
		sx += xx;
		sy += yy;
		sz += zz;
		sxx +=(xx * xx);
		sxy +=(xx * yy);
		sxz +=(xx * zz);
		syy +=(yy * yy);
		syz +=(yy * zz);
	}
	
	d = npoints*sxx*syy + 2*sx*sy*sxy - npoints*sxy*sxy - sx*sx*syy - sy*sy*sxx;
	
	if (d == 0.0) {
		plane_c0 = plane_c1 = plane_c2 = 0.0;
		return;
	}
	
	a = sz*sxx*syy + sx*sxy*syz + sy*sxy*sxz - sz*sxy*sxy - sx*sxz*syy - sy*syz*sxx;
	b = npoints*sxz*syy + sz*sy*sxy + sy*sx*syz - npoints*sxy*syz - sz*sx*syy - sy*sy*sxz;
	c = npoints*sxx*syz + sx*sy*sxz + sz*sx*sxy - npoints*sxy*sxz - sx*sx*syz - sz*sy*sxx;

	plane_c0 = a / d;
	plane_c1 = b / d;
	plane_c2 = c / d;

	for (i = 0; i < npoints; i++) {

		xx = (data[i].x - xmin) * r_xinc;
		yy = (data[i].y - ymin) * r_yinc;
		
		data[i].z -= (float)(plane_c0 + plane_c1 * xx + plane_c2 * yy);
	}

}

void	replace_planar_trend() {
	int	i, j, ij;

	 for (i = 0; i < nx; i++) {
	 	for (j = 0; j < ny; j++) {
	 		ij = ij_sw_corner + i * my + j;
	 		u[ij] = (float)((u[ij] * z_scale) + (plane_c0 + plane_c1 * i + plane_c2 * j));
		}
	}
}

void	throw_away_unusables() {
	/* This is a new routine to eliminate data which will become
		unusable on the final iteration, when grid = 1.
		It assumes grid = 1 and set_grid_parameters has been
		called.  We sort, mark redundant data as OUTSIDE, and
		sort again, chopping off the excess.
		
		Experimental modification 5 Dec 1988 by Smith, as part
		of a new implementation using core memory for b[6]
		coefficients, eliminating calls to temp file.
	*/
	
	int	last_index, n_outside, k;
	
	/* Sort the data  */
	
	qsort ((char *)data, npoints, sizeof (struct SURFACE_DATA), compare_points);
	
	/* If more than one datum is indexed to same node, only the first should be kept.
		Mark the additional ones as OUTSIDE
	*/
	last_index = -1;
	n_outside = 0;
	for (k = 0; k < npoints; k++) {
		if (data[k].index == last_index) {
			data[k].index = OUTSIDE;
			n_outside++;
		}
		else {
			last_index = data[k].index;
		}
	}
	/* Sort again; this time the OUTSIDE points will be thrown away  */
	
	qsort ((char *)data, npoints, sizeof (struct SURFACE_DATA), compare_points);
	npoints -= n_outside;
	data = (struct SURFACE_DATA *) mxCalloc ((size_t)npoints, sizeof(struct SURFACE_DATA));
	if (verbose && (n_outside)) {
		mexPrintf("surface: %d unusable points were supplied; these will be ignored.\n", n_outside);
		mexPrintf("\tYou should have pre-processed the data with blockmean or blockmedian.\n");
	}

}

void	rescale_z_values() {
	int	i;
	double	ssz = 0.0;

	for (i = 0; i < npoints; i++) {
		ssz += (data[i].z * data[i].z);
	}
	
	/* Set z_scale = rms(z):  */
	
	z_scale = sqrt(ssz / npoints);
	r_z_scale = 1.0 / z_scale;

	for (i = 0; i < npoints; i++) {
		data[i].z *= (float)r_z_scale;
	}
}

void load_constraints (char *low, char *high) {
	int i, j, ij;
	/* int n_trimmed;*/
	double yy;
/*	struct GRD_HEADER hdr;*/
	
	/* Load lower/upper limits, verify range, deplane, and rescale */
	
	if (set_low > 0) {
		lower = (float *) mxCalloc ((size_t)(nx * ny), sizeof (float));
		if (set_low < 3)
			for (i = 0; i < nx * ny; i++) lower[i] = (float)low_limit;

			
		for (j = ij = 0; j < ny; j++) {
			yy = ny - j - 1;
			for (i = 0; i < nx; i++, ij++) {
				/*if (GMT_is_fnan ((double)lower[ij])) continue;*/
				lower[ij] -= (float)(plane_c0 + plane_c1 * i + plane_c2 * yy);
				lower[ij] *= (float)r_z_scale;
			}
		}
		constrained = TRUE;
	}
	if (set_high > 0) {
		upper = (float *) mxCalloc ((size_t)(nx * ny), sizeof (float));
		if (set_high < 3)
			for (i = 0; i < nx * ny; i++) upper[i] = (float)high_limit;

		for (j = ij = 0; j < ny; j++) {
			yy = ny - j - 1;
			for (i = 0; i < nx; i++, ij++) {
				/*if (GMT_is_fnan ((double)upper[ij])) continue;*/
				upper[ij] -= (float)(plane_c0 + plane_c1 * i + plane_c2 * yy);
				upper[ij] *= (float)r_z_scale;
			}
		}
		constrained = TRUE;
	}
}

double	guess_surface_time(int nx, int ny) {
	/* Routine to guess a number proportional to the operations
	 * required by surface working on a user-desired grid of
	 * size nx by ny, where nx = (xmax - xmin)/dx, and same for
	 * ny.  (That is, one less than actually used in routine.)
	 *
	 * This is based on the following untested conjecture:
	 * 	The operations are proportional to T = nxg*nyg*L,
	 *	where L is a measure of the distance that data
	 *	constraints must propagate, and nxg, nyg are the
	 * 	current size of the grid.
	 *	For nx,ny relatively prime, we will go through only
	 * 	one grid cycle, L = max(nx,ny), and T = nx*ny*L.
	 *	But for nx,ny whose greatest common divisor is a highly
	 * 	composite number, we will have L equal to the division
	 * 	step made at each new grid cycle, and nxg,nyg will
	 * 	also be smaller than nx,ny.  Thus we can hope to find
	 *	some nx,ny for which the total value of T is small.
	 *
	 * The above is pure speculation and has not been derived
	 * empirically.  In actual practice, the distribution of the
	 * data, both spatially and in terms of their values, will
	 * have a strong effect on convergence.
	 *
	 * W. H. F. Smith, 26 Feb 1992.  */

	int	gcd;		/* Current value of the gcd  */
	int	nxg, nyg;	/* Current value of the grid dimensions  */
	int	nfactors;	/* Number of prime factors of current gcd  */
	int	factor;		/* Currently used factor  */
	/* Doubles are used below, even though the values will be integers,
		because the multiplications might reach sizes of O(n**3)  */
	double	t_sum;		/* Sum of values of T at each grid cycle  */
	double	length;		/* Current propagation distance.  */


	gcd = gcd_euclid(nx, ny);
	if (gcd > 1) {
		nfactors = get_prime_factors(gcd, factors);
		nxg = nx/gcd;
		nyg = ny/gcd;
		if (nxg < 3 || nyg < 3) {
			factor = factors[nfactors - 1];
			nfactors--;
			gcd /= factor;
			nxg *= factor;
			nyg *= factor;
		}
	}
	else {
		nxg = nx;
		nyg = ny;
	}
	length = MAX(nxg, nyg);
	t_sum = nxg * (nyg * length);	/* Make it double at each multiply  */

	/* Are there more grid cycles ?  */
	while (gcd > 1) {
		factor = factors[nfactors - 1];
		nfactors--;
		gcd /= factor;
		nxg *= factor;
		nyg *= factor;
		length = factor;
		t_sum += nxg * (nyg * length);
	}
	return(t_sum);
}


int	get_prime_factors(int n, int f[]) {
	/* Fills the integer array f with the prime factors of n.
	 * Returns the number of locations filled in f, which is
	 * one if n is prime.
	 *
	 * f[] should have been malloc'ed to enough space before
	 * calling prime_factors().  We can be certain that f[32]
	 * is enough space, for if n fits in a long, then n < 2**32,
	 * and so it must have fewer than 32 prime factors.  I think
	 * that in general, ceil(log2((double)n)) is enough storage
	 * space for f[].
	 *
	 * Tries 2,3,5 explicitly; then alternately adds 2 or 4
	 * to the previously tried factor to obtain the next trial
	 * factor.  This is done with the variable two_four_toggle.
	 * With this method we try 7,11,13,17,19,23,25,29,31,35,...
	 * up to a maximum of sqrt(n).  This shortened list results
	 * in 1/3 fewer divisions than if we simply tried all integers
	 * between 5 and sqrt(n).  We can reduce the size of the list
	 * of trials by an additional 20% by removing the multiples
	 * of 5, which are equal to 30m +/- 5, where m >= 1.  Starting
	 * from 25, these are found by alternately adding 10 or 20.
	 * To do this, we use the variable ten_twenty_toggle.
	 *
	 * W. H. F. Smith, 26 Feb 1992, after D.E. Knuth, vol. II  */

	int	current_factor;	/* The factor currently being tried  */
	int	max_factor;	/* Don't try any factors bigger than this  */
	int	n_factors = 0;	/* Returned; one if n is prime  */
	int	two_four_toggle = 0;	/* Used to add 2 or 4 to get next trial factor  */
	int	ten_twenty_toggle = 0;	/* Used to add 10 or 20 to skip_five  */
	int	skip_five = 25;	/* Used to skip multiples of 5 in the list  */
	int	m;	/* Used to keep a working copy of n  */


	/* Initialize m and max_factor  */
	m = abs(n);
	if (m < 2) return(0);
	max_factor = (int)floor(sqrt((double)m));

	/* First find the 2s  */
	current_factor = 2;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 3s  */
	current_factor = 3;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 5s  */
	current_factor = 5;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Now try all the rest  */

	while (m > 1 && current_factor <= max_factor) {

		/* Current factor is either 2 or 4 more than previous value  */

		if (two_four_toggle) {
			current_factor += 4;
			two_four_toggle = 0;
		}
		else {
			current_factor += 2;
			two_four_toggle = 1;
		}

		/* If current factor is a multiple of 5, skip it.  But first,
			set next value of skip_five according to 10/20 toggle:  */

		if (current_factor == skip_five) {
			if (ten_twenty_toggle) {
				skip_five += 20;
				ten_twenty_toggle = 0;
			}
			else {
				skip_five += 10;
				ten_twenty_toggle = 1;
			}
			continue;
		}

		/* Get here when current_factor is not a multiple of 2,3 or 5:  */

		while(!(m%current_factor)) {
			m /= current_factor;
			f[n_factors] = current_factor;
			n_factors++;
		}
	}

	/* Get here when all factors up to floor(sqrt(n)) have been tried.  */

	if (m > 1) {
		/* m is an additional prime factor of n  */
		f[n_factors] = m;
		n_factors++;
	}
	return (n_factors);
}
/* gcd_euclid.c  Greatest common divisor routine  */

#define IABS(i)	(((i) < 0) ? -(i) : (i))

int	gcd_euclid(int a,int b) {
	/* Returns the greatest common divisor of u and v by Euclid's method.
	 * I have experimented also with Stein's method, which involves only
	 * subtraction and left/right shifting; Euclid is faster, both for
	 * integers of size 0 - 1024 and also for random integers of a size
	 * which fits in a long integer.  Stein's algorithm might be better
	 * when the integers are HUGE, but for our purposes, Euclid is fine.
	 *
	 * Walter H. F. Smith, 25 Feb 1992, after D. E. Knuth, vol. II  */

	int	u,v,r;

	u = MAX(IABS(a), IABS(b));
	v = MIN(IABS(a), IABS(b));

	while (v > 0) {
		r = u%v;	/* Knuth notes that u < 2v 40% of the time;  */
		u = v;		/* thus we could have tried a subtraction  */
		v = r;		/* followed by an if test to do r = u%v  */
	}
	return(u);
}



#define GMT_INC_IS_M		1
#define GMT_INC_IS_KM		2
#define GMT_INC_IS_MILES	4
#define GMT_INC_IS_NMILES	8
#define GMT_INC_IS_NNODES	16
#define GMT_INC_IS_EXACT	32
#define GMT_INC_UNITS		15

#define GMT_DEG2MIN_F	60.0
#define GMT_MIN2DEG	(1.0 / GMT_DEG2MIN_F)
#define GMT_DEG2SEC_F	3600.0
#define GMT_SEC2DEG	(1.0 / GMT_DEG2SEC_F)

int GMT_getinc (char *line, double *dx, double *dy)
{	/* Special case of getincn use where n is two. */

	int n;
	double inc[2];

	/* Syntax: -I<xinc>[m|c|e|i|k|n|+|=][/<yinc>][m|c|e|i|k|n|+|=]
	 * Units: m = minutes
	 *	  c = seconds
	 *	  e = meter [Convert to degrees]
	 *	  i = miles [Convert to degrees]
	 *	  k = km [Convert to degrees]
	 *	  n = nautical miles [Convert to degrees]
	 * Flags: = = Adjust -R to fit exact -I [Default modifies -I to fit -R]
	 *	  + = incs are actually nx/ny - convert to get xinc/yinc
	 */

	n = GMT_getincn (line, inc, 2);
	*dx = inc[0] ; *dy = inc[1];
	if (n == 1) {	/* Must copy y info from x */
		*dy = *dx;
		GMT_inc_code[1] = GMT_inc_code[0];	/* Use exact inc codes for both x and y */
	}

	if (GMT_inc_code[0] & GMT_INC_IS_NNODES && GMT_inc_code[0] & GMT_INC_UNITS) {
		mexPrintf ("GMTMBGRID: ERROR: number of x nodes cannot have units\n");
		return (EXIT_FAILURE);
	}
	if (GMT_inc_code[1] & GMT_INC_IS_NNODES && GMT_inc_code[1] & GMT_INC_UNITS) {
		mexPrintf ("GMTMBGRID: ERROR: number of y nodes cannot have units\n");
		return (EXIT_FAILURE);
	}
	return (0);
}


int GMT_getincn (char *line, double inc[], int n) {
	int last, i, pos;
	char p[BUFSIZ];
	double scale = 1.0;

	/* Dechipers dx/dy/dz/dw/du/dv/... increment strings with n items */

	memset ((void *)inc, 0, (size_t)(n * sizeof (double)));

	i = pos = GMT_inc_code[0] = GMT_inc_code[1] = 0;

	while (i < n && (GMT_strtok (line, "/", &pos, p))) {
		last = strlen (p) - 1;
		if (p[last] == '=') {	/* Let -I override -R */
			p[last] = 0;
			if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_EXACT;
			last--;
		}
		else if (p[last] == '+' || p[last] == '!') {	/* Number of nodes given, determine inc from domain (! added since documentation mentioned this once... */
			p[last] = 0;
			if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_NNODES;
			last--;
		}
		switch (p[last]) {
			case 'm':
			case 'M':	/* Gave arc minutes */
				p[last] = 0;
				scale = GMT_MIN2DEG;
				break;
			case 'c':
			case 'C':	/* Gave arc seconds */
				p[last] = 0;
				scale = GMT_SEC2DEG;
				break;
			case 'e':
			case 'E':	/* Gave meters along mid latitude */
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_M;
				break;
			case 'K':	/* Gave km along mid latitude */
			case 'k':
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_KM;
				break;
			case 'I':	/* Gave miles along mid latitude */
			case 'i':
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_MILES;
				break;
			case 'N':	/* Gave nautical miles along mid latitude */
			case 'n':
				p[last] = 0;
				if (i < 2) GMT_inc_code[i] |= GMT_INC_IS_NMILES;
				break;
			default:	/* No special flags or units */
				scale = 1.0;
				break;
		}
		if ( (sscanf(p, "%lf", &inc[i])) != 1) {
			mexPrintf ("GMTMBGRID ERROR: Unable to decode %s as a floating point number\n", p);
			return (EXIT_FAILURE);
		}
		inc[i] *= scale;
		i++;	/* Goto next increment */
	}

	return (i);	/* Returns the number of increments found */
}

int GMT_strtok (const char *string, const char *sep, int *pos, char *token) {
	/* Reentrant replacement for strtok that uses no static variables.
	 * Breaks string into tokens separated by one of more separator
	 * characters (in sep).  Set *pos to 0 before first call.  Unlike
	 * strtok, always pass the original string as first argument.
	 * Returns 1 if it finds a token and 0 if no more tokens left.
	 * pos is updated and token is returned.  char *token must point
	 * to memory of length >= strlen (string).
	 */

	int i, j, string_len;

	string_len = strlen (string);

	/* Wind up *pos to first non-separating character: */
	while (string[*pos] && strchr (sep, (int)string[*pos])) (*pos)++;

	token[0] = 0;	/* Initialize token to NULL in case we are at end */

	if (*pos >= string_len || string_len == 0) return 0;	/* Got NULL string or no more string left to search */

	/* Search for next non-separating character */
	i = *pos; j = 0;
	while (string[i] && !strchr (sep, (int)string[i])) token[j++] = string[i++];
	token[j] = 0;	/* Add terminating \0 */

	/* Wind up *pos to next non-separating character */
	while (string[i] && strchr (sep, (int)string[i])) i++;
	*pos = i;

	return 1;
}

/* -------------------------------------------------------------------- */
int decode_R (char *item, double *w, double *e, double *s, double *n) {
	char *text, string[BUFSIZ];
	
	/* Minimalist code to decode option -R extracted from GMT_get_common_args */
	
	int i, error = 0;
	double *p[4];
	
	p[0] = w;	p[1] = e;	p[2] = s;	p[3] = n;
			
	i = 0;
	strcpy (string, &item[2]);
	text = strtok (string, "/");
	while (text) {
		*p[i] = ddmmss_to_degree (text);
		i++;
		text = strtok (CNULL, "/");
	}
	if (item[strlen(item)-1] == 'r')	/* Rectangular box given, but valid here */
		error++;
	if (i != 4 || check_region (*p[0], *p[1], *p[2], *p[3]))
		error++;
	w = p[0];	e = p[1];
	s = p[2];	n = p[3];
	return (error);
}

/* -------------------------------------------------------------------- */
int check_region (double w, double e, double s, double n) {
	/* If region is given then we must have w < e and s < n */
	return ((w >= e || s >= n));
}

/* -------------------------------------------------------------------- */
double ddmmss_to_degree (char *text) {
	int i, colons = 0, suffix;
	double degree, minute, degfrac, second;

	for (i = 0; text[i]; i++) if (text[i] == ':') colons++;
	suffix = (int)text[i-1];	/* Last character in string */
	if (colons == 2) {	/* dd:mm:ss format */
		sscanf (text, "%lf:%lf:%lf", &degree, &minute, &second);
		degfrac = degree + Loc_copysign (minute / 60.0 + second / 3600.0, degree);
	}
	else if (colons == 1) {	/* dd:mm format */
		sscanf (text, "%lf:%lf", &degree, &minute);
		degfrac = degree + Loc_copysign (minute / 60.0, degree);
	}
	else
		degfrac = atof (text);
	if (suffix == 'W' || suffix == 'w' || suffix == 'S' || suffix == 's') degfrac = -degfrac;	/* Sign was given implicitly */
	return (degfrac);
}

#define EQ_RAD 6371.0087714
#define M_PR_DEG (EQ_RAD * 1000 * M_PI / 180.0)

void GMT_RI_prepare (struct GRD_HEADER *h) {
	int one_or_zero;
	double s = 1.0, f, m_pr_degree;

	/* May have to adjust -R -I depending on how GMT_inc_code was set */

	one_or_zero = !h->node_offset;
	m_pr_degree = M_PR_DEG;
	h->xy_off = 0.5 * h->node_offset;	/* Use to calculate mean location of block */

	/* XINC AND XMIN/XMAX CHECK FIRST */

	if (GMT_inc_code[0] == 0) {	/* Standard -R -I given, just set nx */
		h->nx = GMT_get_n (h->x_min, h->x_max, h->x_inc, h->node_offset);
	}
	else if (GMT_inc_code[0] & GMT_INC_IS_NNODES) {	/* Got nx */
		h->nx = irint (h->x_inc);
		h->x_inc = GMT_get_inc (h->x_min, h->x_max, h->nx, h->node_offset);
	}
	else {	/* Got funny units */
		switch (GMT_inc_code[0] & GMT_INC_UNITS) {
			case GMT_INC_IS_M:	/* Meter */
				s = 1.0;
				break;
			case GMT_INC_IS_KM:	/* km */
				s = 1000.0;
				break;
			case GMT_INC_IS_MILES:	/* miles */
				s = 1609.433;
				break;
			case GMT_INC_IS_NMILES:	/* nmiles */
				s = 1852.0;
				break;
		}
		if (GMT_inc_code[0] & GMT_INC_IS_EXACT && !(GMT_inc_code[0] & GMT_INC_UNITS)) {
			f = m_pr_degree = 1.0;
		}
		else
			f = cosd (0.5 * (h->y_max + h->y_min));	/* Latitude scaling of E-W distances */

		h->x_inc = h->x_inc * s / (m_pr_degree * f);
		h->nx = GMT_get_n (h->x_min, h->x_max, h->x_inc, h->node_offset);
	}
	if (GMT_inc_code[0] & GMT_INC_IS_EXACT) {	/* Want to keep dx exactly as given; adjust x_max accordingly */
		s = (h->x_max - h->x_min) - h->x_inc * (h->nx - one_or_zero);
		if (fabs (s) > 0.0)
			h->x_max -= s;
	}
	else if (!(GMT_inc_code[0] & GMT_INC_IS_NNODES)) {	/* Adjust x_inc to exactly fit west/east */
		s = h->x_max - h->x_min;
		h->nx = irint (s / h->x_inc);
		f = s / h->nx;
		h->nx += one_or_zero;
		if (fabs (f - h->x_inc) > 0.0)
			h->x_inc = f;
	}

	/* YINC AND YMIN/YMAX CHECK SECOND */
	s = 1.0;	/* s was used above with a different purpose */

	if (GMT_inc_code[1] == 0) {	/* Standard -R -I given, just set ny */
		h->ny = GMT_get_n (h->y_min, h->y_max, h->y_inc, h->node_offset);
	}
	else if (GMT_inc_code[1] & GMT_INC_IS_NNODES) {	/* Got ny */
		h->ny = irint (h->y_inc);
		h->y_inc = GMT_get_inc (h->y_min, h->y_max, h->ny, h->node_offset);
		return;
	}
	else {	/* Got funny units */
		switch (GMT_inc_code[1] & GMT_INC_UNITS) {
			case GMT_INC_IS_M:	/* Meter */
				s = 1.0;
				break;
			case GMT_INC_IS_KM:	/* km */
				s = 1000.0;
				break;
			case GMT_INC_IS_MILES:	/* miles */
				s = 1609.433;
				break;
			case GMT_INC_IS_NMILES:	/* nmiles */
				s = 1852.0;
				break;
		}
		if (GMT_inc_code[1] & GMT_INC_IS_EXACT && !(GMT_inc_code[1] & GMT_INC_UNITS))
			m_pr_degree = 1.0;

		else	/* m_pr_degree might have been reset to 1 in the XINC ... case */
			m_pr_degree = M_PR_DEG;

		h->y_inc = (h->y_inc == 0.0) ? h->x_inc : h->y_inc * s / m_pr_degree;
		h->ny = GMT_get_n (h->y_min, h->y_max, h->y_inc, h->node_offset);
	}

	if (GMT_inc_code[1] & GMT_INC_IS_EXACT) {	/* Want to keep dy exactly as given; adjust y_max accordingly */
		s = (h->y_max - h->y_min) - h->y_inc * (h->ny - one_or_zero);
		if (fabs (s) > 0.0)
			h->y_max -= s;
	}
	else if (!(GMT_inc_code[1] & GMT_INC_IS_NNODES)) {	/* Adjust y_inc to exactly fit south/north */
		s = h->y_max - h->y_min;
		h->ny = irint (s / h->y_inc);
		f = s / h->ny;
		h->ny += one_or_zero;
		if (fabs (f - h->y_inc) > 0.0)
			h->y_inc = f;
	}
}

