/*--------------------------------------------------------------------
 *	$Id: surface.c,v 1.18 2004/04/25 09:10:46 pwessel Exp $
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
 * Version 3.0 1 Nov 1988:  New argument switches standardized.
 *	W. H. F. Smith
 *
 * Version 3.2 5 Dec 1988:	Search radius default = 0.0
 *	to skip this step - unnecessary now that planes are
 *	removed, unless grid starts at 1 for diabolical prime lattice.
 *	Ver 3.1 removed best fit plane, undocumented changes.
 *	W. H. F. Smith
 *
 * Version 4.0 -- 5 Dec 1988:	The scratch file has been eliminated,
 *	resulting in a small gain in speed.   It is replaced by the 
 *	in core struct BRIGGS.  Division by xinc or grid_xinc has been
 *	replaced with multiplies by reciprocals, which also made a small
 *	improvement in speed.  The computation of i, j to index a point
 *	was changed from rint(x) to floor(x + 0.5) to make it follow the
 *	convention of blockmean.  This is actually significant at times.
 *	The significance was discovered when the routine "throw_away_unusables"
 *	was written.  This routine and struct BRIGGS are the major differences
 *	between V3.2 and V4.0.  Modifications by w.h.f. smith
 *
 * V 4.1 -- 26 July, 1989:	W.H.F. Smith fixed up the arg loop and usage
 *	message, and added automatic scaling of the range of z values as an
 *	experiment.  Often this improves the accuracy/convergence of numerical
 *	work; simple testing today suggests it makes no difference, but I will
 *	leave it in.
 *
 * V 4.2 -- 4 June, 1991-2000:	P. Wessel upgraded to GMT v2.0 grdfile i/o.
 *	Added feature to constrain solution to be within lower and/or
 *	upper bounds.  Bounds can be min/max input data value, another
 *	value outside the input data range, or be provided by a grdfile
 *	whose values all are outside the input data range.
 *
 * V 4.3 -- 26 February, 1992:	W. H. F. Smith added option to suggest better
 *	dimensions, and fixed bug in -L option to -Lu -Ll so that full paths
 *	to lower/upper bound files may be specified.
 *
 * 21may98 whfs changed warning about prime dimensions to make it use -V
 *
 * 21-JUN-1998 PW: Upgraded to GMT 3.1 with -b option
 *
 * 24-JAN-2000 PW/WHFS: Fixed bug related to use of -L option for "near-node" constraints.
 * 10-JUL-2000 3.3.5 PW: Added -L plain.
 * Version:	4
 *
 *
 */
 /*--------------------------------------------------------------------
 * Mexified version of surface
 * Author:	J. Luis
 * Date: 	24 Oct 2004
 *
 * Usage
 * [Zout,head] = surface_m(x,y,z,'options');
 * 	OR
 * [Zout,head] = surface_m('input_file','options');
 *
 * where	x,y,z are vectors with the points used in the interpolation 
 * 		'input_file' a file with the points used in the interpolation 
 *		and "options" must have at least the '-R...' and '-I<step>'
 *
 * Note: For a unknown reason the verbose state is allwys set (whether we us -V or not)
 *	So I force gmtdefs.verbose = 0 at the begining of arguments parsing.
 *
 *	Also, x,y,z must be doubles for (and again I don't know why) trying to cast the
 *	(for example) plhs[0] like (float *)plhs[0], though compiles well but gives 
 *	"bad things".
 *
 *	Again something in GMT libs did not work. Now was GMT_input, so I had to write a
 *	reading routine. Unfortunately, this also implies that binary or multisegment
 *	files cannot be read.
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */

#include "gmt.h"
#include "mex.h"

#define OUTSIDE 2000000000	/* Index number indicating data is outside usable area */
 
int n_alloc = GMT_CHUNK;
int npoints=0;			/* Number of data points */
int nx=0;			/* Number of nodes in x-dir. */
int ny=0;			/* Number of nodes in y-dir. (Final grid) */
int mx = 0;
int my = 0;
int ij_sw_corner, ij_se_corner, ij_nw_corner, ij_ne_corner;
int block_nx;			/* Number of nodes in x-dir for a given grid factor */
int block_ny;			/* Number of nodes in y-dir for a given grid factor */
int max_iterations=250;		/* Max iter per call to iterate */
int total_iterations = 0;
int grid, old_grid;		/* Node spacings  */
int grid_east;
int n_fact = 0;			/* Number of factors in common (ny-1, nx-1) */
int factors[32];		/* Array of common factors */
int long_verbose = FALSE;
int n_empty;			/* No of unconstrained nodes at initialization  */
int set_low = 0;		/* 0 unconstrained,1 = by min data value, 2 = by user value */
int set_high = 0;		/* 0 unconstrained,1 = by max data value, 2 = by user value */
int constrained = FALSE;	/* TRUE if set_low or set_high is TRUE */
double low_limit, high_limit;	/* Constrains on range of solution */
double x_min, x_max, y_min, y_max;	/* minmax coordinates */
float *lower, *upper;		/* arrays for minmax values, if set */
double xinc, yinc;		/* Size of each grid cell (final size) */
double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
double r_xinc, r_yinc, r_grid_xinc, r_grid_yinc;	/* Reciprocals  */
double converge_limit = 0.0;	/* Convergence limit */
double radius = 0.0;			/* Search radius for initializing grid  */
double	tension = 0.0;		/* Tension parameter on the surface  */
double	boundary_tension = 0.0;
double	interior_tension = 0.0;
double	a0_const_1, a0_const_2;	/* Constants for off grid point equation  */
double	e_2, e_m2, one_plus_e2;
double eps_p2, eps_m2, two_plus_ep2, two_plus_em2;
double	x_edge_const, y_edge_const;
double	l_epsilon = 1.0;
double	z_mean;
double	z_scale = 1.0;		/* Root mean square range of z after removing planar trend  */
double	r_z_scale = 1.0;	/* reciprocal of z_scale  */
double	plane_c0, plane_c1, plane_c2;	/* Coefficients of best fitting plane to data  */
float *u;			/* Pointer to grid array */
float *v2;			/* Pointer to v.2.0 grid array */
char *iu;			/* Pointer to grid info array */
char mode_type[2] = {'I','D'};	/* D means include data points when iterating
				 * I means just interpolate from larger grid */
char format[BUFSIZ];
double	*in0, *in1, *in2;

int	offset[25][12];		/* Indices of 12 nearby points in 25 cases of edge conditions  */
double		coeff[2][12];	/* Coefficients for 12 nearby points, constrained and unconstrained  */

double relax_old, relax_new = 1.4;	/* Coefficients for relaxation factor to speed up convergence */

int compare_points(const void *point_1v, const void *point_2v);

struct DATA {
	float x;
	float y;
	float z;
	int index;
} *data;		/* Data point and index to node it currently constrains  */

struct BRIGGS {
	double b[6];
} *briggs;		/* Coefficients in Taylor series for Laplacian(z) a la I. C. Briggs (1974)  */

struct SUGGESTION {	/* Used to find top ten list of faster grid dimensions  */
	int	nx;
	int	ny;
	double	factor;	/* Speed up by a factor of factor  */
};

FILE	*fp_in = NULL;	/* File pointer  */

int	gcd_euclid(int a, int b);	/* Finds the greatest common divisor  */
int	get_prime_factors(int n, int *f), iterate(int mode);
void set_grid_parameters(void), throw_away_unusables(void), remove_planar_trend(void), rescale_z_values(void);
void load_constraints(char *low, char *high), smart_divide(void), set_offset(void), set_index(void), initialize_grid(void), set_coefficients(void);
void find_nearest_point(void), fill_in_forecast(void), check_errors(void), replace_planar_trend(void);

int to_data(int n_pts), read_data(void);

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	void	suggest_sizes_for_surface(int nx, int ny);
	int	i, j, error = FALSE, size_query = FALSE, n_pts = 0;
	int	argc = 0, n_arg_no_char = 0, index;
	char	modifier, low[128], high[128], **argv;
	float	*pdata;
	double	*info;
	struct GRD_HEADER h;

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
	argv[0] = "surface_m";
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

	GMT_grd_init (&h, argc, argv, FALSE);

	x_min = y_min = 0.0;
	x_max = y_max = 0.0;
	xinc = yinc = 0.0;

	/* New in v4.3:  Default to unconstrained:  */
	set_low = set_high = 0; 

	gmtdefs.verbose = 0;	/* Otherwise it insists in setting it to on all the times */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
              
				/* Common parameters */
                      
				case 'R':
				case 'f':
                                      error += GMT_get_common_args (argv[i], &x_min, &x_max, &y_min, &y_max);
                                      break;
				case ':':
					gmtdefs.xy_toggle[0] = TRUE;
					break;
				case 'V':
					gmtdefs.verbose = 1;
					if (argv[i][2] == 'L' || argv[i][2] == 'l') long_verbose = TRUE;
					break;
				case 'H':
					gmtdefs.n_header_recs = atoi (&argv[i][2]);
					break;

				/* Supplemental parameters */
				case 'A':
					l_epsilon = atof (&argv[i][2]);
					break;
				case 'C':
					converge_limit = atof (&argv[i][2]);
					break;
				case 'I':
					GMT_getinc (&argv[i][2], &xinc, &yinc);
					break;
				case 'L':	/* Set limits */
						/* This is new, to use -Ll and -Lu:  */
					switch (argv[i][2]) {
						case 'l':
							/* Lower limit  */
							if (argv[i][3] == 0) {
								mexPrintf("%s: GMT SYNTAX ERROR -Ll option: No argument given\n", GMT_program);
								error++;
							}
							strcpy (low, &argv[i][3]);
							if (!access (low, R_OK))	/* file exists */
								set_low = 3;
							else if (low[0] == 'd')
								set_low = 1;
							else {
								set_low = 2;
								low_limit = atof (&argv[i][3]);
							}
							break;
						case 'u':
							/* Upper limit  */
							if (argv[i][3] == 0) {
								mexPrintf("%s: GMT SYNTAX ERROR -Lu option: No argument given\n", GMT_program);
								error++;
							}
							strcpy (high, &argv[i][3]);
							if (!access (high, R_OK))	/* file exists */
								set_high = 3;
							else if (high[0] == 'd')
								set_high = 1;
							else {
								set_high = 2;
								high_limit = atof (&argv[i][3]);
							}
							break;
						default:	/* 360-periodicity option */
							GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
							GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
							mexPrintf ("%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);
							break;
					}
					break;
				case 'N':
					max_iterations = atoi (&argv[i][2]);
					break;
				case 'S':
					radius = atof (&argv[i][2]);
					modifier = argv[i][strlen(argv[i])-1];
					if (modifier == 'm' || modifier == 'M') radius /= 60.0;
					break;
				case 'T':
					modifier = argv[i][strlen(argv[i])-1];
					if (modifier == 'b' || modifier == 'B') {
						boundary_tension = atof (&argv[i][2]);
					}
					else if (modifier == 'i' || modifier == 'I') {
						interior_tension = atof (&argv[i][2]);
					}
					else if (modifier >= '0' && modifier <= '9') {
						tension = atof (&argv[i][2]);
					}
					else {
						mexPrintf("%s: GMT SYNTAX ERROR -T option: Unrecognized modifier %c\n", GMT_program, modifier);
						error = TRUE;
					}
					break;
				case 'Q':
					size_query = TRUE;
					break;
				case 'Z':
					relax_new = atof (&argv[i][2]);
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
		else {
			if (n_arg_no_char == 0) {
				if ((fp_in = fopen(argv[i], "r")) == NULL) {
					mexPrintf ("surface: cannot open input data file %s\n", argv[i]);
					mexErrMsgTxt("\n");
				}
			}
		}
	}
	
	if (argc == 1 || GMT_give_synopsis_and_exit || error) {	/* Display usage */
		mexPrintf ("surface - Adjustable tension continuous curvature surface gridding\n\n");
		mexPrintf ("usage: [Zout,head] = surface_m(x,y,z|<xyz-file>, '-I<xinc>[m|c][/<yinc>[m|c]]',\n");
		mexPrintf ("\t'-R<west>/<east>/<south>/<north>', '[-A<aspect_ratio>]', '[-C<convergence_limit>]',\n");
		mexPrintf ("\t'[-Ll<limit>]', '[-Lu<limit>]', '[-N<n_iterations>]', '[-S<search_radius>[m]]', '[-T<tension>[i][b]]',\n");
		mexPrintf ("\t'[-Q]', '[-V[l]]', '[-Z<over_relaxation_parameter>]', '[-f[i|o]<colinfo>]')\n\n");
		
		if (GMT_give_synopsis_and_exit) return;
		
		mexPrintf ("\tsurface will use provided x,y,z vectors or a single <xyz-file>.\n\n");
		mexPrintf ("\tRequired arguments to surface:\n");
		mexPrintf ("\t-I sets the Increment of the grid; enter xinc, optionally xinc/yinc.\n");
		mexPrintf ("\t\tDefault is yinc = xinc.  Append an m [or c] to xinc or yinc to indicate minutes [or seconds]\n");
		mexPrintf ("\t\te.g.  -I10m/5m grids longitude every 10 minutes, latitude every 5 minutes.\n\n");
		/*GMT_explain_option ('R');*/
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-A<aspect_ratio>  = 1.0  by default which gives an isotropic solution.\n");
		mexPrintf ("\t\ti.e. xinc and yinc assumed to give derivatives of equal weight; if not, specify\n");
		mexPrintf ("\t\t<aspect_ratio> such that yinc = xinc / <aspect_ratio>.\n");
		mexPrintf ("\t\te.g. if gridding lon,lat use <aspect_ratio> = cosine(middle of lat range).\n");
		mexPrintf ("\t-C<convergence_limit> iteration stops when max abs change is less than <c.l.>\n");
		mexPrintf ("\t\tdefault will choose 0.001 of the range of your z data (1 ppt precision).\n");
		mexPrintf ("\t\tEnter your own convergence limit in same units as z data.\n");
		mexPrintf ("\t-L constrain the range of output values:\n");
		mexPrintf ("\t\t-Ll<limit> specifies lower limit; forces solution to be >= <limit>.\n");
		mexPrintf ("\t\t-Lu<limit> specifies upper limit; forces solution to be <= <limit>.\n");
		mexPrintf ("\t\t<limit> can be any number, or the letter d for min (or max) input data value,\n");
		mexPrintf ("\t\tor the filename of a grdfile with bounding values.  [Default solution unconstrained].\n");
		mexPrintf ("\t\tExample:  -Ll0 gives a non-negative solution.\n");
		mexPrintf ("\t-N sets max <n_iterations> in each cycle; default = 250.\n");
		mexPrintf ("\t-S sets <search_radius> to initialize grid; default = 0 will skip this step.\n");
		mexPrintf ("\t\tThis step is slow and not needed unless grid dimensions are pathological;\n");
		mexPrintf ("\t\ti.e., have few or no common factors.\n");
		mexPrintf ("\t\tAppend m to give <search_radius> in minutes.\n");
		mexPrintf ("\t-T adds Tension to the gridding equation; use a value between 0 and 1.\n");
		mexPrintf ("\t\tdefault = 0 gives minimum curvature (smoothest; bicubic) solution.\n");
		mexPrintf ("\t\t1 gives a harmonic spline solution (local max/min occur only at data points).\n");
		mexPrintf ("\t\ttypically 0.25 or more is good for potential field (smooth) data;\n");
		mexPrintf ("\t\t0.75 or so for topography.  Experiment.\n");
		mexPrintf ("\t\tAppend B or b to set tension in boundary conditions only;\n");
		mexPrintf ("\t\tAppend I or i to set tension in interior equations only;\n");
		mexPrintf ("\t\tNo appended letter sets tension for both to same value.\n");
		mexPrintf ("\t-Q Query for grid sizes that might run faster than your -R -I give.\n");
		mexPrintf ("\t\tAppend l for long verbose\n");
		mexPrintf ("\t-Z sets <over_relaxation parameter>.  Default = 1.4\n");
		mexPrintf ("\t\tUse a value between 1 and 2.  Larger number accelerates convergence but can be unstable.\n");
		mexPrintf ("\t\tUse 1 if you want to be sure to have (slow) stable convergence.\n\n");
		/*GMT_explain_option ('i');*/
		/*GMT_explain_option ('n');*/
		mexPrintf ("\t\tDefault is 3 input columns.\n\n");
		/*GMT_explain_option ('f');*/
		mexPrintf ("\t(For additional details, see Smith & Wessel, Geophysics, 55, 293-305, 1990.)\n");
		return;
	}

	if (!project_info.region_supplied) {
		mexPrintf ("%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program);
		error++;
	}
	if (xinc <= 0.0 || yinc <= 0.0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program);
		error++;
	}
	if (max_iterations < 1) {
		mexPrintf ("%s: GMT SYNTAX ERROR -N option.  Max iterations must be nonzero\n", GMT_program);
		error++;
	}
	if (relax_new < 1.0 || relax_new > 2.0) {
		mexPrintf ("%s: GMT SYNTAX ERROR -Z option.  Relaxation value must be 1 <= z <= 2\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && gmtdefs.io_header[GMT_IN]) {
		mexPrintf ("%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program);
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
		mexPrintf ("%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 3 columns\n", GMT_program);
		error++;
	}
	
	if (error) mexErrMsgTxt("\n");

	if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("SURFACE ERROR: Must provide one or two outputs.\n");
	if (n_arg_no_char > 0) {
		if (n_arg_no_char != 3)
			mexErrMsgTxt("SURFACE ERROR: Must provide three numeric inputs (x,y,z)\n");

		/* Check that first argument contains at least a mx3 table */
		n_pts = mxGetM (prhs[0]);
		if (!mxIsNumeric(prhs[0]) | !mxIsNumeric(prhs[1]) | !mxIsNumeric(prhs[2]))
			mexErrMsgTxt("SURFACE ERROR: first 3 input args must contain the x,y,z triplets to interpolate.\n");
	}

	if (n_arg_no_char > 0) {		/* Input data was transmited in input*/
		/* Read the input points and convert them to double */
		if (mxIsDouble(prhs[0])) {
			in0 = mxGetPr(prhs[0]);
			in1 = mxGetPr(prhs[1]);
			in2 = mxGetPr(prhs[2]);
		}
		else if (mxIsSingle(prhs[0])) {
			in0 = mxGetData(prhs[0]);
			in1 = mxGetData(prhs[1]);
			in2 = mxGetData(prhs[2]);
		}
	}

	h.x_min = x_min;
	h.x_max = x_max;
	h.y_min = y_min;
	h.y_max = y_max;
	h.x_inc = xinc;
	h.y_inc = yinc;

	/*GMT_grd_RI_verify (&h, 1);*/		/* IF (IVAN == TRUE)  ==> Matlab = BOOM */

	if (tension != 0.0) {
		boundary_tension = tension;
		interior_tension = tension;
	}
	relax_old = 1.0 - relax_new;

	nx = irint ((x_max - x_min)/xinc) + 1;
	ny = irint ((y_max - y_min)/yinc) + 1;
	h.nx = nx;
	h.ny = ny;
	mx = nx + 4;
	my = ny + 4;
	r_xinc = 1.0 / xinc;
	r_yinc = 1.0 / yinc;

	/* New stuff here for v4.3:  Check out the grid dimensions:  */
	grid = gcd_euclid (nx-1, ny-1);

	if (gmtdefs.verbose || size_query) {
		sprintf (format, "W: %s E: %s S: %s N: %s nx: %%d ny: %%d\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		mexPrintf (format, x_min, x_max, y_min, y_max, nx-1, ny-1);
	}
	if (grid == 1 && gmtdefs.verbose) mexPrintf("%s:  WARNING:  Your grid dimensions are mutually prime.\n", GMT_program);
	if (( grid == 1 && gmtdefs.verbose) || size_query) suggest_sizes_for_surface(nx-1, ny-1);
	if (size_query) return;

	/* New idea: set grid = 1, read data, setting index.  Then throw
		away data that can't be used in end game, constraining
		size of briggs->b[6] structure.  */
	
	grid = 1;
	set_grid_parameters();
	if (n_arg_no_char == 0)	{	/* Input data will be read inside next subroutine*/
		if (read_data()) mexErrMsgTxt("\n");
	}
	else {				/* Input data was transmited in input and will be copyied to the data struct*/
		if (to_data(n_pts)) mexErrMsgTxt("\n");
	}

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
	
	briggs = (struct BRIGGS *) GMT_memory (VNULL, (size_t)npoints, sizeof(struct BRIGGS), GMT_program);
	iu = (char *) GMT_memory (VNULL, (size_t)(mx * my), sizeof(char), GMT_program);
	u = (float *) GMT_memory (VNULL, (size_t)(mx * my), sizeof(float), GMT_program);

	if (radius > 0) initialize_grid(); /* Fill in nodes with a weighted avg in a search radius  */

	if (gmtdefs.verbose) mexPrintf("Grid\tMode\tIteration\tMax Change\tConv Limit\tTotal Iterations\n");
	
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
	
	if (gmtdefs.verbose) check_errors ();

	replace_planar_trend();
	
	GMT_free ((void *) data);
	GMT_free ((void *) briggs);
	GMT_free ((void *)iu);
	if (set_low) GMT_free ((void *)lower);
	if (set_high) GMT_free ((void *)upper);

	/*write_output(&h, grdfile);*/

	v2 = (float *) GMT_memory (VNULL, (size_t)(nx * ny), sizeof (float), GMT_program);
	index = ij_sw_corner;
	for(i = 0; i < nx; i++, index += my) for (j = 0; j < ny; j++) v2[j*nx+i] = u[index + ny - j - 1];
	/* Transpose from gmt grd orientation to Matlab orientation */
	for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) u[j*ny+ny-i-1] = v2[i*nx+j];

	plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
	pdata = mxGetData(plhs[0]);
	memcpy(pdata, u, ny*nx * 4);

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

	GMT_free ((void *) u);
	GMT_free ((void *)v2);
	GMT_end (argc, argv);
}

void	set_coefficients(void) {
	double	e_4, loose, a0;
	
	loose = 1.0 - interior_tension;
	e_2 = l_epsilon * l_epsilon;
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

void	set_offset(void) {
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


void fill_in_forecast (void) {

	/* Fills in bilinear estimates into new node locations
	   after grid is divided.   */

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

int compare_points (const void *point_1v, const void *point_2v)
/*struct DATA *point_1, *point_2; { */
{
		/*  Routine for qsort to sort data structure for fast access to data by node location.
		    Sorts on index first, then on radius to node corresponding to index, so that index
		    goes from low to high, and so does radius.
		*/
	int block_i, block_j, index_1, index_2;
	double x0, y0, dist_1, dist_2;
	
struct DATA *point_1, *point_2;
	point_1 = (struct DATA *)point_1v;
	point_2 = (struct DATA *)point_2v;
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
		x0 = x_min + block_i * grid_xinc;
		y0 = y_min + block_j * grid_yinc;
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

void smart_divide (void) {
		/* Divide grid by its largest prime factor */
	grid /= factors[n_fact - 1];
	n_fact--;
}

void set_index (void) {
		/* recomputes data[k].index for new value of grid,
		   sorts data on index and radii, and throws away
		   data which are now outside the usable limits. */
	int i, j, k, k_skipped = 0;

	for (k = 0; k < npoints; k++) {
		i = (int)floor(((data[k].x-x_min)*r_grid_xinc) + 0.5);
		j = (int)floor(((data[k].y-y_min)*r_grid_yinc) + 0.5);
		if (i < 0 || i >= block_nx || j < 0 || j >= block_ny) {
			data[k].index = OUTSIDE;
			k_skipped++;
		}
		else
			data[k].index = i * block_ny + j;
	}
	
	qsort ((void *)data, (size_t)npoints, sizeof (struct DATA), compare_points);
	
	npoints -= k_skipped;
	
}

void find_nearest_point(void) {
	int i, j, ij_v2, k, last_index, block_i, block_j, iu_index, briggs_index;
	double x0, y0, dx, dy, xys, xy1, btemp;
	double b0, b1, b2, b3, b4, b5;
	float z_at_node;
	
	last_index = -1;

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
	 		x0 = x_min + block_i*grid_xinc;
	 		y0 = y_min + block_j*grid_yinc;
	 		dx = (data[k].x - x0)*r_grid_xinc;
	 		dy = (data[k].y - y0)*r_grid_yinc;
	 		if (fabs(dx) < 0.05 && fabs(dy) < 0.05) {	/* Close enough to assign value to node */
	 			iu[iu_index] = 5;
	 			/* v3.3.4: NEW CODE
	 			 * Since point is basically moved from (dx, dy) to (0,0) we must adjust for
	 			 * the small change in the planar trend between the two locations, and then
	 			 * possibly clip the range if constraining surfaces were given.  Note that
	 			 * dx, dy is in -1/1 range normalized by (grid * x|y_inc) so to recover the
	 			 * dx,dy in final grid fractions we must scale by grid */
	 			 
	 			z_at_node = data[k].z + (float) (r_z_scale * grid * (plane_c1 * dx + plane_c2 * dy));
	 			if (constrained) {
					ij_v2 = (ny - block_j * grid - 1) * nx + block_i * grid;
					if (set_low  && !GMT_is_fnan (lower[ij_v2]) && z_at_node < lower[ij_v2])
						z_at_node = lower[ij_v2];
					else if (set_high && !GMT_is_fnan (upper[ij_v2]) && z_at_node > upper[ij_v2])
						z_at_node = upper[ij_v2];
	 			}
	 			u[iu_index] = z_at_node;
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

						
void set_grid_parameters(void)
{			
	block_ny = (ny - 1) / grid + 1;
	block_nx = (nx - 1) / grid + 1;
	grid_xinc = grid * xinc;
	grid_yinc = grid * yinc;
	grid_east = grid * my;
	r_grid_xinc = 1.0 / grid_xinc;
	r_grid_yinc = 1.0 / grid_yinc;
}

void initialize_grid(void)
{	/*
	 * For the initial gridsize, compute weighted averages of data inside the search radius
	 * and assign the values to u[i,j] where i,j are multiples of gridsize.
	 */
	 int	irad, jrad, i, j, imin, imax, jmin, jmax, index_1, index_2, k, ki, kj, k_index;
	 double	r, rfact, sum_w, sum_zw, weight, x0, y0;

	 irad = (int)ceil(radius/grid_xinc);
	 jrad = (int)ceil(radius/grid_yinc);
	 rfact = -4.5/(radius*radius);
	 
	 for (i = 0; i < block_nx; i ++ ) {
	 	x0 = x_min + i*grid_xinc;
	 	for (j = 0; j < block_ny; j ++ ) {
	 		y0 = y_min + j*grid_yinc;
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
	 			sprintf (format, "%%s: Warning: no data inside search radius at: %s %s\n", gmtdefs.d_format, gmtdefs.d_format);
	 			mexPrintf (format, GMT_program, x0, y0);
	 			u[ij_sw_corner + (i * my + j) * grid] = (float)z_mean;
	 		}
	 		else {
	 			u[ij_sw_corner + (i*my+j)*grid] = (float)(sum_zw/sum_w);
	 		}
		}
	}
}


void new_initialize_grid(void) {
	/* For the initial gridsize, load constrained nodes with weighted avg of their data;
	 * and then do something with the unconstrained ones.  */
	 int	k, k_index, u_index, block_i, block_j;
	 double	sum_w, sum_zw, weight, x0, y0, dx, dy, dx_scale, dy_scale;

	dx_scale = 4.0 / grid_xinc;
	dy_scale = 4.0 / grid_yinc;
	n_empty = block_ny * block_nx;
	k = 0;
	while (k < npoints) {
		block_i = data[k].index / block_ny;
		block_j = data[k].index % block_ny;
		x0 = x_min + block_i*grid_xinc;
		y0 = y_min + block_j*grid_yinc;
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

int to_data(n_pts) {
	/* Copy the vectors transmited as inputs into the data structure. For large inputs,
	this implies a significant memory wasting.*/
	int	i, j, k, n, kmax, kmin;
	double	zmin = DBL_MAX, zmax = -DBL_MAX;

	data = (struct DATA *) GMT_memory (VNULL, (size_t)n_pts, sizeof(struct DATA), GMT_program);
	
	/* Read in xyz data and computes index no and store it in a structure */
	
	k = 0;
	z_mean = 0;
			
	for (n = 0; n < n_pts; n++) {
		if (GMT_is_dnan (in2[n])) continue;
		
		i = (int)floor((((float)in0[n]-x_min)*r_grid_xinc) + 0.5);
		if (i < 0 || i >= block_nx) continue;
		j = (int)floor((((float)in1[n]-y_min)*r_grid_yinc) + 0.5);
		if (j < 0 || j >= block_ny) continue;

		data[k].index = i * block_ny + j;
		data[k].x = (float)in0[n];
		data[k].y = (float)in1[n];
		data[k].z = (float)in2[n];
		if (zmin > in2[n]) zmin = in2[n], kmin = k;
		if (zmax < in2[n]) zmax = in2[n], kmax = k;
		k++;
		z_mean += in2[n];
	}
	
	npoints = k;
	
	if (npoints == 0) {
		mexPrintf ("%s:  No datapoints inside region, aborts\n", GMT_program);
		return (1);
	}
	
	z_mean /= k;
	if (gmtdefs.verbose) {
		sprintf(format, "%s %s %s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		mexPrintf("%s: Minimum value of your dataset x,y,z at: ", GMT_program);
		mexPrintf(format, (double)data[kmin].x, (double)data[kmin].y, (double)data[kmin].z);
		mexPrintf("%s: Maximum value of your dataset x,y,z at: ", GMT_program);
		mexPrintf(format, (double)data[kmax].x, (double)data[kmax].y, (double)data[kmax].z);
	}
	data = (struct DATA *) GMT_memory ((void *)data, (size_t)npoints, sizeof(struct DATA), GMT_program);
	
	if (set_low == 1)
		low_limit = data[kmin].z;
	else if (set_low == 2 && low_limit > data[kmin].z)
		mexPrintf ("%s: Warning:  Your lower value is > than min data value.\n", GMT_program);

	if (set_high == 1)
		high_limit = data[kmax].z;
	else if (set_high == 2 && high_limit < data[kmax].z)
		mexPrintf ("%s: Warning:  Your upper value is < than max data value.\n", GMT_program);

	return (0);

}

int read_data(void) {
	/* Again something in GMT libs did not work. Now was GMT_input, so I had to make the
	adaptions in this routine. Unfortunately, this also implies that binary or multisegment
	files cannot be read.*/
	int	i, j, ix, iy, jj, k, kmax, kmin, n_fields, n_expected_fields, n_cols = 0;
	double	*in, zmin = DBL_MAX, zmax = -DBL_MAX;
	char	line[1024], buffer[BUFSIZ], *p;

	data = (struct DATA *) GMT_memory (VNULL, (size_t)n_alloc, sizeof(struct DATA), GMT_program);
	
	/* Read in xyz data and computes index no and store it in a structure */
	
	k = 0;
	z_mean = 0;
	in = (double *) calloc ((size_t)(n_alloc), sizeof(double));
			
	/* Use here the old toggle recipe because things inside GMT_get_common_args are behaving odly */
	if (gmtdefs.xy_toggle[0]) {
		ix = 1;		iy = 0;
	}
	else {
		ix = 0;		iy = 1;
	}
	for (i = 0; i < gmtdefs.n_header_recs; i++) fgets (line, 1024, fp_in);

	while (fgets (line, 1024, fp_in)) {
		if (n_cols == 0) {	/* First time, allocate # of columns */
			strcpy (buffer, line);
			p = (char *)strtok (buffer, " \t\n");
			while (p) {	/* Count # of fields */
				n_cols++;
				p = (char *)strtok ((char *)NULL, " \t\n");
			}
			/* Now we know # of columns */
			if (n_cols < 3) {
				mexPrintf("SURFACE ERROR: input file must have at least 3 columns");
				return(1);
			}
		}
		p = (char *)strtok (line, " \t\n");
		jj = 0;
		while (p && jj < 3) {
			sscanf (p, "%lf", &in[jj]);
			jj++;
			p = (char *)strtok ((char *)NULL, " \t\n");
		}
		if (jj != 3) mexPrintf ("Expected %d but found %d fields in record # %d\n", 3, jj, k);

		if (GMT_is_dnan (in[2])) continue;
		
		i = (int)floor(((in[ix]-x_min)*r_grid_xinc) + 0.5);
		if (i < 0 || i >= block_nx) continue;
		j = (int)floor(((in[iy]-y_min)*r_grid_yinc) + 0.5);
		if (j < 0 || j >= block_ny) continue;

		data[k].index = i * block_ny + j;
		data[k].x = (float)in[ix];
		data[k].y = (float)in[iy];
		data[k].z = (float)in[2];
		if (zmin > in[2]) zmin = in[2], kmin = k;
		if (zmax < in[2]) zmax = in[2], kmax = k;
		k++;
		z_mean += in[2];
		if (k == n_alloc) {
			n_alloc += GMT_CHUNK;
			data = (struct DATA *) GMT_memory ((void *)data, (size_t)n_alloc, sizeof(struct DATA), GMT_program);
		}
	}
	
	fclose (fp_in);

	npoints = k;
	
	if (npoints == 0) {
		mexPrintf ("%s:  No datapoints inside region, aborts\n", GMT_program);
		return (1);
	}
	
	z_mean /= k;
	if (gmtdefs.verbose) {
		sprintf(format, "%s %s %s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		mexPrintf("%s: Minimum value of your dataset x,y,z at: ", GMT_program);
		mexPrintf(format, (double)data[kmin].x, (double)data[kmin].y, (double)data[kmin].z);
		mexPrintf("%s: Maximum value of your dataset x,y,z at: ", GMT_program);
		mexPrintf(format, (double)data[kmax].x, (double)data[kmax].y, (double)data[kmax].z);
	}
	data = (struct DATA *) GMT_memory ((void *)data, (size_t)npoints, sizeof(struct DATA), GMT_program);
	
	if (set_low == 1)
		low_limit = data[kmin].z;
	else if (set_low == 2 && low_limit > data[kmin].z) {
	/*	low_limit = data[kmin].z;	*/
		mexPrintf ("%s: Warning:  Your lower value is > than min data value.\n", GMT_program);
	}
	if (set_high == 1)
		high_limit = data[kmax].z;
	else if (set_high == 2 && high_limit < data[kmax].z) {
	/*	high_limit = data[kmax].z;	*/
		mexPrintf ("%s: Warning:  Your upper value is < than max data value.\n", GMT_program);
	}
	return (0);
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
	double	y_denom = 2 * l_epsilon * (1.0 - boundary_tension) + boundary_tension;
	double	y_0_const = 4 * l_epsilon * (1.0 - boundary_tension) / y_denom;
	double	y_1_const = (boundary_tension - 2 * l_epsilon * (1.0 - boundary_tension) ) / y_denom;

	sprintf(format,"%%4d\t%%c\t%%8d\t%s\t%s\t%%10d\n", gmtdefs.d_format, gmtdefs.d_format);

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
					if (set_low && !GMT_is_fnan (lower[ij_v2]) && sum_ij < lower[ij_v2])
						sum_ij = lower[ij_v2];
					else if (set_high && !GMT_is_fnan (upper[ij_v2]) && sum_ij > upper[ij_v2])
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
		if (long_verbose) mexPrintf (format,
			grid, mode_type[mode], iteration_count, max_change, current_limit, total_iterations);

	} while (max_change > current_limit && iteration_count < max_iterations);
	
	if (gmtdefs.verbose && !long_verbose) mexPrintf(format,
		grid, mode_type[mode], iteration_count, max_change, current_limit, total_iterations);

	return(iteration_count);
}

void check_errors (void) {

	int	i, j, k, ij, n_nodes, move_over[12];	/* move_over = offset[kase][12], but grid = 1 so move_over is easy  */
	
	double	x0, y0, dx, dy, mean_error, mean_squared_error, z_est, z_err, curvature, c;
	double	du_dx, du_dy, d2u_dx2, d2u_dxdy, d2u_dy2, d3u_dx3, d3u_dx2dy, d3u_dxdy2, d3u_dy3;
	
	double	x_0_const = 4.0 * (1.0 - boundary_tension) / (2.0 - boundary_tension);
	double	x_1_const = (3 * boundary_tension - 2.0) / (2.0 - boundary_tension);
	double	y_denom = 2 * l_epsilon * (1.0 - boundary_tension) + boundary_tension;
	double	y_0_const = 4 * l_epsilon * (1.0 - boundary_tension) / y_denom;
	double	y_1_const = (boundary_tension - 2 * l_epsilon * (1.0 - boundary_tension) ) / y_denom;
	
	
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
	 	x0 = x_min + i*xinc;
	 	y0 = y_min + j*yinc;
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

	 mexPrintf("Fit info: N data points  N nodes\tmean error\trms error\tcurvature\n");
	 sprintf (format,"\t%%8d\t%%8d\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
	 mexPrintf (format, npoints, n_nodes, mean_error, mean_squared_error, curvature);
 }

void	remove_planar_trend(void) {
	int	i;
	double	a, b, c, d, xx, yy, zz;
	double	sx, sy, sz, sxx, sxy, sxz, syy, syz;
	
	sx = sy = sz = sxx = sxy = sxz = syy = syz = 0.0;
	
	for (i = 0; i < npoints; i++) {

		xx = (data[i].x - x_min) * r_xinc;
		yy = (data[i].y - y_min) * r_yinc;
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

		xx = (data[i].x - x_min) * r_xinc;
		yy = (data[i].y - y_min) * r_yinc;
		
		data[i].z -= (float)(plane_c0 + plane_c1 * xx + plane_c2 * yy);
	}
}

void	replace_planar_trend(void) {
	int	i, j, ij;

	 for (i = 0; i < nx; i++) {
	 	for (j = 0; j < ny; j++) {
	 		ij = ij_sw_corner + i * my + j;
	 		u[ij] = (float)((u[ij] * z_scale) + (plane_c0 + plane_c1 * i + plane_c2 * j));
		}
	}
}

void	throw_away_unusables(void) {
	/* This is a new routine to eliminate data which will become
		unusable on the final iteration, when grid = 1.
		It assumes grid = 1 and set_grid_parameters has been
		called.  We sort, mark redundant data as OUTSIDE, and
		sort again, chopping off the excess.
		
		Experimental modification 5 Dec 1988 by Smith, as part
		of a new implementation using core memory for b[6]
		coefficients, eliminating calls to temp file.  */
	
	int	last_index, n_outside, k;
	
	/* Sort the data  */
	
	qsort ((void *)data, (size_t)npoints, sizeof (struct DATA), compare_points);
	
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
	
	qsort ((void *)data, (size_t)npoints, sizeof (struct DATA), compare_points);
	npoints -= n_outside;
	data = (struct DATA *) GMT_memory ((void *)data, (size_t)npoints, sizeof(struct DATA), GMT_program);
	if (gmtdefs.verbose && (n_outside)) {
		mexPrintf("%s: %d unusable points were supplied; these will be ignored.\n", GMT_program, n_outside);
		mexPrintf("\tYou should have pre-processed the data with block-mean, -median, or -mode.\n");
	}
}

void	rescale_z_values(void) {
	int	i;
	double	ssz = 0.0;

	for (i = 0; i < npoints; i++) ssz += (data[i].z * data[i].z);
	
	/* Set z_scale = rms(z):  */
	
	z_scale = sqrt (ssz / npoints);
	
	if (z_scale == 0.0) {
		mexPrintf("%s: WARNING: Input data lie exactly on a plane - no solution (aborting).\n", GMT_program);
		r_z_scale = z_scale = 1.0;
	}
	else
		r_z_scale = 1.0 / z_scale;

	for (i = 0; i < npoints; i++) data[i].z *= (float)r_z_scale;

	if (converge_limit == 0.0) converge_limit = 0.001 * z_scale; /* i.e., 1 ppt of L2 scale */
}

void load_constraints (char *low, char *high) {
	int i, j, ij;
	double yy;
	struct GRD_HEADER hdr;
		
	/* Load lower/upper limits, verify range, deplane, and rescale */
	
	if (set_low > 0) {
		lower = (float *) GMT_memory (VNULL, (size_t)(nx * ny), sizeof (float), GMT_program);
		if (set_low < 3)
			for (i = 0; i < nx * ny; i++) lower[i] = (float)low_limit;
		else {
			if (GMT_read_grd_info (low, &hdr)) {
				mexPrintf ("%s: Error opening file %s\n", GMT_program, low);
				exit (EXIT_FAILURE);
			}
			if (hdr.nx != nx || hdr.ny != ny) {
				mexPrintf ("%s: lower limit file not of proper dimension!\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			if (GMT_read_grd (low, &hdr, lower, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
				mexPrintf ("%s: Error reading file %s\n", GMT_program, low);
				exit (EXIT_FAILURE);
			}
/* Comment this out:	n_trimmed = 0;
			for (i = 0; i < nx * ny; i++) if (lower[i] > low_limit) {
				lower[i] = low_limit;
				n_trimmed++;
			}
			if (n_trimmed) fprintf (stderr, "%s: %d lower limit values > min data, reset to min data!\n", GMT_program, n_trimmed);
*/
		}
			
		for (j = ij = 0; j < ny; j++) {
			yy = ny - j - 1;
			for (i = 0; i < nx; i++, ij++) {
				if (GMT_is_fnan (lower[ij])) continue;
				lower[ij] -= (float)(plane_c0 + plane_c1 * i + plane_c2 * yy);
				lower[ij] *= (float)r_z_scale;
			}
		}
		constrained = TRUE;
	}
	if (set_high > 0) {
		upper = (float *) GMT_memory (VNULL, (size_t)(nx * ny), sizeof (float), GMT_program);
		if (set_high < 3)
			for (i = 0; i < nx * ny; i++) upper[i] = (float)high_limit;
		else {
			if (GMT_read_grd_info (high, &hdr)) {
				mexPrintf ("%s: Error opening file %s\n", GMT_program, high);
				exit (EXIT_FAILURE);
			}
			if (hdr.nx != nx || hdr.ny != ny) {
				mexPrintf ("%s: upper limit file not of proper dimension!\n", GMT_program);
				exit (EXIT_FAILURE);
			}
			if (GMT_read_grd (high, &hdr, upper, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE)) {
				mexPrintf ("%s: Error reading file %s\n", GMT_program, high);
				exit (EXIT_FAILURE);
			}
/* Comment this out:	n_trimmed = 0;
			for (i = 0; i < nx * ny; i++) if (upper[i] < high_limit) {
				upper[i] = high_limit;
				n_trimmed++;
			}
			if (n_trimmed) fprintf (stderr, "%s: %d upper limit values < max data, reset to max data!\n", GMT_program, n_trimmed);
*/
		}
		for (j = ij = 0; j < ny; j++) {
			yy = ny - j - 1;
			for (i = 0; i < nx; i++, ij++) {
				if (GMT_is_fnan (upper[ij])) continue;
				upper[ij] -= (float)(plane_c0 + plane_c1 * i + plane_c2 * yy);
				upper[ij] *= (float)r_z_scale;
			}
		}
		constrained = TRUE;
	}
}

double	guess_surface_time(int nx, int ny)
{
	/* Routine to guess a number proportional to the operations
	 * required by surface working on a user-desired grid of
	 * size nx by ny, where nx = (x_max - x_min)/dx, and same for
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


void	suggest_sizes_for_surface(int nx, int ny)
{
	/* Calls guess_surface_time for a variety of trial grid
	 * sizes, where the trials are highly composite numbers
	 * with lots of factors of 2, 3, and 5.  The sizes are
	 * within the range (nx,ny) - (2*nx, 2*ny).  Prints to
	 * stderr the values which are an improvement over the
	 * user's original nx,ny.
	 * Should be called with nx=(x_max-x_min)/dx, and ditto
	 * for ny; that is, one smaller than the lattice used
	 * in surface.c
	 *
	 * W. H. F. Smith, 26 Feb 1992.  */

	double	guess_surface_time(int nx, int ny);
	double	users_time;	/* Time for user's nx, ny  */
	double	current_time;	/* Time for current nxg, nyg  */
	int	i;
	int	nxg, nyg;	/* Guessed by this routine  */
	int	nx2, ny2, nx3, ny3, nx5, ny5;	/* For powers  */
	int	xstop, ystop;	/* Set to 2*nx, 2*ny  */
	int	n_sug = 0;	/* N of suggestions found  */
	int	compare_sugs(const void *point_1, const void *point_2);	/* Sort suggestions decreasing  */
	struct SUGGESTION *sug = NULL;
	
	users_time = guess_surface_time(nx, ny);
	xstop = 2*nx;
	ystop = 2*ny;

	for (nx2 = 2; nx2 <= xstop; nx2 *= 2) {
	  for (nx3 = 1; nx3 <= xstop; nx3 *= 3) {
	    for (nx5 = 1; nx5 <= xstop; nx5 *= 5) {
		nxg = nx2 * nx3 * nx5;
		if (nxg < nx || nxg > xstop) continue;

		for (ny2 = 2; ny2 <= ystop; ny2 *= 2) {
		  for (ny3 = 1; ny3 <= ystop; ny3 *= 3) {
		    for (ny5 = 1; ny5 <= ystop; ny5 *= 5) {
			nyg = ny2 * ny3 * ny5;
			if (nyg < ny || nyg > ystop) continue;

			current_time = guess_surface_time(nxg, nyg);
			if (current_time < users_time) {
				n_sug++;
				sug = (struct SUGGESTION *)GMT_memory ((void *)sug, (size_t)n_sug, sizeof(struct SUGGESTION), GMT_program);
				sug[n_sug-1].nx = nxg;
				sug[n_sug-1].ny = nyg;
				sug[n_sug-1].factor = users_time/current_time;
			}

		    }
		  }
		}

	    }
	  }
	}

	if (n_sug) {
		qsort((void *)sug, (size_t)n_sug, sizeof(struct SUGGESTION), compare_sugs);
		for (i = 0; i < n_sug && i < 10; i++) {
			mexPrintf("%s:  HINT:  Choosing nx = %d, ny = %d might cut run time by a factor of %.8g\n",
				GMT_program, sug[i].nx, sug[i].ny, sug[i].factor);
		}
		GMT_free ((void *)sug);
	}
	else {
		mexPrintf("%s: Cannot suggest any nx,ny better than your -R -I define.\n", GMT_program);
	}
	return;
}

int	compare_sugs(const void *point_1, const void *point_2)
{
	/* Sorts sugs into DESCENDING order!  */
	if ( ((struct SUGGESTION *)point_1)->factor < ((struct SUGGESTION *)point_2)->factor)
		return(1);
	else if ( ((struct SUGGESTION *)point_1)->factor > ((struct SUGGESTION *)point_2)->factor)
		return(-1);
	else
		return(0);
}

int	get_prime_factors(int n, int *f)
{
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

int	gcd_euclid(int a, int b) {
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
