/*--------------------------------------------------------------------
 *	$Id: grdgradient.c,v 1.9 2004/01/02 22:45:13 pwessel Exp $
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
 *  grdgradient.c
 * read a grdfile and compute gradient in azim direction:
 *
 * azim = azimuth clockwise from north in degrees.
 *
 * gradient = -[(dz/dx)sin(azim) + (dz/dy)cos(azim)].
 *
 * the expression in [] is the correct gradient.  We take
 * -[]  in order that data which goes DOWNHILL in the
 * azim direction will give a positive value; this is
 * for image shading purposes.
 *
 *
 * Author:	W.H.F. Smith
 * Date: 	13 Feb 1991
 * Upgraded to v2.0 15-May-1991 Paul Wessel
 *
 * Modified:	1 Mar 94 by WHFS to make -M scale change with j latitude
 *		1 Mar 96 by PW to find gradient direction and magnitude (-S and -D)
 *		13 Mar 96 by WHFS to add exp trans and user-supplied sigma to -N
 *			option, and add optional second azimuth to -A option.
 *		11 Sep 97 by PW now may pass average gradient along with sigma in -N
 *		22 Apr 98 by WHFS to add boundary conditions, switch sense of -S and 
 *			-D, and switch -Da to -Dc, for consistency of args.
 *		6  Sep 05 by J. Luis, added a -E option that allows the Lambertian or
 *			Peucker piecewise linear radiance computations
 * Version:	4
 *--------------------------------------------------------------------*/
/*
 * Mexified version of grdgradient
 * Author:	J. Luis
 * Date: 	17 Aug 2004
 * Modified	12 Jan 2006 - Extended -E (lambertian) option
 *
 * Usage
 * Zout = grdgradient_m(Zin,head,'options');
 *
 * where	Zin is the array containing the input file
 *			head is the header descriptor in the format of grdread/grdwrite mexs
 *			and options may be for example: '-A0',-M','-Ne0.6'
 *
 * The form [Zout,offset,sigma] = grdgradient_m(Zin,head,'options');
 * is also possible and is quite usefull (for Mirone) for it allows ROI image reconstructions
 * whith -N[t][e]1/sigma/offset
 *
 * IMPORTANT NOTE. The data type of Zin is preserved in Zout. That means you can send Zin
 * as a double, single, Int32, Int16, Uint16 or Uint8 and receive Zout in one of those types
 *	 
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *	 
 *		12/07/08 J Luis, Made it standalone (no GMT lib dependency)
 *	 
 *		11/10/11 J Luis, Now works in column major order, which means we can do things faster but
 *		                 above all we save a lot of memory by not making one copy of input array.
 *
 *		12/02/12 J Luis, Add -S<p|d> to compute slope in percentage or degrees. Also -A now multiplies
 *		                 input array by -1, instead of relying that it was done by the caller.
 *
 *		22/02/12 J Luis, The above blind multiplication by -1 if -A was an error. Illumination doesn't
 *		                 want that, but directional deriv does. So add a new -z option to rescale by -1
 *
 *		25/02/12 J Luis, A side effect of the change to column major order was that eastward illuminations
 *		                 and vice-versa. Making the change ang = 360 - ang; solves the issue.
 */
 
#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>

#if HAVE_OPENMP
#include <omp.h>
#endif

#define GMT_SMALL		1.0e-4	/* Needed when results aren't exactly zero but close */

#define	FALSE	0
#define	TRUE	1
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define D2R (M_PI / 180.0)
#define R2D (180.0 / M_PI)
#define cosd(x) cos ((x) * D2R)
#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))

#ifndef M_SQRT2
#define	M_SQRT2		1.41421356237309504880
#endif

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

#define EQ_RAD 6371.0087714
#define M_PR_DEG (EQ_RAD * 1000 * M_PI / 180.0)

/* For floats ONLY */
#define ISNAN_F(x) (((*(int32_T *)&(x) & 0x7f800000L) == 0x7f800000L) && \
                    ((*(int32_T *)&(x) & 0x007fffffL) != 0x00000000L))

struct GRD_HEADER {
/* Do not change the first three items. They are copied verbatim to the native grid header */
	int nx;				/* Number of columns */
	int ny;				/* Number of rows */
	int node_offset;		/* 0 for node grids, 1 for pixel grids */
	/* This section is flexible. It is not copied to any grid header */
	int type;				/* Grid format */
	char name[256];			/* Actual name of the file after any ?<varname> and =<stuff> has been removed */
	char varname[80];		/* NetCDF: variable name */
	int y_order;			/* NetCDF: 1 if S->N, -1 if N->S */
	int z_id;				/* NetCDF: id of z field */
	int ncid;				/* NetCDF: file ID */
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
	double z_scale_factor;	/* grd values must be multiplied by this */
	double z_add_offset;	/* After scaling, add this */
	char x_units[80];		/* units in x-direction */
	char y_units[80];		/* units in y-direction */
	char z_units[80];		/* grid value units */
	char title[80];			/* name of data set */
	char command[320];		/* name of generating command */
	char remark[160];		/* comments re this data set */
}; 

struct GMT_EDGEINFO {
	/* Description below is the final outcome after parse and verify */
	int	nxp;	/* if X periodic, nxp > 0 is the period in pixels  */
	int	nyp;	/* if Y periodic, nxp > 0 is the period in pixels  */
	int	gn;	/* TRUE if top    edge will be set as N pole  */
	int	gs;	/* TRUE if bottom edge will be set as S pole  */
};

double specular(double nx, double ny, double nz, double *s);
void GMT_boundcond_init (struct GMT_EDGEINFO *edgeinfo);
int GMT_boundcond_set (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo, int *pad, float *a);
int GMT_boundcond_param_prep (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo);
int GMT_boundcond_parse (struct GMT_EDGEINFO *edgeinfo, char *edgestring);
void angreloc(double *a);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	i, j, ij, k, n, nm, nx, ny, i2, k1, k2, argc = 0, n_arg_no_char = 0;
	int p[4], n_used = 0, mx, my, nc_h, nr_h, *i_4, *pdata_i4, entry, GMT_pad[4];
	short int *i_2, *pdata_i2;
	unsigned short int *ui_2, *pdata_ui2;
	char	**argv;
	unsigned char *ui_1, *o_ui1, *pdata_ui1;

	int	error = FALSE, map_units = FALSE, normalize = FALSE, atan_trans = FALSE, bad, do_direct_deriv = FALSE;
	int	find_directions = FALSE, do_cartesian = FALSE, do_orientations = FALSE, save_slopes = FALSE;
	int slope_percent = FALSE, slope_deg = FALSE, add_ninety = FALSE, do_change_Zsign = FALSE;
	int	lambertian_s = FALSE, peucker = FALSE, lambertian = FALSE, unknown_nans = TRUE, check_nans = FALSE;
	int	sigma_set = FALSE, offset_set = FALSE, exp_trans = FALSE, two_azims = FALSE;
	int	is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int	is_uint16 = FALSE, is_uint8 = FALSE;
	clock_t tic;
	
	float	*data, *z_4, *pdata_s;
	float	nan = mxGetNaN();
	double	dzdx = 0.0, dzdy = 0.0, ave_gradient = 0., norm_val = 1.0, sigma = 0.0;
	double	azim = 0.0, denom, max_gradient = 0.0, min_gradient = 0.0, rpi, m_pr_degree, lat, azim2;
	double	x_factor2, y_factor2, dzdx2, dzdy2, dzds1, dzds2, offset;
	double	*pdata, *pdata_d, *z_8, *head;
	double	p0, q0, elev, p0q0_cte;
	double	ka = 0.55, kd = 0.6, ks = 0.4, k_ads = 1.55, spread = 10., diffuse, spec;
	double	norm_z, mag, s[3], lim_x, lim_y, lim_z;
	double	dx_grid, dy_grid, x_factor, y_factor, *dx_grid__, *x_factor__, *x_factor2__;
	float	r_min = FLT_MAX, r_max = -FLT_MAX, scale;
	char	input[BUFSIZ], *ptr;
	struct	GRD_HEADER header;
	struct	GMT_EDGEINFO edgeinfo;

#ifdef MIR_TIMEIT
	tic = clock();
#endif

	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}
	argc++;			/* to account for the program's name to be inserted in argv[0] */

	/* get the length of the input string */
	argv=(char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "grdgradient";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	GMT_boundcond_init (&edgeinfo);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				/* Supplemental parameters */
			
				case 'a':	/* NaN value. Use 1 to not change the illuminated NaNs color */
					sscanf(&argv[i][2], "%f", &nan);
					break;
				case 'A':
					do_direct_deriv = TRUE;
					j = sscanf(&argv[i][2], "%lf/%lf", &azim, &azim2);
					two_azims = (j == 2);
					angreloc(&azim);		/* Adjust because of change to column major order */
					if (two_azims) angreloc(&azim2);
					break;
				case 'D':
					find_directions = TRUE;
					j = 2;
					while (argv[i][j]) {
						switch (argv[i][j]) {
							case 'C':
							case 'c':
								do_cartesian = TRUE;
								break;
							case 'O':
							case 'o':
								do_orientations = TRUE;
								break;
							case 'N':
							case 'n':
								add_ninety = TRUE;
								break;
							default:
								mexPrintf("GMT SYNTAX ERROR -S option:  Unrecognized modifier\n");
								error++;
								break;
						}
						j++;
					}
					break;
				case 'E':
					if (argv[i][2] == 'p') {
						peucker = TRUE;
						break;
					}
					else if (argv[i][2] == 's') {	/* "simple" Lmbertian case */
						lambertian_s = TRUE;
						j = sscanf(&argv[i][3], "%lf/%lf", &azim, &elev);
						if (j != 2) {
							mexPrintf("GRDGRADIENT SYNTAX ERROR -Es option: Must give azim & elevation\t=%d\n", j);
							return;
						}
						angreloc(&azim);
						p0 = cos((90 - azim)*D2R) * tan((90 - elev)*D2R);
						q0 = sin((90 - azim)*D2R) * tan((90 - elev)*D2R);
						p0q0_cte = sqrt(1 + p0*p0 + q0*q0);
						break;
					}
					lambertian = TRUE;	/* "full" Lambertian case */
					j = sscanf(&argv[i][2], "%lf/%lf", &azim, &elev);
					if (j < 2) {
						mexPrintf("GRDGRADIENT: SYNTAX ERROR -E option: Must give at least azim & elevation\t=%d\n", j);
						return;
					}
					angreloc(&azim);	/* Adjust because of change to column major order */
					while (azim < 0) azim += 360;
					while (azim > 360) azim -= 360;
					elev = 90 - elev;
					s[0] = sin(azim*D2R) * cos(elev*D2R);
					s[1] = cos(azim*D2R) * cos(elev*D2R);
					s[2] =  sin(elev*D2R);
					strcpy (input, &argv[i][2]);
					ptr = (char *)strtok (input, "/");
					entry = 0;
					while (ptr) {
						switch (entry) {
							case 0:
							case 1:
								break;	/* Cases already processed above */
							case 2:
								if (ptr[0] != '=') ka = atof (ptr);
								break;
							case 3:
								if (ptr[0] != '=') kd = atof (ptr);
								break;
							case 4:
								if (ptr[0] != '=') ks = atof (ptr);
								break;
							case 5:
								if (ptr[0] != '=') spread = atof (ptr);
								break;
							default:
								break;
						}
						ptr = (char *) strtok (NULL, "/");
						entry++;
					}
					k_ads = ka + kd + ks;
					break;
				case 'L':
					error += GMT_boundcond_parse (&edgeinfo, &argv[i][2]);
					break;
				case 'M':
					map_units = TRUE;
					break;
				case 'N':
					normalize = TRUE;
					j = 2;
					if (argv[i][j]) {
						if (argv[i][j] == 't' || argv[i][j] == 'T') {
							atan_trans = TRUE;
							j++;
						}
						else if (argv[i][j] == 'e' || argv[i][j] == 'E') {
							exp_trans = TRUE;
							j++;
						}
						j = sscanf(&argv[i][j], "%lf/%lf/%lf", &norm_val, &sigma, &offset);
						if (j >= 2) sigma_set = TRUE;
						if (j == 3) offset_set = TRUE;
					}
					break;

				case 'S':
					save_slopes = TRUE;
					if (argv[i][2]) {
						if (argv[i][2] == 'p')
							slope_percent = TRUE;
						else if (argv[i][2] == 'd')
							slope_deg = TRUE;
					}
					break;
				case 'z':
					do_change_Zsign = TRUE;		/* Multiply Z by -1 */
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (argc == 1 || error) {
		mexPrintf ("grdgradient - Compute directional gradients from grdfiles\n\n");
		mexPrintf ( "usage: R = grdgradient_m(infile,head,'[-A<azim>[/<azim2>]]', '[-D[a][o][n]]', '[-L<flag>]',\n");
		mexPrintf ( "'[-M]', '[-N[t_or_e][<amp>[/<sigma>[/<offset>]]]]', '[-S<p|d>]', '[-a<nan_val>]')\n\n");
		mexPrintf ("\t<infile> is name of input array\n");
		mexPrintf ("\t<head> is array header descriptor of the form\n");
		mexPrintf ("\t [x_min x_max y_min y_max z_min zmax 0 x_inc y_inc]\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ( "\t-A sets azimuth (0-360 CW from North (+y)) for directional derivatives\n");
		mexPrintf ( "\t  -A<azim>/<azim2> will compute two directions and save the one larger in magnitude.\n");
		mexPrintf ( "\t-D finds the direction of grad z.\n");
		mexPrintf ( "\t   Append c to get cartesian angle (0-360 CCW from East (+x)) [Default:  azimuth]\n");
		mexPrintf ( "\t   Append o to get bidirectional orientations [0-180] rather than directions [0-360]\n");
		mexPrintf ( "\t   Append n to add 90 degrees to the values from c or o\n");
		mexPrintf ( "\t-E Compute Lambertian radiance appropriate to use with grdimage/grdview.\n");
		mexPrintf ( "\t   -E<azim/elev> sets azimuth and elevation of light vector.\n");
		mexPrintf ( "\t   -E<azim/elev/ambient/diffuse/specular/shine> sets azim, elev and\n");
		mexPrintf ( "\t    other parameters that control the reflectance properties of the surface.\n");
		mexPrintf ( "\t    Default values are: 0.55/0.6/0.4/10\n");
		mexPrintf ( "\t    Specify '=' to get the default value (e.g. -E60/30/=/0.5)\n");
		mexPrintf ( "\t   Append s to use a simpler Lambertian algorithm (note that with this form\n");
		mexPrintf ( "\t   you only have to provide the azimuth and elevation parameters)\n");
		mexPrintf ( "\t   Append p to use the Peucker picewise linear aproximation (simpler but faster algorithm)\n");
		mexPrintf ( "\t   Note that in this case the azimuth and elevation are hardwired to 315 and 45 degrees.\n");
		mexPrintf ( "\t   This means that even if you provide other values they will be ignored.\n");
		mexPrintf ( "\t-L sets boundary conditions.  <flag> can be either\n");
		mexPrintf ( "\t   g for geographic boundary conditions\n");
		mexPrintf ( "\t   or one or both of\n");
		mexPrintf ( "\t   x for periodic boundary conditions on x\n");
		mexPrintf ( "\t   y for periodic boundary conditions on y\n");
		mexPrintf ( "\t   [Default:  Natural conditions]\n");
		mexPrintf ( "\t-M to use map units.  In this case, dx,dy of grdfile\n");
		mexPrintf ( "\t   will be converted from degrees lon,lat into meters.\n");
		mexPrintf ( "\t   Default computes gradient in units of data/grid_distance.\n");
		mexPrintf ( "\t-N will normalize gradients so that max |grad| = <amp> [1.0]\n");
		mexPrintf ( "\t  -Nt will make atan transform, then scale to <amp> [1.0]\n");
		mexPrintf ( "\t  -Ne will make exp  transform, then scale to <amp> [1.0]\n");
		mexPrintf ( "\t  -Nt<amp>/<sigma>[/<offset>] or -Ne<amp>/<sigma>[/<offset>] sets sigma\n");
		mexPrintf ( "\t     (and offset) for transform. [sigma, offset estimated from data]\n");
		mexPrintf ( "\t-S output |grad z| (slope)\n");
		mexPrintf ( "\t     append 'p' or 'd' to compute slope in percentage or degrees, respectively.\n");
		mexPrintf ( "\t-a NaN value. Use 1 to not change the illuminated NaNs color [default is NaN]\n");
		return;
	}

	if (!(do_direct_deriv || find_directions || lambertian_s || lambertian || peucker)) {
		mexPrintf ("GMT SYNTAX ERROR:  Must specify -A or -D\n");
		error++;
	}
	if (do_direct_deriv && (azim < 0.0 || azim >= 360.0)) {
		mexPrintf ("GMT SYNTAX ERROR -A option:  Use 0-360 degree range\n");
		error++;
	}
	if (two_azims && (azim2 < 0.0 || azim2 >= 360.0)) {
		mexPrintf ("GMT SYNTAX ERROR -A option:  Use 0-360 degree range\n");
		error++;
	}
	if (norm_val <= 0.0) {
		mexPrintf ("GMT SYNTAX ERROR -N option:  Normalization amplitude must be > 0\n");
		error++;
	}
	if (sigma_set && (sigma <= 0.0) ) {
		mexPrintf ("GMT SYNTAX ERROR -N option:  Sigma must be > 0\n");
		error++;
	}
	if ((lambertian || lambertian_s) && (azim < 0.0 || azim >= 360.0)) {
		mexPrintf ("GRDGRADIENT_M GMT SYNTAX ERROR -E option:  Use 0-360 degree range for azimuth\n");
		error++;
	}
	if ((lambertian || lambertian_s) && (elev < 0.0 || elev > 90.0)) {
		mexPrintf ("GRDGRADIENT_M GMT SYNTAX ERROR -E option:  Use 0-90 degree range for elevation\n");
		error++;
	}
	if ((lambertian || lambertian_s || peucker) && (do_direct_deriv || find_directions || save_slopes)) {
		mexPrintf ("GRDGRADIENT_M WARNING: -E option overrides -A, -D or -S\n");
		do_direct_deriv = find_directions = save_slopes = FALSE;
	}

	if (error) return;

	/* Get non char inputs */
	if (nlhs == 0)
		mexErrMsgTxt("GRDGRADIENT ERROR: Must provide an output.\n");

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
		mexPrintf("GRDGRADIENT ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, In32, In16, UInt16 and Uint8.\n");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2) {
		mexErrMsgTxt("GRDGRADIENT: First non char argument must contain a decent array\n");
	}

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("GRDGRADIENT: Second argument must contain a valid header of the input array\n");
	
	head  = mxGetPr(prhs[1]);		/* Get header info */
	header.x_min = head[0];	header.x_max = head[1];
	header.y_min = head[2];	header.y_max = head[3];
	header.z_min = head[4];	header.z_max = head[5];
	header.x_inc = head[7];	header.y_inc = head[8];
	header.nx = nx;			header.ny = ny;
	header.node_offset = irint(head[6]);
	mx = nx + 4;
	my = ny + 4;
	nm = header.nx * header.ny;

	data = (float *)mxMalloc ((nx+4)*(ny+4) * sizeof (float));

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = i2 = 0; i < nx; i++) {
			k1 = my * (i + 2) + 2;
			for (j = 0; j < ny; j++)
				data[j + k1] = (float)z_8[i2++];
		}
	}
	else if (is_single) {
		for (i = i2 = 0; i < nx; i++) {
			k1 = my * (i + 2) + 2;
			for (j = 0; j < ny; j++)
				data[j + k1] = z_4[i2++];
		}
	}
	else if (is_int32) {
		for (i = i2 = 0; i < nx; i++) {
			k1 = my * (i + 2) + 2;
			for (j = 0; j < ny; j++)
				data[j + k1] = (float)i_4[i2++];
		}
	}
	else if (is_int16) {
		for (i = i2 = 0; i < nx; i++) {
			k1 = my * (i + 2) + 2;
			for (j = 0; j < ny; j++)
				data[j + k1] = (float)i_2[i2++];
		}
	}
	else if (is_uint16) {
		for (i = i2 = 0; i < nx; i++) {
			k1 = my * (i + 2) + 2;
			for (j = 0; j < ny; j++)
				data[j + k1] = (float)ui_2[i2++];
		}
	}
	else if (is_uint8) {
		for (i = i2 = 0; i < nx; i++) {
			k1 = my * (i + 2) + 2;
			for (j = 0; j < ny; j++)
				data[j + k1] = (float)ui_1[i2++];
		}
	}

	if ((is_double || is_single) && unknown_nans) {
		/* Loop over the file and find if we have NaNs. Stop at first NaN occurence. */
		for (i = 2*(my+1); i < mx*(my-2); i++) {
			if (ISNAN_F(data[i])) {
				check_nans = TRUE;
				break;
			}
		}
	}

	/* For direction derivatives we need to revert the sign of Z */
	if (do_change_Zsign)
		for (i = 0; i < (nx+4)*(ny+4); i++)
			data[i] *= -1;

	GMT_boundcond_param_prep (&header, &edgeinfo);

	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
	
	/* set boundary conditions:  */

	/* We must cheat GMT_boundcond_set so it 'thinks' is dealing with a row major array */
	header.nx = ny;      header.ny = nx;
	GMT_boundcond_set (&header, &edgeinfo, GMT_pad, data);
	header.nx = nx;      header.ny = ny;

	/* The dx_grid, x_factor and x_factor2 are computed again later in "if(map_units)" 
	   but it makes the code simpler if we do it this little extra computations
	*/
	if (map_units) {
		/*m_pr_degree = 2.0 * M_PI * gmtdefs.ref_ellipsoid[gmtdefs.ellipsoid].eq_radius / 360.0;*/
		dx_grid = M_PR_DEG * header.x_inc * cosd ((header.y_max + header.y_min) / 2.0);
		dy_grid = M_PR_DEG * header.y_inc;
	}
	else {
		dx_grid = header.x_inc;
		dy_grid = header.y_inc;
	}

	x_factor = -1.0 / (2.0 * dx_grid);
	y_factor = -1.0 / (2.0 * dy_grid);
	if (do_direct_deriv) {
		if (two_azims) {
			azim2 *= (M_PI / 180.0);
			x_factor2 = x_factor * sin(azim2);
			y_factor2 = y_factor * cos(azim2);
		}
		azim *= (M_PI / 180.0);
		x_factor *= sin(azim);
		y_factor *= cos(azim);
	}

	p[0] = 1;	p[1] = -1;	p[2] = my;	p[3] = -my;
	
	min_gradient = DBL_MAX;	max_gradient = -DBL_MAX;

	if (map_units) {
		dx_grid__ = (double *)mxMalloc(header.ny * sizeof(double));
		x_factor__ = (double *)mxMalloc(header.ny * sizeof(double));
		if (do_direct_deriv && two_azims)
			x_factor2__ = (double *)mxMalloc(header.ny * sizeof(double));

		for (j = 0; j < header.ny; j++) {
			lat = (header.node_offset) ? -header.y_inc * (j + 0.5) : -header.y_inc * j;
			lat += header.y_max;
			dx_grid__[j] = M_PR_DEG * header.x_inc * cos (D2R * lat);
			x_factor__[j] = -1.0 / (2.0 * dx_grid__[j]);
			if (do_direct_deriv) {
				if (two_azims)
					x_factor2__[j] = x_factor__[j] * sin(azim2);
				x_factor__[j] *= sin(azim);
			}
		}
	}

	if (lambertian) {
		lim_x = header.x_max - header.x_min;
		lim_y = header.y_max - header.y_min;
		lim_z = header.z_max - header.z_min;
		scale = MAX(lim_z, MAX(lim_x, lim_y));
		lim_x /= scale;
		lim_y /= scale;
		lim_z /= scale;
		dx_grid /= lim_x;
		dy_grid /= lim_y;
		x_factor = -dy_grid / (2 * lim_z);
		y_factor = -dx_grid / (2 * lim_z);
	}

#if HAVE_OPENMP
#pragma omp parallel for  private(i, j, ij, k, n, dzdx, dzdy, dzds1, dzds2, ave_gradient, min_gradient, max_gradient, x_factor, dx_grid) shared(data)
#endif

	for (i = k = 0; i < header.nx; i++) {
		ij = (i + 2) * my + 2;
		for (j = 0; j < header.ny; j++, k++, ij++) {
			if (check_nans) {
				for (n = 0, bad = FALSE; !bad && n < 4; n++) 
					if (ISNAN_F(data[ij+p[n]])) bad = TRUE;
				if (bad) {	/* One of corners = NaN, skip */
					data[k] = nan;
					continue;
				}
			}

			if (map_units) {
				dx_grid = dx_grid__[j];
				x_factor = x_factor__[j];
				if (do_direct_deriv && two_azims)
					x_factor2 = x_factor2__[j];
			}

			dzdy = (data[ij+1] - data[ij-1]) * y_factor;
			if (azim)		/* Otherwise x_factor is zero and dzdx was initialized to zero */
				dzdx = (data[ij-my] - data[ij+my]) * x_factor;
			if (two_azims) {
				dzdy2 = (data[ij+1] - data[ij-1]) * y_factor2;
				dzdx2 = (data[ij-my] - data[ij+my]) * x_factor2;
			}

			/* Write output to unused NW corner */

			if (do_direct_deriv) {	/* Directional derivatives */
				if (two_azims) {
					dzds1 = dzdx + dzdy;
					dzds2 = dzdx2 + dzdy2;
					data[k] = (float)((fabs(dzds1) > fabs(dzds2)) ? dzds1 : dzds2);
				}
				else {
					data[k] = (float)(dzdx + dzdy);
				}
				ave_gradient += data[k];
				min_gradient = MIN (min_gradient, data[k]);
				max_gradient = MAX (max_gradient, data[k]);
			}
			else if (find_directions) {
				azim = (do_cartesian) ? atan2 (-dzdy, -dzdx) * R2D : 90.0 - atan2 (-dzdy, -dzdx) * R2D;
				if (add_ninety) azim += 90.0;
				if (azim < 0.0) azim += 360.0;
				if (azim >= 360.0) azim -= 360.0;
				if (do_orientations && azim >= 180) azim -= 180.0;
				if (!save_slopes)
					data[k] = (float)azim;
				else
					data[k] = (float)hypot (dzdx, dzdy);
			}
			else {
				if (lambertian) {
					norm_z = dx_grid * dy_grid;
					mag = d_sqrt(dzdx*dzdx + dzdy*dzdy + norm_z*norm_z);
					dzdx /= mag;	dzdy /= mag;	norm_z /= mag;
					diffuse = MAX(0,(s[0]*dzdx + s[1]*dzdy + s[2]*norm_z)); 
					spec = specular(dzdx, dzdy, norm_z, s);
					spec = pow(spec, spread);
					data[k] = (float)((ka+kd*diffuse+ks*spec) / k_ads);
				}
				else if (lambertian_s)
					data[k] = (float)( (1 + p0*dzdx + q0*dzdy) / (sqrt(1 + dzdx*dzdx + dzdy*dzdy) * p0q0_cte) );
				else	/* Peucker method */
					data[k] = (float)( -0.4285 * (dzdx - dzdy) - 0.0844 * fabs(dzdx  + dzdy) + 0.6599 );
				r_min = MIN (r_min, data[k]);
				r_max = MAX (r_max, data[k]);
			}
			n_used++;
		}
	}

	if (slope_percent) {
		for (i = 0; i < nm; i++)
			data[i] *= 100;
	}
	else if (slope_deg) {
		for (i = 0; i < nm; i++)
			data[i] = atan(data[i]) * R2D;
	}

	if (map_units) {
		mxFree(dx_grid__);
		mxFree(x_factor__);
		if (do_direct_deriv && two_azims)
			mxFree(x_factor2__);
	}

	if (lambertian || lambertian_s || peucker) {	/* data must be scaled to the [-1,1] interval, but we'll do it into [-.95, .95] to not get too bright */
		scale = (float)(1. / (r_max - r_min));
		for (k = 0; k < nm; k++) {
			if (check_nans && ISNAN_F (data[k])) continue;
			data[k] = (-1. + 2. * ((data[k] - r_min) * scale)) * 0.95;
		}
	}
	
	if (offset_set)
		ave_gradient = offset;
	else
		ave_gradient /= n_used;

	if (do_direct_deriv) {	/* Report some statistics */
	
		if (normalize) {
			if (atan_trans) {
				if (sigma_set) {
					denom = 1.0 / sigma;
				}
				else {
					denom = 0.0;
					if (!check_nans) {
						for (k = 0; k < nm; k++) denom += pow(data[k] - ave_gradient, 2.0);
					}
					else {
						for (k = 0; k < nm; k++) 
							if (!ISNAN_F(data[k])) denom += pow(data[k] - ave_gradient, 2.0);
					}
					denom = sqrt( (n_used - 1) / denom);
					sigma = 1.0 / denom;
				}
				rpi = 2.0 * norm_val / M_PI;
				if (!check_nans) {
					for (k = 0; k < nm; k++) data[k] = (float)(rpi * atan((data[k] - ave_gradient)*denom));
				}
				else {
					for (k = 0; k < nm; k++) 
						if (!ISNAN_F (data[k])) data[k] = (float)(rpi * atan((data[k] - ave_gradient)*denom));
				}
				header.z_max = rpi * atan((max_gradient - ave_gradient)*denom);
				header.z_min = rpi * atan((min_gradient - ave_gradient)*denom);
			}
			else if (exp_trans) {
				if (!sigma_set) {
					sigma = 0.0;
					for (k = 0; k < nm; k++) 
						if (!ISNAN_F (data[k])) sigma += fabs((double)data[k]);
					sigma = M_SQRT2 * sigma / n_used;
				}
				denom = M_SQRT2 / sigma;
				for (k = 0; k < nm; k++) {
					if (check_nans && ISNAN_F (data[k])) continue;
					if (data[k] < ave_gradient) {
						data[k] = (float)(-norm_val * (1.0 - exp((data[k] - ave_gradient)*denom)));
					}
					else {
						data[k] = (float)(norm_val * (1.0 - exp(-(data[k] - ave_gradient)*denom)));
					}
				}
				header.z_max = norm_val * (1.0 - exp(-(max_gradient - ave_gradient)*denom));
				header.z_min = -norm_val * (1.0 - exp((min_gradient - ave_gradient)*denom));
			}
			else {
				if ( (max_gradient - ave_gradient) > (ave_gradient - min_gradient) ) {
					denom = norm_val / (max_gradient - ave_gradient);
				}
				else {
					denom = norm_val / (ave_gradient - min_gradient);
				}
				if (!check_nans) {
					for (k = 0; k < nm; k++) data[k] = (float)((data[k] - ave_gradient) * denom);
				}
				else {
					for (k = 0; k < nm; k++) 
						if (!ISNAN_F (data[k])) data[k] = (float)((data[k] - ave_gradient) * denom);
				}
				header.z_max = (max_gradient - ave_gradient) * denom;
				header.z_min = (min_gradient - ave_gradient) * denom;
			}
		}
	}

	/* Transpose from gmt grd orientation to Matlab orientation */
	/* Because we need to do the transposition and also a type conversion, we need a extra array */
	nx = header.nx;		ny = header.ny;

	if (is_double) {
		plhs[0] = mxCreateDoubleMatrix (0,0, mxREAL);
		mxRealloc((void *)data, nx*ny * sizeof(double));
	}
	else if (is_single) {
		plhs[0] = mxCreateNumericMatrix (0,0,mxSINGLE_CLASS,mxREAL);
		mxRealloc((void *)data, nx*ny * sizeof(float));
	}
	else if (is_int32) {
		plhs[0] = mxCreateNumericMatrix (0,0,mxINT32_CLASS,mxREAL);
		mxRealloc((void *)data, nx*ny * sizeof(float));
	}
	else if (is_int16) {
		plhs[0] = mxCreateNumericMatrix (0,0,mxINT16_CLASS,mxREAL);
		mxRealloc((void *)data, nx*ny * sizeof(short int));
	}
	else if (is_uint16) {
		plhs[0] = mxCreateNumericMatrix (0,0,mxUINT16_CLASS,mxREAL);
		mxRealloc((void *)data, nx*ny * sizeof(short int));
	}
	else if (is_uint8) {
		plhs[0] = mxCreateNumericMatrix (0,0,mxUINT8_CLASS ,mxREAL);
		mxRealloc((void *)data, nx*ny * sizeof(char));
	}

	if (is_double)
		mxSetPr(plhs[0], (void *)data);
	else
		mxSetData(plhs[0], (void *)data);

	mxSetM(plhs[0], ny);
	mxSetN(plhs[0], nx);

	if (nlhs == 2) {
		plhs[1] = mxCreateDoubleMatrix (1,1, mxREAL);
		pdata = mxGetPr(plhs[1]);
		memcpy(pdata, &ave_gradient, 8);
	}
	else if (nlhs == 3) {
		plhs[1] = mxCreateDoubleMatrix (1,1, mxREAL);
		pdata = mxGetPr(plhs[1]);
		memcpy(pdata, &ave_gradient, 8);

		plhs[2] = mxCreateDoubleMatrix (1,1, mxREAL);
		pdata = mxGetPr(plhs[2]);
		if (atan_trans || exp_trans)
			memcpy(pdata, &sigma, 8);
		else {
			sigma = 0;
			memcpy(pdata, &sigma, 8);
		}

	}

#ifdef MIR_TIMEIT
	mexPrintf("GRDGRADIENT_M: CPU ticks = %.3f\tCPS = %d\n", (double)(clock() - tic), CLOCKS_PER_SEC);
#endif

}

void angreloc(double *a) {
	/* Adjustment because of change to column major order.
	   Without this eastward illuminations became westwards and vice-versa. */

	if (*a != 0)
		*a = 360 - *a;
}

double specular(double nx, double ny, double nz, double *s) {
	/* SPECULAR Specular reflectance.
	   R = SPECULAR(Nx,Ny,Nz,S,V) returns the reflectance of a surface with
	   normal vector components [Nx,Ny,Nz].  S and V specify the direction
	   to the light source and to the viewer, respectively. 
	   For the time beeing I'm using V = [azim elev] = [0 90] so the following

	   V[0] =  sin(V[0]*D2R)*cos(V[1]*D2R);
	   V[1] = -cos(V[0]*D2R)*cos(V[1]*D2R);
	   V[2] =  sin(V[1]*D2R);

	   Reduces to V[0] = 0;		V[1] = 0;	V[2] = 1 */

	/*r = MAX(0,2*(s[0]*nx+s[1]*ny+s[2]*nz).*(v[0]*nx+v[1]*ny+v[2]*nz) - (v'*s)*ones(m,n)); */

	return (MAX(0, 2 * (s[0]*nx + s[1]*ny + s[2]*nz) * nz - s[2]));
}


void GMT_boundcond_init (struct GMT_EDGEINFO *edgeinfo) {
	edgeinfo->nxp = 0;
	edgeinfo->nyp = 0;
	edgeinfo->gn = FALSE;
	edgeinfo->gs = FALSE;
	return;
}

int GMT_boundcond_parse (struct GMT_EDGEINFO *edgeinfo, char *edgestring) {
	/* Parse string beginning at argv[i][2] and load user's
		requests in edgeinfo->  Return success or failure.
		Requires that edgeinfo previously initialized to
		zero/FALSE stuff.  Expects g or (x and or y) is
		all that is in string.  */

	int	i, ier;

	i = 0;
	ier = FALSE;
	while (!ier && edgestring[i]) {
		switch (edgestring[i]) {
			case 'g':
			case 'G':
				edgeinfo->gn = TRUE;
				edgeinfo->gs = TRUE;
				break;
			case 'x':
			case 'X':
				edgeinfo->nxp = -1;
				break;
			case 'y':
			case 'Y':
				edgeinfo->nyp = -1;
				break;
			default:
				ier = TRUE;
				break;

		}
		i++;
	}

	if (ier) return (-1);

	return (0);
}

int GMT_boundcond_param_prep (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo) {
	/* Called when edgeinfo holds user's choices.  Sets edgeinfo according to choices and h.  */

	double	xtest;

	/*if (edgeinfo->gn || GMT_grd_is_global(h)) {	I'm not using the is_global test */
	if (edgeinfo->gn) {
		/* User has requested geographical conditions.  */
		if ( (h->x_max - h->x_min) < (360.0 - GMT_SMALL * h->x_inc) ) {
			mexPrintf ("GRDGRADIENT_M Warning: x range too small; g boundary condition ignored.\n");
			edgeinfo->nxp = edgeinfo->nyp = 0;
			edgeinfo->gn  = edgeinfo->gs = FALSE;
			return (0);
		}
		xtest = fmod (180.0, h->x_inc) / h->x_inc;
		/* xtest should be within GMT_SMALL of zero or of one.  */
		if ( xtest > GMT_SMALL && xtest < (1.0 - GMT_SMALL) ) {
			/* Error.  We need it to divide into 180 so we can phase-shift at poles.  */
			mexPrintf ("GRDGRADIENT_M Warning: x_inc does not divide 180; g boundary condition ignored.\n");
			edgeinfo->nxp = edgeinfo->nyp = 0;
			edgeinfo->gn  = edgeinfo->gs = FALSE;
			return (0);
		}
		edgeinfo->nxp = irint(360.0/h->x_inc);
		edgeinfo->nyp = 0;
		edgeinfo->gn = ( (fabs(h->y_max - 90.0) ) < (GMT_SMALL * h->y_inc) );
		edgeinfo->gs = ( (fabs(h->y_min + 90.0) ) < (GMT_SMALL * h->y_inc) );
	}
	else {
		if (edgeinfo->nxp != 0) edgeinfo->nxp = (h->node_offset) ? h->nx : h->nx - 1;
		if (edgeinfo->nyp != 0) edgeinfo->nyp = (h->node_offset) ? h->ny : h->ny - 1;
	}
	return (0);
}

int GMT_boundcond_set (struct GRD_HEADER *h, struct GMT_EDGEINFO *edgeinfo, int *pad, float *a) {
	/* Set two rows of padding (pad[] can be larger) around data according
		to desired boundary condition info in edgeinfo.
		Returns -1 on problem, 0 on success.
		If either x or y is periodic, the padding is entirely set.
		However, if neither is true (this rules out geographical also)
		then all but three corner-most points in each corner are set.

		As written, not ready to use with "surface" for GMT v4, because
		assumes left/right is +/- 1 and down/up is +/- mx.  In "surface"
		the amount to move depends on the current mesh size, a parameter
		not used here.

		This is the revised, two-rows version (WHFS 6 May 1998).
	*/

	int	bok;	/* Counter used to test that things are OK  */
	int	mx;	/* Width of padded array; width as malloc'ed  */
	int	mxnyp;	/* distance to periodic constraint in j direction  */
	int	i, jmx;	/* Current i, j * mx  */
	int	nxp2;	/* 1/2 the xg period (180 degrees) in cells  */
	int	i180;	/* index to 180 degree phase shift  */
	int	iw, iwo1, iwo2, iwi1, ie, ieo1, ieo2, iei1;  /* see below  */
	int	jn, jno1, jno2, jni1, js, jso1, jso2, jsi1;  /* see below  */
	int	jno1k, jno2k, jso1k, jso2k, iwo1k, iwo2k, ieo1k, ieo2k;
	int	j1p, j2p;	/* j_o1 and j_o2 pole constraint rows  */


	/* Check pad  */
	bok = 0;
	for (i = 0; i < 4; i++) {
		if (pad[i] < 2) bok++;
	}
	if (bok > 0) {
		mexPrintf ("GMT BUG:  bad pad for GMT_boundcond_set.\n");
		return (-1);
	}

	/* Check minimum size:  */
	if (h->nx < 2 || h->ny < 2) {
		mexPrintf ("GMT ERROR:  GMT_boundcond_set requires nx,ny at least 2.\n");
		return (-1);
	}

	/* Initialize stuff:  */

	mx = h->nx + pad[0] + pad[1];
	nxp2 = edgeinfo->nxp / 2;	/* Used for 180 phase shift at poles  */

	iw = pad[0];		/* i for west-most data column */
	iwo1 = iw - 1;		/* 1st column outside west  */
	iwo2 = iwo1 - 1;	/* 2nd column outside west  */
	iwi1 = iw + 1;		/* 1st column  inside west  */

	ie = pad[0] + h->nx - 1;	/* i for east-most data column */
	ieo1 = ie + 1;		/* 1st column outside east  */
	ieo2 = ieo1 + 1;	/* 2nd column outside east  */
	iei1 = ie - 1;		/* 1st column  inside east  */

	jn = mx * pad[3];	/* j*mx for north-most data row  */
	jno1 = jn - mx;		/* 1st row outside north  */
	jno2 = jno1 - mx;	/* 2nd row outside north  */
	jni1 = jn + mx;		/* 1st row  inside north  */

	js = mx * (pad[3] + h->ny - 1);	/* j*mx for south-most data row  */
	jso1 = js + mx;		/* 1st row outside south  */
	jso2 = jso1 + mx;	/* 2nd row outside south  */
	jsi1 = js - mx;		/* 1st row  inside south  */

	mxnyp = mx * edgeinfo->nyp;

	jno1k = jno1 + mxnyp;	/* data rows periodic to boundary rows  */
	jno2k = jno2 + mxnyp;
	jso1k = jso1 - mxnyp;
	jso2k = jso2 - mxnyp;

	iwo1k = iwo1 + edgeinfo->nxp;	/* data cols periodic to bndry cols  */
	iwo2k = iwo2 + edgeinfo->nxp;
	ieo1k = ieo1 - edgeinfo->nxp;
	ieo2k = ieo2 - edgeinfo->nxp;

	/* Check poles for grid case.  It would be nice to have done this
		in GMT_boundcond_param_prep() but at that point the data
		array isn't passed into that routine, and may not have been
		read yet.  Also, as coded here, this bombs with error if
		the pole data is wrong.  But there could be an option to
		to change the condition to Natural in that case, with warning.  */

	if (h->node_offset == 0) {
		if (edgeinfo->gn) {
			bok = 0;
			if (ISNAN_F (a[jn + iw])) {
				for (i = iw+1; i <= ie; i++) {
					if (!ISNAN_F (a[jn + i])) bok++;
				}
			}
			else {
				for (i = iw+1; i <= ie; i++) {
					if (a[jn + i] != a[jn + iw]) bok++;
				}
			}
			if (bok > 0) mexPrintf ("GRDGRADIENT_M Warning: Inconsistent grid values at North pole.\n");
		}

		if (edgeinfo->gs) {
			bok = 0;
			if (ISNAN_F (a[js + iw])) {
				for (i = iw+1; i <= ie; i++)
					if (!ISNAN_F (a[js + i])) bok++;
			}
			else {
				for (i = iw+1; i <= ie; i++)
					if (a[js + i] != a[js + iw]) bok++;
			}
			if (bok > 0) mexPrintf ("GRDGRADIENT_M Warning: Inconsistent grid values at South pole.\n");
		}
	}

	/* Start with the case that x is not periodic, because in that
		case we also know that y cannot be polar.  */

	if (edgeinfo->nxp <= 0) {

		/* x is not periodic  */

		if (edgeinfo->nyp > 0) {

			/* y is periodic  */

			for (i = iw; i <= ie; i++) {
				a[jno1 + i] = a[jno1k + i];
				a[jno2 + i] = a[jno2k + i];
				a[jso1 + i] = a[jso1k + i];
				a[jso2 + i] = a[jso2k + i];
			}

			/* periodic Y rows copied.  Now do X naturals.
				This is easy since y's are done; no corner problems.
				Begin with Laplacian = 0, and include 1st outside rows
				in loop, since y's already loaded to 2nd outside.  */

			for (jmx = jno1; jmx <= jso1; jmx += mx) {
				a[jmx + iwo1] = (float)(4.0 * a[jmx + iw])
					- (a[jmx + iw + mx] + a[jmx + iw - mx] + a[jmx + iwi1]);
				a[jmx + ieo1] = (float)(4.0 * a[jmx + ie])
					- (a[jmx + ie + mx] + a[jmx + ie - mx] + a[jmx + iei1]);
			}

			/* Copy that result to 2nd outside row using periodicity.  */
			a[jno2 + iwo1] = a[jno2k + iwo1];
			a[jso2 + iwo1] = a[jso2k + iwo1];
			a[jno2 + ieo1] = a[jno2k + ieo1];
			a[jso2 + ieo1] = a[jso2k + ieo1];

			/* Now set d[laplacian]/dx = 0 on 2nd outside column.  Include
				1st outside rows in loop.  */
			for (jmx = jno1; jmx <= jso1; jmx += mx) {
				a[jmx + iwo2] = (a[jmx + iw - mx] + a[jmx + iw + mx] + a[jmx + iwi1])
					- (a[jmx + iwo1 - mx] + a[jmx + iwo1 + mx])
					+ (float)(5.0 * (a[jmx + iwo1] - a[jmx + iw]));

				a[jmx + ieo2] = (a[jmx + ie - mx] + a[jmx + ie + mx] + a[jmx + iei1])
					- (a[jmx + ieo1 - mx] + a[jmx + ieo1 + mx])
					+ (float)(5.0 * (a[jmx + ieo1] - a[jmx + ie]));
			}

			/* Now copy that result also, for complete periodicity's sake  */
			a[jno2 + iwo2] = a[jno2k + iwo2];
			a[jso2 + iwo2] = a[jso2k + iwo2];
			a[jno2 + ieo2] = a[jno2k + ieo2];
			a[jso2 + ieo2] = a[jso2k + ieo2];

			/* DONE with X not periodic, Y periodic case.  Fully loaded.  */

			return (0);
		}
		else {
			/* Here begins the X not periodic, Y not periodic case  */

			/* First, set corner points.  Need not merely Laplacian(f) = 0
				but explicitly that d2f/dx2 = 0 and d2f/dy2 = 0.
				Also set d2f/dxdy = 0.  Then can set remaining points.  */

	/* d2/dx2 */	a[jn + iwo1]   = (float)(2.0 * a[jn + iw] - a[jn + iwi1]);
	/* d2/dy2 */	a[jno1 + iw]   = (float)(2.0 * a[jn + iw] - a[jni1 + iw]);
	/* d2/dxdy */	a[jno1 + iwo1] = -(a[jno1 + iwi1] - a[jni1 + iwi1] + a[jni1 + iwo1]);


	/* d2/dx2 */	a[jn + ieo1]   = (float)(2.0 * a[jn + ie] - a[jn + iei1]);
	/* d2/dy2 */	a[jno1 + ie]   = (float)(2.0 * a[jn + ie] - a[jni1 + ie]);
	/* d2/dxdy */	a[jno1 + ieo1] = -(a[jno1 + iei1] - a[jni1 + iei1] + a[jni1 + ieo1]);

	/* d2/dx2 */	a[js + iwo1]   = (float)(2.0 * a[js + iw] - a[js + iwi1]);
	/* d2/dy2 */	a[jso1 + iw]   = (float)(2.0 * a[js + iw] - a[jsi1 + iw]);
	/* d2/dxdy */	a[jso1 + iwo1] = -(a[jso1 + iwi1] - a[jsi1 + iwi1] + a[jsi1 + iwo1]);

	/* d2/dx2 */	a[js + ieo1]   = (float)(2.0 * a[js + ie] - a[js + iei1]);
	/* d2/dy2 */	a[jso1 + ie]   = (float)(2.0 * a[js + ie] - a[jsi1 + ie]);
	/* d2/dxdy */	a[jso1 + ieo1] = -(a[jso1 + iei1] - a[jsi1 + iei1] + a[jsi1 + ieo1]);

			/* Now set Laplacian = 0 on interior edge points,
				skipping corners:  */
			for (i = iwi1; i <= iei1; i++) {
				a[jno1 + i] = (float)(4.0 * a[jn + i])
					- (a[jn + i - 1] + a[jn + i + 1]
						+ a[jni1 + i]);

				a[jso1 + i] = (float)(4.0 * a[js + i])
					- (a[js + i - 1] + a[js + i + 1]
						+ a[jsi1 + i]);
			}
			for (jmx = jni1; jmx <= jsi1; jmx += mx) {
				a[iwo1 + jmx] = (float)(4.0 * a[iw + jmx])
					- (a[iw + jmx + mx] + a[iw + jmx - mx]
						+ a[iwi1 + jmx]);
				a[ieo1 + jmx] = (float)(4.0 * a[ie + jmx])
					- (a[ie + jmx + mx] + a[ie + jmx - mx]
						+ a[iei1 + jmx]);
			}

			/* Now set d[Laplacian]/dn = 0 on all edge pts, including
				corners, since the points needed in this are now set.  */
			for (i = iw; i <= ie; i++) {
				a[jno2 + i] = a[jni1 + i]
					+ (float)(5.0 * (a[jno1 + i] - a[jn + i]))
					+ (a[jn + i - 1] - a[jno1 + i - 1])
					+ (a[jn + i + 1] - a[jno1 + i + 1]);
				a[jso2 + i] = a[jsi1 + i]
					+ (float)(5.0 * (a[jso1 + i] - a[js + i]))
					+ (a[js + i - 1] - a[jso1 + i - 1])
					+ (a[js + i + 1] - a[jso1 + i + 1]);
			}
			for (jmx = jn; jmx <= js; jmx += mx) {
				a[iwo2 + jmx] = a[iwi1 + jmx]
					+ (float)(5.0 * (a[iwo1 + jmx] - a[iw + jmx]))
					+ (a[iw + jmx - mx] - a[iwo1 + jmx - mx])
					+ (a[iw + jmx + mx] - a[iwo1 + jmx + mx]);
				a[ieo2 + jmx] = a[iei1 + jmx]
					+ (float)(5.0 * (a[ieo1 + jmx] - a[ie + jmx]))
					+ (a[ie + jmx - mx] - a[ieo1 + jmx - mx])
					+ (a[ie + jmx + mx] - a[ieo1 + jmx + mx]);
			}
			/* DONE with X not periodic, Y not periodic case.
				Loaded all but three cornermost points at each corner.  */

			return (0);
		}
		/* DONE with all X not periodic cases  */
	}
	else {
		/* X is periodic.  Load x cols first, then do Y cases.  */

		for (jmx = jn; jmx <= js; jmx += mx) {
			a[iwo1 + jmx] = a[iwo1k + jmx];
			a[iwo2 + jmx] = a[iwo2k + jmx];
			a[ieo1 + jmx] = a[ieo1k + jmx];
			a[ieo2 + jmx] = a[ieo2k + jmx];
		}

		if (edgeinfo->nyp > 0) {
			/* Y is periodic.  copy all, including boundary cols:  */
			for (i = iwo2; i <= ieo2; i++) {
				a[jno1 + i] = a[jno1k + i];
				a[jno2 + i] = a[jno2k + i];
				a[jso1 + i] = a[jso1k + i];
				a[jso2 + i] = a[jso2k + i];
			}
			/* DONE with X and Y both periodic.  Fully loaded.  */

			return (0);
		}

		/* Do north (top) boundary:  */

		if (edgeinfo->gn) {
			/* Y is at north pole.  Phase-shift all, incl. bndry cols. */
			if (h->node_offset) {
				j1p = jn;	/* constraint for jno1  */
				j2p = jni1;	/* constraint for jno2  */
			}
			else {
				j1p = jni1;	/* constraint for jno1  */
				j2p = jni1 + mx;	/* constraint for jno2  */
			}
			for (i = iwo2; i <= ieo2; i++) {
				i180 = pad[0] + ((i + nxp2)%edgeinfo->nxp);
				a[jno1 + i] = a[j1p + i180];
				a[jno2 + i] = a[j2p + i180];
			}
		}
		else {
			/* Y needs natural conditions.  x bndry cols periodic.
				First do Laplacian.  Start/end loop 1 col outside,
				then use periodicity to set 2nd col outside.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jno1 + i] = (float)(4.0 * a[jn + i])
					- (a[jn + i - 1] + a[jn + i + 1] + a[jni1 + i]);
			}
			a[jno1 + iwo2] = a[jno1 + iwo2 + edgeinfo->nxp];
			a[jno1 + ieo2] = a[jno1 + ieo2 - edgeinfo->nxp];


			/* Now set d[Laplacian]/dn = 0, start/end loop 1 col out,
				use periodicity to set 2nd out col after loop.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jno2 + i] = a[jni1 + i]
					+ (float)(5.0 * (a[jno1 + i] - a[jn + i]))
					+ (a[jn + i - 1] - a[jno1 + i - 1])
					+ (a[jn + i + 1] - a[jno1 + i + 1]);
			}
			a[jno2 + iwo2] = a[jno2 + iwo2 + edgeinfo->nxp];
			a[jno2 + ieo2] = a[jno2 + ieo2 - edgeinfo->nxp];

			/* End of X is periodic, north (top) is Natural.  */

		}

		/* Done with north (top) BC in X is periodic case.  Do south (bottom)  */

		if (edgeinfo->gs) {
			/* Y is at south pole.  Phase-shift all, incl. bndry cols. */
			if (h->node_offset) {
				j1p = js;	/* constraint for jso1  */
				j2p = jsi1;	/* constraint for jso2  */
			}
			else {
				j1p = jsi1;	/* constraint for jso1  */
				j2p = jsi1 - mx;	/* constraint for jso2  */
			}
			for (i = iwo2; i <= ieo2; i++) {
				i180 = pad[0] + ((i + nxp2)%edgeinfo->nxp);
				a[jso1 + i] = a[j1p + i180];
				a[jso2 + i] = a[j2p + i180];
			}
		}
		else {
			/* Y needs natural conditions.  x bndry cols periodic.
				First do Laplacian.  Start/end loop 1 col outside,
				then use periodicity to set 2nd col outside.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jso1 + i] = (float)(4.0 * a[js + i])
					- (a[js + i - 1] + a[js + i + 1] + a[jsi1 + i]);
			}
			a[jso1 + iwo2] = a[jso1 + iwo2 + edgeinfo->nxp];
			a[jso1 + ieo2] = a[jso1 + ieo2 - edgeinfo->nxp];


			/* Now set d[Laplacian]/dn = 0, start/end loop 1 col out,
				use periodicity to set 2nd out col after loop.  */

			for (i = iwo1; i <= ieo1; i++) {
				a[jso2 + i] = a[jsi1 + i]
					+ (float)(5.0 * (a[jso1 + i] - a[js + i]))
					+ (a[js + i - 1] - a[jso1 + i - 1])
					+ (a[js + i + 1] - a[jso1 + i + 1]);
			}
			a[jso2 + iwo2] = a[jso2 + iwo2 + edgeinfo->nxp];
			a[jso2 + ieo2] = a[jso2 + ieo2 - edgeinfo->nxp];

			/* End of X is periodic, south (bottom) is Natural.  */

		}

		/* Done with X is periodic cases.  */

		return (0);
	}
}
