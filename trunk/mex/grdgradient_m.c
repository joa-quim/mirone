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
 * Note: This version differs slightly from the original in C. Namely -S
 * option makes the output contain |grad z| instead of the directional derivatives.
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
 */
 
#include "gmt.h"
#include "mex.h"

double specular(double nx, double ny, double nz, double *s);

BOOLEAN GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	i, j, ij, k, n, nm, nx, ny, i2, argc = 0, n_arg_no_char = 0;
	int p[4], n_used = 0, mx, nc_h, nr_h, *i_4, *o_i4, *pdata_i4, entry;
	short int *i_2, *o_i2, *pdata_i2;
	unsigned short int *ui_2, *o_ui2, *pdata_ui2;
	char	format[BUFSIZ], **argv;
	unsigned char *ui_1, *o_ui1, *pdata_ui1;
	
	BOOLEAN	error = FALSE, map_units = FALSE, normalize = FALSE, atan_trans = FALSE, bad, do_direct_deriv = FALSE;
	BOOLEAN find_directions = FALSE, do_cartesian = FALSE, do_orientations = FALSE, save_slopes = FALSE, add_ninety = FALSE;
	BOOLEAN lambertian_s = FALSE, peucker = FALSE, lambertian = FALSE;
	BOOLEAN sigma_set = FALSE, offset_set = FALSE, exp_trans = FALSE, two_azims = FALSE;
	BOOLEAN is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	BOOLEAN is_uint16 = FALSE, is_uint8 = FALSE;
	
	float	*data, *z_4, *pdata_s, *o_s;
	double	dx_grid, dy_grid, x_factor, y_factor, dzdx, dzdy, ave_gradient, norm_val = 1.0, sigma = 0.0;
	double	azim, denom, max_gradient = 0.0, min_gradient = 0.0, rpi, m_pr_degree, lat, azim2;
	double	x_factor2, y_factor2, dzdx2, dzdy2, dzds1, dzds2, offset;
	double	*pdata, *pdata_d, *z_8, *head, *o_d;
	double	p0, q0, elev, p0q0_cte;
	double	ka = 0.55, kd = 0.6, ks = 0.4, k_ads = 1.55, spread = 10., diffuse, spec;
	double	norm_z, mag, s[3], lim_x, lim_y, lim_z;
	float	r_min = DBL_MAX, r_max = -DBL_MAX, scale;
	char	input[BUFSIZ], *ptr;
	struct	GRD_HEADER header;
	struct	GMT_EDGEINFO edgeinfo;


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

	GMT_lock = FALSE;       /* Override since Matlab would own the lock */

	if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);
	
	GMT_boundcond_init (&edgeinfo);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				/* Supplemental parameters */
			
				case 'A':
					do_direct_deriv = TRUE;
					j = sscanf(&argv[i][2], "%lf/%lf", &azim, &azim2);
					two_azims = (j == 2);
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
							fprintf(stderr,"%s: GMT SYNTAX ERROR -Es option: Must give azim & elevation\t=%d\n", GMT_program, j);
							exit (EXIT_FAILURE);
						}
						p0 = cos((90 - azim)*D2R) * tan((90 - elev)*D2R);
						q0 = sin((90 - azim)*D2R) * tan((90 - elev)*D2R);
						p0q0_cte = sqrt(1 + p0*p0 + q0*q0);
						break;
					}
					lambertian = TRUE;	/* "full" Lambertian case */
					j = sscanf(&argv[i][2], "%lf/%lf", &azim, &elev);
					if (j < 2) {
						fprintf(stderr,"%s: GMT SYNTAX ERROR -E option: Must give at least azim & elevation\t=%d\n", GMT_program, j);
						exit (EXIT_FAILURE);
					}
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
		mexPrintf ( "'[-M]', '[-N[t_or_e][<amp>[/<sigma>[/<offset>]]]]', '[-S]')\n\n");
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
		mexPrintf ( "\t-S output |grad z| instead of directional derivatives; requires -D\n");
		return;
	}

	if (!(do_direct_deriv || find_directions || lambertian_s || lambertian || peucker)) {
		mexPrintf ("GMT SYNTAX ERROR:  Must specify -A or -D\n");
		error++;
	}
	if (save_slopes && !find_directions) {
		mexPrintf ("GMT SYNTAX ERROR -S option:  Must specify -D\n");
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
		mexPrintf ("%s: GMT SYNTAX ERROR -E option:  Use 0-360 degree range for azimuth\n", GMT_program);
		error++;
	}
	if ((lambertian || lambertian_s) && (elev < 0.0 || elev > 90.0)) {
		mexPrintf ("%s: GMT SYNTAX ERROR -E option:  Use 0-90 degree range for elevation\n", GMT_program);
		error++;
	}
	if ((lambertian || lambertian_s || peucker) && (do_direct_deriv || find_directions || save_slopes)) {
		mexPrintf ("%s: WARNING: -E option overrides -A, -D or -S\n", GMT_program);
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
	nm = header.nx * header.ny;

	data = mxCalloc ((nx+4)*(ny+4), sizeof (float));

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--)
			for (j = 0; j < nx; j++) data[i2*mx + j + 2*mx + 2] = (float)z_8[j*ny+i];
	}
	else if (is_single) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + 2*mx + 2] = z_4[j*ny+i];
	}
	else if (is_int32) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + 2*mx + 2] = (float)i_4[j*ny+i];
	}
	else if (is_int16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + 2*mx + 2] = (float)i_2[j*ny+i];
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + 2*mx + 2] = (float)ui_2[j*ny+i];
	}
	else if (is_uint8) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + 2*mx + 2] = (float)ui_1[j*ny+i];
	}

	GMT_boundcond_param_prep (&header, &edgeinfo);

	GMT_grd_init (&header, argc, argv, TRUE);
	GMT_pad[0] = GMT_pad[1] = GMT_pad[2] = GMT_pad[3] = 2;
	
	/* set boundary conditions:  */

	GMT_boundcond_set (&header, &edgeinfo, GMT_pad, data);

	if (map_units) {
		m_pr_degree = 2.0 * M_PI * gmtdefs.ref_ellipsoid[gmtdefs.ellipsoid].eq_radius / 360.0;
		dx_grid = m_pr_degree * header.x_inc * cosd ((header.y_max + header.y_min) / 2.0);
		dy_grid = m_pr_degree * header.y_inc;
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

	p[0] = 1;	p[1] = -1;	p[2] = mx;	p[3] = -mx;
	
	min_gradient = DBL_MAX;	max_gradient = -DBL_MAX;
	ave_gradient = 0.0;
	if (lambertian) {
		lim_x = header.x_max - header.x_min;
		lim_y = header.y_max - header.y_min;
		lim_z = header.z_max - header.z_min;
		scale = MAX(lim_z, MAX(lim_x, lim_y));
		lim_x /= scale;	lim_y /= scale;		lim_z /= scale;
		dx_grid /= lim_x;	dy_grid /= lim_y;
		x_factor = -dy_grid / (2 * lim_z);	y_factor = -dx_grid / (2 * lim_z);
	}
	for (j = k = 0; j < header.ny; j++) {
		if (map_units) {
			lat = (header.node_offset) ? -header.y_inc * (j + 0.5) : -header.y_inc * j;
			lat += header.y_max;
			dx_grid = m_pr_degree * header.x_inc * cos (D2R * lat);
			x_factor = -1.0 / (2.0 * dx_grid);
			if (do_direct_deriv) {
				if (two_azims) {
					x_factor2 = x_factor * sin(azim2);
				}
				x_factor *= sin(azim);
			}
		}
		for (i = 0; i < header.nx; i++, k++) {
			ij = (j + 2) * mx + i + 2;
			for (n = 0, bad = FALSE; !bad && n < 4; n++) if (GMT_is_fnan (data[ij+p[n]])) bad = TRUE;
			if (bad) {	/* One of corners = NaN, skip */
				data[k] = GMT_f_NaN;
				continue;
			}
			
			dzdx = (data[ij+1] - data[ij-1]) * x_factor;
			dzdy = (data[ij-mx] - data[ij+mx]) * y_factor;
			if (two_azims) {
				dzdx2 = (data[ij+1] - data[ij-1]) * x_factor2;
				dzdy2 = (data[ij-mx] - data[ij+mx]) * y_factor2;
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

	if (lambertian || lambertian_s || peucker) {	/* data must be scaled to the [-1,1] interval, but we'll do it into [-.95, .95] to not get too bright */
		scale = 1. / (r_max - r_min);
		for (k = 0; k < nm; k++) {
			if (GMT_is_fnan (data[k])) continue;
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
					for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) denom += pow(data[k] - ave_gradient, 2.0);
					denom = sqrt( (n_used - 1) / denom);
					sigma = 1.0 / denom;
				}
				rpi = 2.0 * norm_val / M_PI;
				for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) data[k] = (float)(rpi * atan((data[k] - ave_gradient)*denom));
				header.z_max = rpi * atan((max_gradient - ave_gradient)*denom);
				header.z_min = rpi * atan((min_gradient - ave_gradient)*denom);
			}
			else if (exp_trans) {
				if (!sigma_set) {
					sigma = 0.0;
					for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) sigma += fabs((double)data[k]);
					sigma = M_SQRT2 * sigma / n_used;
				}
				denom = M_SQRT2 / sigma;
				for (k = 0; k < nm; k++) {
					if (GMT_is_fnan (data[k])) continue;
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
				for (k = 0; k < nm; k++) if (!GMT_is_fnan (data[k])) data[k] = (float)((data[k] - ave_gradient) * denom);
				header.z_max = (max_gradient - ave_gradient) * denom;
				header.z_min = (min_gradient - ave_gradient) * denom;
			}
		}
	}

	/* Now we write out: */
	
	if (do_direct_deriv) {
		if (normalize) {
			strcpy (header.title, "Normalized directional derivative(s)");
		}
		else {
			strcpy (header.title, "Directional derivative(s)");
		}
		sprintf (format, "\t%s\t%s\t%s\t%s\n", gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
	}
	else {
		if (lambertian_s || lambertian)
			strcpy (header.title, "Lambertian radiance");
		else if (peucker)
			strcpy (header.title, "Peucker piecewise linear radiance");
		else
			strcpy (header.title, "Directions of maximum slopes");
	}


	/* Transpose from gmt grd orientation to Matlab orientation */
	/* Because we need to do the transposition and also a type conversion, we need a extra array */
	nx = header.nx;		ny = header.ny;

	if (is_double) {
		o_d = mxCalloc (nx*ny, sizeof (double));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_d[j*ny+ny-i-1] = (double)data[i*nx+j];
		plhs[0] = mxCreateDoubleMatrix (ny,nx, mxREAL);
		pdata_d = mxGetPr(plhs[0]);
		memcpy(pdata_d, o_d, ny*nx * 8);
	}
	else if (is_single) {
		//o_s = mxCalloc (nx*ny, sizeof (float));
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		pdata_s = mxGetData(plhs[0]);
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) pdata_s[j*ny+ny-i-1] = data[i*nx+j];
		//memcpy(pdata_s, o_s, ny*nx * 4);
	}
	else if (is_int32) {
		o_i4 = mxCalloc (nx*ny, sizeof (int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_i4[j*ny+ny-i-1] = irint(data[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT32_CLASS,mxREAL);
		pdata_i4 = mxGetData(plhs[0]);
		memcpy(pdata_i4, o_i4, ny*nx * 4);
	}
	else if (is_int16) {
		o_i2 = mxCalloc (nx*ny, sizeof (short int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_i2[j*ny+ny-i-1] = (short int)irint(data[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxINT16_CLASS,mxREAL);
		pdata_i2 = mxGetData(plhs[0]);
		memcpy(pdata_i2, o_i2, ny*nx * 2);
	}
	else if (is_uint16) {
		o_ui2 = mxCalloc (nx*ny, sizeof (short int));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_ui2[j*ny+ny-i-1] = (unsigned short int)irint(data[i*nx+j]);
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT16_CLASS,mxREAL);
		pdata_ui2 = mxGetData(plhs[0]);
		memcpy(pdata_ui2, o_ui2, ny*nx * 2);
	}
	else if (is_uint8) {
		o_ui1 = mxCalloc (nx*ny, sizeof (char));
		for (i = 0; i < ny; i++) for (j = 0; j < nx; j++) o_ui1[j*ny+ny-i-1] = (unsigned char)data[i*nx+j];
		plhs[0] = mxCreateNumericMatrix (ny,nx,mxUINT8_CLASS ,mxREAL);
		pdata_ui1 = mxGetData(plhs[0]);
		memcpy(pdata_ui1, o_ui1, ny*nx * 1);
	}
	mxFree(data);
	GMT_end_for_mex (argc, argv);

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
