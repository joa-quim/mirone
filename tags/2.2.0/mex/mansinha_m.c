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

/*
 * 
 * Calcula la o que o programa do masinha cacula
 * 
 * Author:	J. Luis
 * Date:	24 April, 2002
 * Revision:	00/00/00 
 */

#include "mex.h"

#include <ctype.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define	FALSE		0
#define	TRUE		1
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif
#define D2R		M_PI / 180.
#define GMT_CONV_LIMIT  1.0e-8  /* Fairly tight convergence limit or "close to zero" limit */ 

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

/* Safe math macros that check arguments */

#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))
#define d_atn(y,x) ((x) == 0.0 && (y) == 0.0 ? 0.2 : atan2 (y, x))

/* In non-Windows this is may not be necessary (or guive conflicts) */
#define copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

struct srf_header {		/* Surfer file header structure */
	char id[4];		/* ASCII Binary identifier (DSAA/DSBB) */
	short int nx;		/* Number of columns */
	short int ny;		/* Number of rows */
	double x_min;		/* Minimum x coordinate */
	double x_max;		/* Maximum x coordinate */
	double y_min;		/* Minimum y coordinate */
	double y_max;		/* Maximum y coordinate */
	double z_min;		/* Minimum z value */
	double z_max;		/* Maximum z value */
};


double uscal(double x1, double x2, double x3, double c, double cc, double dp);
double udcal(double x1, double x2, double x3, double c, double cc, double dp);
void deform (double x_min, double y_min, int i_end, int j_end, float *z, double dx, double dy, double fault_length, double fault_width, double th, double dip, double rake, double d, double top_depth, double xl, double yl);
void meda (double x_min, double y_min, int i_end, int j_end, float *z, double dx, double dy, double side_x, double side_y, double height, double xl, double yl);
void tm (double lon, double lat, double *x, double *y);
void vtm (double lon0, double lat0);
int GMT_getinc (char *line, double *dx, double *dy);
int GMT_getincn (char *line, double inc[], int n);

int	is_geog = FALSE;
double	fault_length = 0.0;		/*  */
double	fault_width = 0.0;		/*  */
double	top_depth = 0.0;	/* Top Fault depth */
double	dip = 0.0;		/* Dip angle */
double	th = 0.0;		/* Strike direction */
double	rake = 0.0;		/* Rake angle */
double	d = 0.0;		/* Dislocation */
double	x_epic = 0.0;		/* x_epicenter coord */
double	y_epic = 0.0;		/* y_epicenter coord */
double	EQ_RAD = 6378137.0;	/* WGS-84 */
double	flattening = 1.0/298.2572235630;
double	ECC2, ECC4, ECC6;
double	one_m_ECC2, i_one_m_ECC2;
double	t_c1, t_c2, t_c3, t_c4, t_e2, t_M0;
double	central_meridian;

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	error = FALSE, is_meda = FALSE;
	int	i, i_end, j_end;
	int	argc = 0, n_arg_no_char = 0, nx, ny;
	float	*z, *pdata;		/* Array with bottom deformation points */
	double  x_min = 0, x_max = 0, y_min = 0, y_max = 0, x_inc = 0, y_inc = 0;
	double	height = 100, side_x = 20, side_y = 20;
	char	**argv;
	
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
	argv[0] = "mansinha";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

        for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
 		
				/* Common parameters */
			
				case '\0':
					error = TRUE;
					break;
				case 'A': /* Fault's angle parameters */
					sscanf (&argv[i][2], "%lf/%lf/%lf/%lf", &dip, &th, &rake, &d);
					break;
				case 'E': /* epicenter */
					sscanf (&argv[i][2], "%lf/%lf", &x_epic, &y_epic);
					break;
				case 'F': /* Fault size  & depth */
					sscanf (&argv[i][2], "%lf/%lf/%lf", &fault_length, &fault_width, &top_depth);
					break;
				case 'I': /* Grid step */
					GMT_getinc (&argv[i][2], &x_inc, &y_inc);
					break;
				case 'M': /* Geographical coords */
					is_geog = TRUE;
					break;
				case 'P': /* Punctual source */
					is_meda = TRUE;
					sscanf (&argv[i][2], "%lf/%lf/%lf", &side_x, &side_y, &height);
					break;
				case 'R':
					sscanf (&argv[i][2], "%lf/%lf/%lf/%lf", &x_min, &x_max, &y_min, &y_max);
					break;
                                default:
                                        error = TRUE;
                                        break;
                        }
                }
        }

	if (argc == 1 || error) {	/* Display usage */
		mexPrintf("mansinha - usage:\n Z = mansinha('-A<dip/azimuth/rake/slip>', '-E<x_epic/y_epic>', '-F<fault_length/fault_width/top_depth>', '-I<grid_inc>', '-R<xmin/xmax/ymin/ymax>', '[-M]', '[-P[x][/y][/h]]')\n\n");
 		
		mexPrintf ("\t-R<limits> Grid coordinate limits where deformation is computed.\n");
		mexPrintf ("\t-A<fault_params> Fault parameters describing Dip/Azimuth/Rake/Slip(m)\n");
		mexPrintf ("\t of focal mechanism. Azimuth is counted cw from North.\n");
		mexPrintf ("\t-E<epic_location> X and Y epicenter coordinates\n");
		mexPrintf ("\t-F<length/height/depth> Fault lenght, height and depth from sea-bottom.\n");
		mexPrintf ("\t They must all be given in km.\n");
		mexPrintf ("\t-I sets the grid spacing for the grid. Append m for minutes, c for seconds.\n");
		mexPrintf ("\t Note that x_max and y_max may be slightly adjusted\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-M create grid in geographical coordinates\n");
		mexPrintf ("\t-P create a exponential source of width x,y in km and heigh h in m (default is 20/20/100)\n");
		mexPrintf ("\n\tEXAMPLE:\n");
		mexPrintf ("\tTo create a geographical defformation grid covering a region between\n");
		mexPrintf ("\t20-5 W; 30-40 N with a 1 minute increment and a pure thrust focal mecanism\n");
		mexPrintf ("\twith a slip of 10 meters of a fault oriented 45N, run:\n");
		mexPrintf ("\n");
		mexPrintf ("\tmansinha -Gteste.grd -R-20/-5/30/40 -I1m -A60/45/90/10 -F60/25/0 -E-15/35 -M\n");
		mexErrMsgTxt("\n");
        }

	if (x_min == 0 && x_max == 0 && y_min == 0 && y_max == 0) {
		mexPrintf ("MANSINHA SYNTAX ERROR: Must specify -R option\n");
		error++;
	}
	if (dip == 0 && th == 0 && rake == 0 && d == 0) {
		mexPrintf ("MANSINHA SYNTAX ERROR: Must specify fault angle parameters (-A option)\n");
		mexPrintf ("                       OR use the -P option\n");
		error++;
	}
	if (fault_length == 0 && fault_width == 0 && top_depth == 0) {
		mexPrintf ("MANSINHA SYNTAX ERROR: Must specify fault size parameters (-F option)\n");
		mexPrintf ("                       OR use the -P option\n");
		error++;
	}
	if (x_inc <= 0 && y_inc <= 0) {
		mexPrintf ("MANSINHA SYNTAX ERROR -I option.  Must specify positive increment(s)\n");
		error++;
	}
	if (x_epic == 0 && y_epic == 0) {
		mexPrintf ("MANSINHA SYNTAX ERROR: Must specify epicenter location (-E option)\n");
		error++;
	}

	if (error) mexErrMsgTxt("\n");

	if (nlhs != 1)
		mexErrMsgTxt("SURFACE ERROR: Must provide one output.\n");

	/* Convert fault dimensions to meters (that's what is used by deform) */
	fault_length *= 1000;
	fault_width *= 1000;
	top_depth *= 1000;

	x_max = x_min + (irint ((x_max - x_min) / x_inc)) * x_inc;
	y_max = y_min + (irint ((y_max - y_min) / y_inc)) * y_inc;

	/* Compute i_end and j_end (that is, last row and last column) */
	i_end = irint ((y_max - y_min) / x_inc) + 1;
	j_end = irint ((x_max - x_min) / y_inc) + 1;

	/* Initialize TM variables. Fault origin will be used as projection's origin. However,
	   this would set it as a singularity point. That's whay it is arbitrarely shifted
	   by a 1/4 of grid step. */ 
	if (is_geog) vtm(x_epic+x_inc/2, y_epic+y_inc/2);

	nx = j_end;	ny = i_end;
	plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
	z = mxGetData(plhs[0]);

	if (is_meda) {	/* Compute a punctual (exponential) source */
		side_x *= 1000;		side_y *= 1000;
		meda (x_min, y_min, i_end, j_end, z, x_inc, y_inc, side_x, side_x, height, x_epic, y_epic);
	}
	else		/* Compute a elastic source */
		deform (x_min, y_min, i_end, j_end, z, x_inc, y_inc, fault_length, fault_width, th, dip, rake, d, top_depth, x_epic, y_epic);

}

void meda (double x_min, double y_min, int i_end, int j_end, float *z, double dx, double dy, double side_x, double side_y, double height, double xl, double yl) {
	/* Build a source with an exponential shape */
	double	c, xx, yy, rx, ry, rx2, ry2, tmp; 
	int	i, j, k = 0;
	c = cos((y_min + y_min + i_end*dx)/2 * D2R);
	side_x = 1. / side_x;
	side_y = 1. / side_y;
	for(i = 0; i < i_end; i++) {
		yy = y_min + dy * i;
		for(j = 0; j < j_end; j++) {
			xx = x_min + dx * j;
			if (is_geog) {
				rx = (xx - xl) * 111317.0 * c;
				ry = (yy - yl) * 111317.0;
			}
			else {
				rx = xx - xl;
				ry = yy - yl;
			}
			rx *= side_x;	ry *= side_y;
			rx2 = rx*rx;	ry2 = ry*ry;
			tmp = height * exp(-(rx2) - (ry2));
			z[k++] = (tmp > 3e-3) ? (float)(tmp) : 0;
		}
	}
}

void deform (double x_min, double y_min, int i_end, int j_end, float *z, double dx, double dy, double fault_length, double fault_width, double th, double dip, double rake, double d, double top_depth, double xl, double yl) {

/*	fault_length - comprimento da falha (m)
	fault_width  - largura da falha (m)
	th - azimute					strike angle
	dip - inclinaç do plano falha em relacao à horiz	dip angle
	rake - deslocação oblíqua da falha		rake or slip angle 
	d - escorregamento (m)				dislocation
	top_depth - profundidade do topo (m)		Top fault depth
*/

	double h1, h2, ds, dd, xx, yy, x1, x2, x3, us, ud, sn_tmp, cs_tmp, tg_tmp;
	double f1, f2, f3, f4, g1, g2, g3, g4, rx, ry;
	int i, j, k = 0;

	h1 = top_depth / sin(D2R * dip);
	h2 = top_depth / sin(D2R * dip) + fault_width;
	ds = -d * cos(D2R * rake);
	dd = d * sin(D2R * rake);
	sn_tmp = sin(D2R*th);	cs_tmp = cos(D2R*th);	tg_tmp = tan(D2R*dip);

	for(j = 0; j < j_end; j++) {
		xx = x_min + dx * j;
		for(i = 0; i < i_end; i++) {
			yy = y_min + dy * i;
			if (is_geog)
				tm(xx, yy, &rx, &ry);	/* Remember that (xl,yl) is already the proj origin */
			else {
				rx = xx - xl;
				ry = yy - yl;
			}
			x1 = rx*sn_tmp + ry*cs_tmp - fault_length/2.0;
			x2 = rx*cs_tmp - ry*sn_tmp + top_depth/tg_tmp;
			x3 = 0.0;
			f1 = uscal(x1, x2, x3, fault_length/2., h2, D2R*dip);
			f2 = uscal(x1, x2, x3, fault_length/2., h1, D2R*dip);
			f3 = uscal(x1, x2, x3, -fault_length/2., h2, D2R*dip);
			f4 = uscal(x1, x2, x3, -fault_length/2., h1, D2R*dip);
			g1 = udcal(x1, x2, x3, fault_length/2., h2, D2R*dip);
			g2 = udcal(x1, x2, x3, fault_length/2., h1, D2R*dip);
			g3 = udcal(x1, x2, x3, -fault_length/2., h2, D2R*dip);
			g4 = udcal(x1, x2, x3, -fault_length/2., h1, D2R*dip);
			us = (f1-f2-f3+f4)*ds / (12.* M_PI);
			ud = (g1-g2-g3+g4)*dd / (12.* M_PI);
			z[k++] = (float)(us + ud);
		}
	}
}

double uscal(double x1, double x2, double x3, double c, double cc, double dp) {
/* Computation of the vertical displacement due to the STRIKE and SLIP component */
	double sn, cs, c1, c2, c3, r, q, r2, r3, q2, q3, h, k, a1, a2, a3, f;
	double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14;

	sn = sin(dp);	cs = cos(dp);
	c1 = c;		c2 = cc * cs;	c3 = cc * sn;
	r = d_sqrt((x1-c1)*(x1-c1) + (x2-c2)*(x2-c2) + (x3-c3)*(x3-c3));
	q = d_sqrt((x1-c1)*(x1-c1) + (x2-c2)*(x2-c2) + (x3+c3)*(x3+c3));
	r2 = x2*sn - x3*cs;	r3 = x2*cs + x3*sn;
	q2 = x2*sn + x3*cs;	q3 = -x2*cs + x3*sn;
	h = d_sqrt(q2*q2 + (q3+cc)*(q3+cc));
	k = d_sqrt(q2*q2 + (x1-c1)*(x1-c1));
	a1 = log(r+r3-cc);	a2 = log(q+q3+cc);	a3 = log(q+x3+c3);
	b1 = 1. + 3. * (tan(dp)*tan(dp));
	b2 = 3. * tan(dp) / cs;
	b3 = 2. * r2 * sn;
	b4 = q2 + x2 * sn;
	b5 = 2. * r2*r2 * cs;
	b6 = r * (r+r3-cc);
	b7 = 4. * q2 * x3 * sn*sn;
	b8 = 2. * (q2+x2*sn) * (x3+q3*sn);
	b9 = q * (q+q3+cc);
	b10 = 4. * q2 * x3 * sn;
	b11 = (x3+c3) - q3 * sn;
	b12 = 4. * q2*q2 * q3 * x3 * cs * sn;
	b13 = 2. * q + q3 + cc;
	b14 = pow(q,3) * pow((q+q3+cc),2);
	f = cs * (a1 + b1*a2 - b2*a3) + b3/r + 2.*sn*b4/q - b5/b6 + (b7-b8)/b9 + b10*b11/(pow(q,3)) - b12*b13/b14;

	return (f);
}


double udcal(double x1, double x2, double x3, double c, double cc, double dp) {
/* Computation of the vertical displacement due to the DIP SLIP component */
	double sn, cs, c1, c2, c3, r, q, r2, r3, q2, q3, h, k, a1, a2;
	double b1, b2, b3, d1, d2, d3, d4, d5, d6, t1, t2, t3, f;

	sn = sin(dp);	cs = cos(dp);
	c1 = c;		c2 = cc * cs;	c3 = cc * sn;
	r = d_sqrt((x1-c1)*(x1-c1) + (x2-c2)*(x2-c2) + (x3-c3)*(x3-c3));
	q = d_sqrt((x1-c1)*(x1-c1) + (x2-c2)*(x2-c2) + (x3+c3)*(x3+c3));
	r2 = x2*sn - x3*cs;	r3 = x2*cs + x3*sn;
	q2 = x2*sn + x3*cs;	q3 = -x2*cs + x3*sn;
	h = d_sqrt(q2*q2 + (q3+cc)*(q3+cc));
	k = d_sqrt(q2*q2 + (x1-c1)*(x1-c1));
	a1 = log(r+x1-c1);	a2 = log(q+x1-c1);
	b1 = q * (q+x1-c1);	b2 = r * (r+x1-c1);	b3 = q * (q+q3+cc);
	d1 = x1 - c1;		d2 = x2 - c2;		d3 = x3 - c3;
	d4 = x3 + c3;		d5 = r3 - cc;		d6 = q3 + cc;
	t1 = d_atn(d1*d2, (h+d4)*(q+h));
	t2 = d_atn(d1*d5, r2*r);
	t3 = d_atn(d1*d6, q2*q);
	f = sn * (d2*(2.*d3/b2 + 4.*d3/b1 - 4.*c3*x3*d4*(2.*q+d1)/(b1*b1*q)) - 6.*t1 + 3.*t2 - 6.*t3) + cs * (a1-a2 - 2.*(d3*d3)/b2 - 4.*(d4*d4 - c3*x3)/b1 - 4.*c3*x3*d4*d4*(2*q+x1-c1)/(b1*b1*q)) + 6.*x3*(cs*sn*(2.*d6/b1 + d1/b3) - q2*(sn*sn - cs*cs)/b1);

	return (f);
}

int	no_sys_mem (char *where, int n) {	
		mexPrintf ("Fatal Error: %s could not allocate memory, n = %d\n", where, n);
		return (-1);
}

void vtm (double lon0, double lat0) {
	/* Set up an TM projection (extract of GMT_vtm)*/
	double lat2, s2, c2;
	
	lat0 *= D2R;
	lat2 = 2.0 * lat0;
	s2 = sin(lat2);
	c2 = cos(lat2);
	ECC2 = 2 * flattening - flattening * flattening;
	ECC4 = ECC2 * ECC2;
	ECC6 = ECC2 * ECC4;
	one_m_ECC2 = 1.0 - ECC2;
	i_one_m_ECC2 = 1.0 / one_m_ECC2;
	t_c1 = 1.0 - (1.0/4.0) * ECC2 - (3.0/64.0) * ECC4 - (5.0/256.0) * ECC6;
	t_c2 = -((3.0/8.0) * ECC2 + (3.0/32.0) * ECC4 + (25.0/768.0) * ECC6);
	t_c3 = (15.0/128.0) * ECC4 + (45.0/512.0) * ECC6;
	t_c4 = -(35.0/768.0) * ECC6;
	t_e2 = ECC2 * i_one_m_ECC2;
	t_M0 = EQ_RAD * (t_c1 * lat0 + s2 * (t_c2 + c2 * (t_c3 + c2 * t_c4)));
	central_meridian = lon0;
}

void tm (double lon, double lat, double *x, double *y) {
	/* Convert lon/lat to TM x/y (adapted from GMT_tm) */
	double N, T, T2, C, A, M, dlon, tan_lat, A2, A3, A5, lat2, s, c, s2, c2;

	if (fabs (fabs (lat) - 90.0) < GMT_CONV_LIMIT) {
		M = EQ_RAD * t_c1 * M_PI_2;
		*x = 0.0;
		*y = M;
	}
	else {
		lat *= D2R;
		lat2 = 2.0 * lat;
		s = sin(lat);	s2 = sin(lat2);
		c = cos(lat);	c2 = cos(lat2);
		tan_lat = s / c;
		M = EQ_RAD * (t_c1 * lat + s2 * (t_c2 + c2 * (t_c3 + c2 * t_c4)));
		dlon = lon - central_meridian;
		if (fabs (dlon) > 360.0) dlon += copysign (360.0, -dlon);
		if (fabs (dlon) > 180.0) dlon = copysign (360.0 - fabs (dlon), -dlon);
		N = EQ_RAD / d_sqrt (1.0 - ECC2 * s * s);
		T = tan_lat * tan_lat;
		T2 = T * T;
		C = t_e2 * c * c;
		A = dlon * D2R * c;
		A2 = A * A;	A3 = A2 * A;	A5 = A3 * A2;
		*x = N * (A + (1.0 - T + C) * (A3 * 0.16666666666666666667) + (5.0 - 18.0 * T + T2 + 72.0 * C - 58.0 * t_e2) * (A5 * 0.00833333333333333333));
		A3 *= A;	A5 *= A;
		*y = (M - t_M0 + N * tan_lat * (0.5 * A2 + (5.0 - T + 9.0 * C + 4.0 * C * C) * (A3 * 0.04166666666666666667) + (61.0 - 58.0 * T + T2 + 600.0 * C - 330.0 * t_e2) * (A5 * 0.00138888888888888889)));
	}
}

int GMT_getinc (char *line, double *dx, double *dy) {
	/* Special case of getincn use where n is two. */
	double inc[2];
	
	*dy = (GMT_getincn (line, inc, 2) == 1) ? inc[0] : inc[1];
	*dx = inc[0];
	return (0);
}

int GMT_getincn (char *line, double inc[], int n) {
	int last, i;
	char tmpstring[512], *p;
	double scale;
	
	/* Dechipers dx/dy/dz/dw/du/dv/... increment strings with n items */
	
	memset ((void *)inc, 0, (size_t)(n * sizeof (double)));
	strcpy (tmpstring, line);	/* Since strtok clobbers the string */
	
	p = (char *)strtok (tmpstring, "/");
	i = 0;
	while (p && i < n) {
		last = strlen (p) - 1;
		if (p[last] == 'm' || p[last] == 'M') {	/* Gave arc minutes */
			p[last] = 0;
			scale = 1.0 / 60.0;
		}
		else if (p[last] == 'c' || p[last] == 'C') {	/* Gave arc seconds */
			p[last] = 0;
			scale = 1.0 / 3600.0;
		}
		else	/*  No units given */
			scale = 1.0;
		if ( (sscanf(p, "%lf", &inc[i])) != 1) {
			mexPrintf ("MANSINHA ERROR: Unable to decode %s as a floating point number\n",p);
			mexErrMsgTxt ("\n");
		}
		inc[i] *= scale;
		p = (char *)strtok (NULL, "/");
		i++;	/* Goto next increment */
	}

	return (i);	/* Returns the number of increments found */
}
