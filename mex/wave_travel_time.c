/*--------------------------------------------------------------------
 *	$Id$
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
 *	Computes the tsunami travel times given a bathymetry grid and a
 *	tsunami location.
 *	This is the mex version of my C code with very little error checking
 *
 *--------------------------------------------------------------------*/
/*
 * Author:      Joaquim Luis
 * Date:        18-Feb-2003
 * Updated:	12-Mar-2011	
 *
 * */

#include "mex.h"
#include <float.h>
#include <math.h>	/* No f.. idea why it must be included to compile on OSX */

/* Macro definition ij_data(i,j) finds the array index to an element
        containing the real data(i,j) in the padded complex array:  */
#define ij_data(i,j) (nx*(j)+(i))
#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif
#define D2R (M_PI / 180.0)
#define rint(x) (floor((x)+0.5))
#define irint(x) ((int)rint(x))
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#define	FALSE		0
#define	TRUE		1

#define EARTH_RAD	6371008.7714	/* GRS-80 sphere */
#define	D_to_M = 2.0 * M_PI * EARTH_RAD / 360.0;        /* GRS-80 sphere m/degree */

struct HEADER {
	int nx;				/* Number of columns */
	int ny;				/* Number of rows */
	double x_min;			/* Minimum x coordinate */
	double x_max;			/* Maximum x coordinate */
	double y_min;			/* Minimum y coordinate */
	double y_max;			/* Maximum y coordinate */
	double x_inc;			/* x increment */
	double y_inc;			/* y increment */
};
struct HEADER h;

int     *rim, *rim1, *rim2, *rim8, *rim16, nx, ny;	/* Still uglly but will stay till a future revision */

int	find_rim_points (double *Z, int i0, int j0, int k);
int	check_in (int i, int j);
void	n_to_ij (int n, int *i, int *j);
void	set_travel_time (double *Z, double *TT, int i0, int j0, int cp, int geo, int *rim0, int k, int ndatac);
void	do_travel_time (double *Z, double *TT, int i_s, int j_s, double t0, int max_range, int ndatac, int geo);
void	bat_to_speed(double *z_8, int ndatac);
double	arc_dist (int i0, int j0, int ic, int jc, int geo);

/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	error = FALSE, stop = FALSE, fill_voids = TRUE;
	int	fonte = FALSE, hours = FALSE, wall = FALSE, geo = TRUE;
	int	i, i2, j, k, n, i_source, j_source, ic, jc, is_double = 0, is_single = 0;
	int	ndatac;
	double	*Z, *z_8, *TT;
	double	lon_source, lat_source, tmp = DBL_MAX, *pdata, *geog;
	double	west = 0.0, east = 0.0, south = 0.0, north = 0.0, x_inc, y_inc, m_NaN, *head, *t_loc;
	float	*z_4, *tout, k_or_m = 1;

	if (nlhs != 1 || nrhs == 0) {
		mexPrintf ("wave_travel_time - Compute the tsunami travel time (ATTENTION z positive up)\n");
		mexPrintf ("usage: tt = wave_travel_time(Z,h_info,t_loc,geog)\n\n");
		mexPrintf ("where: Z contains a bathymetric array (in Matlab orientation) with z in meters positive up\n");
		mexPrintf ("       h_info is a vector with [x_min,x_max,y_min,y_max,x_inc,y_inc]\n");
		mexPrintf ("       t_loc  is a vector with [source_lon,source_lat]\n");
		mexPrintf ("       geog >= 1 if input array is in geogs, or = 0 if it is in crtesian coordinates\n");
		mexPrintf ("NOTE1: Input Z can be single or double and Out tt is single\n");
		mexPrintf ("NOTE2: Input Z array cannot have NaNs\n");
		return;
	}

	if (mxIsDouble(prhs[0])) {
		z_8 = mxGetPr(prhs[0]);
		is_double = 1;
	}
	else if (mxIsSingle(prhs[0])) {
		z_4 = (float *)mxGetData(prhs[0]);
		is_single = 1;
	}
	else {
		mexPrintf("WAVE_TRAVEL_TIME ERROR: Invalid input data type.\n");
		mexErrMsgTxt("Valid types are double or single.\n");
	}
	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	head  = mxGetPr(prhs[1]);		/* Get header info */
	west = head[0];		east = head[1];
	south = head[2];	north = head[3];
	x_inc = head[4];	y_inc = head[5];
	h.x_min = head[0];	h.x_max = head[1];
	h.y_min = head[2];	h.y_max = head[3];
	h.x_inc = head[4];	h.y_inc = head[5];
	h.nx = nx;		h.ny = ny;
	t_loc = mxGetPr(prhs[2]);	/* Where is the tsunami location? */
	lon_source = t_loc[0];
	lat_source = t_loc[1];
	geog = mxGetPr(prhs[3]);	/* Is the grid in geogs? If yes convert from degrees to meters. */
	if (*geog == 0)
		geo = FALSE;

	m_NaN = mxGetNaN();

	ndatac = nx * ny;

	Z = mxCalloc (nx*ny, sizeof (double));
	TT = mxCalloc (nx*ny, sizeof (double));
	rim = (int *) mxCalloc (2*(nx+ny), sizeof(int));
	rim1 = (int *) mxCalloc (2*(nx+ny), sizeof(int));
	rim2 = (int *) mxCalloc (8, sizeof(int));
	rim8 = (int *) mxCalloc (8, sizeof(int));
	rim16 = (int *) mxCalloc (16, sizeof(int));

	/* Transpose from Matlab orientation to gmt grd orientation */
	/* I have to do what comes next because the original code was written
	   for gmt grids. Maybe one day I'll change to a cleaner programing */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--)
			for (j = 0; j < nx; j++) Z[i2*nx+j] = z_8[j*ny+i];
	}
	else {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--)
			for (j = 0; j < nx; j++) Z[i2*nx+j] = z_4[j*ny+i];
	}

	for (j = 0; j < ndatac; j++) TT[j] = 1000000;	/* Initialize TT */


	/* Compute the velocity field from the bathymetry file (and store it on the same variable) */
	bat_to_speed(Z, ndatac);

	i_source = irint((lon_source - west) / x_inc);
	j_source = irint((north - lat_source) / y_inc);

	/* Now the hard work. Compute the wave travel time */
	do_travel_time(Z, TT, i_source, j_source, 0, 0, ndatac, geo);

	if (fill_voids) {
		for (i = 1; i < ndatac; i++) {	/* Search for voids moving W to E */
			n_to_ij(i, &ic, &jc);
			wall = (((float)ic / (float)nx) - ic / nx == 0.) ? 1 : 0;
			if ((TT[i] == 1000000 && TT[i-1] != 1000000) && Z[i] > 0 && !wall)
				do_travel_time (Z, TT, ic, jc, TT[i-1], 10, ndatac, geo);
		}
		for (i = 0, k = ndatac; i < ndatac - nx + 1; k--, i++) { /* Search for voids moving S to N */
			n_to_ij(k, &ic, &jc);
			wall = (((float)ic / (float)ny) - ic / ny == 0.) ? 1: 0;
			if ((TT[k-nx] == 1000000 && TT[k] != 1000000) && Z[k] > 0 && !wall)
				do_travel_time (Z, TT, ic, jc, TT[k], 10, ndatac, geo);
		}
	}

	for (j = 0; j < ndatac; j++) {
		if (TT[j] >= 1000000)
			TT[j] = m_NaN;
		else
			TT[j] /= 3600;		/* output is in hours */
	}

	mxFree(Z);
	mxFree(rim);	mxFree(rim1);	mxFree(rim2);	mxFree(rim8);	mxFree(rim16);

	/* Create and populate matrix for the return array */
	plhs[0] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
	tout = (float *)mxGetData(plhs[0]);
	/* Transpose from gmt grd orientation to Matlab orientation */
	for (i = 0; i < ny; i++)
		for (j = 0; j < nx; j++)
			tout[j*ny+ny-i-1] = (float)TT[i*nx+j];

	mxFree(TT);
}


void bat_to_speed(double *Z, int ndatac) {
	int	j;
	double	g = 9.806199203;	/* Moritz's 1980 IGF value for gravity at 45 degrees latitude */

	for (j = 0; j < ndatac; j++)
		Z[j] = (Z[j] < 0) ? sqrt(-Z[j] * g) : 0;
}

void	do_travel_time (double *Z, double *TT, int i_source, int j_source, double t0, int max_range, int ndatac, int geo) {

	/* t0 != 0 when a recursive use of this funtion is used to fill voids left by the first call
	   In this case, max_range is taken into acount and means the maximum number of grid nodes
	   away from i_source, j_source. This was the "desenrasque" found to fill the voids without
	   extending the computing time to an incomportable value */

	int i_w, i_e;		/* West and East indexes of the working region relative to the source */
	int j_s, j_n;		/* Idem for North-South */
	int k, l, j, np_rim1, np_rim8, np_rim16, i0, j0, nl_max = 0;
	int first = TRUE;

	i0 = i_source;	j0 = j_source;
	TT[ij_data(i0,j0)] = t0; 
	i_e = h.nx - (i0 + 1);	/* (i0 + 1) because i0 starts counting from 0 */
	i_w = h.nx - i_e - 1;
	j_n = (j0 + 1);		/* Idem for j0 */
	j_s = h.ny - j0 - 1;

	nl_max = MAX(i_e, nl_max);	nl_max = MAX(i_w, nl_max);
	nl_max = MAX(j_n, nl_max);	nl_max = MAX(j_s, nl_max);

	if (t0 > 0.01) nl_max = max_range;
	for (k = 1; k < nl_max; k++) {
		/* Find the points of the first rim, 1 grid node away from current point */
		np_rim1 = find_rim_points (Z, i0, j0, k);	/* find the rim points k nodes way from tsun origin */
		for (j = 0; j < np_rim1; j++)
			rim[j] = rim1[j];
		if (first) {
			for (l = 0; l < np_rim1; l++)
				set_travel_time (Z, TT, i0, j0, l, geo, rim1, k, ndatac);
		}
		first = FALSE;

		for (l = 0; l < np_rim1; l++) {			/* Circulate the rim k points */
			n_to_ij(rim[l], &i0, &j0);			/* Make each rim k point the current point */
			if (i0 == 0 || i0 >= h.nx-1 || j0 == 0 || j0 >= h.ny-1) continue; 
			np_rim8 = find_rim_points (Z, i0, j0, 1);		/* find an extra rim8 arround each rim k point */
			for (j = 0; j < np_rim8; j++)
				rim8[j] = rim1[j];				/* Save rim8 vector */
			np_rim16 = find_rim_points (Z, i0, j0, 2);		/* find a rim16 arround each rim k point */
			for (j = 0; j < np_rim16; j++)
				rim16[j] = rim1[j];				/* Save rim16 vector */
			for (j = 0; j < np_rim8; j++) {				/* Circulate arround the local rim8 points */
				set_travel_time (Z, TT, i0, j0, j, geo, rim8, k, ndatac); 
			}
			for (j = 0; j < np_rim16; j++) {			/* Circulate arround the local rim16 points */
				set_travel_time (Z, TT, i0, j0, j, geo, rim16, k, ndatac); 
			}
		}
		i0 = i_source;	j0 = j_source;		/* Reset to true origin to start another big (increasing) loop */
	}
}

void set_travel_time (double *Z, double *TT, int i0, int j0, int cp, int geo, int *rim0, int k, int ndatac) {
	/* Compute the travel time at one of the 16 nearneighbors from the current point
	   If that point has already a time estimate, choose the minimum of the two */
	int ic, jc, ij, i0j0;
	double tmp, v_mean, ds;

	n_to_ij(rim0[cp], &ic, &jc);	/* compute (i,j) from index n of the matrix in vector form */
	ij = ij_data(ic,jc);
	if (Z[ij] > 0) {	/* Point is not on land */
		i0j0 = ij_data(i0,j0);
		if (TT[i0j0] == 1000000) return;
		v_mean = (Z[ij] + Z[i0j0]) / 2;
		ds = arc_dist (i0,j0,ic,jc,geo);
		tmp = ds / v_mean;
		TT[ij] = MIN((TT[i0j0] + tmp), TT[ij]);
	}
}

int	check_in (int i, int j) {
	int in;
	in = (i <= 0 || i >= h.nx-1 || j <= 0 || j >= h.ny-1) ? 0: 1;
	return (in);
}

int	find_rim_points (double *Z, int i0, int j0, int k) {
	int n = 0, ic_min, ic_max, jc_min, jc_max, in_ll, in_ul, in_ur, in_lr, i, in, ij;

	ic_min = i0 - k;	ic_max = i0 + k;
	jc_min = j0 - k;	jc_max = j0 + k;
	in_ll = check_in (ic_min, jc_min);
	in_ul = check_in (ic_min, jc_max);
	in_ur = check_in (ic_max, jc_max);
	in_lr = check_in (ic_max, jc_min);
	if (in_ll && in_ur) { /* ainda estao todos dentro da grelha */
		for (i = ic_min; i <= ic_max; i++) { 	/* north edge */
			ij = ij_data(i,jc_max);
			if (Z[ij] > 0) rim1[n++] = ij;
		}
		for (i = jc_min; i <= jc_max - 1; i++) {	/* east edge */
			ij = ij_data(ic_max,i);
			if (Z[ij] > 0) rim1[n++] = ij;
		}
		for (i = ic_min; i <= ic_max - 1; i++) {	/* south edge */
			ij = ij_data(i,jc_min);
			if (Z[ij] > 0) rim1[n++] = ij;
		}
		for (i = jc_min + 1; i <= jc_max - 1; i++) {	/* west edge */
			ij = ij_data(ic_min,i);
			if (Z[ij] > 0) rim1[n++] = ij;
		}
		return (n);
	}
	else { 			/* saiu, usar forca bruta */
		for (i = ic_min; i <= ic_max; i++) { 		/* chech along the north face */
			in = check_in (i, jc_max); 
			if (in && Z[ij_data(i,jc_max)] > 0) rim1[n++] = ij_data(i,jc_max);
		}
		for (i = jc_min; i <= jc_max - 1; i++) { 	/* chech along the east face */
			in = check_in (ic_max, i); 
			if (in && Z[ij_data(ic_max,i)] > 0) rim1[n++] = ij_data(ic_max,i);
		}
		for (i = ic_min; i <= ic_max - 1; i++) { 	/* check along the south face */
			in = check_in (i, jc_min); 
			if (in && Z[ij_data(i,jc_min)] > 0) rim1[n++] = ij_data(i,jc_min);
		}
		for (i = jc_min + 1; i <= jc_max - 1; i++) { 	/* check along the west face */
			in = check_in (ic_min, i); 
			if (in && Z[ij_data(ic_min,i)] > 0) rim1[n++] = ij_data(ic_min,i);
		}
		return (n);
	}
}

double	arc_dist (int i0, int j0, int ic, int jc, int geo) {
	/* Compute the arc distance in meters between the two points which have
	   the matix indexes of (i0,j0) and (ic,jc) */
	double ds, lon1, lat1, lon2, lat2, tmp;

	/* Due to conversion from C order to ML order and to use the original algo
	   I need to do the following. However, I still suspect there are two errors
	   that compensate to make it right at the end */
	j0 = h.ny - j0;		jc = h.ny - jc;

	lon1 = h.x_min + i0 * h.x_inc;		lon2 = h.x_min + ic * h.x_inc;

	if (geo) {
		lat1 = (h.y_min + j0 * h.y_inc) * D2R;		lat2 = (h.y_min + jc * h.y_inc) * D2R;
		tmp = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos((lon2 - lon1)*D2R);
		ds = fabs (EARTH_RAD * acos(tmp));
	}
	else {
		lat1 = h.y_min + j0 * h.y_inc;		lat2 = h.y_min + jc * h.y_inc;
		ds = sqrt((lon1-lon2)*(lon1-lon2) + (lat1-lat2)*(lat1-lat2));
	}
	return (ds);
}

void	n_to_ij (int n, int *i, int *j) {
	/* compute (i,j) from index n of the matrix in vector form */
	*j = n / h.nx;
	*i = n - *j * h.nx;
}
