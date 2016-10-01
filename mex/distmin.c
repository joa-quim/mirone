/*--------------------------------------------------------------------
 *	$Id$
 *
 *	Copyright (c) 2004-2016 by J. Luis
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
 * Compute the shortest distance between each point in (lon,lat) and the
 * polyline (r_lon,r_lat) SEGLEN holds the segment length of polyline 
 * (r_lon,r_lat) corresponding to elements of DIST.
 *
 * This MEX is 40x faster than the corresponding matlab code (in compute_euler)
 * Equivalent Matlab call
 * [dist, [xy_near, dists, pesos]] = distmin(lon, lat, r_lon, r_lat, lengthsRot)
 *
 * NOTE: This assumes that both lines have the same "growing direction"
 *
 * Author:	Joaquim Luis
 * Date:	22-Dec-2010
 * 
 */

#include "mex.h"
#include <math.h>
#include <string.h>

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef M_PI
#define M_PI          3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2        1.57079632679489661923
#endif

#ifndef D2R
#define D2R (M_PI / 180.0)
#define R2D (180.0 / M_PI)
#endif

#define KMRAD (111.19507973463158 / D2R) 

#define TRUE 1
#define FALSE 0

#if !defined(copysign)
#define copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))
#endif

#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))
#define d_acos(x) (fabs(x) >= 1.0 ? ((x) < 0.0 ? M_PI : 0.0) : acos(x))
#define d_asin(x) (fabs(x) >= 1.0 ? copysign (M_PI_2, (x)) : asin(x))
#define d_atan2(y,x) ((x) == 0.0 && (y) == 0.0 ? 0.0 : atan2(y, x))
#define d_atan2d(y,x) ((x) == 0.0 && (y) == 0.0 ? 0.0 : atan2d(y,x))
#define d_acosd(x) (fabs(x) >= 1.0 ? ((x) < 0.0 ? 180.0 : 0.0) : acosd(x))
#define d_asind(x) (fabs(x) >= 1.0 ? copysign (90.0, (x)) : asind(x))
#define sind(x) sin((x) * D2R)
#define cosd(x) cos((x) * D2R)
#define tand(x) tan((x) * D2R)
#define asind(x) (asin(x) * R2D)
#define acosd(x) (acos(x) * R2D)
#define atan2d(y,x) (atan2(y,x) * R2D)
#define sincosd(x,s,c) sincos((x) * D2R,s,c)

void geo_to_cart (double lat, double lon, double *a, int degrees);
void cart_to_geo (double *lat, double *lon, double *a, int degrees);
void cross3v (double *a, double *b, double *c);
void normalize3v (double *a);
double mag3v (double *a);
double dot3v (double *a, double *b);
double GMT_distance(double x0, double y0, double x1, double y1);
double great_circle_dist (double lon1, double lat1, double lon2, double lat2);
double great_circle_dist2 (double cosa, double sina, double lon1, double lon2, double lat2);
void sincos (double a, double *s, double *c);
double dists_sph (double *lon, double *lat, int n_pt, double *r_lon, double *r_lat, double *lengthsRot, 
                  int n_pt_rot, int nlhs, double *xy_near, double *dists, double *weights, double class_dist[]);
void dists_cart (double *lon, double *lat, double *r_lon, double *r_lat, double *lengthsRot, 
                 int n_pt, int n_pt_rot, double *dist, double *segLen);
int great_circle_intersection (double A[], double B[], double C[], double X[], double *CX_dist);
int near_a_line_spherical (double lon, double lat, double *line_lon, double *line_lat, int n_pt, int row_s,
                           double *dist_min, int *seg_number, double *x_near, double *y_near);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int     n_pt, n_pt_rot, n, k;
	int     regular;	/* Aux variable to decide what the weights for the weighted sum will be */
	double *dist = NULL, *segLen = NULL, *lon, *lat, *r_lon, *r_lat, *lengths, *lengthsRot;
	double *soma = NULL, tmp, *x_near, *y_near, *xy_near = NULL, *pesos = NULL, *dists = NULL;
	double  class_dist[3] = {25, 50, 75};		/* Classes for weigthed sum in dists_sph */

	if (nlhs != 1 && nlhs != 4)
		mexErrMsgTxt("DISTMIN ERROR: number of outputs must be ONE or FOUR.");

	lon = mxGetPr(prhs[0]); 
	lat = mxGetPr(prhs[1]); 
	lengths = mxGetPr(prhs[2]); 
	r_lon   = mxGetPr(prhs[3]); 
	r_lat   = mxGetPr(prhs[4]); 
	lengthsRot = mxGetPr(prhs[5]); 

	n_pt     = mxGetNumberOfElements(prhs[0]);
	n_pt_rot = mxGetNumberOfElements(prhs[3]);

	/* Heuristic to find if data came in radians or in km (* 6371) */
	tmp = MAX (MAX (MAX (fabs(lat[0]), fabs(lat[(int)(n_pt/3)])), fabs(lat[(int)(n_pt*2/3)])), fabs(lat[n_pt-1]));

	plhs[0] = mxCreateDoubleMatrix (1,1, mxREAL);
	soma    = mxGetPr(plhs[0]);

	if (nlhs == 4) {	/* Return also the coordinates of nearest point on line and weight of that segment */
		plhs[1] = mxCreateDoubleMatrix (n_pt,2, mxREAL);
		xy_near = mxGetPr(plhs[1]);
		plhs[2] = mxCreateDoubleMatrix (n_pt,1, mxREAL);
		dists   = mxGetPr(plhs[2]);
		plhs[3] = mxCreateDoubleMatrix (n_pt,1, mxREAL);
		pesos   = mxGetPr(plhs[3]);
	}

	/* Check that there are vertex closer than 75 km apart, otherwise do not do a weighted sum in dist_sph */
	n = mxGetNumberOfElements(prhs[2]);
	regular = 0;
	for (k = 0; k < n; k++) {
		if (lengths[k] < 75) {regular = 1;	break;}
	}
	if (!regular) {
		n = mxGetNumberOfElements(prhs[5]);
		for (k = 0; k < n; k++) {
			if (lengthsRot[k] < 75) {regular = 1;	break;}
		}
	}
	if (!regular) {
		for (k = 0; k < 3; k++) class_dist[k] = 10000;		/* Just a big number in km */
	}
	/* ---------------------------------------------------------------------------------------------------- */

	if (tmp > 10)
		dists_cart (lon, lat, r_lon, r_lat, lengthsRot, n_pt, n_pt_rot, dist, segLen);
	else {
		*soma  = dists_sph (lon, lat, n_pt, r_lon, r_lat, lengthsRot, n_pt_rot, nlhs, xy_near, dists, pesos, class_dist);
		*soma += dists_sph (r_lon, r_lat, n_pt_rot, lon, lat, lengths, n_pt, 1, xy_near, dists, pesos, class_dist);
	}
	*soma /= 2;
}

double dists_sph (double *lon, double *lat, int n_pt, double *r_lon, double *r_lat, double *lengthsRot, 
                  int n_pt_rot, int nlhs, double *xy_near, double *dists, double *weights, double class_dist[]) {
	int k, ind = 0, seg_number;
	double	dist_min, peso, pesos = 0, soma = 0, x_near, y_near;

	for (k = 0; k < n_pt; ++k) {		/* Loop over fixed line vertices */
		if (ind < 0 || ind >= n_pt_rot) ind = 0;	/* Reset this counter */
		x_near = 0;			/* Should not be need but: despair solution to try to avoid a random crash */
		if (near_a_line_spherical (lon[k], lat[k], r_lon, r_lat, n_pt_rot, ind, &dist_min, &seg_number, &x_near, &y_near)) {

			ind = seg_number;

			/* Remember that for 'regular' isochrons -> class_dist[3] = {25, 50, 75}; */
			if (lengthsRot[ind] <= class_dist[0])
				peso = 1;
			else if (lengthsRot[ind] > class_dist[0] && lengthsRot[ind] <= class_dist[1])
				peso = 0.5;
			else if (lengthsRot[ind] > class_dist[1] && lengthsRot[ind] <= class_dist[2])
				peso = 0.25;
			else
				peso = 0;

			soma += dist_min * peso;
			pesos += peso;

			if (nlhs == 4) {		/* If 3 outputs */
				xy_near[k] = x_near;
				xy_near[k + n_pt] = y_near;
				weights[k] = peso;
				dists[k]   = dist_min;
			}
		}
	}
	soma /= pesos;
	return (soma);
}

static int near_a_line_spherical (double lon, double lat, double *line_lon, double *line_lat, int n_pt, int row_s,
	                              double *dist_min, int *seg_number, double *x_near, double *y_near) {
	int row, j0, ind_s, ind_e, row_e, k;
	double d, A[3], B[3], C[3], X[3], xlon, xlat, cx_dist, dist_AB, aux, fraction, DX, DY, coslat, sinlat;

	/* MODIFIED VERSION (OPTIMIZED FOR SPEED UNDER LOCAL CONSTRAINTS) OF THE GMT ORIGINAL FUNCTION */

	/* Return the minimum distance via dist_min. In addition, if > 1:
	   We also return the coordinate of nearest point via x_near, y_near.
	   and segment number in seg_number
	   The function will always return TRUE, except if the point projects onto the extension of the line */

	if (n_pt < 2) return (FALSE);	/* 1-point "line" is a point; skip segment check */

	/* Find nearest point on this line */

	*dist_min = 1e30;	/* Want to find the minimum distance so init to huge */

	/* Make a first, very crude and quick scan using Cartesian dists to find the closest pt to line */
	coslat = cos(lat);	sinlat = sin(lat);
	for (row = row_s; row < n_pt; row++) {	/* loop over nodes on current line */
		DY = (lat-line_lat[row]);	DX = (lon-line_lon[row]) * coslat;
		d = DX*DX + DY*DY;
		if (d < *dist_min) {
			*dist_min = d;
			*seg_number = row;
			*x_near = (double)row;
		}
	}
	*dist_min = 1e30;
	row_s = MAX(0,*seg_number-1);	row_e = MIN(n_pt,*seg_number+2);

	/*for (row = 0; row < n_pt; row++) {*/	/* loop over nodes on current line */
	for (row = row_s; row < row_e; row++) {		/* loop over nodes on current line */
		/* Distance between our pt and row'th node on seg'th line */
		d = great_circle_dist2 (coslat, sinlat, lon, line_lon[row], line_lat[row]);
		if (d < *dist_min) {			/* Update minimum distance */
			*dist_min = d;
			*seg_number = row;			/* Also update pt of nearest pt on the line */
			*x_near = line_lon[row];
			*y_near = line_lat[row];	/* Also update (x,y) of nearest pt on the line */
		}
	}

	if (*seg_number == 0) 				/* Restrict great circle search to the neighbor segments only */
		{ind_s = 0;			ind_e = 1;}
	else if (*seg_number == (n_pt - 1))
		{ind_s = n_pt - 2;		ind_e = n_pt - 1;}
	else
		{ind_s = *seg_number - 1;	ind_e = *seg_number + 1;}

	/* If we get here we must check for intermediate points along the great circle lines between segment nodes.*/

	geo_to_cart (lat, lon, C, FALSE);			/* Our point to test is now C */
	geo_to_cart (line_lat[0], line_lon[0], B, FALSE);	/* 3-D vector of end of last segment */

	/*for (row = 1; row < n_pt; row++) {*/	/* loop over great circle segments on current line */
	for (row = ind_s, k = 0; row <= ind_e; row++) {				/* loop over great circle segments on current line */
		memcpy (A, B, 3 * sizeof(double));						/* End of last segment is start of new segment */
		geo_to_cart (line_lat[row], line_lon[row], B, FALSE);	/* 3-D vector of end of this segment */
		if (great_circle_intersection (A, B, C, X, &cx_dist)) continue;	/* X not between A and B */
		k++;		/* It this increments, point is between A and B */
		/* Get lon, lat of X, calculate distance, and update min_dist if needed */
		cart_to_geo (&xlat, &xlon, X, FALSE);

		/* Distance between our point and closest perpendicular point on seg'th line */
		d = great_circle_dist2 (coslat, sinlat, lon, xlon, xlat);
		if (d < *dist_min) {			/* Update minimum distance */
			*dist_min = d;
			*x_near = xlon;		*y_near = xlat;		/* Update (x,y) of nearest point on the line */

			/* Update pt of nearest point on the line */
			j0 = row - 1;
			dist_AB = great_circle_dist (line_lon[j0], line_lat[j0], line_lon[row], line_lat[row]);
			fraction = (dist_AB > 0.0) ? great_circle_dist (line_lon[j0], line_lat[j0], xlon, xlat) / dist_AB : 0.0;
			aux = (double)j0 + fraction;
		}
	}

	*seg_number = (int)aux;
	if ( (!k && (ind_s == 0)) || (!k && (ind_e == n_pt-1)) )
		return (FALSE);		/* Node outside end points */

	return (TRUE);
}

static int great_circle_intersection (double A[], double B[], double C[], double X[], double *CX_dist) {
	/* A, B, C are 3-D Cartesian unit vectors, i.e., points on the sphere.
	 * Let points A and B define a great circle, and consider a
	 * third point C.  A second great cirle goes through C and
	 * is orthogonal to the first great circle.  Their intersection
	 * X is the point on (A,B) closest to C.  We must test if X is
	 * between A,B or outside.
	 */
	int i;
	double P[3], E[3], M[3], Xneg[3], cos_AB, cos_MX1, cos_MX2, cos_test;

	cross3v (A, B, P);			/* Get pole position of plane through A and B (and origin O) */
	normalize3v (P);			/* Make sure P has unit length */
	cross3v (C, P, E);			/* Get pole E to plane through C (and origin) but normal to A,B (hence going through P) */
	normalize3v (E);			/* Make sure E has unit length */
	cross3v (P, E, X);			/* Intersection between the two planes is oriented line*/
	normalize3v (X);			/* Make sure X has unit length */
	/* The X we want could be +x or -X; must determine which might be closest to A-B midpoint M */
	for (i = 0; i < 3; i++) {
		M[i] = A[i] + B[i];
		Xneg[i] = -X[i];
	}
	normalize3v (M);			/* Make sure M has unit length */
	/* Must first check if X is along the (A,B) segment and not on its extension */

	cos_MX1 = dot3v (M, X);		/* Cos of spherical distance between M and +X */
	cos_MX2 = dot3v (M, Xneg);		/* Cos of spherical distance between M and -X */
	if (cos_MX2 > cos_MX1) memcpy (X, Xneg, 3 * sizeof(double));		/* -X is closest to A-B midpoint */
	cos_AB = fabs (dot3v (A, B));	/* Cos of spherical distance between A,B */
	cos_test = fabs (dot3v (A, X));	/* Cos of spherical distance between A and X */
	if (cos_test < cos_AB) return 1;	/* X must be on the A-B extension if its distance to A exceeds the A-B length */
	cos_test = fabs (dot3v (B, X));	/* Cos of spherical distance between B and X */
	if (cos_test < cos_AB) return 1;	/* X must be on the A-B extension if its distance to B exceeds the A-B length */

	/* X is between A and B.  Now calculate distance between C and X */

	*CX_dist = dot3v (C, X);		/* Cos of spherical distance between C and X */
	return (0);				/* Return zero if intersection is between A and B */
}


static void geo_to_cart (double lat, double lon, double *a, int degrees) {
	/* Convert geographic latitude and longitude (lat, lon)
	   to a 3-vector of unit length (a). If degrees = TRUE,
	   input coordinates are in degrees, otherwise in radian */

	double clat, clon, slon;

	if (degrees) {
		lat *= D2R;
		lon *= D2R;
	}
	sincos (lat, &a[2], &clat);
	sincos (lon, &slon, &clon);
	a[0] = clat * clon;
	a[1] = clat * slon;
}

static void cart_to_geo (double *lat, double *lon, double *a, int degrees) {
	/* Convert a 3-vector (a) of unit length into geographic
	   coordinates (lat, lon). If degrees = TRUE, the output coordinates
	   are in degrees, otherwise in radian. */

	if (degrees) {
		*lat = d_asind (a[2]);
		*lon = d_atan2d (a[1], a[0]);
	}
	else {
		*lat = d_asin (a[2]);
		*lon = d_atan2 (a[1], a[0]);
	}
}

static void cross3v (double *a, double *b, double *c) {
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

static void normalize3v (double *a) {
	double r_length;
	r_length = mag3v(a);
	if (r_length != 0.0) {
		r_length = 1.0 / r_length;
		a[0] *= r_length;
		a[1] *= r_length;
		a[2] *= r_length;
	}
}

static double mag3v (double *a) {
	return (sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
}

static double dot3v (double *a, double *b) {
	return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

double GMT_distance(double x0, double y0, double x1, double y1) {
	return (great_circle_dist (x0, y0, x1, y1) * KMRAD);
}

static double great_circle_dist (double lon1, double lat1, double lon2, double lat2) {
	/* great circle distance on a sphere in degrees */
	double cosa, cosb, sina, sinb, cos_c;

	if (lat1 == lat2 && lon1 == lon2) return (1.0 * KMRAD);

	sincos (lat1, &sina, &cosa);
	sincos (lat2, &sinb, &cosb);

	cos_c = sina*sinb + cosa*cosb*cos(lon1-lon2);
	return (acos (cos_c) * KMRAD);
}

static double great_circle_dist2 (double cosa, double sina, double lon1, double lon2, double lat2) {
	/* great circle distance on a sphere. The difference is that sincos(lat1) are transmiited  */
	double cosb, sinb, cos_c;

	sincos (lat2, &sinb, &cosb);

	cos_c = sina*sinb + cosa*cosb*cos(lon1-lon2);
	return (acos (cos_c) * KMRAD);
}

void sincos (double a, double *s, double *c) {
	*s = sin (a);
	*c = cos (a);
}

void dists_cart (double *lon, double *lat, double *r_lon, double *r_lat, double *lengthsRot, 
                 int n_pt, int n_pt_rot, double *dist, double *segLen) {

	int	i, k, ind = 0;
	double	min = 1e20, Dsts, tmp1, tmp2, Q1[2], Q2[2], Q3[2], DQ[2], D1, D2;

	for (k = 0; k < n_pt; ++k) {			/* Loop over fixed line vertices */
		min = 1e20;	/*ind = 0;*/
		/* Now we only have to look starting at previously found index and ahead */
		for (i = ind; i < n_pt_rot; ++i) {	/* Loop over rotated lined vertices */
			tmp1 = lon[k] - r_lon[i];	tmp2 = lat[k] - r_lat[i];
			Dsts = tmp1*tmp1 + tmp2*tmp2;
			if (Dsts < min) {
				min = Dsts;		ind = i;
			}
		}

		if (ind == 0 || ind == n_pt_rot-1) continue;

		Q1[0] = r_lon[ind-1];	Q1[1] = r_lat[ind-1];
		Q2[0] = r_lon[ind];	Q2[1] = r_lat[ind];
		DQ[0] = Q2[0] - Q1[0];	DQ[1] = Q2[1] - Q1[1];
		D1 = fabs(DQ[0]*(lat[k]-Q1[1]) - DQ[1]*(lon[k]-Q1[0])) / sqrt(DQ[0]*DQ[0] + DQ[1]*DQ[1]);
		Q1[0] = r_lon[ind];	Q1[1] = r_lat[ind];
		Q3[0] = r_lon[ind+1];	Q3[1] = r_lat[ind+1];
		DQ[0] = Q3[0] - Q2[0];	DQ[1] = Q3[1] - Q2[1];
		D2 = fabs(DQ[0]*(lat[k]-Q2[1]) - DQ[1]*(lon[k]-Q2[0])) / sqrt(DQ[0]*DQ[0] + DQ[1]*DQ[1]);

		if (D1 < D2) {
			dist[k] = D1;	segLen[k] = lengthsRot[ind-1];
		}
		else {
			dist[k] = D2;	segLen[k] = lengthsRot[ind];
		}
	}
}
