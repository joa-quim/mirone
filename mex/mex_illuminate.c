/*
 *      Copyright (c) 2004 by J. Luis & P. Wessel
 *      See COPYING file for copying and redistribution conditions.
 *
 *      This program is free software; you can redistribute it and/or modify
 *      it under the terms of the GNU General Public License as published by
 *      the Free Software Foundation; version 2 of the License.
 *
 *      This program is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *      GNU General Public License for more details.
 *
 *      Contact info: w3.ualg.pt/~jluis/mirone
 *--------------------------------------------------------------------*/
/* Program:	mex_illuminate.c
 * Purpose:	matlab callable routine to do illumination on Mirone images
 * Author:	J Luis, but based on much code from GMT
 * Date:	15/02/04
 *
 * Modified:	17/08/04
 *		It now also accepts second input argument in single precision 
*/

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif
#define I_255	(1.0 / 255.0)
#define hsv_max_saturation 0.1
#define hsv_min_saturation 1.0
#define hsv_max_value 1.0
#define hsv_min_value 0.3
/*#define copysign(x,y) _copysign(x,y)
If compiling under unix probably the above line should be removed */
/* In non-Windows this is may not be necessary (or give conflicts) */
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#define mn_data(m,n) (ny*(n)+(m))
#define mnk_data(k,m,n) (k*ny*nx + ny*(n) + m)

#include <math.h>
#include "mex.h"

void GMT_rgb_to_hsv(int rgb[], double *h, double *s, double *v);
void GMT_hsv_to_rgb(int rgb[], double h, double s, double v);
void GMT_illuminate (double intensity, int rgb[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int rgb_pt[3], nx, ny, m, n, nm, k1, k2, ndims[2], nsubs, is_single = 0;
	char *r_out, *g_out, *b_out;
	unsigned char *rgb;
	float  *R_s;
	double intensity, *R_d;

	if (nrhs != 2 || nlhs != 3) {
		mexPrintf ("usage: [r,g,b] = mex_illuminate(rgb,R);\n");
		mexPrintf (" 	where rgb is a [n_row,n_col,3] uint8 matrix\n");
		mexPrintf (" 	R is a double precision matrix with reflectance coeficients\n");
		mexPrintf (" 	and r,g,b are 2D short int matrices with their color changed\n");
		return;
	}

	if (mxIsSingle(prhs[1]))
		is_single = 1;

	/* Check that both matrix have appropriate dimensions */
	nsubs = mxGetNumberOfDimensions(prhs[0]);
	if (nsubs != 3) {
		mexPrintf("mex_illuminate error: First input array must be 3D\n");
		return;
	}
	if (!mxIsUint8(prhs[0])) {
		mexPrintf("mex_illuminate error: Image array MUST be of type byte (or uchar)\n");
		return;
	}
	nsubs = mxGetNumberOfDimensions(prhs[1]);
	if (nsubs != 2) {
		mexPrintf("mex_illuminate error: Second input array must be 2D\n");
		return;
	}
	if ( (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) && (mxGetClassID(prhs[1]) != mxSINGLE_CLASS) ) {
		mexPrintf("mex_illuminate error: Reflectance array MUST be of type single OR double\n");
		return;
	}
	ny = mxGetM(prhs[0]);	nx = mxGetN(prhs[0]);
	nx /= 3;	/* Divide by 3 because Image array should have 3 pages */
	if ( (ny != mxGetM(prhs[1])) || (nx != mxGetN(prhs[1])) ) {
		mexPrintf("mex_illuminate error: Image and Reflectance arrays MUST have the same number of rows & columns\n");
		return;
	}

	/* Create a matrix for the return arrays */
	rgb = (unsigned char *)mxGetData(prhs[0]);
	ny = mxGetM(prhs[0]);
	nx = mxGetN(prhs[0]);
	nx /= 3;
	ndims[0] = ny;		ndims[1] = nx;
	if (is_single)
		R_s = (float *)mxGetData(prhs[1]);
	else
		R_d = mxGetPr(prhs[1]);
	plhs[0] = mxCreateNumericMatrix(ny, nx, mxUINT8_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(ny, nx, mxUINT8_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(ny, nx, mxUINT8_CLASS, mxREAL);
	r_out = (char *)mxGetData(plhs[0]);
	g_out = (char *)mxGetData(plhs[1]);
	b_out = (char *)mxGetData(plhs[2]);
	nm = nx * ny;
	for (m = 0; m < ny; m++){
		for (n = 0; n < nx; n++){
			k1 = mn_data(m,n);
			rgb_pt[0] = (int)rgb[k1];		/* k1 = mnk_data(0,m,n);*/
			rgb_pt[1] = (int)rgb[nm + k1];		/* nm + k1 = mnk_data(1,m,n) */
			rgb_pt[2] = (int)rgb[2*nm + k1];	/* 2*nm + k1 = mnk_data(2,m,n) */ 
			if (is_single)
				intensity = (double)R_s[k1];
			else
				intensity = R_d[k1];
			GMT_illuminate (intensity, rgb_pt);
			r_out[k1] = rgb_pt[0];
			g_out[k1] = rgb_pt[1];
			b_out[k1] = rgb_pt[2];
		}
	}
}

void GMT_rgb_to_hsv (int rgb[], double *h, double *s, double *v) {
	double xr, xg, xb, r_dist, g_dist, b_dist, max_v, min_v, diff, idiff;
	
	xr = rgb[0] * I_255;
	xg = rgb[1] * I_255;
	xb = rgb[2] * I_255;
	max_v = MAX (MAX (xr, xg), xb);
	min_v = MIN (MIN (xr, xg), xb);
	diff = max_v - min_v;
	*v = max_v;
	*s = (max_v == 0.0) ? 0.0 : diff / max_v;
	*h = 0.0;
	if ((*s) == 0.0) return;	/* Hue is undefined */
	idiff = 1.0 / diff;
	/*r_dist = (max_v - xr);
	g_dist = (max_v - xg);
	b_dist = (max_v - xb);*/
	if (xr == max_v)
		/* *h = (b_dist - g_dist) * idiff;*/
		*h = (xg - xb) * idiff;
	else if (xg == max_v)
		/* *h = 2.0 + (r_dist - b_dist) * idiff;*/
		*h = 2.0 + (xb - xr) * idiff;
	else
		/* *h = 4.0 + (g_dist - r_dist) * idiff;*/
		*h = 4.0 + (xr - xg) * idiff;
	(*h) *= 60.0;
	if ((*h) < 0.0) (*h) += 360.0;
}

void GMT_hsv_to_rgb (int rgb[], double h, double s, double v) {
	int i;
	double f, p, q, t, rr, gg, bb;
	
	if (s == 0.0)
		rgb[0] = rgb[1] = rgb[2] = (int) floor (255.999 * v);
	else {
		while (h >= 360.0) h -= 360.0;
		h /= 60.0;
		i = (int)h;
		f = h - i;
		p = v * (1.0 - s);
		q = v * (1.0 - (s * f));
		t = v * (1.0 - (s * (1.0 - f)));
		switch (i) {
			case 0:
				rr = v;	gg = t;	bb = p;
				break;
			case 1:
				rr = q;	gg = v;	bb = p;
				break;
			case 2:
				rr = p;	gg = v;	bb = t;
				break;
			case 3:
				rr = p;	gg = q;	bb = v;
				break;
			case 4:
				rr = t;	gg = p;	bb = v;
				break;
			case 5:
				rr = v;	gg = p;	bb = q;
				break;
		}
		
		rgb[0] = (rr < 0.0) ? 0 : (int) floor (rr * 255.999);
		rgb[1] = (gg < 0.0) ? 0 : (int) floor (gg * 255.999);
		rgb[2] = (bb < 0.0) ? 0 : (int) floor (bb * 255.999);
	}
}

void GMT_illuminate (double intensity, int rgb[])
{
	double h, s, v, di;
	
	/*if (GMT_is_dnan (intensity)) return;*/
	if (intensity == 0.0) return;
	if (fabs (intensity) > 1.0) intensity = Loc_copysign (1.0, intensity);
	
	GMT_rgb_to_hsv (rgb, &h, &s, &v);
	if (intensity > 0.0) {
		di = 1.0 - intensity;
		if (s != 0.0) s = di * s + intensity * hsv_max_saturation;
		v = di * v + intensity * hsv_max_value;
	}
	else {
		di = 1.0 + intensity;
		if (s != 0.0) s = di * s - intensity * hsv_min_saturation;
		v = di * v - intensity * hsv_min_value;
	}
	if (v < 0.0) v = 0.0;
	else if (v > 1.0) v = 1.0;
	if (s < 0.0) s = 0.0;
	else if (s > 1.0) s = 1.0;
	GMT_hsv_to_rgb (rgb, h, s, v);
}
