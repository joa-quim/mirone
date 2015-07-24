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

/* Program:	mex_illuminate.c
 * Purpose:	matlab callable routine to do illumination on Mirone images
 * Author:	J Luis, but based on much code from GMT
 * Date:	15/02/04
 *
 * Modified:	17/08/04
 *		It now also accepts second input argument in single precision
 *
 *		24/02/2012
 *		Use OpenMP when #if HAVE_OPENMP. Can also do inplace illumination
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
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#define mn_data(m,n) (ny*(n)+(m))
#define mnk_data(k,m,n) (k*ny*nx + ny*(n) + m)

#include <math.h>
#include "mex.h"
#include <time.h>

#if HAVE_OPENMP
#include <omp.h>
#endif

void GMT_rgb_to_hsv(int rgb[], double *h, double *s, double *v);
void GMT_hsv_to_rgb(int rgb[], double h, double s, double v);
void GMT_illuminate (double intensity, int rgb[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int rgb_pt[3], nx, ny, m, n, nm, k1, k2, k3, nsubs, is_single = 0, dims[3];
	char *r_out, *g_out, *b_out, *rgb_out;
	unsigned char *rgb;
	float  *R_s;
	double intensity, *R_d;
	clock_t tic;

	if (nrhs != 2 || nlhs == 2 || nlhs > 3) {
		mexPrintf ("usage: rgb_out = mex_illuminate(rgb,R);\n");
		mexPrintf ("	where rgb is a [n_row,n_col,3] uint8 matrix\n");
		mexPrintf ("	R is a single|double precision matrix with reflectance coeficients\n");
		mexPrintf ("	and rgb_out is the illuminated image.\n");
		mexPrintf ("	Optionaly use the syntax:\n");
		mexPrintf ("	[r,g,b] = mex_illuminate(rgb,R);\n");
		mexPrintf (" 	where r,g,b are 2D short int matrices with their color changed\n");
		mexPrintf ("	Finally, without output one can also do inplace illumination\n");
		mexPrintf ("	mex_illuminate(rgb,R);\n");
		return;
	}

#ifdef MIR_TIMEIT
	tic = clock();
#endif

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
	nm = nx * ny;

	/* Create a matrix for the return arrays */
	rgb = (unsigned char *)mxGetData(prhs[0]);
	if (is_single)
		R_s = (float *)mxGetData(prhs[1]);
	else
		R_d = mxGetPr(prhs[1]);

	if (nlhs == 3) {
		plhs[0] = mxCreateNumericMatrix(ny, nx, mxUINT8_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(ny, nx, mxUINT8_CLASS, mxREAL);
		plhs[2] = mxCreateNumericMatrix(ny, nx, mxUINT8_CLASS, mxREAL);
		r_out = (char *)mxGetData(plhs[0]);
		g_out = (char *)mxGetData(plhs[1]);
		b_out = (char *)mxGetData(plhs[2]);
	}
	else if (nlhs == 1) {
		dims[0] = ny;		dims[1] = nx;		dims[2] = 3;
		plhs[0] = mxCreateNumericArray(3, dims, mxUINT8_CLASS, mxREAL);
		rgb_out = (char *)mxGetData(plhs[0]);
	}

#if HAVE_OPENMP
#pragma omp parallel for  private(k1, rgb_pt)
#endif
	for (k1 = 0; k1 < nm; k1++) {
		k2 = k1 + nm;
		k3 = k2 + nm;
		rgb_pt[0] = rgb[k1];
		rgb_pt[1] = rgb[k2];
		rgb_pt[2] = rgb[k3];
		if (is_single)
			intensity = R_s[k1];
		else
			intensity = R_d[k1];
		if (intensity == 0.0) continue;
		GMT_illuminate (intensity, rgb_pt);

		if (nlhs == 0) {
			rgb[k1] = rgb_pt[0];
			rgb[k2] = rgb_pt[1];
			rgb[k3] = rgb_pt[2];
		}
		else if (nlhs == 1) {
			rgb_out[k1] = rgb_pt[0];
			rgb_out[k2] = rgb_pt[1];
			rgb_out[k3] = rgb_pt[2];
		}
		else {
			r_out[k1] = rgb_pt[0];
			g_out[k1] = rgb_pt[1];
			b_out[k1] = rgb_pt[2];
		}
	}


#ifdef MIR_TIMEIT
	mexPrintf("MEX_ILLUMINATE: CPU ticks = %.3f\tCPS = %d\n", (double)(clock() - tic), CLOCKS_PER_SEC);
#endif

}

void GMT_rgb_to_hsv (int rgb[], double *h, double *s, double *v) {
	double xr, xg, xb, r_dist, g_dist, b_dist, max_v, min_v, diff;
	
	xr = rgb[0] * I_255;
	xg = rgb[1] * I_255;
	xb = rgb[2] * I_255;
	max_v = MAX (MAX (xr, xg), xb);
	min_v = MIN (MIN (xr, xg), xb);
	diff = max_v - min_v;
	*h = 0.0;
	*v = max_v;
	*s = (max_v == 0.0) ? 0.0 : diff / max_v;
	if ((*s) == 0.0) return;	/* Hue is undefined */
	if (xr == max_v)
		*h = (xg - xb) / diff;
	else if (xg == max_v)
		*h = 2.0 + (xb - xr) / diff;
	else
		*h = 4.0 + (xr - xg) / diff;
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

void GMT_illuminate (double intensity, int rgb[]) {
	double h, s, v, di;
	
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
