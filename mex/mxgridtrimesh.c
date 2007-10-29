/*  Author:  Willie Brink  [ w.brink@shu.ac.uk ]
             April 2007

	Joaquim Luis
		October 2007

	Modified to increase efficiency and reduce memory consumption
	(long dead to the double ints and other double shits)
	It now runs more than twice fast and consumes < 1/6 of the memory
	Also removed useless NaN manipulations with ints
*/

#include "mex.h"

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
	int	nx, ny, n, m, i, j, k, offset = 0;
	int	east, west, north, south, *F;
	float	*Z, nan;
	double	*X, *Y, *V;
	double	v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z;
	double	w1, w2, w3, z;
	double	minx, maxx, miny, maxy;
	double	tmp1, tmp2, tmp3, tmp4;
	
	if (nrhs == 0) {
		mexPrintf("Z = mxgridtrimesh(F,V,X,Y);\n");
		mexPrintf("where:\n");
		mexPrintf("\tF is a Mx3 integer array with the triangles indices\n");
		mexPrintf("\tV is a Nx3 double array with the triangle vertex\n");
		mexPrintf("\tX and Y are vectors with the grid's coordinates\n\n");
		mexPrintf("Z is a [numel(Y) x numel(X)] SINGLE array\n");
		return;
	}

	if (nrhs == 5)
		offset = 1;		/* Make it 1 if base 1 (Matlab mesh) */ 

	n = mxGetN(prhs[2]);	m = mxGetM(prhs[2]);
	if (n > 1 && m > 1)
		mexErrMsgTxt("X must be a vector.");
	else
		nx = MAX(n,m);

	n = mxGetN(prhs[3]);	m = mxGetM(prhs[3]);
	if (n > 1 && m > 1)
		mexErrMsgTxt("Y must be a vector.");
	else
		ny = MAX(n,m);
	/* dimensions of input data */
	m = mxGetM(prhs[0]);
	n = mxGetM(prhs[1]);
	
	/* pointers to input data */
	F = (int *)mxGetData(prhs[0]);
	V = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	Y = mxGetPr(prhs[3]);
	
	/* create mxArray and point for the output data */
	plhs[0] = mxCreateNumericMatrix(ny,nx,mxSINGLE_CLASS,mxREAL);
	
	/* pointer to the output data */
	Z = (float *)mxGetData(plhs[0]);
	
	/* initialise output */
	nan = mxGetNaN();
	for (i = 0; i < nx*ny; i++) Z[i] = nan;

	/* consider every triangle, projected to the x-y plane and determine whether gridpoints lie inside */
	for (i = 0; i < m; i++) {
		v1x = V[F[i]-offset];		v1y = V[F[i]-offset+n];		v1z = V[F[i]-offset+2*n];
		v2x = V[F[m+i]-offset];		v2y = V[F[m+i]-offset+n];	v2z = V[F[m+i]-offset+2*n];
		v3x = V[F[2*m+i]-offset];	v3y = V[F[2*m+i]-offset+n];	v3z = V[F[2*m+i]-offset+2*n];
		/* we'll use the projected triangle's bounding box: of the form (minx,maxx) x (miny,maxy) */
		minx = v1x;	minx = MIN(MIN(v2x, minx), v3x);
		maxx = v1x;	maxx = MAX(MAX(v2x, maxx), v3x);

		/* find smallest x-grid value > minx, and largest x-grid value < maxx */
		east = west = j = 0;
		while (j < nx && X[j] < minx) j++; 
		if (j < nx) west = j;
		j = nx-1; 
		while (j >= 0 && X[j] > maxx) j--; 
		if (j >= 0) east = j;
		/* if there are gridpoints strictly inside bounds (minx,maxx), continue */        
		if (east-west >= 0) {
			miny = v1y;	miny = MIN(MIN(v2y, miny), v3y);
			maxy = v1y;	maxy = MAX(MAX(v2y, maxy), v3y);

			/* find smallest y-grid value > miny, and largest y-grid value < maxy */
			north = south = j = 0;
			while (j < ny && Y[j] < miny) j++; 
			if (j < ny) north = j;
			j = ny-1; 
			while (j >= 0 && Y[j] > maxy) j--; 
			if (j >= 0) south = j;
			/* if, further, there are gridpoints strictly inside bounds (miny,maxy), continue */
			if (south-north >= 0) {
				/* we now know that there might be gridpoints bounded by (west,east) x (north,south)
				   that lie inside the current triangle, so we'll test each of them */
				tmp1 = v1y*v2x - v1x*v2y;	tmp2 = v2y*v3x - v2x*v3y;
				tmp3 = 1. / (tmp1 - v1y*v3x + tmp2 + v1x*v3y);
				tmp4 = 1. / (v1y*(v2x - v3x) + v2y*v3x - v2x*v3y + v1x*(-v2y + v3y));

				for (j = west; j <= east; j++) {
					int j_ny = j * ny;
					for (k = north; k <= south; k++) {
						/* calculate barycentric coordinates of gridpoint w.r.t. current (projected) triangle */
						w1 = (tmp2 - v2y*X[j] + v3y*X[j] + v2x*Y[k] - v3x*Y[k]) * tmp3;
						w2 = (-(v3y*X[j]) + v1y*(-v3x + X[j]) + v1x*(v3y - Y[k]) + v3x*Y[k]) * tmp4;
						w3 = (tmp1 - v1y*X[j] + v2y*X[j] + v1x*Y[k] - v2x*Y[k]) * tmp3;

						if (w1 > 0 && w2 > 0 && w3 > 0) {
							/* use barycentric coordinates to calculate z-value */
							z = (float)(w1*v1z + w2*v2z + w3*v3z);
							if (mxIsNaN(Z[j_ny+k]) || z > Z[j_ny+k]) Z[j_ny+k] = z;
						}
					}
				}
			}
		}

	}
}
