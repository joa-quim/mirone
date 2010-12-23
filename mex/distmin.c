/*
 * Compute the shortest distance between each point in (lon,lat) and the
 * polyline (r_lon,r_lat) SEGLEN holds the segment length of polyline 
 * (r_lon,r_lat) corresponding to elements of DIST
 *
 * This MEX is 20x faster than the corresponding matlab code (in compute_euler)
 * Equivalent Matlab call
 * [dist, segLen] = distmin(lon, lat, r_lon, r_lat, lengthsRot)
 *
 * Author:	Joaquim Luis
 * Date:	22-Dec-2010
 * 
 */

#include "mex.h"
#include <math.h>

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	i, k, n_pt, n_pt_rot, ind;
	double	min, Dsts, tmp1, tmp2, Q1[2], Q2[2], DQ[2], D1, D2;
	double	*dist, *segLen, *lon, *lat, *r_lon, *r_lat, *lengthsRot;
	double	*lon0, *lat0, *r_lon0, *r_lat0;

	lon0 = mxGetPr(prhs[0]); 
	lat0 = mxGetPr(prhs[1]); 
	r_lon0 = mxGetPr(prhs[2]); 
	r_lat0 = mxGetPr(prhs[3]); 
	lengthsRot = mxGetPr(prhs[4]); 

	n_pt = mxGetNumberOfElements(prhs[0]);
	n_pt_rot = mxGetNumberOfElements(prhs[2]);

	lon = (double *)mxMalloc(n_pt * sizeof(double));
	lat = (double *)mxMalloc(n_pt * sizeof(double));
	r_lon = (double *)mxMalloc(n_pt_rot * sizeof(double));
	r_lat = (double *)mxMalloc(n_pt_rot * sizeof(double));

	for (i = 0; i < n_pt; ++i) {
		lon[i] = lon0[i] * cos(lat0[i]) * 6371;
		lat[i] = lat0[i] * 6371;
	}
	for (i = 0; i < n_pt_rot; ++i) {
		r_lon[i] = r_lon0[i] * cos(r_lat0[i]) * 6371;
		r_lat[i] = r_lat0[i] * 6371;
	}

	plhs[0] = mxCreateDoubleMatrix (n_pt,1, mxREAL);
	dist = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix (n_pt,1, mxREAL);
	segLen = mxGetPr(plhs[1]);

	for (k = 0; k < n_pt; ++k) {
		min = 1e20;	ind = -1;
		for (i = 0; i < n_pt_rot; ++i) {
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
		D1 = fabs(DQ[0]*(lat[k]-Q1[1]) - DQ[1]*(lon[k]-Q1[0])) / 
			sqrt(DQ[0]*DQ[0] + DQ[1]*DQ[1]);
		Q1[0] = r_lon[ind];	Q1[1] = r_lat[ind];
		Q2[0] = r_lon[ind+1];	Q2[1] = r_lat[ind+1];
		DQ[0] = Q2[0] - Q1[0];	DQ[1] = Q2[1] - Q1[1];
		D2 = fabs(DQ[0]*(lat[k]-Q1[1]) - DQ[1]*(lon[k]-Q1[0])) / 
			sqrt(DQ[0]*DQ[0] + DQ[1]*DQ[1]);

		if (D1 < D2) {
			dist[k] = D1;	segLen[k] = lengthsRot[ind-1];
		}
		else {
			dist[k] = D2;	segLen[k] = lengthsRot[ind];
		}
	}

	mxFree(lon); 	mxFree(lat); 	mxFree(r_lon); 	mxFree(r_lat); 
}
