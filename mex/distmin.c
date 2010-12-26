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
/*#ifdef _OPENMP
#include <omp.h>
#endif*/

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	i, k, n_pt, n_pt_rot, ind;
	double	min, Dsts, tmp1, tmp2, Q1[2], Q2[2], DQ[2], D1, D2;
	double	*dist, *segLen, *lon, *lat, *r_lon, *r_lat, *lengthsRot;

	lon = mxGetPr(prhs[0]); 
	lat = mxGetPr(prhs[1]); 
	r_lon = mxGetPr(prhs[2]); 
	r_lat = mxGetPr(prhs[3]); 
	lengthsRot = mxGetPr(prhs[4]); 

	n_pt = mxGetNumberOfElements(prhs[0]);
	n_pt_rot = mxGetNumberOfElements(prhs[2]);

	plhs[0] = mxCreateDoubleMatrix (n_pt,1, mxREAL);
	dist = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix (n_pt,1, mxREAL);
	segLen = mxGetPr(plhs[1]);

/*#pragma omp parallel for*/
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

	if (nlhs == 3) {	/* Compute a weighted sum of dist */
		double peso, pesos, *soma;

		plhs[2] = mxCreateDoubleMatrix (1,1, mxREAL);
		soma = mxGetPr(plhs[2]);

		for (k = 0; k < n_pt; ++k) {
			if (dist[k] == 0) continue;	/* Pts outside the fixed line */

			if (segLen[k] <= 50)
				peso = 1;
			else if (segLen[k] > 50 && segLen[k] < 80)
				peso = 0.25;
			else
				peso = 0;

			*soma += dist[k] * peso;
			pesos += peso;
		}
		*soma /= pesos;
	}
}














