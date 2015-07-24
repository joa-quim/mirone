/*--------------------------------------------------------------------
 *	$Id: $
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
 * Tool to compute moving window block operations. NaNs tolerant and respecteful.
 *
 * Currently implemented:
 *
 * (From definitions in Wilson et al 2007, Marine Geodesy 30:3-35)
 * - TRI: Terrain Ruggedness Index 
 * - TPI: Topographic Position Index
 * - Roughness
 * 
 * - mean  
 * - minimum  
 * - maximum  
 * 
 * 
 * NOTE: Since I don't care about double arrays this program works only with singles.
 *
 * This program does not need a padded array. It does a kinda of In-situ border conditions
 * by mirroring. With this we don't need to make a copy of input array thus leading to
 * important memory savings. Some troubles at array's corners still persist though, which
 * are desguised with pixel corner replacements using closest non perturbed neighbor value.
 *
 * Author:	Joaquim Luis
 * Date:	11-OCT-2009
 * Revised:	
 * 
 */

#include "mex.h"
#include <float.h>
#include <math.h>
#include <string.h>

#define TRUE	1
#define FALSE	0

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define D2R		M_PI / 180.
#define R2D		180. / M_PI
#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#define MIRBLOCK_ALGO_TRI	0
#define MIRBLOCK_ALGO_TPI	1
#define MIRBLOCK_ALGO_ROUGH	2
#define MIRBLOCK_ALGO_AVG	3
#define MIRBLOCK_ALGO_MIN	4
#define MIRBLOCK_ALGO_MAX	5
#define MIRBLOCK_ALGO_SLOPE	6
#define MIRBLOCK_ALGO_ASPECT	7
#define MIRBLOCK_ALGO_RMS	8
#define MIRBLOCK_ALGO_TREND	9
#define MIRBLOCK_ALGO_RESIDUE	10
#define MIRBLOCK_ALGO_RES_RMS	11
#define MIRBLOCK_ALGO_AGC_FAMP	12
#define MIRBLOCK_ALGO_AGC_LAMP	13
#define MIRBLOCK_N_ALGOS	14

typedef void (*PFV) ();		/* PFV declares a pointer to a function returning void */

void TPI        (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void TRI        (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void roughness  (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void average    (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void block_min  (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void block_max  (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void callAlgo   (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void surface_fit(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void agc_fullAmp(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void block_rms  (int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans);
void load_pstuff(double *pstuff, int n_model, double x, double y, int newx, int newy, int basis);
void load_gtg_and_gtd (float *data, int nx, int ny, double *xval, double *yval, double *pstuff, double *gtg, double *gtd, int n_model);
void GMT_gauss (double *a, double *vec, int n_in, int nstore_in, double test, int *ierror, int itriag, int *line, int *isub);
void GMT_cheb_to_pol (double *c, int n, double a, double b, double *d, double *dd);
void prepVars(float *out, int n_win, int nx, int ny, int *n_win2, int *nHalfWin, int *m_stop, 
	int *n_stop, int *n_off, int *inc);

PFV MIR_block_funs[MIRBLOCK_N_ALGOS];

static double globalStore[4];

/* For floats ONLY */
#define ISNAN_F(x) (((*(int32_T *)&(x) & 0x7f800000L) == 0x7f800000L) && \
                    ((*(int32_T *)&(x) & 0x007fffffL) != 0x00000000L))
/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	n_win = 3, i, j;
	int	algo = 1, nm, nx, ny, nfound = 0, argc = 0, n_arg_no_char = 0;
	int	error = FALSE, unknown_nans = TRUE, check_nans = FALSE, subsample = FALSE, overlap = FALSE;
	char	**argv;
	float	*zdata, *out;
	double	*hdr = NULL;

	globalStore[0] = globalStore[1] = globalStore[2] = globalStore[3] = 0;
	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}
	argc++;			/* to account for the program's name to be inserted in argv[0] */

	if ( (n_arg_no_char == 2) && mxIsNumeric(prhs[1]) ) {
		if (mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 9) {
			mexErrMsgTxt("MIRBLOCK: header array must be a 1x9 row vector.\n");
		}
		else {
			hdr = (double *)mxCalloc(10, sizeof(double));
			memcpy((void *)hdr, mxGetData(prhs[1]), 10 * sizeof(double));
			/* The extra field will contain the isGeog [9] info */
		}
	}

	/* get the length of the input string */
	argv = (char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "mirblock";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'A':
					algo = atoi(&argv[i][2]);
					break;
				case 'G':
					if (hdr) hdr[9] = 1;	/* Otherwise, ignore this option that is useless */
					break;
				case 'H':
					overlap = TRUE;
					break;
				case 'N':
					j = atoi(&argv[i][2]);
					if (j && j == 0) {		/* Info is that there are no NaNs here */
						unknown_nans = FALSE;
						check_nans = FALSE;
					}
					else if (j && j == 1) {		/* Info is that input has NaNs */
						unknown_nans = FALSE;
						check_nans = TRUE;
					}
					break;
				case 'S':
					subsample = TRUE;
					break;
				case 'W':
					n_win = atoi(&argv[i][2]);
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (n_arg_no_char == 0 || error) {
		mexPrintf ("mirblock - Compute morphological quantities from DEMs single precision arrays\n\n");
		mexPrintf ("usage: out = mirblock(input, ['-A<0|...|13>'], [-N<0|1>], [-S], [-W<winsize>], [-G])\n");
		
		mexPrintf ("\t<input> is name of input array (singles only)\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t   ---------------------------------\n");
		mexPrintf ("\t-A select algorithm (see Wilson et al 2007, Marine Geodesy 30:3-35):\n");
		mexPrintf ("\t   0 -> Terrain Ruggedness Index\n");
		mexPrintf ("\t   1 -> Topographic Position Index\n");
		mexPrintf ("\t   2 -> Roughness\n");
		mexPrintf ("\t   3 -> Mean\n");
		mexPrintf ("\t   4 -> Minimum\n");
		mexPrintf ("\t   5 -> Mmaximum\n");
		mexPrintf ("\t   6 -> Slope\n"); 
		mexPrintf ("\t   7 -> Aspect\n"); 
		mexPrintf ("\t   8 -> RMS\n"); 
		mexPrintf ("\t   9 -> Trend\n"); 
		mexPrintf ("\t  10 -> Residue\n"); 
		mexPrintf ("\t  11 -> RMS of Residue\n"); 
		mexPrintf ("\t  12 -> Automatic Gain Control (Full Amplitude)\n"); 
		mexPrintf ("\t  13 -> Automatic Gain Control (Local Amplitude)\n"); 
		mexPrintf ("\t-N Inform if input has NaNs (1) or not (0) thus avoiding wasting time with repeated test.\n");
		mexPrintf ("\t-S Subsample option. Means that output will be at the center the window set by -W.\n");
		mexPrintf ("\t-W select the rectangular window size [default is 3, which means 3x3].\n");
		mexPrintf ("\t-G Convert dx,dy in degrees of longitude,latitude into meters (ignored when not needed).\n");
		mexErrMsgTxt("\n");
	}
	
	if (nlhs == 0) {
		mexPrintf("MIRBLOCK ERROR: Must provide an output.\n");
		return;
	}

	if (!mxIsSingle(prhs[0]))
		mexErrMsgTxt("MIRBLOCK ERROR: Invalid input data type. Only valid type is: Single.\n");

	if (n_win % 2 == 0)
		mexErrMsgTxt("MIRBLOCK -W ERROR: Window size must be an odd number.\n");

	if (algo < 0 || algo > MIRBLOCK_N_ALGOS)
		mexErrMsgTxt("MIRBLOCK -A ERROR: unknown algorithm selection.\n");

	nx = mxGetN(prhs[0]);
	ny = mxGetM(prhs[0]);
	nm = nx * ny;

	if (!mxIsNumeric(prhs[0]) || ny < 3 || nx < 3)
		mexErrMsgTxt("MIRBLOCK ERROR: First argument must contain a decent array\n");

	zdata = (float *)mxGetData(prhs[0]);

	if (subsample)
		plhs[0] = mxCreateNumericMatrix (ny / n_win, nx / n_win, mxSINGLE_CLASS, mxREAL);

	else if (overlap) {
		int nnx, nny;
		nnx = nx / n_win;	nnx += nnx / 2;		/* ex: nx = 26; n_win = 5;	nnx = (26 / 5) = 5; nnx = 5 + 5 / 2 = 7; */
		nny = ny / n_win;	nny += nny / 2;
		plhs[0] = mxCreateNumericMatrix (nny, nnx, mxSINGLE_CLASS, mxREAL);
	}
	else		/* Most common use case */
		plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), mxGetDimensions(prhs[0]), 
						mxSINGLE_CLASS, mxREAL);

	out = (float *)(mxGetData(plhs[0]));
	if (subsample) out[0] = 1;	/* Flag to subsample in the algorithm functions */
	if (overlap)   out[0] = 2;	/*  */

	if (unknown_nans) {
		/* Loop over the file and find if we have NaNs. Stop at first NaN occurence. */
		for (i = 0; i < nm; i++) {
			if (ISNAN_F(zdata[i])) {
				check_nans = TRUE;
				break;
			}
		}
	}

	MIR_block_funs[MIRBLOCK_ALGO_TRI]     = (PFV) TRI;
	MIR_block_funs[MIRBLOCK_ALGO_TPI]     = (PFV) TPI;
	MIR_block_funs[MIRBLOCK_ALGO_ROUGH]   = (PFV) roughness;
	MIR_block_funs[MIRBLOCK_ALGO_AVG]     = (PFV) average;
	MIR_block_funs[MIRBLOCK_ALGO_MIN]     = (PFV) block_min;
	MIR_block_funs[MIRBLOCK_ALGO_MAX]     = (PFV) block_max;
	MIR_block_funs[MIRBLOCK_ALGO_SLOPE]   = (PFV) surface_fit;
	MIR_block_funs[MIRBLOCK_ALGO_ASPECT]  = (PFV) surface_fit;
	MIR_block_funs[MIRBLOCK_ALGO_TREND]   = (PFV) surface_fit;
	MIR_block_funs[MIRBLOCK_ALGO_RESIDUE] = (PFV) surface_fit;
	MIR_block_funs[MIRBLOCK_ALGO_RES_RMS] = (PFV) surface_fit;
	MIR_block_funs[MIRBLOCK_ALGO_RMS]     = (PFV) block_rms;
	MIR_block_funs[MIRBLOCK_ALGO_AGC_FAMP]= (PFV) agc_fullAmp;
	MIR_block_funs[MIRBLOCK_ALGO_AGC_LAMP]= (PFV) surface_fit;

	callAlgo(algo, hdr, zdata, out, n_win, nx, ny, check_nans);
}

void callAlgo(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	int	j, n, m, nHalfWin, n_extra, nWinExtended, nHalfWinExtended, nny, idSave = -1;
	float	*pad_cols, *out_cols;

	if (id == MIRBLOCK_ALGO_AGC_LAMP) {id = MIRBLOCK_ALGO_RES_RMS;	idSave = MIRBLOCK_ALGO_AGC_LAMP;}

	MIR_block_funs[id](id, hdr, in, out, n_win, nx, ny, check_nans);

	if (out[0]) return;	/* When subsampling no need to do BC, so we are done */

	nHalfWin = n_win / 2;
	n_extra  = nHalfWin - 1;		/* Will be 0 for W=3, 1 for W=5, 2 for W=7, etc ... */
	nWinExtended     = n_win + n_extra;
	nHalfWinExtended = nHalfWin + n_extra;
	nny = ny + 2 * nHalfWin;

	/* ------------------------- Do the W & E frontier conditions ------------------------ */

	/* Temp arrays for grid padding */
	pad_cols = (float *)mxCalloc((size_t)(nny * nWinExtended), sizeof(float));
	out_cols = (float *)mxCalloc((size_t)(nny * nWinExtended), sizeof(float));

	/*------  West ------ */
	for (n = nHalfWin, m = 0; n > 0; n--, m++)		 /* Mirror nHalfWin columns */
		memcpy((void *)&pad_cols[m * nny + nHalfWin],(void *)&in[n * ny], ny * sizeof(float));
	for (n = 0; n <= nHalfWinExtended; n++, m++)	 	/* Copy first (nHalfWin + n_extra + 1) columns */
		memcpy((void *)&pad_cols[m * nny + nHalfWin],(void *)&in[n * ny], ny * sizeof(float));

	/* Replicate the first/last values of pad_cols into the nHalfWin padding zone */
	for (n = 0; n < n_win; n++) {
		for (m = 0; m < nHalfWin; m++)      pad_cols[m + n*nny] = pad_cols[nHalfWin + n*nny];
		for (m = ny+nHalfWin; m < nny; m++) pad_cols[m + n*nny] = pad_cols[ny+nHalfWin-1 + n*nny];
	}

	MIR_block_funs[id](id, hdr, pad_cols, out_cols, n_win, nWinExtended, nny, check_nans);

	for (n = 0, m = nHalfWin - 1; n < nHalfWin; n++, m--)		/* Put the result in the first nHalfWin columns */
		memcpy((void *)&out[m * ny],(void *)&out_cols[(n+nHalfWin) * nny + nHalfWin], ny * sizeof(float));

	/* ------ East ------ */
	for (n = nx - nHalfWinExtended - 1, m = 0; n < nx; n++, m++)	/* Copy last (nHalfWin + n_extra + 1) columns */
		memcpy((void *)&pad_cols[m * nny + nHalfWin],(void *)&in[n * ny], ny * sizeof(float));
	for (n = nx - 2; n >= nx - nHalfWin - 1; n--, m++)		/* Mirror nHalfWin columns */
		memcpy((void *)&pad_cols[m * nny + nHalfWin],(void *)&in[n * ny], ny * sizeof(float));


	/* Replicate the first/last values of pad_cols into the nHalfWin padding zone */
	for (n = 0; n < n_win; n++) {
		for (m = 0; m < nHalfWin; m++) pad_cols[m + n*nny]      = pad_cols[nHalfWin + n*nny];
		for (m = ny+nHalfWin; m < nny; m++) pad_cols[m + n*nny] = pad_cols[ny+nHalfWin-1 + n*nny];
	}

	MIR_block_funs[id](id, hdr, pad_cols, out_cols, n_win, nWinExtended, nny, check_nans);

	for (n = nx - 1, m = nHalfWin; n > nx - nHalfWin -1; n--, m++)	/* Put the result in the last nHalfWin columns */
		memcpy((void *)&out[n * ny],(void *)&out_cols[m * nny + nHalfWin], ny * sizeof(float));

	mxFree((void *)pad_cols);
	mxFree((void *)out_cols);

	/* ----------------------- Now the N & S frontier conditions ------------------------ */

	pad_cols = (float *)mxCalloc((size_t)(nx * nWinExtended), sizeof(float));
	out_cols = (float *)mxCalloc((size_t)(nx * nWinExtended), sizeof(float));

	/* Starting rows (if this is N or S that depends on grid's orientation) */
	for (n = nHalfWin, m = 0; n > 0; n--, m++)			/* Mirror nHalfWin rows */
		for (j = 0; j < nx; j++)
			pad_cols[j * nWinExtended + m] = in[j*ny + n];
	for (n = 0; n <= nHalfWinExtended; n++, m++)			/* Copy first (nHalfWin + n_extra + 1) rows */
		for (j = 0; j < nx; j++)
			pad_cols[j * nWinExtended + m] = in[j*ny + n];

	MIR_block_funs[id](id, hdr, pad_cols, out_cols, n_win, nx, nWinExtended, check_nans);

	for (n = nHalfWin-1, m = nHalfWin; n >= 0 ; n--, m++)		/* Put the result in the first nHalfWin rows*/
		for (j = nHalfWin; j < nx - nHalfWin; j++)
			out[j*ny + n] = out_cols[j * nWinExtended + m];

	/* Ending rows */
	for (n = ny - nHalfWinExtended - 1, m = 0; n < ny; n++, m++)	/* Copy last (nHalfWin + n_extra + 1) rows */
		for (j = 0; j < nx; j++)
			pad_cols[j * nWinExtended + m] = in[j*ny + n];
	for (n = ny - 2; n >= ny - (nHalfWin + 1); n--, m++)		/* Mirror nHalfWin rows */
		for (j = 0; j < nx; j++)
			pad_cols[j * nWinExtended + m] = in[j*ny + n];

	MIR_block_funs[id](id, hdr, pad_cols, out_cols, n_win, nx, nWinExtended, check_nans);

	for (n = ny - nHalfWin, m = 2*nHalfWin-1; n < ny; n++, m--)	/* Put the result in the last nHalfWin rows */
		for (j = nHalfWin; j < nx - nHalfWin; j++)
			out[j*ny + n] = out_cols[j * nWinExtended + m];

	mxFree((void *)pad_cols);
	mxFree((void *)out_cols);

	/* Since this bloody thing still f on the corners I give up and replace them by their next diagonal neighbor */
	for (n = 0; n < nHalfWin; n++) {
		for (m = 0; m < nHalfWin; m++)					/* First_rows/West */
			out[m + n*ny] = out[nHalfWin * ny + nHalfWin];
		for (m = ny - nHalfWin; m < ny; m++)				/* Last_rows/West */
			out[m + n*ny] = out[(nHalfWin+1) * ny - (nHalfWin+1)];
	}
	for (n = nx - nHalfWin; n < nx; n++) {
		for (m = 0; m < nHalfWin; m++)					/* First_rows/East */
			out[m + n*ny] = out[(nx - nHalfWin - 1) * ny + nHalfWin];
		for (m = ny - nHalfWin; m < ny; m++)
			out[m + n*ny] = out[(nx - nHalfWin) * ny - (nHalfWin+1)];
	}

	/* ================================================================================================ */
	/*	In some cases we might still have unfinished work					
	/* ================================================================================================ */
	if (id == MIRBLOCK_ALGO_AGC_FAMP) {
		float fac, rmsMax = (float)globalStore[1];		/* We stored it there in fullAmp() */
		for (j = 0; j < nx * ny; j++) {
			if (check_nans && ISNAN_F(out[j])) continue;
			fac = (float)MIN(rmsMax / out[j], 10.);		/* Limit amplification to 10 */
			out[j] = in[j] * fac;
		}
	}
	else if (idSave == MIRBLOCK_ALGO_AGC_LAMP) {
		/* Here the "out" array contains the local RMS of data-local_trend at each point */
		float	*trend, fac, rmsMax = 0;

		trend = (float *) mxCalloc ((size_t)(nx * ny), sizeof (float));
		callAlgo(MIRBLOCK_ALGO_TREND, hdr, in, trend, n_win, nx, ny, check_nans);

		for (j = 0; j < nx * ny; j++) {				/* Compute maximum RMS */
			if (check_nans && ISNAN_F(out[j])) continue;
			if (out[j] > rmsMax) rmsMax = out[j];
		}

		for (j = 0; j < nx * ny; j++) {
			if (check_nans && ISNAN_F(in[j]))
				continue;
			else {
				fac = (float)MIN(rmsMax / out[j], 20.);	/* Limit amplification to 20 */
				out[j] = (in[j] - trend[j]) * fac + trend[j];
			}
		}
		mxFree((void *)trend);
	}

}

void TPI(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes TPI. TPI stands for Topographic Position Index, which is defined as the 
	   difference between a central pixel and the mean of its surrounding cells 
	   (Wilson et al 2007, Marine Geodesy 30:3-35). */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, nNaNs, no = 0, inc = 1, o = -1;
	float	*p;
	double	mean = 0;

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */

	/* Option -S is not yet finished. When resume it we should have to have something like
	 * if (!-S) no = nm;
	 * and at the end
	 * if (-S) no++;
	 * and than use out[no] instead of out[nm]
	 */
	for (n = nHalfWin; n < n_stop; n += inc) {
		k = (n - nHalfWin) * ny - inc;		/* Index of window's UL corner */
		for (m = nHalfWin; m < m_stop; m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			mean = 0;	nNaNs = 0;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && ISNAN_F(*p))	/* Ignore this value */
						nNaNs++;
					else
						mean += *p;
					p++;
				}
			}
			out[o] = in[nm] - (float) (mean / (n_win2 - nNaNs));
		}
	}
}

void TRI(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes TRI. TRI stands for Terrain Ruggedness Index, which is defined as the 
	   mean difference between a central pixel and its surrounding cells 
	   (Wilson et al 2007, Marine Geodesy 30:3-35). */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, nNaNs, no = 0, inc = 1, o = -1;
	float	*p;
	double	mean = 0, tmp;

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */

	for (n = nHalfWin; n < n_stop; n += inc) {
		k = (n - nHalfWin) * ny - inc;
		for (m = nHalfWin; m < m_stop; m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			tmp = in[nm];
			mean = 0;	nNaNs = 0;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && ISNAN_F(*p))	/* Ignore this value */
						nNaNs++;
					else
						mean += fabs(tmp - *p);
					p++;
				}
			}
			out[o] = (float) (mean / (n_win2 - nNaNs));
		}
	}
}

void roughness(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes Roughness defined as the the largest inter-cell difference of a central pixel 
	   and its surrounding cell (Wilson et al 2007, Marine Geodesy 30:3-35). */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, no = 0, inc = 1, o = -1;
	float	*p, tmp, min, max;

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */

	for (n = nHalfWin; n < n_stop; n += inc) {
		k = (n - nHalfWin) * ny - inc;
		for (m = nHalfWin; m < m_stop; m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			tmp = in[nm];
			min = max = tmp;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && ISNAN_F(*p)) continue;	/* Ignore this value */
					max = MAX(max, *p);
					min = MIN(min, *p);
					p++;
				}
			}
			out[o] = (max - min);
		}
	}
}

void average(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes the mean of cells  inside rectangular window */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, nNaNs, no = 0, inc = 1, o = -1;
	float	*p;
	double	mean = 0;

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */

	for (n = nHalfWin; n < n_stop; n += inc) {
		k = (n - nHalfWin) * ny - inc;		/* Index of window's UL corner */
		for (m = nHalfWin; m < m_stop; m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			mean = 0;	nNaNs = 0;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && ISNAN_F(*p))	/* Ignore this value */
						nNaNs++;
					else
						mean += *p;
					p++;
				}
			}
			out[o] = (float) (mean / (n_win2 - nNaNs));
		}
	}
}

void block_min(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes the minimum of cells  inside rectangular window */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, no = 0, inc = 1, o = -1;
	float	*p, min;

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */

	for (n = nHalfWin; n < n_stop; n += inc) {
		k = (n - nHalfWin) * ny - inc;
		for (m = nHalfWin; m < m_stop; m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			min = in[nm];
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && ISNAN_F(*p)) continue;	/* Ignore this value */
					min = MIN(min, *p);
					p++;
				}
			}
			out[o] = min;
		}
	}
}

void block_max(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes the minimum of cells  inside rectangular window */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, no = 0, inc = 1, o = -1;
	float	*p, max;

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */

	for (n = nHalfWin; n < n_stop; n += inc) {
		k = (n - nHalfWin) * ny - inc;
		for (m = nHalfWin; m < m_stop; m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			max = in[nm];
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && ISNAN_F(*p)) continue;	/* Ignore this value */
					max = MAX(max, *p);
					p++;
				}
			}
			out[o] = max;
		}
	}
}

void agc_fullAmp(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Maximum Gain Correction to apply at any position (Full Amplitude) */
	/* The subtility here is that the algo has to be called 5 times. One for the whole region and 4
	   others for the boundary conditions. Since the AGC algo comprises two passes we need to wait till
	   all 5 rounds are done before we can finish the job. We use a global variable to control that. */
	int	i;
	float	rmsMax = 0;

	block_rms(id, hdr, in, out, n_win, nx, ny, check_nans);

	for (i = 0; i < nx * ny; i++) {			/* Compute maximum RMS */
		if (check_nans && ISNAN_F(out[i])) continue;
		if (out[i] > rmsMax) rmsMax = out[i];
	}

	globalStore[1] = MAX(rmsMax, globalStore[1]);

}

void block_rms(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes the RMS of cells inside the rectangular window */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, nNaNs, ngood, no = 0, inc = 1, o = -1;
	float	*p;
	double	mean = 0, rms;

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */

	for (n = nHalfWin; n < n_stop; n += inc) {
		k = (n - nHalfWin) * ny - inc;			/* Index of window's UL corner */
		for (m = nHalfWin; m < m_stop; m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			mean = rms = 0.;	nNaNs = 0;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && ISNAN_F(*p))	/* Ignore this value */
						nNaNs++;
					else {
						mean += *p;
						rms  += (*p) * (*p);
					}
					p++;
				}
			}
			ngood = n_win2 - nNaNs;
			if (ngood == 0) {out[nm] = 0;	continue;}
			mean /= (double)ngood; 
			rms /= (double)ngood;
			out[o] = (float) (sqrt(rms - mean*mean));
		}
	}
}

void surface_fit(int id, double *hdr, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, ij, k, l, o, n_off, ierror = 0, n_model = 3;
	int	*line, *isub, n_error = 0, no = 0, inc = 1;
	float	*p, aspect;
	double	zero_test = 1.0e-08, dv, dy, dx, m_per_deg = 1., *co, c[3], *d, *dd;
	double	trend, mean, rms, *diff;
	double	*xval;		/* Pointer for array of change of variable:  x[i]  */
	double	*yval;		/* Pointer for array of change of variable:  y[j]  */
	double	*gtg;		/* Pointer for array for matrix G'G normal equations  */
	double	*gtd;		/* Pointer for array for vector G'd normal equations  */
	double	*pstuff;	/* Pointer for array for Legendre polynomials of x[i],y[j]  */
	float	*data;		/* Pointer for array with a copy of 'in' inside 'n_win' */

	prepVars(out, n_win, nx, ny, &n_win2, &nHalfWin, &m_stop, &n_stop, &n_off, &inc); /* Set vals for vars in pointers */
	co = (double *)mxMalloc((m_stop - nHalfWin) * sizeof(double));

	if (id == MIRBLOCK_ALGO_SLOPE || id == MIRBLOCK_ALGO_ASPECT) {
		if (hdr[9]) {		/* Geog */
			m_per_deg = 111195.01524;		/* Spherical approx (Authalic radius 6371005.076 m) */
			for (m = nHalfWin, k = 0; m < m_stop; m++, k++)
				co[k] = cos((hdr[2] + m * hdr[8]) * D2R);
		}
		else
			for (m = nHalfWin, l = 0; m < m_stop; m++, l++)
				co[l] = 1;

		dy = (n_win - 1) * hdr[8] * m_per_deg;
		dx = (n_win - 1) * hdr[7] * m_per_deg;	/* In Geog case DX must be recomputed inside the loop below */
	}

	xval = (double *) mxCalloc ((size_t)n_win, sizeof (double));
	yval = (double *) mxCalloc ((size_t)n_win, sizeof (double));
	gtg  = (double *) mxCalloc ((size_t)(n_model*n_model), sizeof (double));
	gtd  = (double *) mxCalloc ((size_t)n_model, sizeof (double));
	pstuff = (double *) mxCalloc ((size_t)n_model, sizeof (double));
	data = (float *) mxCalloc ((size_t)n_win2, sizeof (float));

	/* Set up xval and yval lookup tables:  */
	dv = 2. / (double)(n_win - 1);
	for (j = 0; j < n_win-1; j++) xval[j] = -1 + j * dv;
	for (j = 0; j < n_win-1; j++) yval[j] = -1 + j * dv;
	xval[n_win - 1] = yval[n_win - 1] = 1;

	/* Work arrays */
	line =  (int *) mxCalloc ((size_t)n_model, sizeof (int));
	isub =  (int *) mxCalloc ((size_t)n_model, sizeof (int));
	d  = (double *) mxCalloc ((size_t)n_model, sizeof (double));
	dd = (double *) mxCalloc ((size_t)n_model, sizeof (double));
	if (id == MIRBLOCK_ALGO_RES_RMS)
		diff = (double *) mxCalloc ((size_t)n_win2, sizeof (double));

	for (n = nHalfWin, o = -1, m = 0; n < n_stop; n += inc) {		/* Loop over columns */
		k = (n - nHalfWin) * ny - inc;
		for (m = nHalfWin, l = 0; m < m_stop; l++, m += inc) {
			k += inc;
			nm = k + n_off;
			o = (inc == 1) ? nm : (++o);
			if (check_nans && ISNAN_F(in[nm])) {out[o] = in[nm];	continue;}
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					data[i * n_win + j] = *p;
					p++;
				}
			}

			load_gtg_and_gtd(data, n_win, n_win, xval, yval, pstuff, gtg, gtd, n_model);
			GMT_gauss (gtg, gtd, n_model, n_model, zero_test, &ierror, 1, line, isub);
			if (ierror) n_error++;
			/* In the following, remember that since we are working by columns (ML is column major)
			   the X & Y are swapped with respect to GMT (i.e. row major) logic */
			if (n_model == 3) {
				if (id == MIRBLOCK_ALGO_TREND) {	/* Try these first because no need to cheb_to_pol */
					out[o] = (float)gtd[0];	/* gtd[0] = val at window center */
					continue;
				}
				else if (id == MIRBLOCK_ALGO_RESIDUE) {
					out[o] = data[n_win2 / 2] - (float)gtd[0];
					continue;
				}
				else if (id == MIRBLOCK_ALGO_RES_RMS) {
					for (i = ij = 0; i < n_win; i++) {
						for (j = 0; j < n_win; j++, ij++) {
							trend  = gtd[0];
							trend += xval[i] * gtd[2];
							trend += yval[j] * gtd[1];
							diff[ij] = data[ij] - trend;
						}
					}
					/* Now compute RMS of diff */
					mean = rms = 0.;
					for (ij = 0; ij < n_win2; ij++) {
						mean += diff[ij];
						rms  += (diff[ij] * diff[ij]);
					}
					mean /= (double)n_win2; 
					rms  /= (double)n_win2;
					out[o] = (float) (sqrt(rms - mean*mean));
					continue;
				}
				c[0] = gtd[0];	c[1] = gtd[2]; 			/* For X we must make a copy */
				GMT_cheb_to_pol (c, 2, 0, dx * co[l], d, dd);
				gtd[2] = c[1];
				GMT_cheb_to_pol (gtd, 2, 0, dy, d, dd);		/* Here we can use gtd directly */
				if (id == MIRBLOCK_ALGO_SLOPE)
					out[o] = (float)(atan(sqrt(gtd[1]*gtd[1] + gtd[2]*gtd[2])) * R2D);
				else {
					aspect = (float)(-(90 + atan2(gtd[1], gtd[2]) * R2D));
					if (aspect < 0) aspect += 360;
					out[o] = aspect;
				}
			}
			else {		/* This will (eventually) hold the n_model > 3 cases */
				/* Not correct. We still need to find out what to do with the cross term */
				c[0] = gtd[0];	c[1] = gtd[2]; 	c[2] = gtd[4];
				GMT_cheb_to_pol (c, 3, 0, dx * co[l], d, dd);
				gtd[2] = c[1];	gtd[4] = c[2];
				c[1] = gtd[1]; 	c[2] = gtd[3];
				GMT_cheb_to_pol (c, 3, 0, dy, d, dd);
				gtd[1] = c[1];	gtd[3] = c[2];
			}
		}
	}

	if (n_error == (m-1))
		mexPrintf("Gauss returned error codes for all nodes %d\n", ierror);

	mxFree((void *) xval);		mxFree((void *) yval);
	mxFree((void *) gtg);		mxFree((void *) gtd);
	mxFree((void *) pstuff);	mxFree((void *) data);
	mxFree((void *)isub);		mxFree ((void *)line);
	mxFree((void *) d);		mxFree((void *) dd);
	mxFree((void *)co);
	if (id == MIRBLOCK_ALGO_RES_RMS) mxFree((void *)diff);
}

void load_pstuff (double *pstuff, int n_model, double x, double y, int newx, int newy, int basis) {
	/* If either x or y has changed, compute new Legendre polynomials as needed  */


	if (basis == 0) {		/* Tchebyshev polynomials */
		if (newx) {
			pstuff[1] = x;
			if (n_model >= 5) pstuff[4] = 2.0*pstuff[1]*pstuff[1] - 1.0;
		}
		if (newy) {
			pstuff[2] = y;
			if (n_model >= 6) pstuff[5] = 2.0*pstuff[2]*pstuff[2] - 1.0;
		}
	}
	else if (basis == 1) {		/* Legendre polynomials */
		if (newx) {
			pstuff[1] = x;
			if (n_model >= 5) pstuff[4] = 0.5*(3.0*pstuff[1]*pstuff[1] - 1.0);
		}
		if (newy) {
			pstuff[2] = y;
			if (n_model >= 6) pstuff[5] = 0.5*(3.0*pstuff[2]*pstuff[2] - 1.0);
		}
	}
	else {
		if (newx) {
			pstuff[1] = x;
			if (n_model >= 5) pstuff[4] = pstuff[1]*pstuff[1];
		}
		if (newy) {
			pstuff[2] = y;
			if (n_model >= 6) pstuff[5] = pstuff[2]*pstuff[2];
		}
	}

	/* In either case, refresh cross term:  */

	if (n_model >= 4) pstuff[3] = pstuff[1]*pstuff[2];

	return;
}

void load_gtg_and_gtd (float *data, int nx, int ny, double *xval, double *yval, double *pstuff, double *gtg, double *gtd, int n_model) {
	/* Routine to load the matrix G'G (gtg) and vector G'd (gtd)
	for the normal equations.  Routine uses indices i,j to refer
	to the grid file of data, and k,l to refer to the k_row, l_col
	of the normal equations matrix.  We need sums of
	data and model functions in gtg and gtd.  We save time by
	loading only lower triangular part of gtg and then filling
	by symmetry after i,j loop.  */

	int	i, j, k, l, n_used, ij;

	/*	First zero things out to start:  */
	n_used = 0;
	memset(gtd, 0, n_model * sizeof(double));
	memset(gtg, 0, n_model*n_model * sizeof(double));

	/*  Now get going.  Have to load_pstuff separately in i and j,
	because it is possible that we skip data when i = 0.
	Loop over all data:  */

	for (ij = 0, j = 0; j < ny; j++ ) {
		load_pstuff(pstuff, n_model, xval[0], yval[j], 0, 1, 0);
		for (i = 0; i < nx; i++, ij++) {

			if (ISNAN_F (data[ij])) continue;

			n_used++;
			load_pstuff(pstuff, n_model, xval[i], yval[j], 1, 0, 0);

			/* Loop over all gtg and gtd elements:  */
			gtd[0] += data[ij];
			for (k = 1; k < n_model; k++) {
				gtd[k] += (data[ij] * pstuff[k]);
				gtg[k] += pstuff[k];
				for (l = k; l < n_model; l++)
					gtg[k + l*n_model] += (pstuff[k]*pstuff[l]);
			}
		}
	}

	/* Now use more accurate sum for gtg[0], and set symmetry:  */

	gtg[0] = n_used;

	for (k = 0; k < n_model; k++) {
		for (l = 0; l < k; l++)
			gtg[l + k*n_model] = gtg[k + l*n_model];
	}

	return;
}

void GMT_gauss (double *a, double *vec, int n_in, int nstore_in, double test, int *ierror, int itriag, int *line, int *isub) {
 
/* subroutine gauss, by william menke */
/* july 1978 (modified feb 1983, nov 85) */
 
/* a subroutine to solve a system of n linear equations in n unknowns*/
/* gaussian reduction with partial pivoting is used */
/*      a               (sent, destroyed)       n by n matrix           */
/*      vec             (sent, overwritten)     n vector, replaced w/ solution*/
/*      nstore          (sent)                  dimension of a  */
/*      test            (sent)                  div by zero check number*/
/*      ierror          (returned)              zero on no error*/
/*      itriag          (sent)                  matrix triangularized only*/
/*                                               on TRUE useful when solving*/
/*                                               multiple systems with same a */
/*      line            (sent)                  Work array of size n_in. Transmited because this function may be called many times */
/*      isub            (sent)                  		"" */

        static int l1;
        int i = 0, j, k, l, j2, n, nstore;
	int iet, ieb;
        double big, testa, b, sum;

        iet=0;  /* initial error flags, one for triagularization*/
        ieb=0;  /* one for backsolving */
	n = n_in;
	nstore = nstore_in;
	memset(line, 0, n * sizeof(int));
	memset(isub, 0, n * sizeof(int));

	/* triangularize the matrix a*/
	/* replacing the zero elements of the triangularized matrix */
	/* with the coefficients needed to transform the vector vec */

	if (itriag) {   /* triangularize matrix */
 
		for( j = 0; j < n; j++ ) {      /*line is an array of flags*/
                        line[j] = 0; 
                        /* elements of a are not moved during pivoting*/
                        /* line=0 flags unused lines */
		}
                        
                for( j = 0; j < n-1; j++ ) {
                        /*  triangularize matrix by partial pivoting */
			big = 0.0; /* find biggest element in j-th column*/
                                  /* of unused portion of matrix*/
			for( l1 = 0; l1 < n; l1++ ) {
				if( line[l1] == 0 ) {
					testa = (double) fabs((double) (*(a+l1*nstore+j)));
					if (testa > big) {
						i = l1;
						big = testa;
					}
				}
			}
			if( big <= test)	/* test for div by 0 */
				iet = 1;
 
			line[i] = 1;  /* selected unused line becomes used line */
			isub[j] = i;  /* isub points to j-th row of tri. matrix */
 
			sum = 1.0/(*(a+i*nstore+j)); 
                                /*reduce matrix towards triangle */
			for( k = 0; k < n; k++ ) {
				if( line[k] == 0 ) {
					b = (*(a+k*nstore+j))*sum;
					for( l = j+1; l < n; l++ )
						*(a+k*nstore+l) = (*(a+k*nstore+l)) -b*(*(a+i*nstore+l));

					*(a+k*nstore+j) = b;
				}
			}
		}
 
		for( j = 0; j < n; j++ ) {
                        /*find last unused row and set its pointer*/
                        /*  this row contains the apex of the triangle*/
			if( line[j] == 0) {
				l1 = j;		/*apex of triangle*/
				isub[n-1] = j;
				break;
			}
		}
	}
                
        /*start backsolving*/
        
	for( i = 0; i < n; i++ ) {  /* invert pointers. line(i) now gives*/
                                /* row no in triang matrix of i-th row*/
                                /* of actual matrix */
                line[isub[i]] = i;
	}
 
	for( j = 0; j < n-1; j++) {		/*transform the vector to match triang. matrix*/
               b = vec[isub[j]];
               for( k = 0; k < n; k++ ) {
                      if (line[k] > j)		/* skip elements outside of triangle*/
                                vec[k] = vec[k] - (*(a+k*nstore+j))*b;
		}
	}
 
	b = *(a+l1*nstore+(n-1));		/*apex of triangle*/
	if( (fabs( (double) b)) <= test) 	/*check for div by zero in backsolving*/
		ieb = 2;

	vec[isub[n-1]] = vec[isub[n-1]]/b;
 
	for( j = n-2; j >= 0; j-- ) {		/* backsolve rest of triangle*/
                sum=vec[isub[j]];
                for( j2 = j+1; j2 < n; j2++ )
			sum = sum - vec[isub[j2]] * (*(a+isub[j]*nstore+j2));

		b = *(a+isub[j]*nstore+j);
		if( (fabs((double)b))<=test) 	/* test for div by 0 in backsolving */
			ieb = 2;
                vec[isub[j]] = sum / b;		/*solution returned in vec*/
	}

	/*put the solution vector into the proper order*/

	for( i = 0; i < n; i++ ) {		/* reorder solution */
		for( k = i; k < n; k++ ) {	/* search for i-th solution element */
			if( line[k]==i ) {
                                j = k;
                                break;
			}
		}
		b = vec[j];			/* swap solution and pointer elements*/
		vec[j] = vec[i];
		vec[i] = b;
		line[j] = line[i];
	}
 
	*ierror = iet + ieb;			/* set final error flag*/
}

void GMT_cheb_to_pol (double *c, int n, double a, double b, double *d, double *dd) {
	/* Convert from Chebyshev coefficients used on a t =  [-1,+1] interval
	 * to polynomial coefficients on the original x = [a b] interval.
	 * d & dd are work arrays of the same size as c
	 * Modified from Numerical Miracles, ...eh Recipes (Paul Wessel) */
 
	int j, k;
	double sv, cnst, fac;
 
	/* First we generate coefficients for a polynomial in t */
 
	d[0] = c[n-1];
	 
	for (j = n - 2; j >= 1; j--) {
	 	for (k = n - j; k >= 1; k--) {
			sv = d[k];
			d[k] = 2.0 * d[k-1] - dd[k];
			dd[k] = sv;
		}
		sv = d[0];
		d[0] = -dd[0] + c[j];
		dd[0] = sv;
	}
	for (j = n - 1; j >= 1; j--) d[j] = d[j-1] - dd[j];
	/* d[0] = -dd[0] + 0.5 * c[0]; */	/* This is what Num. Rec. says, but we do not do the approx with 0.5 * c[0] */
	d[0] = -dd[0] + c[0];

	/* Next step is to undo the scaling so we can use coefficients with x */

	cnst = fac = 2.0 / (b - a);
	for (j = 1; j < n; j++) {
		d[j] *= fac;
		fac *= cnst;
	}
	cnst = 0.5 * (a + b);
	for (j = 0; j <= n - 2; j++) for (k = n - 2; k >= j; k--) d[k] -= cnst * d[k+1];

	/* Return the new coefficients via c */

	memcpy ((void *)c, (void *)d, (size_t)(n * sizeof (double)));
}

void prepVars(float *out, int n_win, int nx, int ny, int *n_win2, int *nHalfWin, int *m_stop, int *n_stop, int *n_off, int *inc) {
	/* This is common to all algo functions */

	*nHalfWin = n_win / 2;
	*n_win2   = n_win * n_win;	*n_off  = *nHalfWin * (ny + 1);
	*n_stop   = nx - *nHalfWin;	*m_stop = ny - *nHalfWin;
	if (out[0] == 1) *inc = n_win;		/* The Subsample case */
	if (out[0] == 2) *inc = *nHalfWin;	/* The half-overlap case */
}
