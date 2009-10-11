/*
 * Tool to compute moving window block operations
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

#define TRUE	1
#define FALSE	0

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

typedef void (*PFV) ();		/* PFV declares a pointer to a function returning void */

void TPI(float *in, float *out, int n_win, int nx, int ny, int check_nans);
void TRI(float *in, float *out, int n_win, int nx, int ny, int check_nans);
void roughness(float *in, float *out, int n_win, int nx, int ny, int check_nans);
void average(float *in, float *out, int n_win, int nx, int ny, int check_nans);
void block_min(float *in, float *out, int n_win, int nx, int ny, int check_nans);
void block_max(float *in, float *out, int n_win, int nx, int ny, int check_nans);
void callAlgo(int id, float *in, float *out, int n_win, int nx, int ny, int check_nans);

PFV MIR_block_funs[6];

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	n_win = 3, i, j;
	int	algo = 1, nm, nx, ny, nfound = 0, argc = 0, n_arg_no_char = 0;
	int	error = FALSE, unknown_nans = TRUE, check_nans = FALSE;
	char	**argv;
	float	*zdata, *out;

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
	argv[0] = "mirdem";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'A':
					algo = atoi(&argv[i][2]);
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
		mexPrintf ("mirdem - Compute morphological quantities from DEMs single precision arrays\n\n");
		mexPrintf ("usage: out = mirdem(input, ['-A<0|...|5>'], [-N<0|1>], [-W<winsize>]\n");
		
		mexPrintf ("\t<input> is name of input array (singles only)\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t   ---------------------------------\n");
		mexPrintf ("\t-A select algorithm (see Wilson et al 2007, Marine Geodesy 30:3-35):\n");
		mexPrintf ("\t   0 -> Terrain Ruggedness Index\n");
		mexPrintf ("\t   1 -> Topographic Position Index\n");
		mexPrintf ("\t   2 -> Roughness\n");
		mexPrintf ("\t   3 -> Mean\n");
		mexPrintf ("\t   3 -> Minimum\n");
		mexPrintf ("\t   5 -> Mmaximum\n");
		mexPrintf ("\t-N Inform if input has NaNs (1) or not (0) thus avoiding wasting time with repeated test.\n");
		mexPrintf ("\t-W select the rectangular window size [default is 3, which means 3x3].\n");
		mexErrMsgTxt("\n");
	}
	
	if (nlhs == 0) {
		mexPrintf("MIRDEM ERROR: Must provide an output.\n");
		return;
	}

	if (!mxIsSingle(prhs[0]))
		mexErrMsgTxt("MIRDEM ERROR: Invalid input data type. Only valid type is: Single.\n");

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	nm = nx * ny;

	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2)
		mexErrMsgTxt("MIRDEM ERROR: First argument must contain a decent array\n");

	zdata = (float *)mxGetData(prhs[0]);

	plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), mxGetDimensions(prhs[0]), mxSINGLE_CLASS, mxREAL);
	out = (float *)(mxGetData(plhs[0]));

	if (unknown_nans) {
		/* Loop over the file and find if we have NaNs. Stop at first NaN occurence. */
		for (i = 0; i < nm; i++) {
			if (mxIsNaN(zdata[i])) {
				check_nans = TRUE;
				break;
			}
		}
	}

	MIR_block_funs[0] = (PFV) TRI;
	MIR_block_funs[1] = (PFV) TPI;
	MIR_block_funs[2] = (PFV) roughness;
	MIR_block_funs[3] = (PFV) average;
	MIR_block_funs[4] = (PFV) block_min;
	MIR_block_funs[5] = (PFV) block_max;

	callAlgo(algo, zdata, out, n_win, nx, ny, check_nans);
}

void callAlgo(int id, float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	int	j, n, m, nHalfWin, n_extra, nWinExtended, nHalfWinExtended, nny;
	float	*pad_cols, *out_cols;

	MIR_block_funs[id](in, out, n_win, nx, ny, check_nans);

	nHalfWin = n_win / 2;
	n_extra = nHalfWin / 2;		/* Will be 0 for W=3, 1 for W=5, 2 for W=7, etc ... */
	nWinExtended = n_win + n_extra;
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
		for (m = 0; m < nHalfWin; m++) pad_cols[m + n*nny]      = pad_cols[nHalfWin + n*nny];
		for (m = ny+nHalfWin; m < nny; m++) pad_cols[m + n*nny] = pad_cols[ny+nHalfWin-1 + n*nny];
	}
			
	MIR_block_funs[id](pad_cols, out_cols, n_win, nWinExtended, nny, check_nans);

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

	MIR_block_funs[id](pad_cols, out_cols, n_win, nWinExtended, nny, check_nans);

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
	for (n = 0; n <= nHalfWin + n_extra; n++, m++)			/* Copy first (nHalfWin + n_extra + 1) rows */
		for (j = 0; j < nx; j++)
			pad_cols[j * nWinExtended + m] = in[j*ny + n];

	MIR_block_funs[id](pad_cols, out_cols, n_win, nx, nWinExtended, check_nans);

	for (n = nHalfWin-1, m = nHalfWin; n >= 0 ; n--, m++)		/* Put the result in the first nHalfWin rows*/
		for (j = nHalfWin; j < nx - nHalfWin; j++)
			out[j*ny + n] = out_cols[j * nWinExtended + m];

	/* Ending rows */
	for (n = ny - nHalfWinExtended - 1, m = 0; n < ny; n++, m++)	/* Copy first (nHalfWin + n_extra + 1) rows */
		for (j = 0; j < nx; j++)
			pad_cols[j * nWinExtended + m] = in[j*ny + n];
	for (n = ny - 2; n >= ny - (nHalfWinExtended + 1); n--, m++)	/* Mirror nHalfWin rows */
		for (j = 0; j < nx; j++)
			pad_cols[j * nWinExtended + m] = in[j*ny + n];

	MIR_block_funs[id](pad_cols, out_cols, n_win, nx, nWinExtended, check_nans);

	for (n = ny - nHalfWin, m = 2*nHalfWin-1; n < ny; n++, m--)		/* Put the result in the last nHalfWin rows */
		for (j = nHalfWin; j < nx - nHalfWin; j++)
			out[j*ny + n] = out_cols[j * nWinExtended + m];

	mxFree((void *)pad_cols);
	mxFree((void *)out_cols);

	/* Since this bloody thing still f.. on the corners I give up and replace them by their next diagonal neighbor */
	for (n = 0; n < nHalfWin; n++) {
		for (m = 0; m < nHalfWin; m++)           out[m + n*ny] = out[nHalfWin * ny + nHalfWin];	/* First_rows/West */
		for (m = ny - nHalfWin - 1; m < ny; m++) out[m + n*ny] = out[(nHalfWin+1) * ny - nHalfWin];
	}
	for (n = nx - nHalfWin; n < nx; n++) {
		for (m = 0; m < nHalfWin; m++)           out[m + n*ny] = out[(nx - nHalfWin - 1) * ny + nHalfWin];/* First_rows/East */
		for (m = ny - nHalfWin - 1; m < ny; m++) out[m + n*ny] = out[(nx - nHalfWin) * ny - nHalfWin];
	}

}

void TPI(float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes TPI. TPI stands for Topographic Position Index, which is defined as the 
	   difference between a central pixel and the mean of its surrounding cells 
	   (Wilson et al 2007, Marine Geodesy 30:3-35). */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, nNaNs;
	float	*p;
	double	mean = 0;

	nHalfWin = n_win / 2;
	n_win2  = n_win * n_win;	n_off   = nHalfWin * (ny + 1);
	n_stop = nx - nHalfWin;		m_stop = ny - nHalfWin;

	for (n = nHalfWin; n < n_stop; n++) {
		k = (n - nHalfWin) * ny - 1;		/* Index of window's UL corner */
		for (m = nHalfWin; m < m_stop; m++) {
			k++;
			nm = k + n_off;
			if (check_nans && mxIsNaN(in[nm])) {out[nm] = in[nm];	continue;}
			mean = 0;	nNaNs = 0;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && mxIsNaN(*p))	/* Ignore this value */
						nNaNs++;
					else
						mean += *p;
					p++;
				}
			}
			out[nm] = in[nm] - (float) (mean / (n_win2 - nNaNs));
		}
	}
}

void TRI(float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes TRI. TRI stands for Terrain Ruggedness Index, which is defined as the 
	   mean difference between a central pixel and its surrounding cells 
	   (Wilson et al 2007, Marine Geodesy 30:3-35). */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, nNaNs;
	float	*p;
	double	mean = 0, tmp;

	nHalfWin = n_win / 2;
	n_win2  = n_win * n_win;	n_off   = nHalfWin * (ny + 1);
	n_stop = nx - nHalfWin;		m_stop = ny - nHalfWin;

	for (n = nHalfWin; n < n_stop; n++) {
		k = (n - nHalfWin) * ny - 1;
		for (m = nHalfWin; m < m_stop; m++) {
			k++;
			nm = k + n_off;
			if (check_nans && mxIsNaN(in[nm])) {out[nm] = in[nm];	continue;}
			tmp = in[nm];
			mean = 0;	nNaNs = 0;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && mxIsNaN(*p))	/* Ignore this value */
						nNaNs++;
					else
						mean += fabs(tmp - *p);
					p++;
				}
			}
			out[nm] = (float) (mean / (n_win2 - nNaNs));
		}
	}
}

void roughness(float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes Roughness defined as the the largest inter-cell difference of a central pixel and its surrounding cell
	   (Wilson et al 2007, Marine Geodesy 30:3-35). */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off;
	float	*p, tmp, min, max;

	nHalfWin = n_win / 2;
	n_win2  = n_win * n_win;	n_off  = nHalfWin * (ny + 1);
	n_stop = nx - nHalfWin;		m_stop = ny - nHalfWin;

	for (n = nHalfWin; n < n_stop; n++) {
		k = (n - nHalfWin) * ny - 1;
		for (m = nHalfWin; m < m_stop; m++) {
			k++;
			nm = k + n_off;
			if (check_nans && mxIsNaN(in[nm])) {out[nm] = in[nm];	continue;}
			tmp = in[nm];
			min = max = tmp;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && mxIsNaN(*p)) continue;	/* Ignore this value */
					max = MAX(max, *p);
					min = MIN(min, *p);
					p++;
				}
			}
			out[nm] = (max - min);
		}
	}
}

void average(float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes the mean of cells  inside rectangular window */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off, nNaNs;
	float	*p;
	double	mean = 0;

	nHalfWin = n_win / 2;
	n_win2  = n_win * n_win;	n_off   = nHalfWin * (ny + 1);
	n_stop = nx - nHalfWin;		m_stop = ny - nHalfWin;

	for (n = nHalfWin; n < n_stop; n++) {
		k = (n - nHalfWin) * ny - 1;		/* Index of window's UL corner */
		for (m = nHalfWin; m < m_stop; m++) {
			k++;
			nm = k + n_off;
			if (check_nans && mxIsNaN(in[nm])) {out[nm] = in[nm];	continue;}
			mean = 0;	nNaNs = 0;
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && mxIsNaN(*p))	/* Ignore this value */
						nNaNs++;
					else
						mean += *p;
					p++;
				}
			}
			out[nm] = (float) (mean / (n_win2 - nNaNs));
		}
	}
}

void block_min(float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes the minimum of cells  inside rectangular window */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off;
	float	*p, min;

	nHalfWin = n_win / 2;
	n_win2  = n_win * n_win;	n_off  = nHalfWin * (ny + 1);
	n_stop = nx - nHalfWin;		m_stop = ny - nHalfWin;

	for (n = nHalfWin; n < n_stop; n++) {
		k = (n - nHalfWin) * ny - 1;
		for (m = nHalfWin; m < m_stop; m++) {
			k++;
			nm = k + n_off;
			if (check_nans && mxIsNaN(in[nm])) {out[nm] = in[nm];	continue;}
			min = in[nm];
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && mxIsNaN(*p)) continue;	/* Ignore this value */
					min = MIN(min, *p);
					p++;
				}
			}
			out[nm] = min;
		}
	}
}

void block_max(float *in, float *out, int n_win, int nx, int ny, int check_nans) {
	/* Computes the minimum of cells  inside rectangular window */
	int	n_win2, nHalfWin, n_stop, m_stop, n, m, nm, i, j, k, n_off;
	float	*p, max;

	nHalfWin = n_win / 2;
	n_win2  = n_win * n_win;	n_off  = nHalfWin * (ny + 1);
	n_stop = nx - nHalfWin;		m_stop = ny - nHalfWin;

	for (n = nHalfWin; n < n_stop; n++) {
		k = (n - nHalfWin) * ny - 1;
		for (m = nHalfWin; m < m_stop; m++) {
			k++;
			nm = k + n_off;
			if (check_nans && mxIsNaN(in[nm])) {out[nm] = in[nm];	continue;}
			max = in[nm];
			for (i = 0; i < n_win; i++) {		/* Loop columns inside window */
				p = &in[k + i * ny];
				for (j = 0; j < n_win; j++) {
					if (check_nans && mxIsNaN(*p)) continue;	/* Ignore this value */
					max = MAX(max, *p);
					p++;
				}
			}
			out[nm] = max;
		}
	}
}
