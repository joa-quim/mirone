/* ********************************************************************* */
/* Program for automatic extraction of ridge and valley axes from the */
/* digital elevation data set. The main steps are: */
/* step 1: data and parameters input by sub._input */
/* step 2: TarGet ReCognition and CONnection by sub._TGRCON */
/* step 3: SEGment checK-Out by sub._SEGKO */
/* step 4: LiNe SMooth and OutPut by sub._LNSMOP
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
/* ********************************************************************* */

#include "gmt.h"
#include "mex.h"

#ifndef rint
#define rint(x) (floor((x)+0.5))
#endif
#ifndef irint
#define irint(x) ((int)rint(x))
#endif

#define ij(i,j) ((i) + ((j)-1)*n_cols -1)
#define ijk(i,j,k) ((i) + ((j)-1)*n_cols + ((k)-1)*n_cols*n_rows -1)

/* Table of constant values */
static int c__0 = 0;
static int c__1 = 1;
static int c__2 = 2;

int segko(), input(), tgrcon();
int smooth(int i, int j, float *x, float *y, float *w);
int kst(int i, int j, int k, int jc);
int neb(int i);
int no_sys_mem (char *where, int n);
void con(int i, int j, int k, int *n, int jc, float *s);
void grd_FLIPUD (float data[], int nx, int ny);

int verbose = FALSE, first = TRUE;
int	ic = 1, ib = 3, n_cols, n_rows;
int	neigh_x[8] = {0, 1, 1, 1, 0, -1, -1, -1};
int	neigh_y[8] = {1, 1, 0, -1, -1, -1, 0, 1};
float	x_inc, y_inc, x_min, y_min, z_min, z_max, z_scale, *data, *w;
char	*v;
/*int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int error = FALSE, global = FALSE;
	int is_double = FALSE, is_single = FALSE, is_int32 = FALSE, is_int16 = FALSE;
	int is_uint16 = FALSE;
	int	argc = 0, n_arg_no_char = 0, nc_h, nr_h, i2, *i_4;
        int	i, j, k, n, ic, nx, ny,  mx;
	int	nx_new, ny_new, one_or_zero, ndatac, p_alloc, bytes_to_copy;
	short int *i_2;
	unsigned short int *ui_2;
	char    *infile = NULL, **argv;
	float	*z_4, ww, x, y, v1, v2, x1, y1, heig, x_save, y_save;
	double	w_new = 0.0, e_new = 0.0, s_new = 0.0, n_new = 0.0, nan, *p_out, *pdata, *z_8, *head;
	struct	GRD_HEADER h;

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
	argv[0] = "grdppa";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);*/
	argc = GMT_begin (argc, argv);

        for (i = 1; !error && i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
				case 'R':
                                case '\0':
					error += GMT_get_common_args (argv[i], &w_new, &e_new, &s_new, &n_new);
                                        break;
				case 'L':
					sscanf (&argv[i][2], "%d", &ib);
					break;
				case 'T':
					sscanf (&argv[i][2], "%d", &ic);
					break;
				case 'V':
					verbose = TRUE;
					break;
				default:
                                        error = TRUE;
					GMT_default_error (argv[i][1]);
                                        break;
                        }
                }
                else
                        infile = argv[i];
        }

	if (n_arg_no_char < 2 || error) {
		mexPrintf ("grdppa %s - automatic extraction of ridge or valley axes\n\n", GMT_VERSION);
		mexPrintf ("usage: out = grdppa(infile, head, '[-L<npoints>]', '[-T<1|2>]', '[-Rw/e/s/n]', '[-V]')\n\n");
		
		mexPrintf ("\t<infile> is name of input array\n");
		mexPrintf ("\t<head> is array header descriptor of the form\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-L npoints for polygon recognition [default = 3]\n");
		mexPrintf ("\t-T 1 [default] for ridges, 2 for valleys\n");
		return;
	}
	if (error) return;

	if (nlhs == 0) {
		mexPrintf("GRD2CDF ERROR: Must provide an output.\n");
		return;
	}

	/* Find out in which data type was given the input array */
	if (mxIsDouble(prhs[0])) {
		z_8 = mxGetPr(prhs[0]);
		is_double = TRUE;
	}
	else if (mxIsSingle(prhs[0])) {
		z_4 = mxGetData(prhs[0]);
		is_single = TRUE;
	}
	else if (mxIsInt32(prhs[0])) {
		i_4 = mxGetData(prhs[0]);
		is_int32 = TRUE;
	}
	else if (mxIsInt16(prhs[0])) {
		i_2 = mxGetData(prhs[0]);
		is_int16 = TRUE;
	}
	else if (mxIsUint16(prhs[0])) {
		ui_2 = mxGetData(prhs[0]);
		is_uint16 = TRUE;
	}
	else {
		mexPrintf("GRDPPA ERROR: Unknown input data type.\n");
		mexErrMsgTxt("Valid types are:double, single, Int32, Int16 and UInt16.\n");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 2 || nx < 2)
		mexErrMsgTxt("GRDPPA ERROR: First argument must contain a decent array\n");

	nc_h = mxGetN (prhs[1]);
	nr_h = mxGetM (prhs[1]);
	if (!mxIsNumeric(prhs[1]) || nr_h > 1 || nc_h < 9)
		mexErrMsgTxt("GRDPPA ERROR: Second argument must contain a valid header of the input array.\n");

	head  = mxGetPr(prhs[1]);		/* Get header info */

	h.x_min = head[0];	h.x_max = head[1];
	h.y_min = head[2];	h.y_max = head[3];
	h.z_min = head[4];	h.z_max = head[5];
	h.x_inc = head[7];	h.y_inc = head[8];
	h.nx = nx;		h.ny = ny;
	h.node_offset = irint(head[6]);
	mx = nx + 2;
	
	/* If -R was not used */
	if (!project_info.region_supplied) {
		w_new = h.x_min;	e_new = h.x_max;
		s_new = h.y_min;	n_new = h.y_max;
	}

	if (s_new < h.y_min || s_new > h.y_max) error = TRUE;
	if (n_new < h.y_min || n_new > h.y_max) error = TRUE;
	global = (fabs (h.x_max - h.x_min) == 360.0);
		
	if ( !global && ((w_new < h.x_min) || (e_new > h.x_max)) ) error = TRUE;
	if (error) {
		mexPrintf ("%s: Subset exceeds data domain!\n", GMT_program);
		return;
	}

	/* Check if new wesn differs from old wesn by integer dx/dy */

	if (GMT_minmaxinc_verify (h.x_min, w_new, h.x_inc, GMT_SMALL) == 1) {
		mexPrintf ("%s: Old and new x_min do not differ by N * dx\n", GMT_program);
		return;
	}
	if (GMT_minmaxinc_verify (e_new, h.x_max, h.x_inc, GMT_SMALL) == 1) {
		mexPrintf ("%s: Old and new x_max do not differ by N * dx\n", GMT_program);
		return;
	}
	if (GMT_minmaxinc_verify (h.y_min, s_new, h.y_inc, GMT_SMALL) == 1) {
		mexPrintf ("%s: Old and new y_min do not differ by N * dy\n", GMT_program);
		return;
	}
	if (GMT_minmaxinc_verify (n_new, h.y_max, h.y_inc, GMT_SMALL) == 1) {
		mexPrintf ("%s: Old and new y_max do not differ by N * dy\n", GMT_program);
		return;
	}
	
       	/*GMT_grd_init (&h, argc, argv, TRUE);*/

	data = mxCalloc ((nx+2)*(ny+2), sizeof (float));

	/* Transpose from Matlab orientation to gmt grd orientation */
	if (is_double) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + mx + 1] = (float)z_8[j*ny+i];
	}
	else if (is_single) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + mx + 1] = z_4[j*ny+i];
	}
	else if (is_int32) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + mx + 1] = (float)i_4[j*ny+i];
	}
	else if (is_int16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + mx + 1] = (float)i_2[j*ny+i];
	}
	else if (is_uint16) {
		for (i = 0, i2 = ny - 1; i < ny; i++, i2--) 
			for (j = 0; j < nx; j++) data[i2*mx + j + mx + 1] = (float)ui_2[j*ny+i];
	}

	one_or_zero = (h.node_offset) ? 0 : 1;
	nx_new = irint ((e_new - w_new) / h.x_inc) + one_or_zero;
	ny_new = irint ((n_new - s_new) / h.y_inc) + one_or_zero;
	
        ndatac = (nx_new + 2) * (ny_new + 2);
	n_rows = ny_new + 2;		n_cols = nx_new + 2;
	x_inc = (float)h.x_inc;		y_inc = (float)h.y_inc;
	x_min = (float)w_new;		y_min = (float)s_new;
	z_min = (float)h.z_min;		z_max = (float)h.z_max;
	z_scale = 499 / (z_max - z_min);

	/* Regarding the original fortran algo we need to flipud the grid */
	grd_FLIPUD(data,n_cols,n_rows);

	if ((w = (float *) calloc ((size_t)(n_rows*n_cols*4), sizeof(float)) ) == NULL) 
		{if (no_sys_mem("ppa --> (w)", n_rows*n_cols*4))	return;}
	if ((v = (char *) calloc ((size_t)(n_rows*n_cols*4), sizeof(char)) ) == NULL) 
		{if (no_sys_mem("ppa --> (v)", n_rows*n_cols*4))	return;} 

	tgrcon();
	segko();

	/* This was the lnsmop() routine */
	p_alloc = 5000;
	p_out = mxCalloc (p_alloc, sizeof (double));
	n = 0;
	nan = mxGetNaN();
	for (i = 2; i <= n_cols - 1; ++i) {
		for (j = 2; j <= n_rows - 1; ++j) {
			for (k = 1; k <= 4; ++k) {
				if (kst(i, j, k, 1) == 0) continue;
				smooth(i, j, &x, &y, &v1);
				smooth(i + neigh_x[k-1], j + neigh_y[k-1], &x1, &y1, &v2);
				ww = (v1 + v2) / 2;
				if (ic == 1) heig = z_min + (ww - 1) / z_scale;
				else	heig = z_max - (ww - 1) / z_scale;

				if (!first) {	/* Atempt to reduce the number of segments */
					if (fabs(x - x_save) < 1e-5 && fabs(y - y_save) < 1e-5) {
						p_out[n] = x1;	p_out[n+1] = y1; n += 2;
					}
					else if (fabs(x1 - x_save) < 1e-5 && fabs(y1 - y_save) < 1e-5) {
						p_out[n] = x;	p_out[n+1] = y; n += 2;
					}
					else {
						p_out[n] = nan;	p_out[n+1] = nan; n += 2;
						p_out[n] = x;	p_out[n+1] = y; n += 2;
						p_out[n] = x1;	p_out[n+1] = y1; n += 2;
					}
					x_save = x1;	y_save = y1;
				}
				else {
					p_out[n] = nan;	p_out[n+1] = nan; n += 2;
					p_out[n] = x;	p_out[n+1] = y; n += 2;
					p_out[n] = x1;	p_out[n+1] = y1; n += 2;
					x_save = x1;	y_save = y1;
					first = FALSE;
				}
				if (p_alloc <= (n + 6)) {
					p_alloc += 5000;
					if ((p_out = mxRealloc(p_out, p_alloc * sizeof(double))) == 0) {
						mexErrMsgTxt("GRDPPA ERROR: Could not reallocate memory\n");
					}
				}
			}
		}
	}

	/* Create and populate matrices for the return array(s) */
	plhs[0] = mxCreateDoubleMatrix (2,(int)(n/2), mxREAL);
	pdata = mxGetPr(plhs[0]);
	bytes_to_copy = n * 8;
	memcpy(pdata, p_out, bytes_to_copy);

	free((void *) w);
	free((void *) v);
	GMT_end (argc, argv);
}


int tgrcon() {
	/* ********************************************************************* */
	/* Subroutine for target recognition and connection */
	/* ********************************************************************* */
	int *b, i, j, k, n[8], mc, ii, jj, ml, mm, itg;

	if ((b = (int *) calloc ((size_t)(n_cols*n_rows), sizeof(int)) ) == NULL) 
		{no_sys_mem("ppa --> (b)", 400*400);} 

	/* # recognize the targets along profiles in four direction */
	itg = 0;
	for (i = 3; i <= n_cols - 2; ++i) {
		for (j = 2; j <= n_rows - 2; ++j) {
			for (k = 1; k <= 8; ++k) {
				n[k-1] = 0;
				for (ml = 1; ml <= ib / 2; ++ml) {
					ii = i + ml * neigh_x[k-1];
					jj = j + ml * neigh_y[k-1];
					if (ii < 2 || ii > n_cols - 1 || jj < 2 || jj > n_rows - 1) continue;
					if (data[ij(ii,jj)] < data[ij(i,j)]) n[k-1] = 1;
				}
			}
			for (k = 1; k <= 4; ++k) {
				if (n[k-1] + n[k+3] > 1) b[ij(i,j)] = 1;
			}
			if (b[ij(i,j)] == 1) ++itg;
		}
	}
	if (verbose) mexPrintf(" %d\ttargets found\n", itg);
	/* # target connection */
	for (i = 2; i <= n_cols - 1; ++i) {
		for (j = 2; j <= n_rows - 1; ++j) {
			for (k = 1; k <= 4; ++k) {
				mm = kst(i, j, k, 4);
				if (b[ij(i,j)] + b[ij(i + neigh_x[k-1],(j + neigh_y[k-1]))] == 2) {
					mc = kst(i, j, k, 2);
				}
			}
		}
	}
	free((void *) b);
	return 0;
}

int segko() {
	/* ********************************************************************* */
	/* check-out improper segments by polygon breaking and branch reduction */
	/* ********************************************************************* */
	float v, z, wn;
	int id, ii, jj, kk, in, jn, mm, nv, i, j, k, l, m, n, *b;

	if ((b = (int *) calloc ((size_t)(n_cols*n_rows), sizeof(int)) ) == NULL) 
		{no_sys_mem("ppa --> (b)", 400*400);} 

	if (verbose) mexPrintf (" polygon breaking ...\n");

	m = 0;
	/* # pick the weakest segment */
L1:
	++m;
	wn = 1001.;
	for (i = 2; i <= n_cols - 1; ++i) {
		for (j = 2; j <= n_rows - 1; ++j) {
			if (b[ij(i,j)] == 1)  continue;
			nv = 0;
			for (k = 1; k <= 4; ++k) {
				con(i, j, k, &n, 1, &v);
				if (v == 2e3) ++nv;
				if (v >= wn) continue;
				/* # skip the end-segment */
				if (kst(i,j,k,-8) == 1 || kst(i + neigh_x[k-1], j + neigh_y[k-1], k, -8) == 1) {
					mm = kst(i, j, k, -4);
					continue;
				}
				wn = v;
				ii = i;
				jj = j;
				kk = k;
			}
			if (nv == 4) b[ij(i,j)] = 1;
		}
	}
	if (wn == 1001.) goto L4;
	if (verbose && m % 100 == 1) {
		if (ic == 1) {
			z = z_min + (wn / 2 - 1) / z_scale;
			mexPrintf (" segments below z = %.3f\tchecked %f\t%d\r", z,wn,m);
		} else {
			z = z_max - (wn / 2 - 1) / z_scale;
			mexPrintf (" segments above z = %.3f\tchecked\r", z);
		}
	}
	/* # polygon tracing */
	mm = kst(ii, jj, kk, -4);
	if (kst(ii, jj, kk, -8) == 0) goto L1;
	in = ii + neigh_x[kk-1];
	jn = jj + neigh_y[kk-1];
	if (kst(in, jn, kk, -8) == 0) goto L1;
	id = 1;
L6:
	i = in;		j = jn;		k = kk;
L8:
	k = neb(k+4);
L11:
	k = neb(k+id);
	if (kst(i, j, k, -1) == 0) goto L11;
	i += neigh_x[k-1];
	j += neigh_y[k-1];
	if (i == ii && j == jj) {
		mm = kst(ii, jj, kk, 4);
		goto L1;
	}
	if (i == in && j == jn) {
		if (id == -1) goto L1;
		id = -1;
		goto L6;
	}
	goto L8;
	/* # branch reduction */
L4:
	for (l = 1; l <= ib / 2; ++l) {
		for (i = 2; i <= n_cols - 1; ++i) {
			for (j = 2; j <= n_rows - 1; ++j) {
				b[ij(i,j)] = 0;
				if (kst(i, j, k, 8) != 1) continue;
				for (k = 1; k <= 8; ++k) {
					if (kst(i, j, k, 1) == 1) b[ij(i,j)] = k;
				}
			}
		}
		for (i = 2; i <= n_cols - 1; ++i) {
			for (j = 2; j <= n_rows - 1; ++j) {
				mm = kst(i, j, b[ij(i,j)], 4);
			}
		}
	}
	free((void *) b);
	return 0;
}

int smooth(int i, int j, float *x, float *y, float *w) {
/* ********************************************************************* */
/* subroutine to smooth target position according to the weights of */
/* connected neighbors */
/* ********************************************************************* */
	static float f;
	static int k;

	*w = data[ij(i,j)];
	*x = 0.;
	*y = 0.;
	for (k = 1; k <= 8; ++k) {
		f = kst(i, j, k, 1);
		*w += data[ij(i + neigh_x[k-1], j + neigh_y[k-1])] * f;
		*x += neigh_x[k-1] * data[ij(i + neigh_x[k-1], j + neigh_y[k-1])] * f;
		*y += neigh_y[k-1] * data[ij(i + neigh_x[k-1], j + neigh_y[k-1])] * f;
	}
	*x = x_min + ((float)(i - 2) + *x / *w) * x_inc;
	*y = y_min + ((float)(j - 2) + *y / *w) * y_inc;
	*w /= (kst(i, j, k, 8) + 1);
	return 0;
}

int kst(int i, int j, int k, int jc) {
	/* ********************************************************************* */
	/* function to handle connection status and tracing route tables */
	/* jc=1: check the connection status; jc=-1: check the route status */
	/* jc=2: target connection, crossed segment will be checked out */
	/* jc=4: break the connection,        jc=-4: break the route */
	/* jc=8: check number of connections, jc=-8: check number of routes */
	/* ********************************************************************* */
	int ret_val = 0, m, n;
	static float v, w1, w2;

	if (jc == 2) {
		con(i, j, k, &c__2, 2, &v);
		if (k % 2 == 1) return ret_val;
		con(i + neigh_x[k-2], j + neigh_y[k-2], neb(k+2), &n, 1, &w2);
		if (n == 0) return ret_val;
		con(i, j, k, &n, 1, &w1);
		if (w2 > w1)
			con(i, j, k, &c__0, 2, &v);
		else
			con(i + neigh_x[k-2], j + neigh_y[k-2], neb(k+2), &c__0, 2, &v);
		return ret_val;
	}
	ret_val = 0;
	if (abs(jc) == 8) {
		for (m = 1; m <= 8; ++m) {
			con(i, j, m, &n, 1, &v);
			if (jc == 8 && n != 0) ++ret_val;
			if (jc == -8 && n == 2) ++ret_val;
		}
		return ret_val;
	}
	con(i, j, k, &n, 1, &v);
	if (jc == 1 && n != 0) ret_val = 1;
	if (jc == -1 && n == 2) ret_val = 1;
	if (jc == 4) con(i, j, k, &c__0, 2, &v);
	if (jc == -4) con(i, j, k, &c__1, 2, &v);
	return ret_val;
}

void con(int i, int j, int k, int *n, int jc, float *s) {
	/* ********************************************************************* */
	/* check or change the status of connection or route */
	/* ********************************************************************* */
	int ii, jj, kk;

	if (k <= 0) return;
	ii = i;		jj = j;		kk = k;
	if (k > 4) {
		ii = i + neigh_x[k-1];
		jj = j + neigh_y[k-1];
		kk += -4;
	}
	if (jc == 1) {
		*n = v[ijk(ii, jj, kk)];
		*s = w[ijk(ii, jj, kk)];
		return;
	}
	v[ijk(ii, jj, kk)] = (char) (*n);
	w[ijk(ii, jj, kk)] = 2e3;
	if (*n == 2) w[ijk(ii,jj,kk)] = data[ij(i,j)] + data[ij(i + neigh_x[k-1], (j + neigh_y[k-1]))];
	return;
}

int neb(int i) {
	/* ********************************************************************* */
	/* function to scale the cycled neighbor order */
	/* ********************************************************************* */
	int ret_val;
	ret_val = i;
	if (ret_val > 8) ret_val += -8;
	if (ret_val < 1) ret_val += 8;
	return ret_val;
}

int	no_sys_mem (char *where, int n) {	
	mexPrintf ("Fatal Error: %s could not allocate memory, n = %d\n", where, n);
	return (1);
}


void grd_FLIPUD (float data[], int nx, int ny) {
	int i, j, k, ny1, ny_half;

	/* Reverse order of all columns */
	ny_half = ny / 2;
	ny1 = ny - 1;
	for (i = 0; i < nx; i++) {
		for (j = 0, k = ny1; j < ny_half; j++, k--) {	/* Do this to all rows */
			f_swap (data[j*nx+i], data[k*nx+i]);
		}
	}
}
