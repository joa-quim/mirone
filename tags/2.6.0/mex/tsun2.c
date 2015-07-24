/*   ************** TUNAMI-N2 ****************
*      Tohoku University's Numerical-Analysis Model
*         for Investigation of tsunami
*           Near-field Tsunami version

*            WITH SHALLOW WATER THEORY
*     		INCLUDING EFFECTS OF RUNUP

*                1991.1.25
*                   BY

*         F.IMAMURA, TOHOKU UNIV., JAPAN */

/*       Z; WAVE SURFACE LEVEL    M,N; WATER DISCHARGE
 *   HM,HN; STILL WATER DEPTH AT POINT OF WATER DISCHARGE
 *   DM,DN; TOTAL WATER DEPTH AT POINT OF WATER DISCHARGE
 *      HZ; STILL WATER DEPTH      ZD; TOTAL WATER DEPTH
 *       G; GRAVITATIONAL ACCES    KL; TOTAL TIME STEP
 *      ZM; MAXIMUM WATER LEVEL */

/* 
 * Tradutor:	J. Luis
 * Date:	10 May, 2002
 * Revision:	16/04/2004 
 *          	28/06/2004 	Numero de maregs estava contado em mais 1 na func bnc
 *		23/08/2004	Retuched the indexes in bnc_n function.
 *				Finally programed the North border option
 *
 *		12/01/2006	Updated nlmass function with (hopefully) corrected end indexes
 *				Needs to finish the option to output maregraphs as well.
 *		08/10/2007	grids are numbered from time in sconds and not cycle number as before
 *				Friction coefficient may be passed as argument (-F option)
 *				Accepts output maregraphs via numeric argin and -O option.
 *
 *	version WITH waitbar
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "mex.h"

#define	FALSE		0
#define	TRUE		1
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif
#define	CHUNK 		2000

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef rint
#define rint(x) (floor((x)+0.5))
#endif
#ifndef irint
#define irint(x) ((int)rint(x))
#endif

/* Safe math macros that check arguments */

#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))
#define d_acos(x) (fabs (x) >= 1.0 ? ((x) < 0.0 ? M_PI : 0.0) : acos (x))
#define d_asin(x) (fabs (x) >= 1.0 ? copysign (M_PI_2, (x)) : asin (x))
#define d_atan2(y,x) ((x) == 0.0 && (y) == 0.0 ? 0.0 : atan2 (y, x))
#define d_atn(y,x) ((x) == 0.0 && (y) == 0.0 ? 0.2 : atan2 (y, x))

#define ij2(i,j,k) (((i) + (j)*nx) + (k) * nx * ny)
#define ij0(i,j) ((i) + (j)*nx)
#define ijb(i,j) ((i)*n_mareg + (j))

#define ij(i,j,k) (((i)*11 + (j)) + (k)*10*11)

struct srf_header {		/* Surfer file header structure */
	char id[4];		/* ASCII Binary identifier (DSAA/DSBB) */
	short int nx;		/* Number of columns */
	short int ny;		/* Number of rows */
	double x_min;		/* Minimum x coordinate */
	double x_max;		/* Maximum x coordinate */
	double y_min;		/* Minimum y coordinate */
	double y_max;		/* Maximum y coordinate */
	double z_min;		/* Minimum z value */
	double z_max;		/* Maximum z value */
};

void no_sys_mem (char *where, int n);
void intl(int nx, int ny, double *z, double *m, double *n, double *d, double *h); 
void hmn(int nx, int ny, double *hz, double *hm, double *hn);
void nlmass(int nx, int ny, double *z, double *m, double *n, double *dz, double *hz, double r, int kk); 
void nlmmt(int nx, int ny, double *z, double *m, double *n, double *dz, double *dm, double *dn, double *hz, double *hm, double *hn, double r, double dt, double fm); 
void change(int nx, int ny, double *z, double *m, double *n, double *d); 
int data(char *alt_ond, int *n_tstep_in, double *dt, float time_jump);
void out1(int nx, int ny, double *z, double *m, double *n, double *dm, double *dn, double *h); 
void max_z (int nx, int ny, double *z, double *zm); 
int write_grd_bin (int kk, double x_min, double y_min, double dx, int nx, int ny, float *work);
int intp_lin (double *x, double *y, int n, int m, double *u, double *v, int mode);
int count_col (char *line);
int read_index(char *file, int *nl_w, int *nl_e, int *nc_s, int *nc_n);
void bnc_n(int nx, int ny, double *z, int kk, int nl_w, int nl_e, int nc_s, int nc_n);


int  *ip_s, *ip_n, *jp_w, *jp_e, n_mareg;
double GX = 1e-10;
double GY = 1e-10;
double GG = 9.8;
double *b;			/* pointer to maregraph heights */
double *x_w, *y_w, *u_w, *v_w;	/* working variables for using in intp_lin */
char	stem[80];
struct	srf_header hdr;


/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	step_t = 60;	/* Time step at which grids are writen */
	int	n_of_cycles = 0;	/* Total time step */
	int	argc, n_arg_no_char = 0, nx, ny, ny_tmp, dims[3], n_frames = 0;
	int	nl_w, nl_e, nc_s, nc_n;	/* Number of possible elements on border index */
	int	i, i2, k, kk, j, ncl, n_tstep_in;
	int	error = FALSE, surf_level = TRUE, max_level = FALSE;
	int	maregs_out = FALSE, water_depth = FALSE;
	int	n_mareg_out, *lcum_p;
	float	*work, *ptr_mov_32, *mov_32, time_jump = 0.;
	double	dt = 1.;	/* Time step increment */
	double	r, *head, *b_tmp, *outLoc, dx, time_h = 0.0, amp_fact = 1.;
	double	*z, *m, *n, *hz, *hz1, *dz, *dm, *dn, *hm, *hn, *zm;
	double	*ptr, *h_bar, tmp_ptr[1];	/* Pointers to be used in the waitbar */
	double	manning_coeff = 0.025;

	char	*bathy = NULL;  /* Name pointer for bathymetry file */
	char 	*wave_height = NULL;	/* Name pointer for wave heights file */
	char 	*mareg_index = NULL;	/* Name pointer for maregraph index file */
	char 	*maregs_out_hgt = NULL;	/* Name pointer for output maregraph wave height file */
	char	**argv, w_bar_title[] = "Aguenta ai";
	int	params_in_input = FALSE, bat_in_input = FALSE;
	int	write_grids = FALSE, movie = FALSE, movie_char = FALSE, movie_float = FALSE;
	int	maregs_in_input = FALSE;
	mxArray *rhs[2], *lhs[1];

	int ip[7] = {0,20,40,60,80,100,120};
	int jp[7] = {0,0,0,0,0,0,0};
	double bt[10] = {5.26, 6.7, 8., 5.5, 6., 10., 10., 10., 10., 10.};
	double pz[7];
	FILE	*fp, *fpOutMaregs;

	wave_height = "alt_onda";
	mareg_index = "tsun2.par";

	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}

	if (!(n_arg_no_char > 0 && n_arg_no_char <= 4)) {		/* Not all combinations are tested */
		mexErrMsgTxt("Wrong number of numeric inputs.");
	}

	/* 1st is the bat array */ 
	hz1 = mxGetPr(prhs[0]);
	ny = mxGetM (prhs[0]);
	nx = mxGetN (prhs[0]);
	if (ny < 50 || nx < 50) 	/* This still might fail for maregraphs as 1 arg */
		mexErrMsgTxt("First non char argument must contain a decent bathymetry array\n");


	/* 2ndt is the bathymetry header info array */ 
	head  = mxGetPr(prhs[1]);
	ny_tmp = mxGetM (prhs[1]);
	if (ny_tmp > 1)
		mexErrMsgTxt("Second non char argument must contain the header of the bathymetry array\n");

	hdr.x_min = head[0];		hdr.x_max = head[1];
	hdr.y_min = head[2];		hdr.y_max = head[3];
	hdr.z_min = head[4];		hdr.z_max = head[5];
	hdr.nx = nx;			hdr.ny = ny;
	dx = head[8];		/* Square grid cells are assumed */
	bat_in_input = TRUE;

	/*  */ 
	if (n_arg_no_char >= 3) {	/* The 3th argument must be the maregraph array */
		b_tmp = mxGetPr(prhs[2]);
		n_tstep_in = mxGetM (prhs[2]);
		n_mareg = mxGetN (prhs[2]);
		if (n_tstep_in < 10)
			mexErrMsgTxt("Maregraph series too short. It can't be true.\n");

		/* Transpose from Matlab orientation to scanline orientation */
		b = (double *)mxCalloc(n_mareg*n_tstep_in, sizeof(double));
		for (j = 0; j < n_tstep_in; j++)
			for (i = 0; i < n_mareg; i++) b[j*n_mareg + i] = b_tmp[i*n_tstep_in + j];

		maregs_in_input = TRUE;

	}

	if (n_arg_no_char >= 4) {	/* The 4th arg must be the optional output maregragh array */
		double	x, y, x_tmp, y_tmp, x_inc, y_inc;
		int	ix, jy, k;
		n_mareg_out = mxGetM (prhs[3]);
		lcum_p = (int *) mxCalloc ((size_t)(n_mareg_out), sizeof(int));
		outLoc = mxGetPr(prhs[3]);

		x_inc = (hdr.x_max - hdr.x_min) / (hdr.nx - 1);
		y_inc = (hdr.y_max - hdr.y_min) / (hdr.ny - 1);

		mexPrintf("Adjusted maregraph positions\n");
		mexPrintf("old x\t\told y\t\tnew x\t\tnew y\t\tdepth\n");
		for (i = 0; i < n_mareg_out; i++) {
			x = outLoc[i];		y = outLoc[i + n_mareg_out];
			ix = irint((x - hdr.x_min) / x_inc);
			jy = irint((y - hdr.y_min) / y_inc); 
			lcum_p[i] = jy * hdr.nx + ix; 
			x_tmp = hdr.x_min + x_inc * ix;	/* Adjusted x maregraph pos */
			y_tmp = hdr.y_min + y_inc * jy;	/* Adjusted y maregraph pos */
			k = ix*hdr.ny + jy;		/* The column-major row-major story */
			mexPrintf("%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",x,y,x_tmp,y_tmp,hz1[k]);

		}
		maregs_out = TRUE;
	}

	if (n_arg_no_char == 5) {	/* The 4th argument must be the parameter array */
		/* Very complicated to decode. Not used yet. */
	}

	/* get the length of the input string */
	argv = (char **)mxCalloc(argc, sizeof(char *));
	for (i = 0; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char]);

	for (i = 0; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'A':	/* Amplify wave heights by this factor */
					sscanf (&argv[i][2], "%lf", &amp_fact);
					break;
				case 'f':	/* Movie */
					movie = TRUE;
					movie_float = TRUE;
					break;
				case 'm':	/* Movie */
					movie = TRUE;
					break;
				case 'D':
					water_depth = TRUE;
					surf_level = FALSE;
					max_level = FALSE;
					break;
				case 'F':	/* friction coefficient */
					manning_coeff = atof(&argv[i][2]);
					break;
				case 'I':
					sscanf (&argv[i][2], "%lf/%lf", &dx, &dt);
					break;
				case 'J':
					sscanf (&argv[i][2], "%f", &time_jump);
					break;
				case 'G':	/* Write grids at step_t (see -T) intervals */
					write_grids = TRUE;
					strcpy (stem, &argv[i][2]);
					break;
				case 'M':
					max_level = TRUE;
					water_depth = FALSE;
					surf_level = FALSE;
					break;
				case 'N':	/* Number of cycles to compute */
					n_of_cycles = atoi(&argv[i][2]);
					break;
				case 'O':	/* File name for output maregraph data */
					maregs_out_hgt = &argv[i][2];
					break;
				case 'P':	/* Name of params file */
					mareg_index = &argv[i][2];
					break;
				case 'T':	/* Time step for writing grids */
					step_t = atoi(&argv[i][2]);
					break;
				case 'W':	/* Maregraph file */
					wave_height = &argv[i][2];
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (argc == 0 || error) {
		mexPrintf ("tsun2 %s - Computes the ...\n\n", argv[0]);
		mexPrintf ("usage: tsun2(bat. head, maregs, [out_maregs], [-G<stem>], [-T<grid_step>], [-F<friction>], [-I<dx/dt>], [-J<jmp>], [-M], [-D], [-O<optmareg>], [m])\n");
		mexPrintf ("\t-m outputs a 3D grid used to do a movie\n");
		mexPrintf ("\t-D write grids with the total water depth\n");
		mexPrintf ("\t-G<stem> write grids at the step_t intervals. Append file prefix. Files will be called <stem>#.grd\n");
		mexPrintf ("\t-F<friction> Manning coefficient [Default = 0.025]\n");
		mexPrintf ("\t-I space (read from grid) and time increments [Default -I<dx/1>] \n");
		mexPrintf ("\t-J Times in virtual maregraphs < to <jump_time> will not be loaded.\n");
		mexPrintf ("\t   Use this option to start the simulation at a time close to the time\n");
		mexPrintf ("\t   that the main waves arrive to grid borders. [Default load them all]\n");
		mexPrintf ("\t-M write the last grid of max water level [Default wave surface level]\n");
		mexPrintf ("\t-N number of cycles [Default 1010].\n");
		mexPrintf ("\t-O name of optional output cumulative hight file (default maregs_out_heights.dat).\n");
		mexPrintf ("\t-P<file> name of params file\n");
		mexPrintf ("\t-T time step for writing grids [Default 60 s].\n");
	}

	if (error) return;

	if (!movie && !write_grids && !maregs_out)
		mexErrMsgTxt("Nothing selected for output (grids, movie or out maregs), exiting\n");

	if (!bat_in_input)
		mexErrMsgTxt("Bathymetry was not transmited, exiting\n");
 
	if (maregs_out && maregs_out_hgt == NULL) {
		if ((fpOutMaregs = fopen ("maregs_out_heights.dat", "w")) == NULL)
			mexErrMsgTxt("TSUN2: Unable to open default file name - exiting");
	}
	else if (maregs_out) {
		if ((fpOutMaregs = fopen (maregs_out_hgt, "w")) == NULL) {
			mexPrintf("TSUN2: Unable to open file %s - exiting\n", maregs_out_hgt);
			mexErrMsgTxt("");
		}
	}

	/* Take into account the dt value so that step_t is in seconds */
	step_t = irint(step_t / dt);

	if (water_depth) {
		surf_level = FALSE;
		max_level = FALSE;
	}

	if (time_jump > 0.)	/* If we will jump "time_jump" times, update time_h */
		time_h += time_jump;

	ncl = nx * ny;

	/* Allocate memory	*/
	if ((z = (double *) mxCalloc ((size_t)(ncl * 2), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (z)", ncl * 2);	return;} 
	if ((m = (double *) mxCalloc ((size_t)(ncl * 2), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (m)", ncl * 2);	return;} 
	if ((n = (double *) mxCalloc ((size_t)(ncl * 2), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (n)", ncl * 2);	return;} 
	if ((dz = (double *) mxCalloc ((size_t)(ncl * 2), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (dz)", ncl * 2);	return;} 
	if ((dm = (double *) mxCalloc ((size_t)(ncl * 2), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (dm)", ncl * 2);	return;}
	if ((dn = (double *) mxCalloc ((size_t)(ncl * 2), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (dn)", ncl * 2);	return;}
	if ((hz = (double *) mxCalloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (hz)", ncl);	return;}
	if ((hm = (double *) mxCalloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (hm)", ncl);	return;}
	if ((hn = (double *) mxCalloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (hn)", ncl);	return;}
	if ((zm = (double *) mxCalloc ((size_t)(ncl), sizeof(double)) ) == NULL) 
		{no_sys_mem("tsun2 --> (zm)", ncl);	return;}
	if ((work = (float *) mxCalloc ((size_t)(ncl), sizeof(float)) ) == NULL) 
		{no_sys_mem("tsun2 --> (work)", ncl);	return;}

	/*  ********* INITIAL CONDITION ************* */

	if (!maregs_in_input) {		/* If maregraphs where not given as argument, load them */
		if (n_mareg = data (wave_height, &n_tstep_in, &dt, time_jump) < 0)
			return;
	}
	hmn (nx, ny, hz, hm, hn);
	intl (nx, ny, z, m, n, dz, hz);			/* Set inicial condition */ 
	if (read_index(mareg_index, &nl_w, &nl_e, &nc_s, &nc_n) < 0)	/* Decode the file with maregraph indexes */ 
		return;

	if (amp_fact != 1.) {
		for (i = 0; i < n_tstep_in; i++) {
			for (j = 1; j < n_mareg; j++)
				b[ijb(i,j)] = b[ijb(i,j)] * amp_fact;
		}
	}

	/* Transpose from Matlab orientation to scanline orientation */
	for (i = 0; i < ny; i++)
		for (j = 0; j < nx; j++) hz[i*nx+j] = -hz1[j*ny+i];

	/* ------------------------------------------------------------------------------ */

	if (n_of_cycles == 0) n_of_cycles = n_tstep_in;
	if (n_tstep_in < n_of_cycles) {
		mexPrintf("WARNING: number of time steps (%d) in mareg file is lower than\n", n_tstep_in);
		mexPrintf("\twhat is asked here (%d). Current value is changed accordingly\n", n_of_cycles);
		n_of_cycles = n_tstep_in;
	}

	i2 = MAX(nx,ny);	/* This is for excess, but it's much easier */
	if ((x_w = (double *) mxCalloc ((size_t)(i2), sizeof(double))) == NULL) {
		no_sys_mem("tsun2 --> (x_w)", i2);	return;}
	if ((y_w = (double *) mxCalloc ((size_t)(i2), sizeof(double))) == NULL) {
		no_sys_mem("tsun2 --> (y_w)", i2);	return;}
	if ((u_w = (double *) mxCalloc ((size_t)(i2), sizeof(double))) == NULL) {
		no_sys_mem("tsun2 --> (u_w)", i2);	return;}
	if ((v_w = (double *) mxCalloc ((size_t)(i2), sizeof(double))) == NULL) {
		no_sys_mem("tsun2 --> (v_w)", i2);	return;}

	/* Declarations for the (if) movie option */
	if (movie) {
		mov_32 = (float *)mxCalloc(nx*ny, sizeof(float));
		dims[0] = ny;	dims[1] = nx;	dims[2] = (int)(n_of_cycles / step_t + 1);
		if ((plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL)) == NULL) {
			mexPrintf("TSUN2 ERROR: Could not reallocate memory\n");
			mxFree(mov_32);
			return;
		}
		ptr_mov_32 = (float *)mxGetData(plhs[0]);
	}

	/*  ********* MAIN CALCULATION ************* */

	lhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	rhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	ptr = mxGetPr(rhs[0]);
	tmp_ptr[0] = 0.;				/* Start the waitbar with zero length */
	memcpy(ptr, tmp_ptr, 8);
	rhs[1] = mxCreateString(w_bar_title);		/* Waitbar message */
	mexCallMATLAB(1,lhs,2,rhs,"waitbar");
	h_bar = mxGetPr(lhs[0]);			/* Save the waitbar handle */

	/*  ********* MAIN CALCULATION ************* */

	r  = dt / dx;
	n_of_cycles--;		/* decrease 1 to not having to -1 at each loop step testing */

	for (k = 0; k <= n_of_cycles; k++) {
		if (k % 5 == 0) {
			tmp_ptr[0] = (double)k/(double)(n_of_cycles + 1);
			memcpy(ptr, tmp_ptr, 8);
			mexCallMATLAB(0,NULL,1,rhs,"waitbar");
		}
		kk = k;
		nlmass (nx, ny, z, m, n, dz, hz, r, kk);	/* Change z e dz */
		bnc_n (nx, ny, z, kk, nl_w, nl_e, nc_s, nc_n);
		nlmmt (nx, ny, z, m, n, dz, dm, dn, hz,hm,hn,r,dt,manning_coeff);
		max_z  (nx, ny, z, zm);
		if (maregs_out) {			/* Want time series at maregraph positions */
			time_h += dt;
			fprintf (fpOutMaregs, "%.3f", time_h);
			for (i = 0; i < n_mareg_out; i++) {
				if (surf_level)
					fprintf (fpOutMaregs, "\t%.3f", z[lcum_p[i]-1 + ncl]);
				else if (max_level)
					fprintf (fpOutMaregs, "\t%.3f", zm[lcum_p[i]-1]);
				else
					fprintf (fpOutMaregs, "\t%.3f", dz[lcum_p[i]-1 + ncl]);
			}
			fprintf (fpOutMaregs,"\n");
			fflush(fpOutMaregs);
		}
		if ((k % step_t) == 0 || k == n_of_cycles) {
			if (surf_level) {
				for (j = 0; j < ny; j++) {
					for (i = 0; i < nx; i++)
						work[ij0(i,j)] = (float) z[ij2(i,j,1)];
				}
			}
			else if (max_level) {	/* Max surface level */
				for (j = 0; j < ny; j++) {
					for (i = 0; i < nx; i++)
						work[ij0(i,j)] = (float) zm[ij0(i,j)];
				}
			}
			else if (water_depth) {
				for (j = 0; j < ny; j++) {
					for (i = 0; i < nx; i++)
						work[ij0(i,j)] = (float) dz[ij2(i,j,1)];
				}
			}
			else
				mexPrintf ("Shit! I should't pass here\n");

			if (write_grids && (surf_level || water_depth) )
				write_grd_bin ((int)(k*dt+time_jump), hdr.x_min, hdr.y_min, dx, nx, ny, work);
			else if (write_grids && max_level && (k == n_of_cycles) )	/* Write only the last grid */
				write_grd_bin (n_of_cycles, hdr.x_min, hdr.y_min, dx, nx, ny, work);

			if (movie) {		/* Output the "movie" in a 3D float matrix */
				/* Transpose to matlab ordering */
				if (n_frames < dims[2]) {
					for (i = 0, i2 = ny-1; i < ny; i++, i2--)
					for (j = 0; j < nx; j++) mov_32[j*ny+i2] = work[i*nx+j];

					memcpy(ptr_mov_32, mov_32, nx*ny*4);
					ptr_mov_32 += nx*ny;
				}
			}
			n_frames++;
		}
		change (nx, ny, z, m, n, dz);
	}

	mxSetPr(rhs[0],h_bar);	mxSetPr(rhs[1],NULL);	mxSetPr(lhs[0],NULL);
	mexCallMATLAB(0,lhs,1,rhs,"close");
	mxDestroyArray(rhs[0]); mxDestroyArray(rhs[1]);
	mxDestroyArray(lhs[0]);

	if (movie && movie_float) mxFree(mov_32);

	if (maregs_out) fclose(fpOutMaregs);

	mxFree ((void *) b);
	mxFree ((void *) z);
	mxFree ((void *) m);
	mxFree ((void *) n);
	mxFree ((void *) dz);
	mxFree ((void *) dm);
	mxFree ((void *) dn);
	mxFree ((void *) hz);
	mxFree ((void *) hm);
	mxFree ((void *) hn);
	mxFree ((void *) zm);
	mxFree ((void *) work);
	mxFree ((void *) x_w);
	mxFree ((void *) y_w);
	mxFree ((void *) u_w);
	mxFree ((void *) v_w);
	mxFree ((void *) jp_w);		mxFree ((void *) jp_e);
	mxFree ((void *) ip_s);		mxFree ((void *) ip_n);
	if (maregs_out) mxFree ((void *) lcum_p);	 
}

/* ------------------------------------------------------------------ */
void intl(int nx, int ny, double *z, double *m, double *n, double *dz, double *hz) {
/*	Set inicial condition	*/
	int i, j, k;
	
	for (k = 0; k < 2; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				dz[ij2(i,j,k)] = hz[ij0(i,j)];
				if(hz[ij0(i,j)] <= 0.) {
					dz[ij2(i,j,k)] = 0.;
					z[ij2(i,j,k)] = -hz[ij0(i,j)];
				}
			}
		}
	}
}

/* ------------------------------------------------------------------ */
void hmn (int nx, int ny, double *hz, double *hm, double *hn) {
/*	CAL. OF WATER DEPTH AT POINT OF DISCHARGE	*/
	int i, j;

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			if (i == (nx-1))
				hm[ij0(i,j)] = hz[ij0(i,j)];
			else
				hm[ij0(i,j)] = 0.5*(hz[ij0(i,j)] + hz[ij0(i+1,j)]);
			if (j == (ny-1))
				hn[ij0(i,j)] = hz[ij0(i,j)];
			else
				hn[ij0(i,j)] = 0.5*(hz[ij0(i,j)] + hz[ij0(i,j+1)]);
		}
	}	
}

/* ------------------------------------------------------------------ */
void nlmass (int nx, int ny, double *z, double *m, double *n, double *dz, double *hz, double r, int kk) {
/*	MASS CONSERVATION	*/
	double xm, xn, zzz, dd;
	int i, j;

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx-0; i++) {
			if (hz[ij0(i,j)] < -19.9) {
				dz[ij2(i,j,1)] = 0.;
				z[ij2(i,j,1)] = - hz[ij0(i,j)];
				continue;
			}
			xm = xn = 0.;
			if (i != 0) xm = m[ij2(i-1,j,0)];  /* STUPID test, why? */
			if (j != 0) xn = n[ij2(i,j-1,0)];
			zzz = z[ij2(i,j,0)] - r*(m[ij2(i,j,0)] - xm + n[ij2(i,j,0)] - xn);
			if (fabs(zzz) < GY) zzz = 0.;
			dd = zzz + hz[ij0(i,j)];
			if (dd < GX) {
				dz[ij2(i,j,1)] = 0.;
				z[ij2(i,j,1)] = - hz[ij0(i,j)];
			}
			else {
				dz[ij2(i,j,1)] = dd;
				z[ij2(i,j,1)] = zzz;
			}

		}
	}
}

/* ------------------------------------------------------------------ */
void nlmmt(int nx, int ny, double *z, double *m, double *n, double *dz, double *dm, double *dn, double *hz, double *hm, double *hn, double r, double dt, double fm) {
/*	MOMENTUM CONSERVATION	*/
	/* This funtion is impossible to translate to decent c */
	double dm1, dm2, dn1, dn2, df, zzz, dd, fn, ff, xnn, xne;
	double xn, xm, xmm, xme, z1, z2, zx, seven_three;
	int i, j;

/*    ------ CAL. OF TOTAL DEPTH AT POINT OF DISCHARGE ------- */
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx-1; i++) {
			dm2 = 0.5*(dz[ij2(i,j,1)] + dz[ij2(i+1,j,1)]);
			dm1 = 0.25*(dz[ij2(i,j,1)] + dz[ij2(i,j,0)] + dz[ij2(i+1,j,1)] + dz[ij2(i+1,j,0)]);
	    		if (dm1 < GX) dm1 = 0.;
	    		if (dm2 < GX) dm2 = 0.;
	    		dm[ij2(i,j,0)] = dm1;
	    		dm[ij2(i,j,1)] = dm2;
		}
	}
	for (j = 0; j < ny-1; j++) {
		for (i = 0; i < nx; i++) {
			dn2 = 0.5*(dz[ij2(i,j,1)] + dz[ij2(i,j+1,1)]);
			dn1 = 0.25*(dz[ij2(i,j,1)] + dz[ij2(i,j,0)] + dz[ij2(i,j+1,1)] + dz[ij2(i,j+1,0)]);
	    		if (dn1 < GX) dn1 = 0.;
	    		if (dn2 < GX) dn2 = 0.;
	    		dn[ij2(i,j,0)] = dn1;
	    		dn[ij2(i,j,1)] = dn2;
		}
	}

/*   ------- CAL. OF LINEAR TERMS (X-DIRECTION) ------- */
	fn = 0.5 * dt * GG * fm*fm; 
	seven_three = 7. / 3.;
	for (j = 1; j < ny; j++) {
		for (i = 1; i < nx-1; i++) {
			if (hz[ij0(i,j)] < -20.) continue;
			if (hm[ij0(i,j)] < -20.) goto L30;
			if (dz[ij2(i,j,1)] <= 0.) goto L31;
			else goto L32;
L31:
			if (dz[ij2(i+1,j,1)] <= 0.) goto L30;
			else goto L34;
L32:
			if (dz[ij2(i+1,j,1)] <= 0.) goto L35;
			else goto L36;
L34:
			if (z[ij2(i+1,j,1)] + hz[ij0(i,j)] <= 0.) goto L30;
			else goto L37;
L35:
			if (z[ij2(i,j,1)] + hz[ij0(i+1,j)] <= 0.) goto L30;
			else goto L38;
L36:
			dd = dm[ij2(i,j,1)];
			goto L39;
L37:
			dd = z[ij2(i+1,j,1)] + hz[ij0(i,j)];
			goto L39;
L38:
			dd = z[ij2(i,j,1)] + hz[ij0(i+1,j)];
L39:
			xnn = 0.25 * (n[ij2(i,j,0)] + n[ij2(i+1,j,0)] + n[ij2(i,j-1,0)] + n[ij2(i+1,j-1,0)]);
			df = dd;
			if (df < .01) df = .01;
			ff = fn * d_sqrt(pow(m[ij2(i,j,0)],2) + xnn*xnn) / (pow(df, seven_three));
			if (dd < GX) goto L30;
	    		xm = (1. - ff)*m[ij2(i,j,0)] - GG*r*dd*(z[ij2(i+1,j,1)] - z[ij2(i,j,1)]);

/*  ----- CAL. OF NON-LINEAR TERMS (CONVECTION TERMS) ------ */

	    		if (i <= 1 || j <= 2) goto L40;
	    		if (i+2 > nx - 1 || j+2 > ny - 1) goto L40;
	    		if (dm[ij2(i,j,0)] < GX) goto L40;
	    		if (m[ij2(i,j,0)] <= 0.) goto L41;
	    		else goto L42;
L41:
	    		if (dm[ij2(i+1,j,0)] < GX) goto L40;
	    		if (dz[ij2(i+2,j,1)] < GX) goto L40;
	    		if (dz[ij2(i+1,j,1)] < GX) goto L40;
	    		xm -= r*(pow(m[ij2(i+1,j,0)],2)/dm[ij2(i+1,j,0)] - pow(m[ij2(i,j,0)],2)/dm[ij2(i,j,0)]);
	    		goto L43;
L42:
			if (dm[ij2(i-1,j,0)] < GX) goto L40;
			if (dz[ij2(i-1,j,1)] < GX) goto L40;
			if (dz[ij2(i,j,1)] < GX) goto L40;
	    		xm -= r*(pow(m[ij2(i,j,0)],2)/dm[ij2(i,j,0)] - pow(m[ij2(i-1,j,0)],2)/dm[ij2(i-1,j,0)]);
L43:
			if (xnn <= 0.) goto L44;
			else goto L45;
L44:
			xne = 0.25*(n[ij2(i,j+1,0)] + n[ij2(i+1,j+1,0)] + n[ij2(i,j,0)] + n[ij2(i+1,j,0)]);
	    		if (dm[ij2(i,j+1,0)] < GX) goto L40;
	    		if (dz[ij2(i,j+1,1)] < GX) goto L40;
	    		if (dz[ij2(i,j+2,1)] < GX) goto L40;
	    		if (dz[ij2(i+1,j+1,1)] < GX) goto L40;
	    		if (dz[ij2(i+1,j+2,1)] < GX) goto L40;
	    		xm -= r*(m[ij2(i,j+1,0)]*xne/dm[ij2(i,j+1,0)] - m[ij2(i,j,0)]*xnn/dm[ij2(i,j,0)]);
			goto L40;
L45:
			xne = 0.25*(n[ij2(i,j-1,0)] + n[ij2(i+1,j-1,0)] + n[ij2(i,j-2,0)] + n[ij2(i+1,j-2,0)]);
	    		if (dm[ij2(i,j-1,0)] < GX) goto L40;
	    		if (dz[ij2(i,j-2,1)] < GX) goto L40;
	    		if (dz[ij2(i,j-1,1)] < GX) goto L40;
	    		if (dz[ij2(i+1,j-1,1)] < GX) goto L40;
	    		if (dz[ij2(i+1,j-2,1)] < GX) goto L40;
	    		xm -= r*(m[ij2(i,j,0)]*xnn/dm[ij2(i,j,0)] - m[ij2(i,j-1,0)]*xne/dm[ij2(i,j-1,0)]);
L40:
			xm /= ff + 1.;
			if (fabs(xm) < GX) xm = 0.;

/*   ------ LIMITING OF DISCHARGE ------- */

			if (fabs(xm) < GX) xm = 0.;
			if (xm > dd * 10.) xm = dd * 10.;
			if (xm < dd * -10.) xm = dd * -10.;
	    		m[ij2(i,j,1)] = xm;
	    		continue;
L30:
	    		m[ij2(i,j,1)] = 0.0;
		}
	}

/*   ------- CAL. OF LINEAR TERMS (Y-DIRECTION) ------- */

	for (j = 1; j < ny-1; j++) {
		for (i = 1; i < nx; i++) {
			if (hz[ij0(i,j)] < -20.) continue;
			if (hn[ij0(i,j)] < -20.) goto L130;
			if (dz[ij2(i,j,1)] <= 0.) goto L131;
			else goto L132;
L131:
			if (dz[ij2(i,j+1,0)] <= 0.) goto L130;
			else goto L134;
L132:
			if (dz[ij2(i,j+1,0)] <= 0.) goto L135;
			else goto L136;
L134:
			if (z[ij2(i,j+1,1)] + hz[ij0(i,j)] <= 0.) goto L130;
			else goto L137;
L135:
			if (z[ij2(i,j,1)] + hz[ij0(i,j)] <= 0.) goto L130;
			else goto L138;
L136:
	    		dd = dn[ij2(i,j,1)];
	    		goto L139;
L137:
	    		dd = z[ij2(i,j+1,1)] + hz[ij0(i,j)];
	    		goto L139;
L138:
	    		dd = z[ij2(i,j,1)] + hz[ij0(i,j+1)];
L139:
	    		xmm = 0.25*(m[ij2(i,j,0)] + m[ij2(i,j+1,0)] + m[ij2(i-1,j,0)] + m[ij2(i-1,j+1,0)]);
	    		df = dd;
	    		if (df < .01) df = .01;
	    		ff = fn * d_sqrt(pow(n[ij2(i,j,0)],2) + xmm*xmm) / (pow(df, seven_three));
	    		if (dd < GX) goto L130;
	    		xn = (1. - ff)*n[ij2(i,j,0)] - GG*r*dd*(z[ij2(i,j+1,1)] - z[ij2(i,j,1)]);

/*  ----- CAL. OF NON-LINEAR TERMS (CONVECTION TERMS) ------ */

	    		if (i <= 1 || j <= 1) goto L140;
	    		if (i+2 > nx - 1 || j+2 > ny - 1) goto L140;
	    		if (dn[ij2(i,j,0)] < GX) goto L140;
	    		if (n[ij2(i,j,0)] <= 0.) goto L141;
	    		else goto L142;
L141:
	    		if (dn[ij2(i,j+1,0)] < GX) goto L140;
	    		if (dz[ij2(i,j+2,1)] < GX) goto L140;
	    		if (dz[ij2(i,j+1,1)] < GX) goto L140;
	    		xn -= r*(pow(n[ij2(i,j+1,0)],2)/dn[ij2(i,j+1,0)] - pow(n[ij2(i,j,0)],2)/dn[ij2(i,j,0)]);
	    		goto L143;
L142:
	    		if (dn[ij2(i,j-1,0)] < GX) goto L140;
	    		if (dz[ij2(i,j-1,1)] < GX) goto L140;
	    		if (dz[ij2(i,j,1)] < GX) goto L140;
	    		xn -= r*(pow(n[ij2(i,j,0)],2)/dn[ij2(i,j,0)] - pow(n[ij2(i,j-1,0)],2)/dn[ij2(i,j-1,0)]);
	    		goto L143;
L143:
	    		if (xmm <= 0.) goto L144;
	    		else goto L145;
L144:
			xme = 0.25*(m[ij2(i+1,j,0)] + m[ij2(i+1,j+1,0)] + m[ij2(i,j,0)] + m[ij2(i,j+1,0)]);

	    		if (dn[ij2(i+1,j,0)] < GX) goto L140;
	    		if (dz[ij2(i+1,j,1)] < GX) goto L140;
	    		if (dz[ij2(i+2,j,1)] < GX) goto L140;
	    		if (dz[ij2(i+1,j+1,1)] < GX) goto L140;
	    		if (dz[ij2(i+2,j+1,1)] < GX) goto L140;
	    		xn -= r*(n[ij2(i+1,j,0)]*xme/dn[ij2(i+1,j,0)] - n[ij2(i,j,0)]*xmm/dn[ij2(i,j,0)]);
	    		goto L140;
L145:
			xme = 0.25*(m[ij2(i-1,j,0)] + m[ij2(i-1,j+1,0)] + m[ij2(i-2,j,0)] + m[ij2(i-2,j+1,0)]);
	    		if (dn[ij2(i-1,j,0)] < GX) goto L140;
	    		if (dz[ij2(i-2,j,1)] < GX) goto L140;
	    		if (dz[ij2(i-2,j+1,1)] < GX) goto L140;
	    		if (dz[ij2(i-1,j,1)] < GX) goto L140;
	    		if (dz[ij2(i-1,j+1,1)] < GX) goto L140;
	    		xn -= r*(n[ij2(i,j,0)]*xmm/dn[ij2(i,j,0)] - n[ij2(i-1,j,0)]*xme/dn[ij2(i-1,j,0)]);
	    		goto L140;
L140:
	    		xn /= (ff + 1.);
	    		if (fabs(xn) < GY) xn = 0.;

/*   ------ LIMITING OF DISCHARGE ------- */

	    		if (xn > dd * 20.) xn = dd * 20.;
	    		if (xn < dd * -20.) xn = dd * -20.;
	    		n[ij2(i,j,1)] = xn;
	    		continue;
L130:
	    		n[ij2(i,j,1)] = 0.;
		}
	}
}

/* ***** EXCHANGE FOR LAST STEP DATA TO NEXT STEP DATA ***** */

void change(int nx, int ny, double *z, double *m, double *n, double *d) {
	int i, j;

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			z[ij2(i,j,0)] = z[ij2(i,j,1)];
			m[ij2(i,j,0)] = m[ij2(i,j,1)];
			n[ij2(i,j,0)] = n[ij2(i,j,1)];
			d[ij2(i,j,0)] = d[ij2(i,j,1)];
		}
	}
}


/*  ******* INPUT OF BOUNDARY DATA (WATER LEVEL) ******* */
int data(char *alt_ond, int *n_tstep_in, double *dt, float time_jump) {
	/* Read wave heights */
	int j, n_col = 0, k = 0, n_alloc = CHUNK, ndata = 0;
	char buffer[1024], line[1024], *p;
	FILE *fp;

	if ((fp = fopen (alt_ond, "r")) == NULL) {
		mexPrintf ("%s: Unable to open file %s - exiting\n", "tsun2_j", alt_ond);
		return (-1);
	}

	if ((b = (double *) mxCalloc ((size_t)(n_alloc), sizeof(double)) ) == NULL) {no_sys_mem("tsun2 --> (b)", n_alloc);} 

	while (fgets (line, 1024, fp)) {
		if (n_col == 0) {	/* First time, allocate # of columns */
			strcpy (buffer, line);
			p = (char *)strtok (buffer, " \t\n");
			while (p) {	/* Count # of fields */
				n_col++;
				p = (char *)strtok ((char *)NULL, " \t\n");
			}
			/* Now we know # of columns */
		}
		p = (char *)strtok (line, " \t\n");
		j = 0;
		while (p && j < n_col) {
			sscanf (p, "%lf", &b[k*n_col+j]);
			j++;
			ndata = ijb(k,j);
			p = (char *)strtok ((char *)NULL, " \t\n");
		}
		if (b[k*n_col] < time_jump) continue;	/* Jump this time by not incrementing k */
		if (j != n_col) mexPrintf ("%s %s:  Expected %d but found %d fields in record # %d\n", "DATA", alt_ond, n_col, j, k);
		k++;
		if (ndata + j >= n_alloc) {
			n_alloc += CHUNK;
			if ((b = (double *) mxRealloc ((void *)b, (size_t)(n_alloc * sizeof(double)))) == NULL) {no_sys_mem("tsun2 --> (b)", n_alloc);}
		}
	}
	*n_tstep_in = k;
	*dt = b[j] - b[0];		/* Simulation time step */
	fclose (fp);
	return (n_col-1);		/* Return number of maregraphs (-1 because fist col has time) */
}

/*   *************  OUTPUT OF DATA ****************** */
int write_grd_bin (int kk, double x_min, double y_min, double dx, int nx, int ny, float *work) {

	/* Writes a grid in the Surfer binary format */
	int i, j;
	double x_max, y_max;
	float work_min = FLT_MAX, work_max = -FLT_MAX;
	char name[80];
	struct srf_header h;
	FILE *fp;

	if (stem[0] == 0)
		sprintf (name,"%d%s", kk,".grd");
	else
		sprintf (name, "%s%d%s", stem, kk,".grd");

	if ((fp = fopen (name, "wb")) == NULL) {
		mexPrintf ("Fatal Error: Could not create file %s!\n", name);
		return (-1);
	}

	x_max = x_min + (nx - 1) * dx;
	y_max = y_min + (ny - 1) * dx;

	/* Find zmin/zmax */
	for (i = 0; i < nx*ny; i++) {
		work_max = MAX(work[i], work_max);
		work_min = MIN(work[i], work_min);
	}

	/* store header information and array */
	strcpy (h.id,"DSBB");
	h.nx = nx;	 	h.ny = ny;
	h.x_min = x_min;	h.x_max = x_max;
	h.y_min = y_min;	h.y_max = y_max;
	h.z_min = (double)work_min;	h.z_max = (double)work_max;

	if (fwrite ((void *)&h, sizeof (struct srf_header), (size_t)1, fp) != 1) {
		mexPrintf ("Fatal Error: Error writing file %s!\n", name);
		return (-1);
	}
	for(j = 0; j < ny; j++) {
		for(i = 0; i < nx; i++) {
			fwrite ((void *)&work[ij0(i,j)], sizeof(float), (size_t)1, fp);
		}
	}

	fclose(fp);
	return (0);
}

/*     ********* CHECK OF MAXIMUM VALUE *********** */
void max_z (int nx, int ny, double *z, double *zm) {
	int i, j;

	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			if(zm[ij0(i,j)] < z[ij2(i,j,1)]) zm[ij0(i,j)] = z[ij2(i,j,1)];
		}
	}
}

/* ------------------------------------------------------------------ */
int read_index(char *file, int *nl_w, int *nl_e, int *nc_s, int *nc_n) {
	/* Decode the file with maregraph indexes	*/
	int i, j, k, n_col, n_border, n_node[3], n_alloc1 = 1500, side[4];
	int n_alloc2 = 500, nb, *i_tmp;
	char *p, line[512], buffer[512];
	FILE *fp;

	*nl_w = 0; *nl_e = 0; *nc_s = 0; *nc_n = 0;
	n_node[0] = n_node[1] = n_node[2] = 0;
	side[0] = side[1] = side[2] = side[3] = 0;

	if ((i_tmp = (int *) mxCalloc ((size_t)(n_alloc1), sizeof(int)) ) == NULL) 
		{no_sys_mem("tsun2 --> (i_tmp)", n_alloc1);} 

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("%s: Unable to open file %s - exiting\n", "tsun2_j", file);
		return (-1);
	}

	fgets (line, 512, fp);
	while (line[0] == '#')	/* Jump comment lines */
		fgets (line, 512, fp);
	strcpy (buffer, line);
	n_border = count_col (buffer);	/* Count number of borders */

	if (n_border > 3)
		mexPrintf ("WARNING: number of borders must be <= 3. Something bad will hapen\n");

	for (i = 0; line[i]; i++) {
		switch (line[i]) {
			case ' ':
				break;
			case 'w':
			case 'W':
				side[0] = 1;
				jp_w = (int *) mxCalloc ((size_t)(n_alloc2), sizeof(int));
				break;
			case 's':
			case 'S':
				side[1] = 1;
				ip_s = (int *) mxCalloc ((size_t)(n_alloc2), sizeof(int));
				break;
			case 'e':
			case 'E':
				side[2] = 1;
				jp_e = (int *) mxCalloc ((size_t)(n_alloc2), sizeof(int));
				break;
			case 'n':
			case 'N':
				side[3] = 1;
				ip_n = (int *) mxCalloc ((size_t)(n_alloc2), sizeof(int));
				break;
			case '\n':
			case '\r':
				break;
			default:
				mexPrintf ("ERROR: border code (%c) unknown. Must be one the following: W S E N\n", line[i]);
				return(-1);
		}
	}

	for (i = k = 0; i < n_border; i++) {
		fgets (line, 512, fp);
		strcpy (buffer, line);
		n_col = count_col (buffer); /* conta numero de nodos da fronteira */
		p = (char *)strtok (line, " \t\n");
		j = 0;
		while (p && j < n_col) {
			sscanf (p, "%d", &i_tmp[j+k]);
			j++;
			p = (char *)strtok ((char *)NULL, " \t\n");
		}
		n_node[i] = n_col;
		k += n_node[i];
	}

	nb = k = 0;
	if (side[0] == 1) {		/* W border */
		for (j = 0; j < n_node[nb]; j++) jp_w[j] = i_tmp[j+k];
		k += n_node[nb];	nb++;	*nl_w = j;
	}
	if (side[1] == 1) {	/* S border */
		for (j = 0; j < n_node[nb]; j++) ip_s[j] = i_tmp[j+k];
		k += n_node[nb];	nb++;	*nc_s = j;
	}
	if (side[2] == 1) {	/* E border */
		for (j = 0; j < n_node[nb]; j++) jp_e[j] = i_tmp[j+k];
		k += n_node[nb];	nb++;	*nl_e = j;
	}
	if (side[3] == 1) {	/* N border */
		for (j = 0; j < n_node[nb]; j++) ip_n[j] = i_tmp[j+k];
		k += n_node[nb];	nb++;	*nc_n = j;
	}

	fclose (fp);
	free ((void *) i_tmp);
	return(0);
}

/* ------------------------------------------------------------------ */
int count_col (char *line) {
	/* Count # of fields contained in line */
	int n_col = 0;
	char *p;

	p = (char *)strtok (line, " \t\n\r");
	while (p) {	/* Count # of fields */
		n_col++;
		p = (char *)strtok ((char *)NULL, " \t\n\r");
	}
	return (n_col);
}

/* ------------------------------------------------------------------ */
void bnc_n(int nx, int ny, double *z, int kk, int nl_w, int nl_e, int nc_s, int nc_n) {

	int i, j;

	if (nl_w != 0) {  /* Interpolate along the West border of the smaller grid */
		for (i = 0; i < nl_w; i++) {
			x_w[i] = jp_w[i];
			y_w[i] = b[ijb(kk,i)];
		}
		for (j = 0; j < ny; j++) u_w[j] = j;
		intp_lin (x_w, y_w, nl_w, ny, u_w, v_w, 0);
		for (j = 0; j < ny; j++)
			z[ij2(1,j,1)] = v_w[j];
	}

	if (nc_s != 0) {  /* Interpolate along the South border of the smaller grid */
		for (i = 0; i < nc_s; i++) {
			x_w[i] = ip_s[i];
			y_w[i] = b[ijb(kk, i + nl_w)];
		}
		for (i = 0; i < nx; i++) u_w[i] = i;
		intp_lin (x_w, y_w, nc_s, nx, u_w, v_w, 0);
		for (i = 0; i < nx; i++)
			z[ij2(i,1,1)] = v_w[i];
	}

	if (nl_e != 0) {  /* Interpolate along the East border of the smaller grid */
		for (j = 0; j < nl_e; j++) {
			x_w[j] = jp_e[j];
			y_w[j] = b[ijb(kk, j + nl_w + nc_s)];
		}
		for (j = 0; j < ny; j++) u_w[j] = j;
		intp_lin (x_w, y_w, nl_e, ny, u_w, v_w, 0);
		for (j = 0; j < ny; j++)
			z[ij2(nx-2,j,1)] = v_w[j];
	}
	if (nc_n != 0) {  /* Interpolate along the North border of the smaller grid */
		for (i = 0; i < nc_n; i++) {
			x_w[i] = ip_n[i];
			y_w[i] = b[ijb(kk, i + nl_w + nc_s + nl_e)];
		}
		for (i = 0; i < nx; i++) u_w[i] = i;
		intp_lin (x_w, y_w, nc_n, nx, u_w, v_w, 0);
		for (i = 0; i < nx; i++)
			z[ij2(i,ny-2,1)] = v_w[i];
	}
}

/* ------------------------------------------------------------------ */
void	no_sys_mem (char *where, int n) {	
		mexPrintf ("Fatal Error: %s could not allocate memory, n = %d\n", where, n);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * INTP_LIN will interpolate from the dataset <x,y> onto a new set <u,v>
 * where <x,y> and <u> is supplied by the user. <v> is returned. The 
 * parameter mode governs what interpolation scheme that will be used.
 *
 * input:  x = x-values of input data
 *	   y = y-values "    "     "
 *	   n = number of data pairs
 *	   m = number of desired interpolated values
 *	   u = x-values of these points
 *   	  mode = 0 : Linear interpolation
 * output: v = y-values at interpolated points
 *
 * Programmer:	Paul Wessel
 * Date:	16-JAN-1987
 * Ver:		v.2 --> cripled version for linear interpolation (J.L.)
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */

int intp_lin (double *x, double *y, int n, int m, double *u, double *v, int mode) {
	int i, j, err_flag = 0;
	int down = FALSE;
	double dx;
	
	/* Check to see if x-values are monotonically increasing/decreasing */

	dx = x[1] - x[0];
	if (dx > 0.0) {
		for (i = 2; i < n && err_flag == 0; i++)
			if ((x[i] - x[i-1]) < 0.0) err_flag = i;
	}
	else {
		down = TRUE;
		for (i = 2; i < n && err_flag == 0; i++)
			if ((x[i] - x[i-1]) > 0.0) err_flag = i;
	}

	if (err_flag) {
		mexPrintf ("%s: Fatal Error: x-values are not monotonically increasing/decreasing!\n", "intp_lin");
		return (err_flag);
	}
	
	if (down) {	/* Must flip directions temporarily */
		for (i = 0; i < n; i++) x[i] = -x[i];
		for (i = 0; i < m; i++) u[i] = -u[i];
	}

	/* Compute the interpolated values from the coefficients */
	
	j = 0;
	for (i = 0; i < m; i++) {
		while (x[j] > u[i] && j > 0) j--;	/* In case u is not sorted */
		while (j < n && x[j] <= u[i]) j++;
		if (j == n) j--;
		if (j > 0) j--;

		dx = u[i] - x[j];
		v[i] = (y[j+1]-y[j])*dx/(x[j+1]-x[j]) + y[j];
	}

	if (down) {	/* Must reverse directions */
		for (i = 0; i < n; i++) x[i] = -x[i];
		for (i = 0; i < m; i++) u[i] = -u[i];
	}

	return (0);
}
