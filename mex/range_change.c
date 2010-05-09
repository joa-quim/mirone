/* rngchn_change.c - executable for Matlab (Mex file) RNGCHN program. */
/* Kurt Feigl */

/* input parameter list, all real*8 */
/* variable                               dimen- suggested */
/* name   meaning                         sion   units */

/* x      patch easting coordinates       nseg    km */
/* y      patch northing coordinates      nseg    km */
/* d      fault depth (lower left)        nseg    km */
/* strike fault azimuth                   nseg    deg CW */
/* delta  fault dip                       nseg    deg */
/* u1     fault up-dip slip               nseg    mm */
/* u2     fault left-lateral strike-slip  nseg    mm */
/* u3     fault tensile slip              nseg    mm */
/* l      fault length                    nseg    km */
/* w      fault width                     nseg    km */
/* sx     ground point easting coords.    npts    km */
/* sy     ground point northing coords.   npts    km */
/* plook  unit vector from sat to ground  3       - */

/* nseg = number of fault segments (patches) */
/* npts = number of ground points */

/* For an arbitrary list of points, then sx and sy are both */
/* column vectors of dimension npts. */
/* For an array (image), then set sx to be a row vector, and r will be an */
/* array of dimesions dim(sx) by dim(sy), starting with r(1,1) at x=sx(1), */
/* y=sy(1). */

/*
 *	Translated to C & mexified (+ improvments) By
 *	Joaquim Luis - 2005
 *
 * Revision 17/10/2007 JL	Condensed terms on chin() routine
 *
 */

#include "mex.h"
#include <math.h>
#include <string.h>

#define	FALSE	0
#define	TRUE	1
#ifndef TWO_PI
#define TWO_PI  6.28318530717958647692
#endif
#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif
#define D2R	0.01745329251994329577
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
/*#define D2R	M_PI / 180.*/
#define SUMOL	1.0e-8
#define d_sqrt(x) ((x) < 0.0 ? 0.0 : sqrt (x))
/* In non-Windows this is may not be necessary (or guive conflicts) */
#define copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

/* Table of constant values */

static double c_b235 = 0.;

void rngchn_comp(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, int);
int form_all_part3(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int);
int okada_patch(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int chin(double *, double *, double * , double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru_delta(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int ech_(double *, double *), okada_patch(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru_d(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru_l(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru_u(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru_w( double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru_x( double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru1_x(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru_p(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru1_delta(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
int deru1_d(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
void tm (double lon, double lat, double *x, double *y);
void vtm (double lon0, double lat0);

int	is_geog = FALSE;
double	EQ_RAD = 6378137.0;	/* WGS-84 */
double	flattening = 1.0/298.2572235630;
double	ECC2, ECC4, ECC6;
double	one_m_ECC2, i_one_m_ECC2;
double	t_c1, t_c2, t_c3, t_c4, t_e2, t_M0;
double	central_meridian;

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	i, i1, i2, j, k, m, n, nseg, npts, nnl, mmx, mmy, nnx, nny;
	int	argc = 0, n_arg_no_char = 0;
	int	isarr = FALSE;
	char	**argv;
	double *ptrsx, *ptrsy;
	double *sx, *sy, *tx, *ty, *faultdelta, *partials, *faultlambda, *faultstrike, *plook;
	double *faultd, *ranges, *patchx, *patchy, *faultl, *faultw;
	double *ptrrng, *ptrprt;
	double *faultu1, *faultu2, *faultu3, *faultmu, rx, ry, d1, d2;

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
	argv[0] = "range_change";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	/* Maybe this test will grow, but for now it only goes for detecting the -M option */
        for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
				case 'M': /* Geographical coords */
					is_geog = TRUE;
					break;
                        }
                }
        }

	/* ----------------------------------------------------------------------- */
	/*     fault patch parameters */
	/*     ground point coordinates */
	/*     look vector */
	/*     ouptut range change */
	/*     ouptut partial derivatives */
	/* ----------------------------------------------------------------------- */
	/*     all these have dimensions which we do not know in advance. */
	/*     Choose a maximum possible value here. */

	/* CHECK FOR PROPER NUMBER OF ARGUMENTS */
	if (nrhs < 13) {
		mexPrintf("variable				dimen-	\n");
		mexPrintf("name		meaning			sion	units-\n");
		mexPrintf("\n");
		mexPrintf("x      patch easting coordinates       nseg  km\n");
		mexPrintf("y      patch northing coordinates      nseg  km\n");
		mexPrintf("d      fault depth (lower left)        nseg  km\n");
		mexPrintf("strike fault azimuth                   nseg  deg\n");
		mexPrintf("delta  fault dip                       nseg  deg\n");
		mexPrintf("u1     fault up-dip slip               nseg  mm\n");
		mexPrintf("u2     fault left-lateral strike-slip  nseg  mm\n");
		mexPrintf("u3     fault tensile slip              nseg  mm\n");
		mexPrintf("l      fault length                    nseg  km\n");
		mexPrintf("w      fault width                     nseg  km\n");
		mexPrintf("e      ground point easting coords.    npts  km\n");
		mexPrintf("n      ground point northing coords.   npts  km\n");
		mexPrintf("s      unit vector from ground to sat  3     - \n");
		mexPrintf("\n");
		mexPrintf(" nseg = number of fault segments\n");
		mexPrintf(" npts = number of ground points \n");
		mexPrintf("\n");
		mexPrintf(" For an arbitrary list of points, then e and n\nare both column vectors of dimension npts.\n");
		mexPrintf("\n");
		mexPrintf(" For an array (image), then set sx to be a row\n");
		mexPrintf(" vector, and r will be an array of dimesions dim(e) by dim(n)\n");
		mexPrintf(" starting with r(1,1) at east=e(1), north=n(1).\n");
		mexPrintf("\n");
		mexPrintf(" Example use:\n");
		mexPrintf(" r=range_change(x,y,strike,depth,dip,u1,u2,u3,l,w,e,n,[0.,0.,1.])\n");
		mexPrintf("\n");
		mexPrintf(" To calculate partials in an array a:\n");
		mexPrintf(" [r,a]=range_change(x,y,strike,depth,dip,u1,u2,u3,l,w,e,n,[0.,0.,1.])\n");
		mexPrintf("\n");
		mexPrintf(" One row of A will contain, in order:\n");
		mexPrintf("\n");
		mexPrintf("1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n");
		mexPrintf("dr/dX\tdr/dY\tdr/dstrike\tdr/ddepth\tdr/ddip\tdr/dU1\tdr/dU2\tdr/dU3\tdr/dL\tdr/dW\n");
		mexPrintf("mm/km\tmm/km\tmm/deg\tmm/km\tmm/deg\tmm/mm\tmm/mm\tmm/mm\tmm/km\tmm/km\n");
		mexPrintf("\n");
		mexPrintf(" NOTE THE FOLLOWING ASSUMPTIONS:\n");
		mexPrintf("    Poisson solid:  lambda = mu = 1\n");
		mexPrintf("    Finite fault (f = 1)\n");
		mexPrintf("    Only one fault if partials are requested\n");
		mexPrintf("\n");
		mexErrMsgTxt("RNGCHN requires 13 input arguments\n");
	}
	else if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("RNGCHN requires 1 or 2 output arguments\n");

	/* CHECK THE DIMENSIONS OF FAULT PARAMETERS */

	nseg = mxGetM(prhs[0]);
	for (i = 1; i <= 9; ++i) {
		m = mxGetM(prhs[i]);
		if (m != nseg)
			mexErrMsgTxt("RNGCHN: dimension error in input fault parameters\n");
	}

	/* CHECK THE DIMENSIONS OF GROUND POINT ARRAYS */

	nnx = mxGetN(prhs[10]);		mmx = mxGetM(prhs[10]);
	nny = mxGetN(prhs[11]);		mmy = mxGetM(prhs[11]);
	if (mmx == mmy && nnx == 1 && nny == 1)
		npts = mmx;
	else if (mmx == 1 && nny == 1) {
		isarr = TRUE;
		npts = nnx * mmy;
	}
	/*     make sure that input fault variables are column vectors */
	for (i = 0; i <= 9; ++i) {
		n = mxGetN(prhs[i]);
		if (n != 1) {
 			mexPrintf("Bad parameter number: %d\n", i);
 			mexErrMsgTxt("RNGCHN: inputs should be column vectors (:,1)\n");
		}
	}

	/* CHECK THE DIMENSIONS OF LOOK UNIT VECTOR (E,N,U) */

	/* Computing MAX */
	i1 = mxGetM(prhs[12]), i2 = mxGetN(prhs[12]);
	nnl = MAX(i1,i2);
	if (nnl != 3) {
		mexErrMsgTxt("RNGCHN: look unit vector needs dimension 3\n");
	}

	if (isarr)
		plhs[0] = mxCreateDoubleMatrix (mmy,nnx, mxREAL);
	else
		plhs[0] = mxCreateDoubleMatrix (npts,1, mxREAL);
	ptrrng = mxGetPr(plhs[0]);
	/*     partial derivatives for one fault only! */
	if (nlhs == 2) {
		if (nseg == 1) {
			plhs[1] = mxCreateDoubleMatrix (npts,10, mxREAL);
			ptrprt = mxGetPr(plhs[1]);
			partials = mxCalloc (npts*10, sizeof (double));
		}
		else
			mexErrMsgTxt("RNGCHN: Partials available for only 1 fault at a time!\n");
	}

	faultlambda = mxCalloc (nseg+1, sizeof (double));
	faultmu = mxCalloc (nseg+1, sizeof (double));
	sx = mxCalloc (npts, sizeof (double));
	sy = mxCalloc (npts, sizeof (double));
	tx = mxCalloc (npts, sizeof (double));
	ty = mxCalloc (npts, sizeof (double));
	ranges = mxCalloc (npts, sizeof (double));

	/* ASSIGN POINTERS TO THE VARIOUS PARAMETERS */
	/* order is x,y,strike,depth,dip,u1,u2,u3,l,w,sx,sy */
	patchx = mxGetPr(prhs[0]);		patchy = mxGetPr(prhs[1]);
	faultstrike = mxGetPr(prhs[2]);		faultd = mxGetPr(prhs[3]);
	faultdelta = mxGetPr(prhs[4]);		faultu1 = mxGetPr(prhs[5]);
	faultu2 = mxGetPr(prhs[6]);		faultu3 = mxGetPr(prhs[7]);
	faultl = mxGetPr(prhs[8]);		faultw = mxGetPr(prhs[9]);
	ptrsx = mxGetPr(prhs[10]);		ptrsy = mxGetPr(prhs[11]);
	plook = mxGetPr(prhs[12]);

	/* COPY RIGHT HAND ARGUMENTS TO LOCAL ARRAYS OR VARIABLES */
	memcpy(sx, ptrsx, nnx*8);
	memcpy(sy, ptrsy, mmy*8);

	/* ASSUME POISSON SOLID */
	for (i = 0; i < nseg; i++)
		faultlambda[i] = faultmu[i] = 1.;

	if (isarr) {	/* make array */ 
		k = 0;
		for (i = 0; i < nnx; i++) {
			for (j = 0; j < mmy; j++) {
				tx[k] = sx[i];
				ty[k] = sy[j];
				k++;
			}
		}
		if (is_geog) {
			vtm(patchx[0]+0.01, patchy[0]+0.01);	/*  Set patch(x|y)[0] as the projection origin*/
			for (k = 0; k < npts; k++) {
				tm(tx[k], ty[k], &rx, &ry);
				sx[k] = rx/1000.;	sy[k] = ry/1000.;
			}
			for (k = 0; k < nseg; k++) {
				tm(patchx[k], patchy[k], &rx, &ry);
				patchx[k] = rx/1000.;	patchy[k] = ry/1000.;
			}
		}
		else {
			for (k = 0; k < npts; k++) {
				sx[k] = tx[k];
				sy[k] = ty[k];
			}
		}
	}
	rngchn_comp(ranges, sx, sy, patchx, patchy, faultl, faultw, faultd, faultdelta, faultstrike, faultlambda, faultmu, faultu1, faultu2, faultu3, plook, npts, nseg);
	if (nlhs == 2)
		form_all_part3(partials, sx, sy, patchx, patchy, faultl, faultw, faultd, faultdelta, faultstrike, faultlambda, faultmu, faultu1, faultu2, faultu3, plook, npts);

	/* COPY OUTPUT WHICH IS STORED IN LOCAL ARRAY TO MATRIX OUTPUT */
	memcpy(ptrrng, ranges, npts*8);
	if (nlhs == 2) {
		n = npts * 10;
		memcpy(ptrprt, partials, n*8);
	}
	if (nlhs == 2) mxFree(partials);
	mxFree(faultlambda);	mxFree(faultmu);	mxFree(ranges);
	mxFree(sx);		mxFree(sy);		mxFree(tx);		mxFree(ty);
}
/* cccc -------------- END OF MEX GATEWAY ---------------------- */

void rngchn_comp(double *ranges, double *sx, double *sy, double *patchx, double *patchy, double *faultl,
	double *faultw, double *faultd, double *faultdelta, double *faultstrike, double *faultlambda, 
	double *faultmu, double *faultu1, double *faultu2, double *faultu3, double *plook, int npts, int nseg) {
	/*     number of x,y ground points */
	/*     number of faults */
	/*     range change */
	/*     strike-slip, dip-slip, and tensile slip on surface, in fault geometry */
	/*     strike-slip, dip-slip, and tensile slip on surface, in map geometry */
	/*     displacements due to 1 fault in map geometry */

	int i, j;
	double uf[9], um[9], xp, yp, *ca1, *ca2, *sa1, *sa2, range, fdism[3], faultx, faulty;

	sa1 = mxCalloc (npts, sizeof (double));
	sa2 = mxCalloc (npts, sizeof (double));
	ca1 = mxCalloc (npts, sizeof (double));
	ca2 = mxCalloc (npts, sizeof (double));
	for (i = 0; i < nseg; i++) {
		/* For rotation from map to fault geometry. Direction sine and cosine for rotation */
		sa1[i] = sin((90. - faultstrike[i]) * D2R);
		ca1[i] = cos((90. - faultstrike[i]) * D2R);
		/* For rotation from fault to map geometry. Direction sine and cosine for rotation */
		sa2[i] = sin((faultstrike[i] - 90.) * D2R);
		ca2[i] = cos((faultstrike[i] - 90.) * D2R);
	}
	/* loop over all ground points */
	for (j = 0; j < npts; j++) {
		fdism[0] = fdism[1] = fdism[2] = 0.;
		/*        loop over all fault segments */
		for (i = 0; i < nseg; i++) {
			/*  rotate from map to fault geometry */
			xp = sx[j] - patchx[i];
			yp = sy[j] - patchy[i];
			faultx =  xp * ca1[i] + yp * sa1[i];
			faulty = -xp * sa1[i] + yp * ca1[i];
			okada_patch(&faultu1[i], &faultu2[i], &faultu3[i], &faultx, &faulty, &faultl[i], 
				    &faultw[i], &faultd[i], &faultdelta[i], &faultlambda[i], &faultmu[i], uf);
			/*  rotate from fault to map */
			um[0] =  uf[0] * ca2[i] + uf[1] * sa2[i];
			um[1] = -uf[0] * sa2[i] + uf[1] * ca2[i];
			um[2] =  uf[3] * ca2[i] + uf[4] * sa2[i];
			um[3] = -uf[3] * sa2[i] + uf[4] * ca2[i];
			um[4] =  uf[6] * ca2[i] + uf[7] * sa2[i];
			um[5] = -uf[6] * sa2[i] + uf[7] * ca2[i];
			um[6] = uf[2];
			um[7] = uf[5];
			um[8] = uf[8];
			/*  strike-slip, dip-slip, tensile parts for e,n,v components of displacement */
			fdism[0] = fdism[0] + um[0] + um[2] + um[4];
			fdism[1] = fdism[1] + um[1] + um[3] + um[5];
			fdism[2] = fdism[2] + um[6] + um[7] + um[8];
		}
		/* OK. The vector plook should contain a unit vector from the ground point to the satellite.
        	   Now take its scalar product with the motion vector, fdism */
		range = 0.;
		for (i = 0; i < 3; i++) range += plook[i] * fdism[i];
		ranges[j] = range;
	}
	mxFree(sa1);	mxFree(sa2);	mxFree(ca1);	mxFree(ca2);
}

int okada_patch(double *fltu1, double *fltu2, double *fltu3, double *fltx, double *flty, double * fltl, double *fltw, double *fltd, double *fltdelta, double *fltlambda, double *fltmu, double *uf) {
	/*     calculate slip on a finite rectangular patch using Okada's expressions. 
	/*     input parameters in flt geometry
	/*     see Figure 1 of paper for definition of these variables
	/*     strike-slip, dip-slip, and tensile slip components on flt
	/*     these are input.
	/*     strike-slip, dip-slip, and tensile slip on surface, in flt geometry
	/*     these are output
	/*     local variables
	/*     dip sine and cosine, useful term
	/*     Okada variables */

	int i;
	double p, cd, sd, xi, up[9], eta, fmat, small;

	small = 1e-8;
	if (fabs(fabs(*fltdelta) - 90.) > small)
		cd = cos(*fltdelta * D2R);
	else
		cd = 0.;

	sd = sin(*fltdelta * D2R);
	fmat = *fltmu / (*fltlambda + *fltmu);
	p = *flty * cd + *fltd * sd;
	for (i = 0; i < 9; i++) uf[i] = 0.;
	xi = *fltx;
	eta = p;
	chin(fltx, flty, &sd, &cd, &fmat, &xi, &eta, fltd, fltu1, fltu2, fltu3, up);
	for (i = 0; i < 9; i++) uf[i] += up[i];
	xi = *fltx;
	eta = p - *fltw;
	chin(fltx, flty, &sd, &cd, &fmat, &xi, &eta, fltd, fltu1, fltu2, fltu3, up);
	for (i = 0; i < 9; i++) uf[i] -= up[i];
	xi = *fltx - *fltl;
	eta = p;
	chin(fltx, flty, &sd, &cd, &fmat, &xi, &eta, fltd, fltu1, fltu2, fltu3, up);
	for (i = 0; i < 9; i++) uf[i] -= up[i];
	xi = *fltx - *fltl;
	eta = p - *fltw;
	chin(fltx, flty, &sd, &cd, &fmat, &xi, &eta, fltd, fltu1, fltu2, fltu3, up);
	for (i = 0; i < 9; i++) uf[i] += up[i];
	for (i = 0; i < 9; i++)
		if (fabs(uf[i]) < small) uf[i] = 0.;

	return 0;
}

int chin(double *fltx, double *flty, double *sd, double *cd, double *fmat, double *xi, double *eta, 
	double *fltd, double *fltu1, double *fltu2, double *fltu3, double *up) {

	double	p, q, i1, i2, i3, i4, i5, t1, t2, t3, lg, at, dt, rrxi, tmp2, sd2, sd_cd;
	double	rr, yt, xx, xi2, xi3, rr3, jnk, rpe, rdt, tmp, jnk2, rdt2, frpe, fmat2, small;
	double	q_rpe, q_rrxi;

	/*     perform Chinnery double bar substition */
	/* Parameter adjustments */
	--up;

	small = 1e-8;
	p = *flty * *cd + *fltd * *sd;
	q = *flty * *sd - *fltd * *cd;
	yt = *eta * *cd + q * *sd;
	dt = *eta * *sd - q * *cd;
	tmp = *xi * *xi + *eta * *eta + q * q;
	if (fabs(tmp) > small)
		rr = sqrt(tmp);
	else
		rr = 0.;
	tmp = *xi * *xi + q * q;
	if (fabs(tmp) > small)
		xx = sqrt(tmp);
	else
		xx = 0.;
	fmat2 = *fmat / 2.;
	t1 = -(*fltu1) / TWO_PI;
	t2 = -(*fltu2) / TWO_PI;
	t3 = *fltu3 / TWO_PI;
	rdt = rr + dt;
	rdt2 = rdt * rdt;
	xi2 = *xi * *xi;
	xi3 = xi2 * *xi;
	rr3 = rr * rr * rr;
	rpe = rr * (rr + *eta);
	rrxi = rr * (rr + *xi);
	sd2 = *sd * *sd;
	sd_cd = *sd * *cd;
	if (fabs(rr + *eta) < small)
		frpe = 0.;
	else
		frpe = rr + *eta;
	/*     find displacement parts */
	jnk2 = *xi * *eta / (q * rr);
	if (fabs(q) < small)
		at = 0.;
	else
		at = atan(jnk2);
	if (fabs(frpe) < small)
		lg = -log(rr - *eta);
	else
		lg = log(rr + *eta);
	if (fabs(*cd) < small) {
		i1 = -fmat2 * *xi * q / rdt2;
		i3 = fmat2 * (*eta / rdt + yt * q / rdt2 - lg);
		i4 = -(*fmat) * q / rdt;
		/*        Andrea does not handle this case */
		if (fabs(*xi) < small)
	    		i5 = 0.;
		else
	    		i5 = -(*fmat) * *xi * *sd / (rr + dt);
		i2 = *fmat * (-lg) - i3;
	} else {
		jnk = (*eta * (xx + q * *cd) + xx * (rr + xx) * *sd) / (*xi * (rr + xx) * *cd);
		if (fabs(*xi) < small)
			i5 = 0.;
		else
			i5 = *fmat * 2. / *cd * atan(jnk);
		i4 = *fmat * (1. / *cd) * (log(rdt) - *sd * lg);
		i3 = *fmat * (1. / *cd * yt / rdt - lg) + *sd * i4 / *cd;
		i2 = *fmat * (-lg) - i3;
		i1 = *fmat * (-1. / *cd * (*xi / rdt)) - *sd * i5 / *cd;
	}
	q_rrxi = q / rrxi;
	if (fabs(frpe) < small) {
		up[1] = t1 * (at + i1 * *sd);
		up[2] = t1 * (i2 * *sd);
		up[3] = t1 * (i4 * *sd);
		up[4] = t2 * (q / rr - i3 * sd_cd);
		up[5] = t2 * (yt * q_rrxi + *cd * at - i1 * sd_cd);
		up[6] = t2 * (dt * q_rrxi + *sd * at - i5 * sd_cd);
		up[7] = t3 * (-i3 * sd2);
		up[8] = t3 * (-dt * q_rrxi - *sd * (-at) - i1 * sd2);
		up[9] = t3 * ( yt * q_rrxi + *cd * (-at) - i5 * sd2);
	} else {
		tmp2 = (*xi * q / rpe - at);
		q_rpe = q / rpe;
		up[1] = t1 * (*xi * q_rpe + at + i1 * *sd);
		up[2] = t1 * (yt * q_rpe + q * *cd / (rr + *eta) + i2 * *sd);
		up[3] = t1 * (dt * q_rpe + q * *sd / (rr + *eta) + i4 * *sd);
		up[4] = t2 * (q / rr - i3 * sd_cd);
		up[5] = t2 * (yt * q_rrxi + *cd * at - i1 * sd_cd);
		up[6] = t2 * (dt * q_rrxi + *sd * at - i5 * sd_cd);
		up[7] = t3 * ( q * q_rpe - i3 * sd2);
		up[8] = t3 * (-dt * q_rrxi - *sd * tmp2 - i1 * sd2);
		up[9] = t3 * ( yt * q_rrxi + *cd * tmp2 - i5 * sd2);
	}
	return 0;
}

int form_all_part3(double *partials, double *sx, double *sy, double *patchx, double *patchy, double *faultl,
	double *faultw, double *faultd, double *faultdelta, double *faultstrike, double *faultlambda, 
	double *faultmu, double *faultu1, double *faultu2, double *faultu3, double *plook, int npts) {

	int i, j, k, l, n;
	double a, d1, vf[9], vm[9], xp, yp, ca1, ca2, faultx, faulty;
	double sa1, sa2, vf1[9], vf2[9], range, gdism[3], partl[10];

	/*     This is the version called by the matlab gateway. */
	/*     No common blocks! */
	/*     EXACT PARTIAL DERIVATIVES */
	/*     strike-slip, dip-slip, and tensile slip on surface, in fault geometry */
	/*     strike-slip, dip-slip, and tensile slip on surface, in map geometry */
	/*     partial derivatives in same units and order as parameters */
	/*     displacements due to 1 fault in map geometry */
	/* Parameter adjustments */
	--plook;
	--sy;
	--sx;
	--partials;
	--faultu3;
	--faultu2;
	--faultu1;
	--faultmu;
	--faultlambda;
	--faultstrike;
	--faultdelta;
	--faultd;
	--faultw;
	--faultl;
	--patchy;
	--patchx;

	/*  Only 1 fault segment at a time! */
	i = 1;
	/*        for rotation from map to fault geometry */
	/*        direction sine and cosine for rotation */
	sa1 = sin((90. - faultstrike[i]) * D2R);
	ca1 = cos((90. - faultstrike[i]) * D2R);
	sa2 = sin((faultstrike[i] - 90.) * D2R);
	ca2 = cos((faultstrike[i] - 90.) * D2R);
	/*     loop over all points */
	for (j = 1; j <= npts; ++j) {
		range = 0.;
		gdism[0] = gdism[1] = gdism[2] = 0.;
		/*        FROM HERE DOWN, THE CODE IS THE SAME AS FORM_ALL_PART2.F */
		/*        rotate from map to fault geometry */
		xp = sx[j] - patchx[i];
		yp = sy[j] - patchy[i];
		faultx = xp * ca1 + yp * sa1;
		faulty = -xp * sa1 + yp * ca1;
		for (k = 1; k <= 10; ++k) {
			if (k == 1 || k == 2 || k == 3) {
				/*              DERIVATIVES WITH RESPECT TO U1,U2,U3 */
				deru_u(&faultx, &faulty, &faultdelta[i], &faultd[i], &faultu1[i], &faultu2[i], &faultu3[i], &faultl[i], &faultw[i], &faultlambda[i], &faultmu[i], vf, vf1, vf2);
				if (k == 2) ech_(vf1, vf);
				if (k == 3) ech_(vf2, vf);
	    		}
			/*         The following subroutines need delta in radian */
			if (k == 4 || k == 5 || k == 10) {
				/*            DERIVATIVES WITH RESPECT TO X , Y  AND ALPHA */
				/*            X is patchx, Y is patchy */
				d1 = faultdelta[i] * D2R;
				deru_x(&faultx, &faulty, &d1, &faultd[i], &faultu1[i], &faultu2[i], &faultu3[i], &faultl[i], &faultw[i], &faultlambda[i], &faultmu[i], vf, vf1);
				/*            Transformation to obtain the derivatives with respect to patchx and patchy */
				for (n = 1; n <= 9; ++n) {
					a = vf[n - 1];
					vf2[n-1] = (-faulty * a + faultx * vf1[n - 1]) * D2R;
					vf[n-1] = -ca1 * a + sa1 * vf1[n - 1];
					vf1[n-1] = -ca1 * vf1[n - 1] - sa1 * a;
				}
				if (k == 5) ech_(vf1, vf);
				if (k == 10) ech_(vf2, vf);
	    		}
			if (k == 6) { /*            DERIVATIVES WITH RESPECT TO L -- length */
				d1 = faultdelta[i] * D2R;
				deru_l(&faultx, &faulty, &d1, &faultd[i], &faultu1[i], &faultu2[i], &faultu3[i], &faultl[i], &faultw[i], &faultlambda[i], &faultmu[i], vf);
			}
			if (k == 7) { /*      DERIVATIVES WITH RESPECT TO W --  width */
				d1 = faultdelta[i] * D2R;
				deru_w(&faultx, &faulty, &d1, &faultd[i], &faultu1[i], &faultu2[i], &faultu3[i], &faultl[i], &faultw[i], &faultlambda[i], &faultmu[i], vf);
	    		}
			if (k == 8) { /*      DERIVATIVES WITH RESPECT TO depth */
				d1 = faultdelta[i] * D2R;
				deru_d(&faultx, &faulty, &d1, &faultd[i], &faultu1[i], &faultu2[i], &faultu3[i], &faultl[i], &faultw[i], &faultlambda[i], &faultmu[i], vf);
			}
			if (k == 9) { /*      DERIVATIVES WITH RESPECT TO DELTA --  dip */
				d1 = faultdelta[i] * D2R;
				deru_delta(&faultx, &faulty, &d1, &faultd[i], &faultu1[ i], &faultu2[i], &faultu3[i], &faultl[i], & faultw[i], &faultlambda[i], &faultmu[i], vf);
				for (n = 1; n <= 9; ++n) vf[n-1] *= D2R;
			}
			if (k == 10) {
				okada_patch(&faultu1[i], &faultu2[i], &faultu3[i], &faultx, &faulty, &faultl[i], &faultw[i], &faultd[i], &faultdelta[i], &faultlambda[i], &faultmu[i], vf1);
				/*             Rotate from fault to map with respect to alpha */
				vm[0] = (vf[0] + vf1[1] * D2R) * ca2 + (vf[1] - vf1[0] * D2R) * sa2;
				vm[1] = -sa2 * (vf1[1] * D2R + vf[0]) + (vf[1] - vf1[0] * D2R) * ca2;
				vm[2] = ca2 * (vf[3] + vf1[4] * D2R) + sa2 * (vf[4] - vf1[3] * D2R);
				vm[3] = ca2 * (vf[4] - vf1[3] * D2R) - sa2 * (vf[3] + vf1[4] * D2R);
				vm[4] = ca2 * (vf[6] + vf1[7] * D2R) + sa2 * (vf[7] - vf1[6] * D2R);
				vm[5] = ca2 * (vf[7] - vf1[6] * D2R) - sa2 * (vf[6] + vf1[7] * D2R);
				vm[6] = vf[2];
				vm[7] = vf[5];
				vm[8] = vf[8];
	    		} else {
				/*             rotate from fault to map */
				vm[0] = vf[0] * ca2 + vf[1] * sa2;
				vm[1] = -vf[0] * sa2 + vf[1] * ca2;
				vm[2] = vf[3] * ca2 + vf[4] * sa2;
				vm[3] = -vf[3] * sa2 + vf[4] * ca2;
				vm[4] = vf[6] * ca2 + vf[7] * sa2;
				vm[5] = -vf[6] * sa2 + vf[7] * ca2;
				vm[6] = vf[2];
				vm[7] = vf[5];
				vm[8] = vf[8];
			}
			/*           strike-slip, dip-slip, tensile parts for */
			/*           e,n,v components of displacement */
			gdism[0] = vm[0] + vm[2] + vm[4];
			gdism[1] = vm[1] + vm[3] + vm[5];
			gdism[2] = vm[6] + vm[7] + vm[8];
			/*           scalar product with s vector : plook. */
			range = 0.;
			for (l = 1; l <= 3; ++l) range += plook[l] * gdism[l-1];
			partl[k-1] = range;
		}
		partials[j] = partl[3]; 	/*     dr/dX */
		partials[j + npts] = partl[4]; 	/*     dr/dY */
		partials[j + 2*npts] = partl[9]; /*     dr/dstrike */
		partials[j + 3*npts] = partl[7]; /*     dr/ddepth */
		partials[j + 4*npts] = partl[8]; /*     dr/ddip */
		partials[j + 5*npts] = partl[0]; /*     dr/dU1 */
		partials[j + 6*npts] = partl[1]; /*     dr/dU2 */
		partials[j + 7*npts] = partl[2]; /*     dr/dU3 */
		partials[j + 8*npts] = partl[5]; /*     dr/dL */
		partials[j + 9*npts] = partl[6]; /*     dr/dW */
	}
	return 0;
}

/* ***************************************************************************** */
int deru_x(double *x, double *y, double *delta, double *d, double *u1, double *u2, double *u3, double *l, double *w, double *lambda, double *mu, double *dux, double *duy) {

	int i;
	double k, p, q, et, dz, duux[9], duuy[9];

	/* Parameter adjustments */
	--duy;
	--dux;

	k = *mu / (*lambda + *mu);
	p = *y * cos(*delta) + *d * sin(*delta);
	q = *y * sin(*delta) - *d * cos(*delta);
	/*      initialisation of dux and duy */
	for (i = 1; i <= 9; ++i) {
		duy[i] = 0.;
		dux[i] = 0.;
	}
	dz = *x;
	et = p;
	deru1_x(&dz, &et, delta, d, &p, &q, &k, u1, u2, u3, duux, duuy);
	for (i = 1; i <= 9; ++i) {
		duy[i] = duuy[i-1] + duy[i];
		dux[i] = duux[i-1] + dux[i];
	}
	dz = *x;
	et = p - *w;
	deru1_x(&dz, &et, delta, d, &p, &q, &k, u1, u2, u3, duux, duuy);
	for (i = 1; i <= 9; ++i) {
		duy[i] = -duuy[i-1] + duy[i];
		dux[i] = -duux[i-1] + dux[i];
	}
	dz = *x - *l;
	et = p;
	deru1_x(&dz, &et, delta, d, &p, &q, &k, u1, u2, u3, duux, duuy);
	for (i = 1; i <= 9; ++i) {
		duy[i] = -duuy[i-1] + duy[i];
		dux[i] = -duux[i-1] + dux[i];
	}
	dz = *x - *l;
	et = p - *w;
	deru1_x(&dz, &et, delta, d, &p, &q, &k, u1, u2, u3, duux, duuy);
	for (i = 1; i <= 9; ++i) {
		duy[i] = duuy[i-1] + duy[i];
		dux[i] = duux[i-1] + dux[i];
	}
	return 0;
}

int deru1_x(double *xi, double *et, double * delta, double *d, double *p, double *q, double *k, double *u1, double *u2, double *u3, double *duux, double *duuy) {

	double d1, c, r, s, t, j1, k1, j2, k3, j3, j4, k2, h5, q3, r3, bn, dt, rx, yt, xx, g2p, i1p, et3, h5p, i5p, xi3, axi, aeta, small;

	/* Parameter adjustments */
	--duuy;
	--duux;

	small = 1e-6;
	c = cos(*delta);
	s = sin(*delta);
	t = tan(*delta);
	r = sqrt(*xi * *xi + *et * *et + *q * *q);
	xx = sqrt(*xi * *xi + *q * *q);
	rx = r + xx;
	yt = *et * c + *q * s;
	dt = *et * s - *q * c;
	r3 = r * r * r;
	/*     special case of cos(delta) = 0 */
	if (fabs(c) > small) {
		if (fabs(r + *et) > small) {
			k1 = *k * *xi * (1. / (r * (r + dt)) - s / (r * (r + *et) )) / c;
			k3 = *k * (*q / (r * (r + *et)) - yt / (r * (r + dt))) / c;
		} else {
			k1 = *k * *xi / (r * (r + dt)) / c;
			k3 = *k * (-yt / (r * (r + dt))) / c;
		}
	} else {
		k1 = *k * *xi * *q / (r * (r + dt) * (r + dt));
		k3 = *k * (s / (r + dt)) * (*xi * *xi / (r * (r + dt)) - 1.);
	}
	if (fabs(r + *et) > small)
		k2 = *k * (-s / r + *q * c / (r * (r + *et))) - k3;
	else
		k2 = *k * (-s / r) - k3;
	/*     special case of cos(delta) = 0 */
	if (fabs(c) > small) {
		/* Computing 2nd power */
		d1 = r + dt;
		j1 = *k * (*xi * *xi / (r * (d1 * d1)) - 1. / (r + dt)) / c - t * k3;
		d1 = r + dt;
		j2 = *k * *xi * yt / (r * c * (d1 * d1)) - t * k1;
	} else {
		/* Computing 2nd power */
		d1 = r + dt;
		j1 = *k * .5 * *q / (d1 * d1) * (*xi * 2 / (r * (r + dt)) - 1.);
		d1 = r + dt;
		j2 = *k * .5 * *xi * s / (d1 * d1) * (*q * 2 / (r * (r + dt)) - 1.);
	}
	if (fabs(r + *et) > small) {
		j3 = *k * (-(*xi) / (r * (r + *et))) - j2;
		j4 = *k * (-c / r - *q * s / (r * (r + *et))) - j1;
	} else {
		j3 = j2 * -1.;
		j4 = *k * (-c / r) - j1;
	}
	h5 = (*et * (xx + *q * c) + xx * rx * s) / (*xi * rx * c);
	h5p = (*xi * *xi - xx * (xx + *q * c) * (*xi * *xi / (r * xx) + 1)) / rx;
	h5p = (*et * h5p / c - *q * *q * t) / (xx * *xi * *xi);
	i5p = *k * 2. * h5p / (h5 * h5 + 1) / c;
	g2p = *q / ((r + *et) * (r - *et));
	/* Computing 2nd power */
	d1 = r + dt;
	i1p = -(*k) * (1. / (r + dt) - *xi * *xi / (r * (d1 * d1))) / c - t * i5p;
	/*     (36a) */
	if (fabs(r + *et) > small) {
		/* Computing 2nd power */
		d1 = r + *et;
		aeta = (r * 2. + *et) / (r3 * (d1 * d1));
	} else
		aeta = 0.;
	/*     (36b) */
	/* Computing 2nd power */
	d1 = r + *xi;
	axi = (r * 2. + *xi) / (r3 * (d1 * d1));
	d1 = r + *xi;
	bn = (*xi * 2. * r + r * r + *xi * *xi) / (r3 * (d1 * d1));
	xi3 = *xi * *xi * *xi;
	et3 = *et * *et * *et;
	q3 = *q * *q * *q;
	/*     DERIVATIVES WITH RESPECT TO X */
	/*     Strike slip */
	/*     duxdx (31a) */
	duux[1] = *u1 * (*xi * *xi * *q * aeta - j1 * s) / TWO_PI;
	/*     duydx (31c) */
	duux[2] = *u1 * (*xi * *q * c / r3 + (*xi * *q * *q * aeta - j2) * s) / TWO_PI;
	/*     duzdx (37a) */
	duux[3] = *u1 * (-(*xi) * *q * *q * aeta * c + (*xi * *q / r3 - k1) * s) / TWO_PI;
	/*     Dip slip */
	/*     duxdx (32a) */
	duux[4] = *u2 * (*xi * *q / r3 + j3 * s * c) / TWO_PI;
	if (fabs(r + *et) > small) {
		/*         duydx (32c) */
		duux[5] = *u2 * (yt * *q / r3 + *q * c / (r * (r + *et)) + j1 * s * c) / TWO_PI;
		/*         duzdx (38a) */
		duux[6] = *u2 * (dt * *q / r3 + *q * s / (r * (r + *et)) + k3 * s * c) / TWO_PI;
	} else {
		duux[5] = *u2 * (yt * *q / r3 + j1 * s * c) / TWO_PI;
		duux[6] = *u2 * (dt * *q / r3 + *q * s / (r * (r + *et)) + k3 * s * c) / TWO_PI;
	}
	/*     Tensile */
	/*     duxdx (33a) */
	duux[7] = -(*u3) * (*xi * *q * *q * aeta + j3 * s * s) / TWO_PI;
	/*     duydx (33c) */
	duux[8] = -(*u3) * (*q * *q * c / r3 + q3 * aeta * s + j1 * s * s) / TWO_PI;
	/*     duzdx (39a) */
	duux[9] = -(*u3) * (*q * *q * s / r3 - q3 * aeta * c + k3 * s * s) / TWO_PI;
	/*     DERIVATVES WITH RESPECT TO Y */
	/*     Strike-slip */
	/*     duxdy (31b) */
	duuy[1] = *u1 * (xi3 * dt / (r3 * (*et * *et + *q * *q)) - (xi3 * aeta + j2) * s) / TWO_PI;
	/*     duydy (31d) */
	if (fabs(r + *et) > small) {
		duuy[2] = *u1 * (s * (q3 * aeta * s - *q * 2. * s / (r * (r + * et)) - (*xi * *xi + *et * *et) * c / r3 - j4) + yt * *q * c / r3) / TWO_PI;
	} else {
		duuy[2] = *u1 * (s * (q3 * aeta * s - (*xi * *xi + *et * *et) * c / r3 - j4) + yt * *q * c / r3) / TWO_PI;
	}
	/*     duzdy  (37b) */
	duuy[3] = *u1 * (s * (*xi * *xi * *q * aeta * c - s / r + yt * *q / r3 - k2) + dt * *q * c / r3) / TWO_PI;
	/*     Dip slip */
	/*     duxdy (32b) */
	duuy[4] = *u2 * (yt * *q / r3 - s / r + j1 * s * c) / TWO_PI;
	/*     duydy (32d) */
	/*     duydz (38b) */
	if (fabs(r + *et) > small) {
		duuy[5] = *u2 * (-s * (yt * 2. / (r * (r + *xi)) + *xi * c / ( r * (r + *et))) + yt * yt * *q * axi + j2 * s * c) / TWO_PI;
		duuy[6] = *u2 * (-s * (dt * 2. / (r * (r + *xi)) + *xi * s / ( r * (r + *et))) + yt * dt * *q * axi + k1 * s * c) / TWO_PI;
	} else {
		duuy[5] = *u2 * (-s * 2. * yt / (r * (r + *xi)) + yt * yt * *q * axi + j2 * s * c) / TWO_PI;
		duuy[6] = *u2 * (-s * 2. * dt / (r * (r + *xi)) + yt * dt * *q * axi + k1 * s * c) / TWO_PI;
	}
	/*     Tensile */
	/*     duydx (33b) */
	duuy[7] = -(*u3) * (-dt * *q / r3 - *xi * *xi * *q * aeta * s + j1 * s * s) / TWO_PI;
	/*     duydy (33d) */
	duuy[8] = -(*u3) * ((yt * c - dt * s) * *q * *q * axi - s * 2 * c * *q / (r * (r + *xi)) - (*xi * *q * *q * aeta - j2) * s * s) / TWO_PI;
	/*     duydz (39b) */
	duuy[9] = -(*u3) * ((yt * s + dt * c) * *q * *q * axi + *xi * *q * *q *aeta * s * c - (*q * 2. / (r * (r + *xi)) - k1) * s * s) / TWO_PI;
	return 0;
}

/* ***************************************************************************** */
int deru_w(double *x, double *y, double *delta, double *d, double *u1, double *u2, double *u3, double *l, double *w, double *lambda, double *mu, double *duw) {

	int i;
	double d1, k, p, q, dz, duup[9];

	/* Parameter adjustments */
	--duw;

	k = *mu / (*lambda + *mu);
	p = *y * cos(*delta) + *d * sin(*delta);
	q = *y * sin(*delta) - *d * cos(*delta);
	/*      initialisation of  duw */
	for (i = 1; i <= 9; ++i) duw[i] = 0.;
	dz = *x;
	d1 = p - *w;
	deru_p(&dz, delta, d, &d1, &q, &k, u1, u2, u3, duup);
	for (i = 1; i <= 9; ++i) duw[i] = duup[i-1] + duw[i];
	dz = *x - *l;
	d1 = p - *w;
	deru_p(&dz, delta, d, &d1, &q, &k, u1, u2, u3, duup);
	for (i = 1; i <= 9; ++i) duw[i] = -duup[i-1] + duw[i];
	return 0;
}

/* ***************************************************************************** */
int deru_l(double *x, double *y, double *delta, double *d, double *u1, double *u2, double *u3, double *l, double *w, double *lambda, double *mu, double *dul) {

	int i;
	double k, p, q, et, dz, duux[9], duuy[9];

	/* Parameter adjustments */
	--dul;

	k = *mu / (*lambda + *mu);
	p = *y * cos(*delta) + *d * sin(*delta);
	q = *y * sin(*delta) - *d * cos(*delta);
	/*      initialisation of dul */
	for (i = 1; i <= 9; ++i) dul[i] = 0.;
	dz = *x - *l;
	et = p;
	deru1_x(&dz, &et, delta, d, &p, &q, &k, u1, u2, u3, duux, duuy);
	for (i = 1; i <= 9; ++i) dul[i] = duux[i - 1] + dul[i];
	dz = *x - *l;
	et = p - *w;
	deru1_x(&dz, &et, delta, d, &p, &q, &k, u1, u2, u3, duux, duuy);
	for (i = 1; i <= 9; ++i) dul[i] = -duux[i-1] + dul[i];
	return 0;
}

int deru_delta(double *x, double *y, double * delta, double *d, double *u1, double *u2, double * u3, double *l, double *w, double *lambda, double *mu, double *dud) {

	int i;
	double d1, k, p, q, duud[9];

	/* dud : vector of partial derivatives */
	/* intermediate variables (cf OKADA) */
	/* Parameter adjustments */
	--dud;

	k = *mu / (*lambda + *mu);
	p = *y * cos(*delta) + *d * sin(*delta);
	q = *y * sin(*delta) - *d * cos(*delta);
	/*      initialisation of dud */
	for (i = 1; i <= 9; ++i) dud[i] = 0.;
	deru1_delta(x, y, delta, d, &p, &q, &c_b235, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = duud[i-1] / TWO_PI + dud[i];
	deru1_delta(x, y, delta, d, &p, &q, w, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = -duud[i-1] / TWO_PI + dud[i];
	d1 = *x - *l;
	deru1_delta(&d1, y, delta, d, &p, &q, &c_b235, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = -duud[i-1] / TWO_PI + dud[i];
	d1 = *x - *l;
	deru1_delta(&d1, y, delta, d, &p, &q, w, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = duud[i-1] / TWO_PI + dud[i];
	return 0;
}

int deru1_delta(double *x, double *y, double *delta, double *d, double *p, double *q, double *w, double *k, double *u1, double *u2, double *u3, double *duud) {

	int i;
	double c, d1, d2, r, s, t, g1, g2, g3, g4, i1, g6, i2, i3, i4, i5;
	double h5, r3, pp, rx, ww, xx, g1p, g2p, g3p, g4p, g5p, g6p, i1p, i2p, i3p, i4p, i5p, h5p;

	/* Parameter adjustments */
	--duud;

	r = sqrt(*x * *x + (*p - *w) * (*p - *w) + *q * *q);
	xx = sqrt(*x * *x + *q * *q);
	r3 = r * r * r;
	c = cos(*delta);
	s = sin(*delta);
	t = tan(*delta);
	pp = *p - *w;
	rx = r + xx;
	ww = r + *d - *w * s;
	/*      definition des fonctions intermediaires */
	g1 = *q / (r * (r + pp));
	g2 = *x * pp / (*q * r);
	g3 = 1 / (c * ww);
	g4 = *q * c / (r + pp);
	g6 = *q / (r * (r + *x));
	/* Computing 2nd power */
	d1 = r + pp;
	g1p = ((*p * r * r - *q * *q * *w) * (pp + r) - *q * *q * r * (*w - r)) / (r3 * (d1 * d1));
	g2p = *x * ((-(*p) * pp - *q * *q) / (*q * *q * r) + *w * (*w - *p) / r3);
	g3p = s / (c * c * ww) - *w * (*q / r - c) / (c * ww * ww);
	/* Computing 2nd power */
	d1 = r + pp;
	g4p = ((*p * c - *q * s) * (r + pp) - *q * *q * c * (*w / r - 1.) ) / (d1 * d1);
	/* Computing 2nd power */
	d1 = r + pp;
	g5p = (*p * s + *q * c) / (r + pp) - *q * *q * s * (*w - r) / (r * (d1 * d1));
	/* Computing 2nd power */
	d1 = r + *x;
	g6p = *p / (r * (r + *x)) - *q * *q * *w * (r * 2. + *x) / (r3 * ( d1 * d1));
	h5 = (pp * (xx + *q * c) + xx * rx * s) / (*x * rx * c);
	i5 = *k * 2. * atan(h5) / c;
	i4 = *k * (log(ww) - s * log(r + pp)) / c;
	i3 = *k * ((*y - *w * c) / (c * ww) - log(r + pp)) + t * i4;
	i2 = -(*k) * log(r + pp) - i3;
	i1 = *k * (-(*x) * g3) - t * i5;
	/*      Expression des derivees des fonctions intermediaires */
	/* Computing 2nd power */
	d1 = rx;
	/* Computing 2nd power */
	d2 = rx;
	h5p = -(*q) * (xx * xx - *p * pp) / (xx * rx * c) - pp * *q * (*p + *w * xx / r) / (d1 * d1 * c) - *q * *q * pp * (*p + xx * *w / r) / (xx * (d2 * d2)) + *q * *p * t / xx + xx / (c * c) + xx * pp * s / (c * c * rx) + (*p * pp - *q * *q) / rx;
	h5p /= *x;
	i5p = *k * 2. * (h5p / (h5 * h5 + 1) + t * atan(h5)) / c;
	i4p = *k * (t * log(ww) / c - log(r + pp) * (t * t + 1) + 1 / c * (* w * (*q / r - c) / ww - *q * s * (*w / r - 1) / (r + pp))) ;
	i3p = *k * ((*y - *w * c) * g3p + g3 * *w * s - *q * (*w - r) / (r * (r + pp))) + t * i4p + i4 / (c * c);
	i2p = *k * (-(*q) * (*w - r) / (r * (r + pp))) - i3p;
	i1p = -(*k) * *x * g3p - t * i5p - i5 / (c * c);
	/*      Strike-slip */
	duud[1] = *x * g1p + g2p / (g2 * g2 + 1) + i1p * s + i1 * c;
	duud[2] = *w * s * g1 + (*y - *w * c) * g1p + g4p + s * i2p + c * i2;
	duud[3] = (*d - *w * s) * g1p - *w * c * g1 + g5p + i4 * c + i4p * s;
	for (i = 1; i <= 3; ++i) duud[i] *= -(*u1);
	/*      Dip-slip */
	duud[4] = *p / r - *q * *q * *w / r3 - i3 * cos(*delta * 2.) - i3p * s * c;
	duud[5] = g6p * (*y - *w * c) + *w * s * g6 + c * g2p / (g2 * g2 + 1) - s * atan(g2) - i1p * s * c - cos(*delta * 2.) * i1;
	duud[6] = -(*w) * c * g6 + (*d - *w * s) * g6p + s * g2p / (g2 * g2 + 1) + c * atan(g2) - i5p * s * c - i5 * cos(*delta * 2.);
	for (i = 4; i <= 6; ++i) duud[i] *= -(*u2);
	/*     Tensil-fault */
	duud[7] = *q * g1p + *p * g1 - i3p * s - i3 * 2. * s * c - i3p * s * s;
	duud[8] = -(*d - *w * s) * g6p + *w * c * g6 - s * (*x * g1p - g2p / ( g2 * g2 + 1)) - c * (*x * g1 - atan(g2)) - i1p * s * s - i1 * 2 * s * c;
	duud[9] = (*y - *w * c) * g6p + *w * s * g6 + c * (*x * g1p - g2p / ( g2 * g2 + 1)) - s * (*x * g1 - atan(g2)) - i5p * s * s - s * 2 * c * i5;
	for (i = 7; i <= 9; ++i) duud[i] *= *u3;
	return 0;
}


int deru1_d(double *x, double *y, double *delta, double *d, double *p, double *q, double *w, double *k, double *u1, double *u2, double *u3, double *duud) {

	int i;
	double d1, d2, c, r, s, t, g1, g2, g3, g4, i1, g6, i2, i3, i4, i5, h5, g7, r3, dt, pp, rx, yt, ww, xx, g1p, g2p, g3p, g4p, i1p, g6p, i2p, i3p, i4p, i5p, h5p;

	/* Parameter adjustments */
	--duud;

	r = sqrt(*x * *x + (*p - *w) * (*p - *w) + *q * *q);
	xx = sqrt(*x * *x + *q * *q);
	r3 = r * r * r;
	c = cos(*delta);
	s = sin(*delta);
	t = tan(*delta);
	pp = *p - *w;
	rx = r + xx;
	ww = r + *d - *w * s;
	dt = *d - *w * s;
	yt = *y - *w * c;
	/*     Intermediate functions */
	g1 = *q / (r * (r + pp));
	g2 = *x * pp / (*q * r);
	g3 = 1 / (c * ww);
	g4 = *q / (r + pp);
	g6 = *q / (r * (r + *x));
	g7 = -(c * r * r + dt * *q) / r3;
	/* Computing 2nd power */
	d1 = r + pp;
	g1p = -c / (r * (r + pp)) - *q * (r * 2. * dt + r * r * s + dt * pp) / (r3 * (d1 * d1));
	/* Computing 2nd power */
	d1 = *q * r;
	/* Computing 2nd power */
	d2 = *x * pp;
	g2p = *x * (r * r * yt - pp * *q * dt) / (r * (d1 * d1 + d2 * d2));
	g3p = 1. / (c * r * (r + dt));
	/* Computing 2nd power */
	d1 = r + pp;
	g4p = -c / (r + pp) - *q * (dt + r * s) / (r * (d1 * d1));
	/* Computing 2nd power */
	d1 = r + *x;
	g6p = -c / (r * (r + *x)) - *q * dt * (r * 2. + *x) / (r3 * ( d1 * d1));
	h5 = (pp * (xx + *q * c) + xx * rx * s) / (*x * rx * c);
	i5 = *k * 2. * atan(h5) / c;
	i4 = *k * (log(ww) - s * log(r + pp)) / c;
	i3 = *k * ((*y - *w * c) / (c * ww) - log(r + pp)) + t * i4;
	i2 = -(*k) * log(r + pp) - i3;
	i1 = *k * (-(*x) * g3) - t * i5;
	/* Computing 2nd power */
	d1 = rx;
	h5p = (-c * pp * (c + *q / xx) + s * (xx + *q * c)) / rx - pp * (xx + *q * c) * (dt - *q * r * c / xx) / (r * (d1 * d1));
	h5p = h5p / (*x * c) - s * *q / (xx * *x);
	i5p = *k * 2. * (h5p / (h5 * h5 + 1)) / c;
	i4p = *k * (1 - s * (dt + r * s) / (r + pp)) / (r * c);
	i3p = -(*k) * (yt / (c * (r + dt)) + (dt + r * s) / (r + pp)) / r + t * i4p;
	i2p = -(*k) * (dt + r * s) / (r * (r + pp)) - i3p;
	i1p = *k * *x * g3p - t * i5p;
	/*      Strike-slip */
	duud[1] = *x * g1p + g2p + i1p * s;
	duud[2] = yt * g1p + g4p * c + i2p * s;
	duud[3] = dt * g1p + g1 + s * g4p + i4p * s;
	for (i = 1; i <= 3; ++i) duud[i] *= -(*u1);
	/*      Dip-slip */
	duud[4] = g7 - i3p * s * c;
	duud[5] = g6p * yt + c * g2p - i1p * s * c;
	duud[6] = dt * g6p + g6 + s * g2p - i5p * s * c;
	for (i = 4; i <= 6; ++i) duud[i] *= -(*u2);
	/*      Tensil-fault */
	duud[7] = *q * g1p - g1 * c - i3p * s * s;
	duud[8] = -dt * g6p - g6 - s * (*x * g1p - g2p) - i1p * s * s;
	duud[9] = yt * g6p + c * (*x * g1p - g2p) - i5p * s * s;
	for (i = 7; i <= 9; ++i) duud[i] *= *u3;
	return 0;
}

/* ***************************************************************************** */
int deru_d(double *x, double *y, double *delta, double *d, double *u1, double *u2, double *u3, double *l, double *w, double *lambda, double *mu, double *dud) {

	int i;
	double k, p, q, duud[9], d1;

	/* Parameter adjustments */
	--dud;

	k = *mu / (*lambda + *mu);
	p = *y * cos(*delta) + *d * sin(*delta);
	q = *y * sin(*delta) - *d * cos(*delta);
	/*   Initialisation of vector dud */
	for (i = 1; i <= 9; ++i) dud[i] = 0.;
	deru1_d(x, y, delta, d, &p, &q, &c_b235, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = duud[i-1] / TWO_PI + dud[i];
	deru1_d(x, y, delta, d, &p, &q, w, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = -duud[i-1] / TWO_PI + dud[i];
	d1 = *x - *l;
	deru1_d(&d1, y, delta, d, &p, &q, &c_b235, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = -duud[i-1] / TWO_PI + dud[i];
	d1 = *x - *l;
	deru1_d(&d1, y, delta, d, &p, &q, w, &k, u1, u2, u3, duud);
	for (i = 1; i <= 9; ++i) dud[i] = duud[i-1] / TWO_PI + dud[i];
	return 0;
}

/* **************************************************************************** */
int deru_u(double *x, double *y, double *delta, double *d, double *u1, double *u2, double *u3, double *l, double *w, double *lambda, double *mu, double *du1, double *du2, double *du3) {
	int i;
	double uf[9], a;

	/* Parameter adjustments */
	--du3;
	--du2;
	--du1;

	a = 1.;
	okada_patch(&a, &a, &a, x, y, l, w, d, delta, lambda, mu, uf);
	for (i = 1; i <= 3; ++i) {
		du1[i] = uf[i-1] / a;
		du2[i] = 0.;
		du3[i] = 0.;
	}
	for (i = 4; i <= 6; ++i) {
		du1[i] = 0.;
		du2[i] = uf[i-1] / a;
		du3[i] = 0.;
	}
	for (i = 7; i <= 9; ++i) {
		du1[i] = 0.;
		du2[i] = 0.;
		du3[i] = uf[i-1] / a;
	}
	return 0;
}

int deru_p(double *dz, double *delta, double *d, double *p, double *q, double *k, double *u1, double *u2, double *u3, double *duup) {
	double d1, d2, c, r, s, t, g1, g2, g6, h5, r3, dt, rx, yt, xx, g1p, g2p, i1p, i2p, i3p, h5p, i5p, i4p, g3p, g6p;

	/* Parameter adjustments */
	--duup;

	c = cos(*delta);
	s = sin(*delta);
	t = tan(*delta);
	r = sqrt(*dz * *dz + *p * *p + *q * *q);
	xx = sqrt(*dz * *dz + *q * *q);
	rx = r + xx;
	yt = *p * c + *q * s;
	dt = *p * s - *q * c;
	r3 = r * r * r;
	g1 = *q / (r * (r + *p));
	g2 = atan(*dz * *p / (*q * r));
	g6 = *q / (r * (r + *dz));
	g1p = -(*q) / r3;
	/* Computing 2nd power */
	d1 = *q * r;
	/* Computing 2nd power */
	d2 = *dz * *p;
	g2p = *q * *dz * (r * r - *p * *p) / (r * (d1 * d1 + d2 * d2));
	g3p = -(*q) / (r * (r + *p));
	/* Computing 2nd power */
	d1 = r + *dz;
	g6p = -(*q) * *p * (r * 2 + *dz) / (r3 * (d1 * d1));
	h5 = (*p * (xx + *q * c) + xx * rx * s) / (*dz * rx * c);
	/* Computing 2nd power */
	d1 = rx;
	h5p = (xx + *q * c) * (*dz * *dz + *q * *q + r * xx) / (*dz * c * r * (d1 * d1));
	i5p = *k * 2 * h5p / (h5 * h5 + 1) / c;
	/* Computing 2nd power */
	d1 = r + dt;
	i1p = *k * *dz * (r * s + *p) / (r * c * (d1 * d1)) - t * i5p;
	i4p = *k * (*p - dt * s) / (c * r * (r + dt));
	i3p = (1. - yt * (*p + r * s) / (r * (r + dt) * c)) / (r + dt) - 1. / r;
	i3p = *k * i3p + t * i4p;
	i2p = -(*k) / r - i3p;
	duup[1] = -(*u1) * (*dz * g1p + g2p + i1p * s) / TWO_PI;
	duup[2] = -(*u1) * (yt * g1p + c * g1 + c * g3p + i2p * s) / TWO_PI;
	duup[3] = -(*u1) * (dt * g1p + s * g1 + s * g3p + i4p * s) / TWO_PI;
	duup[4] = -(*u2) * (-(*q) * *p / r3 - i3p * s * c) / TWO_PI;
	duup[5] = -(*u2) * (yt * g6p + c * g6 + c * g2p - i1p * s * c) / TWO_PI;
	duup[6] = -(*u2) * (dt * g6p + s * g6 + s * g2p - i5p * s * c) / TWO_PI;
	duup[7] = *u3 * (*q * g1p - i3p * s * s) / TWO_PI;
	duup[8] = *u3 * (-dt * g6p - s * (g6 + *dz * g1p - g2p) - i1p * s * s) / TWO_PI;
	duup[9] = *u3 * (yt * g6p + g6 * c + c * (g1p * *dz - g2p) - i5p * s * s) / TWO_PI;
	return 0;
}

int ech_(double *u, double *v) {
	int i;
	/* Parameter adjustments */
	--v; --u;

	for (i = 1; i <= 9; ++i) v[i] = u[i];
	return 0;
}

void vtm (double lon0, double lat0) {
	/* Set up an TM projection (extract of GMT_vtm)*/
	double lat2, s2, c2;
	
	lat0 *= D2R;
	lat2 = 2.0 * lat0;
	s2 = sin(lat2);
	c2 = cos(lat2);
	ECC2 = 2 * flattening - flattening * flattening;
	ECC4 = ECC2 * ECC2;
	ECC6 = ECC2 * ECC4;
	one_m_ECC2 = 1.0 - ECC2;
	i_one_m_ECC2 = 1.0 / one_m_ECC2;
	t_c1 = 1.0 - (1.0/4.0) * ECC2 - (3.0/64.0) * ECC4 - (5.0/256.0) * ECC6;
	t_c2 = -((3.0/8.0) * ECC2 + (3.0/32.0) * ECC4 + (25.0/768.0) * ECC6);
	t_c3 = (15.0/128.0) * ECC4 + (45.0/512.0) * ECC6;
	t_c4 = -(35.0/768.0) * ECC6;
	t_e2 = ECC2 * i_one_m_ECC2;
	t_M0 = EQ_RAD * (t_c1 * lat0 + s2 * (t_c2 + c2 * (t_c3 + c2 * t_c4)));
	central_meridian = lon0;
}

void tm (double lon, double lat, double *x, double *y) {
	/* Convert lon/lat to TM x/y (adapted from GMT_tm) */
	double N, T, T2, C, A, M, dlon, tan_lat, A2, A3, A5, lat2, s, c, s2, c2;

	if (fabs (fabs (lat) - 90.0) < SUMOL) {
		M = EQ_RAD * t_c1 * M_PI_2;
		*x = 0.0;
		*y = M;
	}
	else {
		lat *= D2R;
		lat2 = 2.0 * lat;
		s = sin(lat);	s2 = sin(lat2);
		c = cos(lat);	c2 = cos(lat2);
		tan_lat = s / c;
		M = EQ_RAD * (t_c1 * lat + s2 * (t_c2 + c2 * (t_c3 + c2 * t_c4)));
		dlon = lon - central_meridian;
		if (fabs (dlon) > 360.0) dlon += copysign (360.0, -dlon);
		if (fabs (dlon) > 180.0) dlon = copysign (360.0 - fabs (dlon), -dlon);
		N = EQ_RAD / d_sqrt (1.0 - ECC2 * s * s);
		T = tan_lat * tan_lat;
		T2 = T * T;
		C = t_e2 * c * c;
		A = dlon * D2R * c;
		A2 = A * A;	A3 = A2 * A;	A5 = A3 * A2;
		*x = N * (A + (1.0 - T + C) * (A3 * 0.16666666666666666667) + (5.0 - 18.0 * T + T2 + 72.0 * C - 58.0 * t_e2) * (A5 * 0.00833333333333333333));
		A3 *= A;	A5 *= A;
		*y = (M - t_M0 + N * tan_lat * (0.5 * A2 + (5.0 - T + 9.0 * C + 4.0 * C * C) * (A3 * 0.04166666666666666667) + (61.0 - 58.0 * T + T2 + 600.0 * C - 330.0 * t_e2) * (A5 * 0.00138888888888888889)));
	}
}
