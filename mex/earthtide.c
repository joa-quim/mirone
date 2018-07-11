/* SOLID EARTH TIDE
   Dennis Milbert, 2017
   http://geodesyworld.github.io/SOFTS/solid.htm

   Converted to C with f2c, cleaned and MEXified by Joaquim Luis
*/

#include <math.h>
#include <time.h>

#define I_AM_MEX        /* Build as a MEX */
//#define TIMEIT

#ifdef I_AM_C           /* Build as a stand-alone exe */
#	undef I_AM_MEX
#endif

#ifdef I_AM_MEX
#	include "mex.h"
#	include "mwsize.h"
#else
#	include <stdio.h>
#	include <stdlib.h>
#	include <stdint.h>
#	define mxCalloc calloc
#	define mxMalloc malloc
#	define mxRealloc realloc
#	define mxFree free
#	define mexPrintf(...) fprintf(stderr, __VA_ARGS__);

#	ifndef NAN
#		ifdef _MSC_VER
#			include <ymath.h>
#			define NAN _Nan._Double
#		else
			static const double _NAN = 20;//(HUGE_VAL-HUGE_VAL);
#			define NAN _NAN
#		endif
#	endif
#	define mxGetNaN() (NAN)

#endif
#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2  1.57079632679489661923
#endif
#ifndef TWO_PI
#define TWO_PI  6.28318530717958647692
#endif
#ifndef D2R
#define D2R (M_PI / 180.0)
#endif
#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#define EARTH_RAD 6378137.0		// GRS80
#define ECC2 0.00669438002290341574957

void solid_grd(struct solid_time_series *TS, struct solid_grid *GRD);
void detide(double *, int, double, double *, double *, double *, bool *);
static void civmjd(int, int, int, int, int, double, int *, double*);
static void setjd0(int, int, int);
static void rge(double, double, double *, double *, double *, double, double, double);
static void solid_ts(struct solid_time_series *TS);
static void mjdciv(int, double, int *, int *, int *, int *, int *, double *);
static void geoxyz(double, double, double, double *, double *, double *);
static void sunxyz(int, double, double *, bool *);
static void moonxyz(int, double, double *, bool *);
static void rot1(double, double, double, double, double *, double *, double *);
static void rot3(double, double, double, double, double *, double *, double *);
static void getghar(int, double, double * );
static void st1isem(double *, double *, double *, double, double, double *);
static void st1idiu(double *, double *, double *, double, double, double *);
static void st1l1(double *, double *, double *, double, double, double *);
static void sprod(double *, double *, double *, double *, double *);
static void step2diu(double *, double, double, double *);
static void step2lon(double *, double, double *);
static double enorm8(double *);
static double utc2ttt(double tutc, bool *leapflag);
static double gps2ttt(double tgps);
static double utc2tai(double tutc, bool *leapflag);
static double getutcmtai(double tsec, bool *leapflag);
static double gps2tai(double tgps);
static double tai2tt(double ttai);

struct solid_time_series {	/* To hold the IO when computing a time series */
	int year, month, day, hour, min;
	int length;                 /* Number of minutes of the computed time series */
	double sec;
	double t_inc;
	double lon, lat;
	double *t, *x, *y, *z;		/* To hold up the result data */
};

struct solid_grid {			/* To hold the IO when computing one to three grids */
	int n_rows, n_cols;
	int do_north, do_east, do_up;
	double x_min, x_max, y_min, y_max, x_inc, y_inc;
	float *north, *east, *up;
};

struct {
	int mjd0;
} mjdoff_;

#define mjdoff_1 mjdoff_

double d_mod(const double x, const double y) {
	double quotient;
	if ((quotient = x / y) >= 0)
		quotient = floor(quotient);
	else
		quotient = -floor(-quotient);
	return (x - y * quotient );
}

int Return(int code) {		/* To handle return codes between MEX and standalone code */
#ifdef I_AM_MEX
	mexErrMsgTxt("\n");		/* Most of the cases (no_sys_mem) this instruction was executed already */
#endif
	return(code);
}

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

#ifdef I_AM_MEX
#	define Return(code) {mexErrMsgTxt("\n");}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
#else
#	define Return(code) {return(code);}
int main(int argc, char **argv) {
#endif

	int i;
	struct solid_time_series TS;
	struct solid_grid GRD;
#ifdef TIMEIT
	clock_t tic = clock();
#endif

#ifdef I_AM_MEX
	int dims[2], n_minutes, n_grids = 0;
	char *str;
	static const char *fieldnames[3] = {"north", "east", "up"};
	double *head, *ti, *tf, *ptr;
	mxArray *grd_north = NULL, *grd_east = NULL, *grd_up = NULL, *ts = NULL;

	if (nrhs == 0) {	/* Usage */
		mexPrintf("EarthTide - A Solid tides calculator\n\n");
		mexPrintf("\tOUT = earthtide(time_start, time_end, hdr, '-G[<neu>]|-T[n_minutes]\n\n");
		mexPrintf("'time_start' is a 1x6 vector with intial time in the form of [year month day hour min sec]\n");
		mexPrintf("'time_end' is either a 1x6 vector similar to 'time_start' with end time or an empty array [].\n");
		mexPrintf("Empty end time arrays are used when computing grids (see option -G) or when the -Tn_minutes is used\n");
		mexPrintf("'hdr' is either 1x2 or a 1x6 array with. First case is used when computing time series and the two\n");
		mexPrintf("elements contain the [lon lat] coordinates where to compute the time series. Second form is used\n");
		mexPrintf("when computing grid(s) and must contain [x_min x_max y_min y_max x_inc y_inc] all in geogs.\n\n");
		mexPrintf("Fourth argument is is either a string with -T[n_minutes] to compute time series or -G[neu]\n");
		mexPrintf("The n_minutes arg is the number of minutes to compute the tide starting at 'time_start' and\n");
		mexPrintf("takes precedence over 'time_end' (in case it was provided instead of passing the [] array).\n");
		mexPrintf("For time series the OUT array contains a Mx4 array whe columns 1 to 4 contain [minuts north east up]\n\n");
		mexPrintf("The -G[neu] option tells the program to compute 1 to 3 grids with the n(orth), e(ast), u(p)\n");
		mexPrintf("component(s). If only one component is selected (or none, which defaults to UP) the OUT variable\n");
		mexPrintf("is a MxN array. If more than one component is requested, the OUT var is a struct with fields\n");
		mexPrintf("'north', 'east', 'up'. The non requested components have corresponding field set to empty array.\n");
		return;
	}

	if (nrhs != 4 || !mxIsChar(prhs[3]))
		mexErrMsgTxt("Fourth argument must be a -T[n] or -G<neu> string option\n");

	ti = mxGetPr(prhs[0]);
	if (mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 6)
		mexErrMsgTxt("Initial time array (first arg) must be a 1x6 array\n");

	TS.year = (int)ti[0];		TS.month = (int)ti[1];		TS.day = (int)ti[2];
	TS.hour = (int)ti[3];		TS.min = (int)ti[4];		TS.sec = ti[5];

	str = (char *)mxArrayToString(prhs[3]);
	if (str[1] == 'T') {				/* Time-series */
		n_minutes = atoi(&str[2]);
		if (n_minutes)
			TS.length = n_minutes;
		else {		/* Must read final time in prhs[2] */
			char d1[32] = {""}, d2[32] = {""};
			double sd1, sd2;
			mxArray *lhs[1], *rhs[2], *rhs2[2];
			tf = mxGetPr(prhs[1]);
			if (mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 6)
				mexErrMsgTxt("Final time array (Second arg) must be a 1x6 array\n");

			/* OK, now we must find the difference in minutes between the two dates.
			   Not simple so we resort to call a MEX, but ... */
			tf = mxGetPr(prhs[1]);
			lhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
			sprintf(d1, "%d-%0.2d-%0.2d %0.2d:%0.2d:%0.2d", 
					(int)ti[0], (int)ti[1], (int)ti[2], (int)ti[3], (int)ti[4], (int)ti[5]);
			rhs[0] = mxCreateString(d1);
			rhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
			ptr = mxGetPr(rhs[1]);		ptr[0] = 31;
			mexCallMATLAB(1,lhs,2,rhs,"DateStr2Num");
			ptr = mxGetPr(lhs[0]);
			sd1 = ptr[0];
			sprintf(d2, "%d-%0.2d-%0.2d %0.2d:%0.2d:%0.2d", 
					(int)tf[0], (int)tf[1], (int)tf[2], (int)tf[3], (int)tf[4], (int)tf[5]);
			rhs2[0] = mxCreateString(d2);
			rhs2[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
			ptr = mxGetPr(rhs2[1]);		ptr[0] = 31;
			mexCallMATLAB(1,lhs,2,rhs2,"DateStr2Num");
			ptr = mxGetPr(lhs[0]);
			sd2 = ptr[0];
			TS.length = (int)ceil((sd2 - sd1) * 24 * 60);
		}
		head = mxGetPr(prhs[2]);
		TS.lon = head[0];	TS.lat = head[1];
	}
	else if (str[1] == 'G') {			/* Grid(s) */
		int k = 2;
		GRD.do_north = GRD.do_east = GRD.do_up = 0;		/* Need to initialize these */
		while (str[k]) {
			if (str[k] == 'n')
				{GRD.do_north = 1;	n_grids++;}
			else if (str[k] == 'e')
				{GRD.do_east = 1;	n_grids++;}
			else if (str[k] == 'u')
				{GRD.do_up = 1;		n_grids++;}
			k++;
		}
		if (GRD.do_north == 0 && GRD.do_east == 0)	/* If -G was used then at least UP component is computed */
			{GRD.do_up = 1;	n_grids = 1;}

		head = mxGetPr(prhs[2]);
		if (mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != 6)
			mexErrMsgTxt("Header array (Third arg) must be a 1x6 array\n");

		GRD.x_min = head[0];		GRD.x_max = head[1];
		GRD.y_min = head[2];		GRD.y_max = head[3];
		GRD.x_inc = head[4];		GRD.y_inc = head[5];
		GRD.n_cols = (int)(GRD.x_max - GRD.x_min) / GRD.x_inc + 1;		/* Use grid registration */
		GRD.n_rows = (int)(GRD.y_max - GRD.y_min) / GRD.y_inc + 1;
	}
	else
		mexErrMsgTxt("Option method (Fourth argument must be -T or -G\n");

	if (TS.year < 1901 || TS.year > 2099)
		mexErrMsgTxt("Dates must be within the [1901 2099] period.\n");

	dims[0] = GRD.n_rows;	dims[1] = GRD.n_cols;
	if (GRD.do_north) {
		grd_north = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		GRD.north = (float *)mxGetData(grd_north);
	}
	if (GRD.do_east) {
		grd_east = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		GRD.east = (float *)mxGetData(grd_east);
	}
	if (GRD.do_up) {
		grd_up = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
		GRD.up = (float *)mxGetData(grd_up);
	}

	if (n_grids == 1) {
		solid_grd(&TS, &GRD);		/* DO THE WORK */
		if (grd_north)
			plhs[0] = grd_north;
		else if (grd_east)
			plhs[0] = grd_east;
		else
			plhs[0] = grd_up;
	}
	else if (n_grids > 1) {
		solid_grd(&TS, &GRD);		/* DO THE WORK */
		/* Output in a structure. Non requested component will be an empty field */
		plhs[0] = mxCreateStructMatrix(1, 1, 3, fieldnames);
		mxSetField(plhs[0], 0, "north", grd_north);
		mxSetField(plhs[0], 0, "east", grd_east);
		mxSetField(plhs[0], 0, "up", grd_up);
	}
	else {			/* Time series */
		dims[0] = TS.length;	dims[1] = 4;
		ts = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
		ptr = (double *)mxGetPr(ts);
		TS.t = ptr;					/* Because this mem is column-major we can split the columns like this */
		TS.x = &ptr[TS.length];
		TS.y = &ptr[2*TS.length];
		TS.z = &ptr[3*TS.length];
		plhs[0] = ts;
		solid_ts(&TS);
	}
#else
	TS.year = 2017;
	TS.month = 1;
	TS.day = 1;
	TS.hour = 0;
	TS.min = 0;
	TS.sec = 0;
	TS.length = 60*24;		/* minutes */
	TS.lon = -7;
	TS.lat = 37;

	TS.t = (double *)malloc(TS.length * sizeof(double));
	TS.x = (double *)malloc(TS.length * sizeof(double));
	TS.y = (double *)malloc(TS.length * sizeof(double));
	TS.z = (double *)malloc(TS.length * sizeof(double));
	if (TS.t == NULL) {
		mexPrintf("Deu merda a alocar\n");
		Return(0);
	}

	GRD.x_min = 0;		GRD.x_max = 360;
	GRD.y_min = -90;	GRD.y_max = 90;
	GRD.x_inc = 1;		GRD.y_inc = 1;
	GRD.n_rows = 181;	GRD.n_cols = 361;
	GRD.do_north = 0;	GRD.do_east = 0;	GRD.do_up = 1;

	if (GRD.do_north)
		GRD.north = (float *)malloc((size_t)GRD.n_rows * GRD.n_cols * sizeof(float));

	if (GRD.do_east)
		GRD.east = (float *)malloc((size_t)GRD.n_rows * GRD.n_cols * sizeof(float));

	if (GRD.do_up)
		GRD.up = (float *)malloc((size_t)GRD.n_rows * GRD.n_cols * sizeof(float));

	solid_ts(&TS);
	//solid_grd(&TS, &GRD);
#endif

#ifdef TIMEIT
	mexPrintf("CPU ticks = %.3f\tCPS = %d\n", (double)(clock() - tic), CLOCKS_PER_SEC);
#endif

#ifndef I_AM_MEX
	//for (i = 0; i < TS.length; i++)
	for (i = 0; i < 5; i++)
		fprintf(stdout, "%f\t%f\t%f\t%f\n", TS.t[i], TS.x[i], TS.y[i], TS.z[i]);

	return 0;
#endif
}

/* ----------------------------------------------------------------------- */
void solid_grd(struct solid_time_series *TS, struct solid_grid *GRD) {
	int row, col, mjd;
	bool leapflag;
	size_t ij_n = 0, ij_e = 0, ij_u = 0, n_inc = 0, e_inc = 0, u_inc = 0;
	float *grd_n, *grd_e, *grd_u;
	double fmjd, xsta[3], rsun[3], tdel2, etide[3], rmoon[3];
	double lon, lat, ut, vt, wt;

	/* Select which indices to increment based on user selection */
	/* Use the trick of not incrementing the indices of unwanted arrays to avoid IF branches inside the loops */
	if (GRD->do_north) {
		grd_n = GRD->north;
		n_inc = 1;
	}
	else
		grd_n = (float *)malloc(1 * sizeof(float));

	if (GRD->do_east) {
		grd_e = GRD->east;
		e_inc = 1;
	}
	else
		grd_e = (float *)malloc(1 * sizeof(float));

	if (GRD->do_up) {
		grd_u = GRD->up;
		u_inc = 1;
	}
	else
		grd_u = (float *)malloc(1 * sizeof(float));

	leapflag = false;                       /* false means flag not raised */
	for (col = 0; col < GRD->n_cols; col++) {
		lon = (GRD->x_min + col * GRD->x_inc) * D2R;
		for (row = 0; row < GRD->n_rows; row++) {
			lat = (GRD->y_min + row * GRD->y_inc) * D2R;
			geoxyz(lat, lon, 0, &xsta[0], &xsta[1], &xsta[2]);
			civmjd(TS->year, TS->month, TS->day, TS->hour, TS->min, TS->sec, &mjd, &fmjd);
			mjdciv(mjd, fmjd, &TS->year, &TS->month, &TS->day, &TS->hour, &TS->min, &TS->sec);	/* normalize civil time */
			setjd0(TS->year, TS->month, TS->day);

			sunxyz(mjd, fmjd, rsun, &leapflag);
			moonxyz(mjd, fmjd, rmoon, &leapflag);
			detide(xsta, mjd, fmjd, rsun, rmoon, etide, &leapflag);
			/* determine local geodetic horizon components (topocentric) */
			rge(lat, lon, &ut, &vt, &wt, etide[0], etide[1], etide[2]);		/* tide vect */
			grd_n[ij_n] = (float)ut;
			grd_e[ij_e] = (float)vt;
			grd_u[ij_u] = (float)wt;
			ij_n += n_inc;
			ij_e += e_inc;
			ij_u += u_inc;
		}
	}
	if (leapflag)
		mexPrintf("time crossed leap seconds table boundaries. Boundary edge used instead.");

	/* Free these that were never used anyway */
	if (!GRD->do_north) free(grd_n);
	if (!GRD->do_east) free(grd_e);
	if (!GRD->do_up) free(grd_u);
}

/* ----------------------------------------------------------------------- */
void solid_ts(struct solid_time_series *TS) {
	/* iyr	year    [1980-2018] */
	/* imo	month number [1-12] */
	/* idy	day          [1-31] */
	/* lat	Lat. (pos N.) [- 90, +90] */
	/* lon	Lon. (pos E.) [-360,+360] */
	int k, mjd;
	bool leapflag;
	double d__2, ut, vt, wt;
	double fmjd, xsta[3], rsun[3], tdel2, etide[3], rmoon[3];

	/* position of observing point (positive East) */
	if (TS->lon < 0)    TS->lon += 360;
	if (TS->lon >= 360) TS->lon += -360;
 
	TS->lat *= D2R;
	TS->lon *= D2R;
	geoxyz(TS->lat, TS->lon, 0, &xsta[0], &xsta[1], &xsta[2]);

	/* ** here comes the sun  (and the moon)  (go, tide!) */
	civmjd(TS->year, TS->month, TS->day, TS->hour, TS->min, TS->sec, &mjd, &fmjd);
	mjdciv(mjd, fmjd, &TS->year, &TS->month, &TS->day, &TS->hour, &TS->min, &TS->sec);	/* normalize civil time */
	setjd0(TS->year, TS->month, TS->day);
	tdel2 = 1.0 / 24 / 60;			/* 1 minute steps */
	for (k = 0; k < TS->length; k++) {
		leapflag = false;                       /* false means flag not raised */
		sunxyz(mjd, fmjd, rsun, &leapflag);      /* mjd/fmjd in UTC */
		moonxyz(mjd, fmjd, rmoon, &leapflag);
		detide(xsta, mjd, fmjd, rsun, rmoon, etide, &leapflag);
		/* determine local geodetic horizon components (topocentric) */
		rge(TS->lat, TS->lon, &ut, &vt, &wt, etide[0], etide[1], etide[2]);		/* tide vect */
		d__2 = TS->sec - 0.001;
		mjdciv(mjd, fmjd + 1.1574074074074074e-8, &TS->year, &TS->month, &TS->day, &TS->hour, &TS->min, &d__2);
		TS->t[k] = TS->hour * 3600 + TS->min * 60 + TS->sec;
		TS->x[k] = ut;
		TS->y[k] = vt;
		TS->z[k] = wt;
		fmjd += tdel2;
		fmjd = (int)(round(fmjd * 86400)) / 86400.0;	/* force 1 sec. granularity */
	}

	if (leapflag)
		mexPrintf("time crossed leap seconds table boundaries. Boundary edge used instead.");
}

/* ----------------------------------------------------------------------- */
static double getutcmtai(double tsec, bool *leapflag) {

	int mjd0t;
	double ttsec, tai_utc = 0;

	/*  get utc - tai (s) */
	/*  "Julian Date Converter" */
	/*  http://aa.usno.navy.mil/data/docs/JulianDate.php */
	/*  parameter(MJDUPPER=58299)    !*** upper limit, leap second table, 2018jun30 */
	/* upper limit, leap second table, */
	/* lower limit, leap second table, */
	/* leap second table limit flag */
	/* leap second table limit flag */
	/* clone for tests (and do any rollover) */
	ttsec = tsec;
	mjd0t = mjdoff_1.mjd0;

	while (ttsec >= 86400.) {
		ttsec += -86400.;
		mjd0t++;
	}

	while (ttsec < 0) {
		ttsec += 86400.;
		mjd0t--;
	}

	/*  test upper table limit (upper limit set by bulletin C memos) */
	if (mjd0t > 58664) {
		*leapflag = true;		/* true means flag *IS* raised */
		return -37;				/* return the upper table valu */
	}
	else if (mjd0t < 41317) {	/*  test lower table limit */
		*leapflag = true;		/* true means flag *IS* raised */
		return -10;				/* return the lower table valu */
	}

	/*  http://maia.usno.navy.mil/ser7/tai-utc.dat */
	/* 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0s */
	/* 1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0s */
	/* 1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0s */
	/* 1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0s */
	/* 1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0s */
	/* 1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0s */
	/* 1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0s */
	/* 1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0s */
	/* 1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0s */
	/* 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0s */
	/* 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0s */
	/* 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0s */
	/* 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0s */
	/* 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0s */
	/* 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0s */
	/* 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0s */
	/* 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0s */
	/* 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0s */
	/* 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0s */
	/* 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0s */
	/* 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0s */
	/* 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0s */
	/* 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0s */
	/* 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0s */
	/* 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0s */
	/* 2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0s */
	/* 2015 JUL  1 =JD 2457204.5  TAI-UTC=  36.0s */
	/* 2017 JAN  1 =JD 2457754.5  TAI-UTC=  37.0s */
	/*  other leap second references at: */
	/*  http://hpiers.obspm.fr/eoppc/bul/bulc/Leap_Second_History.dat */
	/*  http://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.dat */
	/* test against newest leaps first */
	if (mjd0t >= 57754)			/* 2017 JAN 1 = 57754 */
		tai_utc = 37.;
	else if (mjd0t >= 57204)    /* 2015 JUL 1 = 57204 */
		tai_utc = 36.;
	else if (mjd0t >= 56109)    /* 2012 JUL 1 = 56109 */
		tai_utc = 35.;
	else if (mjd0t >= 54832)    /* 2009 JAN 1 = 54832 */
		tai_utc = 34.;
	else if (mjd0t >= 53736)    /* 2006 JAN 1 = 53736 */
		tai_utc = 33.;
	else if (mjd0t >= 51179)    /* 1999 JAN 1 = 51179 */
		tai_utc = 32.;
	else if (mjd0t >= 50630)    /* 1997 JUL 1 = 50630 */
		tai_utc = 31.;
	else if (mjd0t >= 50083)    /* 1996 JAN 1 = 50083 */
		tai_utc = 30.;
	else if (mjd0t >= 49534)    /* 1994 JUL 1 = 49534 */
		tai_utc = 29.;
	else if (mjd0t >= 49169)    /* 1993 JUL 1 = 49169 */
		tai_utc = 28.;
	else if (mjd0t >= 48804)    /* 1992 JUL 1 = 48804 */
		tai_utc = 27.;
	else if (mjd0t >= 48257)    /* 1991 JAN 1 = 48257 */
		tai_utc = 26.;
	else if (mjd0t >= 47892)    /* 1990 JAN 1 = 47892 */
		tai_utc = 25.;
	else if (mjd0t >= 47161)    /* 1988 JAN 1 = 47161 */
		tai_utc = 24.;
	else if (mjd0t >= 46247)    /* 1985 JUL 1 = 46247 */
		tai_utc = 23.;
	else if (mjd0t >= 45516)    /* 1983 JUL 1 = 45516 */
		tai_utc = 22.;
	else if (mjd0t >= 45151)    /* 1982 JUL 1 = 45151 */
		tai_utc = 21.;
	else if (mjd0t >= 44786)    /* 1981 JUL 1 = 44786 */
		tai_utc = 20.;
	else if (mjd0t >= 44239)    /* 1980 JAN 1 = 44239 */
		tai_utc = 19.;
	else if (mjd0t >= 43874)    /* 1979 JAN 1 = 43874 */
		tai_utc = 18.;
	else if (mjd0t >= 43509)    /* 1978 JAN 1 = 43509 */
		tai_utc = 17.;
	else if (mjd0t >= 43144)    /* 1977 JAN 1 = 43144 */
		tai_utc = 16.;
	else if (mjd0t >= 42778)    /* 1976 JAN 1 = 42778 */
		tai_utc = 15.;
	else if (mjd0t >= 42413)    /* 1975 JAN 1 = 42413 */
		tai_utc = 14.;
	else if (mjd0t >= 42048)    /* 1974 JAN 1 = 42048 */
		tai_utc = 13.;
	else if (mjd0t >= 41683)    /* 1973 JAN 1 = 41683 */
		tai_utc = 12.;
	else if (mjd0t >= 41499)    /* 1972 JUL 1 = 41499 */
		tai_utc = 11.;
	else if (mjd0t >= 41317)    /* 1972 JAN 1 = 41317 */
		tai_utc = 10.;
	/* return utc - tai (in seconds) */
	return -tai_utc;
}

/* *********************************************************************** */
/* ** new supplemental time functions ************************************ */
/* *********************************************************************** */
/* ----------------------------------------------------------------------- */
static double tai2tt(double ttai) {
	/*  convert tai (sec) to terrestrial time (sec) */
	/*  http://tycho.usno.navy.mil/systime.html */
	return ttai + 32.184;
}

/* ----------------------------------------------------------------------- */
static double utc2tai(double tutc, bool *leapflag) {
	/* convert utc (sec) to tai (sec) */
	return tutc - getutcmtai(tutc, leapflag);
}

/* ----------------------------------------------------------------------- */
static double utc2ttt(double tutc, bool *leapflag) {
	double ttai;

	/* convert utc (sec) to terrestrial time (sec) */
	ttai = utc2tai(tutc, leapflag);
	return tai2tt(ttai);
}

/* ----------------------------------------------------------------------- */
static void sprod(double *x, double *y, double *scal, double *r1, double *r2) {
	/*  computation of the scalar-product of two vectors and their norms */
	/*  input:   x(i),i=1,2,3  -- components of vector x */
	/*           y(i),i=1,2,3  -- components of vector y */
	/*  output:  scal          -- scalar product of x and y */
	/*           r1,r2         -- lengths of the two vectors x and y */

	*r1 = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	*r2 = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
	*scal =    x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

/* ----------------------------------------------------------------------- */
static double enorm8(double *a) {
	/* compute euclidian norm of a vector (of length 3) */
	return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/* ----------------------------------------------------------------------- */
static  void st1idiu(double *xsta, double *xsun, double *xmon, double fac2sun, double fac2mon, double *xcorsta) {
	/* this subroutine gives the out-of-phase corrections induced by */
	/* mantle inelasticity in the diurnal band */
	/*  input: xsta,xsun,xmon,fac2sun,fac2mon */
	/* output: xcorsta */
	double dhi = -0.0025;
	double dli = -7e-4;
	double de, dn, dr, rsta, rmon, rsun, cosla, demon, sinla, dnmon, desun, drmon, dnsun, drsun;
	double cosphi, sinphi, cos2phi, inv_rsun2, inv_rmon2;

	rsta = enorm8(&xsta[0]);
	sinphi = xsta[2] / rsta;
	cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
	cos2phi = cosphi * cosphi - sinphi * sinphi;
	sinla = xsta[1] / cosphi / rsta;
	cosla = xsta[0] / cosphi / rsta;
	rmon = enorm8(&xmon[0]);
	rsun = enorm8(&xsun[0]);
	inv_rsun2 = 1 / (rsun * rsun);
	inv_rmon2 = 1 / (rmon * rmon);
	drsun = dhi * -3 * sinphi  * cosphi  * fac2sun * xsun[2]  * (xsun[0] * sinla - xsun[1] * cosla) * inv_rsun2;
	drmon = dhi * -3 * sinphi  * cosphi  * fac2mon * xmon[2]  * (xmon[0] * sinla - xmon[1] * cosla) * inv_rmon2;
	dnsun = dli * -3 * cos2phi * fac2sun * xsun[2] * (xsun[0] * sinla - xsun[1] * cosla) * inv_rsun2;
	dnmon = dli * -3 * cos2phi * fac2mon * xmon[2] * (xmon[0] * sinla - xmon[1] * cosla) * inv_rmon2;
	desun = dli * -3 * sinphi  * fac2sun * xsun[2] * (xsun[0] * cosla + xsun[1] * sinla) * inv_rsun2;
	demon = dli * -3 * sinphi  * fac2mon * xmon[2] * (xmon[0] * cosla + xmon[1] * sinla) * inv_rmon2;
	dr = drsun + drmon;
	dn = dnsun + dnmon;
	de = desun + demon;
	xcorsta[0] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
	xcorsta[1] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
	xcorsta[2] = dr * sinphi + dn * cosphi;
}

/* ----------------------------------------------------------------------- */
static void st1isem(double *xsta, double *xsun, double *xmon, double fac2sun, double fac2mon, double *xcorsta) {
	/* this subroutine gives the out-of-phase corrections induced by */
	/* mantle inelasticity in the diurnal band */
	/*  input: xsta,xsun,xmon,fac2sun,fac2mon */
	/* output: xcorsta */

	const double dhi = -0.0022;
	const double dli = -7e-4;
	double costwola, sintwola, de, dn, dr, rsta, rmon, rsun, cosla, demon, sinla, dnmon, desun, drmon, dnsun, drsun;
	double cosphi, sinphi, cosphi2, inv_rsun2, inv_rmon2, dif_xsun2, dif_xmon2, t;

	rsta = enorm8(&xsta[0]);
	sinphi = xsta[2] / rsta;
	cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
	cosphi2 = cosphi * cosphi;
	sinla = xsta[1] / cosphi / rsta;
	cosla = xsta[0] / cosphi / rsta;
	costwola = cosla * cosla - sinla * sinla;
	sintwola = cosla * 2 * sinla;
	rmon = enorm8(&xmon[0]);
	rsun = enorm8(&xsun[0]);
	inv_rsun2 = 1 / (rsun * rsun);
	inv_rmon2 = 1 / (rmon * rmon);
	dif_xsun2 = xsun[0] * xsun[0] - xsun[1] * xsun[1];
	dif_xmon2 = xmon[0] * xmon[0] - xmon[1] * xmon[1];
	t = dhi * -0.75 * cosphi2;
	drsun = t * fac2sun * (dif_xsun2 * sintwola - xsun[0] * 2 * xsun[1] * costwola) * inv_rsun2;
	drmon = t * fac2mon * (dif_xmon2 * sintwola - xmon[0] * 2 * xmon[1] * costwola) * inv_rmon2;
	t = dli * 1.5 * sinphi;
	dnsun = t * cosphi * fac2sun * (dif_xsun2 * sintwola - xsun[0] * 2 * xsun[1] * costwola) * inv_rsun2;
	dnmon = t * cosphi * fac2mon * (dif_xmon2 * sintwola - xmon[0] * 2 * xmon[1] * costwola) * inv_rmon2;
	t = dli * -1.5 * cosphi;
	desun = t * fac2sun * (dif_xsun2 * costwola + xsun[0] * 2 * xsun[1] * sintwola) * inv_rsun2;
	demon = t * fac2mon * (dif_xmon2 * costwola + xmon[0] * 2 * xmon[1] * sintwola) * inv_rmon2;
	dr = drsun + drmon;
	dn = dnsun + dnmon;
	de = desun + demon;
	xcorsta[0] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
	xcorsta[1] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
	xcorsta[2] = dr * sinphi + dn * cosphi;
}

/* ----------------------------------------------------------------------- */
static void st1l1(double *xsta, double *xsun, double *xmon, double fac2sun, double fac2mon, double *xcorsta) {
	/* this subroutine gives the corrections induced by the latitude dependence */
	/* given by l^(1) in mahtews et al (1991) */
	/*  input: xsta,xsun,xmon,fac3sun,fac3mon */
	/* output: xcorsta */

	double l1d = .0012;
	double l1sd = .0024;
	double costwola, sintwola, l1, de, dn, rsta, rmon, rsun, cosla, demon, sinla, dnmon, desun, dnsun;
	double cosphi, sinphi, cosphi2, sinphi2, inv_rsun2, inv_rmon2, dif_xsun2, dif_xmon2, t;

	rsta = enorm8(&xsta[0]);
	sinphi = xsta[2] / rsta;
	cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
	sinla = xsta[1] / cosphi / rsta;
	cosla = xsta[0] / cosphi / rsta;
	rmon = enorm8(&xmon[0]);
	rsun = enorm8(&xsun[0]);
	/* ** for the diurnal band */
	l1 = l1d;
	sinphi2 = sinphi * sinphi;
	cosphi2 = cosphi * cosphi;
	inv_rsun2 = 1 / (rsun * rsun);
	inv_rmon2 = 1 / (rmon * rmon);
	dnsun = -l1 * sinphi2 * fac2sun * xsun[2] * (xsun[0] * cosla + xsun[1] * sinla) * inv_rsun2;
	dnmon = -l1 * sinphi2 * fac2mon * xmon[2] * (xmon[0] * cosla + xmon[1] * sinla) * inv_rmon2;
	t = l1 * sinphi * (cosphi2 - sinphi2);
	desun = t * fac2sun * xsun[2] * (xsun[0] * sinla - xsun[1] * cosla) * inv_rsun2;
	demon = t * fac2mon * xmon[2] * (xmon[0] * sinla - xmon[1] * cosla) * inv_rmon2;
	de = (desun + demon) * 3.;
	dn = (dnsun + dnmon) * 3.;
	xcorsta[0] = -de * sinla - dn * sinphi * cosla;
	xcorsta[1] =  de * cosla - dn * sinphi * sinla;
	xcorsta[2] =  dn * cosphi;
	/* ** for the semi-diurnal band */
	l1 = l1sd;
	costwola = cosla * cosla - sinla * sinla;
	sintwola = cosla * 2 * sinla;
	dif_xsun2 = xsun[0] * xsun[0] - xsun[1] * xsun[1];
	dif_xmon2 = xmon[0] * xmon[0] - xmon[1] * xmon[1];
	t = -l1 / 2 * sinphi  * cosphi;
	dnsun = t * fac2sun * (dif_xsun2 * costwola + xsun[0] * 2 * xsun[1] * sintwola) * inv_rsun2;
	dnmon = t * fac2mon * (dif_xmon2 * costwola + xmon[0] * 2 * xmon[1] * sintwola) * inv_rmon2;
	t = -l1 / 2 * sinphi2 * cosphi;
	desun = t * fac2sun * (dif_xsun2 * sintwola - xsun[0] * 2 * xsun[1] * costwola) * inv_rsun2;
	demon = t * fac2mon * (dif_xmon2 * sintwola - xmon[0] * 2 * xmon[1] * costwola) * inv_rmon2;
	de = (desun + demon) * 3;
	dn = (dnsun + dnmon) * 3;
	xcorsta[0] = xcorsta[0] - de * sinla - dn * sinphi * cosla;
	xcorsta[1] = xcorsta[1] + de * cosla - dn * sinphi * sinla;
	xcorsta[2] += dn * cosphi;
}

/* ----------------------------------------------------------------------- */
static void step2diu(double *xsta, double fhr, double t, double *xcorsta) {
	/* last change:  vd   17 may 00   1:20 pm */
	/* these are the subroutines for the step2 of the tidal corrections. */
	/* they are called to account for the frequency dependence */
	/* of the love numbers. */

	static double datdi[279]	/* was [9][31] */ = { -3.,0.,2.,0.,0.,
		-.01,-.01,0.,0.,-3.,2.,0.,0.,0.,-.01,-.01,0.,0.,-2.,0.,1.,-1.,0.,
		-.02,-.01,0.,0.,-2.,0.,1.,0.,0.,-.08,0.,.01,.01,-2.,2.,-1.,0.,0.,
		-.02,-.01,0.,0.,-1.,0.,0.,-1.,0.,-.1,0.,0.,0.,-1.,0.,0.,0.,0.,
		-.51,0.,-.02,.03,-1.,2.,0.,0.,0.,.01,0.,0.,0.,0.,-2.,1.,0.,0.,.01,
		0.,0.,0.,0.,0.,-1.,0.,0.,.02,.01,0.,0.,0.,0.,1.,0.,0.,.06,0.,0.,
		0.,0.,0.,1.,1.,0.,.01,0.,0.,0.,0.,2.,-1.,0.,0.,.01,0.,0.,0.,1.,
		-3.,0.,0.,1.,-.06,0.,0.,0.,1.,-2.,0.,1.,0.,.01,0.,0.,0.,1.,-2.,0.,
		0.,0.,-1.23,-.07,.06,.01,1.,-1.,0.,0.,-1.,.02,0.,0.,0.,1.,-1.,0.,
		0.,1.,.04,0.,0.,0.,1.,0.,0.,-1.,0.,-.22,.01,.01,0.,1.,0.,0.,0.,0.,
		12.,-.78,-.67,-.03,1.,0.,0.,1.,0.,1.73,-.12,-.1,0.,1.,0.,0.,2.,0.,
		-.04,0.,0.,0.,1.,1.,0.,0.,-1.,-.5,-.01,.03,0.,1.,1.,0.,0.,1.,.01,
		0.,0.,0.,1.,1.,0.,1.,-1.,-.01,0.,0.,0.,1.,2.,-2.,0.,0.,-.01,0.,0.,
		0.,1.,2.,0.,0.,0.,-.11,.01,.01,0.,2.,-2.,1.,0.,0.,-.01,0.,0.,0.,
		2.,0.,-1.,0.,0.,-.02,.02,0.,.01,3.,0.,0.,0.,0.,0.,.01,0.,.01,3.,
		0.,0.,1.,0.,0.,.01,0.,0. };

	int i, j;
	double h, t2, t3, t4, cosphi2, sinphi2, sin_tf, cos_tf;
	double p, s, de, dn, dr, pr, ps, zla, tau, zns, rsta, cosla, sinla, thetaf, cosphi, sinphi;

	/* ** note, following table is derived from dehanttideinelMJD.f (2000oct30 16:10) */
	/* ** has minor differences from that of dehanttideinel.f (2000apr17 14:10) */
	/* ** D.M. edited to strictly follow published table 7.5a (2006aug08 13:46) */
	/* ** cf. table 7.5a of IERS conventions 2003 (TN.32, pg.82) */
	/* ** columns are s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op) */
	/* ** units of mm */
	/* ****----------------------------------------------------------------------- */
	/* ***** -2., 0., 1., 0., 0.,-0.08,-0.05, 0.01,-0.02,      !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ****----------------------------------------------------------------------- */
	/* ***** -1., 0., 0.,-1., 0.,-0.10,-0.05, 0.0 ,-0.02,      !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ***** -1., 0., 0., 0., 0.,-0.51,-0.26,-0.02,-0.12,      !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ****----------------------------------------------------------------------- */
	/* *****  0., 0., 1., 0., 0., 0.06, 0.02, 0.0 , 0.01,      !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ****----------------------------------------------------------------------- */
	/* *****  1.,-2., 0., 0., 0.,-1.23,-0.05, 0.06,-0.06,      !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ****----------------------------------------------------------------------- */
	/* *****  1., 0., 0., 0., 0.,12.02,-0.45,-0.66, 0.17,      !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* *****  1., 0., 0., 1., 0., 1.73,-0.07,-0.10, 0.02,      !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ****----------------------------------------------------------------------- */
	/* *****  1., 1., 0., 0.,-1.,-0.50, 0.0 , 0.03, 0.0,       !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ****----------------------------------------------------------------------- */
	/* *****  0., 1., 0., 1.,-1.,-0.01, 0.0 , 0.0 , 0.0,       !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* ****----------------------------------------------------------------------- */
	/* *****  1., 2., 0., 0., 0.,-0.12, 0.01, 0.01, 0.0,       !*** original entry */
	/* ****----------------------------------------------------------------------- */
	/* *** table 7.5a */
	/* *** v.dehant 2 */
	t2 = t * t;
	t3 = t * t2;
	t4 = t2 * t2;
	s = 218.31664563 + 481267.88194 * t - .0014663889 * t2 + 1.85139e-6 * t3;
	tau = fhr * 15. + 280.4606184 + t * 36000.7700536 + t * 3.8793e-4 * t - t3 * 2.58e-8 - s;
	pr = t * 1.396971278f + t * 3.08889e-4f * t + t3 * 2.1e-8f + t2 * 7e-9f;
	s += pr;
	h = t * 36000.7697489 + 280.46645 + t * 3.0322222e-4 * t + t3 * 2e-8f - t2 * 6.54e-9f;
	p = t * 4069.01363525 + 83.35324312 - t * .01032172222 * t - t3 * 1.24991e-5 + t2 * 5.263e-8;
	zns = t * 1934.13626197 + 234.95544499 - t * .00207561111 * t - t3 * 2.13944e-6 + t2 * 1.65e-8;
	ps = t * 1.71945766667 + 282.93734098 + t * 4.5688889e-4 * t - t3 * 1.778e-8 - t2 * 3.34e-9;
	/* ** reduce angles to between 0 and 360 */
	s = d_mod(s, 360.);
	tau = d_mod(tau, 360.);
	h = d_mod(h, 360.);
	p = d_mod(p, 360.);
	zns = d_mod(zns, 360.);
	ps = d_mod(ps, 360.);
	rsta = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1] + xsta[2] * xsta[2]);
	sinphi = xsta[2] / rsta;
	cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
	sinphi2 = sinphi * sinphi;
	cosphi2 = cosphi * cosphi;
	cosla = xsta[0] / cosphi / rsta;
	sinla = xsta[1] / cosphi / rsta;
	zla = atan2(xsta[1], xsta[0]);
	for (i = 0; i < 3; i++) xcorsta[i] = 0;

	for (j = 1; j <= 31; j++) {
		thetaf = (tau + datdi[j * 9 - 9] * s + datdi[j * 9 - 8] * h + datdi[j * 9 - 7] * p + datdi[j * 9 - 6] * zns + datdi[j * 9 - 5] * ps) * D2R + zla;
		sin_tf = sin(thetaf);		cos_tf = cos(thetaf);
		dr     = datdi[j * 9 - 4] * 2 * sinphi * cosphi * sin_tf + datdi[j * 9 - 3] * 2 * sinphi * cosphi * cos_tf;
		dn     = datdi[j * 9 - 2] * (cosphi2 - sinphi2) * sin_tf + datdi[j * 9 - 1] * (cosphi2 - sinphi2) * cos_tf;
		/* following correction by V.Dehant to match eq.16b, p.81, 2003 Conventions */
		/* de=datdi(8,j)*sinphi*cos(thetaf+zla)+ */
		de = datdi[j * 9 - 2] * sinphi * cos_tf - datdi[j * 9 - 1] * sinphi * sin_tf;
		xcorsta[0] += (dr * cosla * cosphi - de * sinla - dn * sinphi * cosla);
		xcorsta[1] += (dr * sinla * cosphi + de * cosla - dn * sinphi * sinla);
		xcorsta[2] += (dr * sinphi + dn * cosphi);
	}

	for (i = 0; i < 3; ++i) xcorsta[i] /= 1e3;
}

/* ----------------------------------------------------------------------- */
static void step2lon(double *xsta, double t, double *xcorsta) {
	/* cf. table 7.5b of IERS conventions 2003 (TN.32, pg.82) */
	/* columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op) */
	/* IERS cols.= s,h,p,N',ps, dR(ip),dR(op),dT(ip),dT(op) */
	/* units of mm */

	static double datdi[45]	/* was [9][5] */ = { 0.,0.,0.,1.,0.,.47,.23,
		.16,.07,0.,2.,0.,0.,0.,-.2,-.12,-.11,-.05,1.,0.,-1.,0.,0.,-.11,
		-.08,-.09,-.04,2.,0.,0.,0.,0.,-.13,-.11,-.15,-.07,2.,0.,0.,1.,0.,
		-.05,-.05,-.06,-.03 };

	int i, j;
	double h, t2, t3, t4;
	double p, s, de, dn, dr, pr, ps, zns, rsta, cosla, sinla, thetaf, cosphi, sinphi;

	t2 = t * t;		t3 = t * t2;		t4 = t2 * t2;
	s = 218.31664563 + 481267.88194 * t - .0014663889 * t2 + 1.85139e-6 * t3;
	pr = t * 1.396971278f + t2 * 3.08889e-4f + t3 * 2.1e-8f + t4 * 7e-9f;
	s += pr;
	h   = t * 36000.7697489 + 280.46645    + t2 * 3.0322222e-4 + t3 * 2e-8f      - t4 * 6.54e-9f;
	p   = t * 4069.01363525 + 83.35324312  - t2 * .01032172222 - t3 * 1.24991e-5 + t4 * 5.263e-8;
	zns = t * 1934.13626197 + 234.95544499 - t2 * .00207561111 - t3 * 2.13944e-6 + t4 * 1.65e-8;
	ps = t * 1.71945766667 + 282.93734098  + t2 * 4.5688889e-4 - t3 * 1.778e-8   - t4 * 3.34e-9;
	rsta = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1] + xsta[2] * xsta[2]);
	sinphi = xsta[2] / rsta;
	cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
	cosla = xsta[0] / cosphi / rsta;
	sinla = xsta[1] / cosphi / rsta;
	/* ** reduce angles to between 0 and 360 */
	s = d_mod(s, 360.);
	/* **** tau=dmod(tau,360.d0)       !*** tau not used here--09jul28 */
	h = d_mod(h, 360.);
	p = d_mod(p, 360.);
	zns = d_mod(zns, 360.);
	ps = d_mod(ps, 360.);
	for (i = 0; i < 3; i++) xcorsta[i] /= 1e3;

	/* **             1 2 3 4   5   6      7      8      9 */
	/* ** columns are s,h,p,N',ps, dR(ip),dT(ip),dR(op),dT(op) */
	for (j = 1; j <= 5; j++) {
		thetaf = (datdi[j * 9 - 9] * s + datdi[j * 9 - 8] * h + datdi[j * 9 - 7] * p + datdi[j * 9 - 6] * zns + datdi[j * 9 - 5] * ps) * D2R;
		dr = datdi[j * 9 - 4] * (sinphi * sinphi * 3 - 1) / 2. * cos(thetaf) + datdi[j * 9 - 2] * (sinphi * sinphi * 3 - 1) / 2. * sin(thetaf);
		dn = datdi[j * 9 - 3] * (cosphi * sinphi * 2) * cos(thetaf) + datdi[j * 9 - 1] * (cosphi * sinphi * 2) * sin(thetaf);
		de = 0.;
		xcorsta[0] += dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
		xcorsta[1] += dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
		xcorsta[2] += dr * sinphi + dn * cosphi;
	}
	for (i = 0; i < 3; i++) xcorsta[i] /= 1e3;
}

/* ------------------------------------------------------------------------------- */
static void detide(double *xsta, int mjd, double fmjd, double *xsun, double *xmon, double *dxtide, bool *leapflag) {
	/* Computation of tidal corrections of station displacements caused
	 * by lunar and solar gravitational attraction.  UTC version.
	 * step 1 (here general degree 2 and 3 corrections +
	 *         call st1idiu + call st1isem + call st1l1)
	 *   + step 2 (call step2diu + call step2lon + call step2idiu)
	 * It has been decided that the step 3 un-correction for permanent tide
	 * would *not* be applied in order to avoid jump in the reference frame
	 * (this step 3 must added in order to get the mean tide station position
	 * and to be conformed with the iag resolution.)
	 * inputs:
	 *   xsta(i),i=1,2,3   -- geocentric position of the station (ITRF/ECEF)
	 *   xsun(i),i=1,2,3   -- geoc. position of the sun (ECEF)
	 *   xmon(i),i=1,2,3   -- geoc. position of the moon (ECEF)
	 *   mjd,fmjd          -- modified julian day (and fraction) (in GPS time)
	 * ***old calling sequence*****************************************************
	 *   dmjd               -- time in mean julian date (including day fraction)
	 *   fhr=hr+zmin/60.+sec/3600.   -- hr in the day
	 * outputs:
	 *   dxtide(i),i=1,2,3           -- displacement vector (ITRF)
	 *   flag              -- leap second table limit flag, false:flag not raised
	 * Author iers 1996 :  V. Dehant, S. Mathews and J. Gipson
	 *    (test between two subroutines)
	 * Author iers 2000 :  V. Dehant, C. Bruyninx and S. Mathews
	 *    (test in the bernese program by C. Bruyninx)
	 * Created:  96/03/23 (see above)
	 * Modified from dehanttideinelMJD.f by Dennis Milbert 2006sep10
	 * Bug fix regarding fhr (changed calling sequence, too)
	 * Modified to reflect table 7.5a and b IERS Conventions 2003
	 * Modified to use TT time system to call step 2 functions
	 * Sign correction by V.Dehant to match eq.16b, p.81, Conventions
	 * applied by Dennis Milbert 2007may05
	 * UTC version by Dennis Milbert 2018june01
	 */

	int i;
	double h20 = .6078;
	double l20 = .0847;
	double h3 = .292;
	double l3 = .015;
	double re_over_rsun, re_over_rmon, re;
	double mass_ratio_moon, mass_ratio_sun, t, t2, h2, l2, fhr, scm, scs;
	double rsta, rmon, rsun, p2mon, p3mon, x2mon, x3mon, p2sun, p3sun, x2sun, x3sun, scmon;
	double scsun, cosphi, dmjdtt, fmjdtt, tsectt, fac2mon, fac3mon, fac2sun, fac3sun;
	double tsecutc, xcorsta[3], inv_rsta;

	/* nominal second degree and third degree love numbers and shida numbers */

	/* internal support for new calling sequence */
	/* first, convert UTC time into TT time (and, bring leapflag into variable) */
	tsecutc = fmjd * 86400.;			/* UTC time (sec of day) */
	tsectt = utc2ttt(tsecutc, leapflag);			/* TT  time (sec of day) */
	fmjdtt = tsectt / 86400.;			/* TT  time (fract. day) */
	dmjdtt = mjd + fmjdtt;
	/*  commented line was live code in dehanttideinelMJD.f */
	/*  changed on the suggestion of Dr. Don Kim, UNB -- 09mar21 */
	/*  Julian date for 2000 January 1 00:00:00.0 UT is  JD 2451544.5 */
	/*  MJD         for 2000 January 1 00:00:00.0 UT is MJD   51544.0 */
	/*  t=(dmjdtt-51545.d0)/36525.d0                !*** days to centuries, TT */
	/*  float MJD in TT */
	t = (dmjdtt - 51544.) / 36525.;			/* days to centuries, TT */
	fhr = (dmjdtt - (int) dmjdtt) * 24.;	/* hours in the day, TT */
	/* ** scalar product of station vector with sun/moon vector */
	sprod(&xsta[0], &xsun[0], &scs, &rsta, &rsun);
	sprod(&xsta[0], &xmon[0], &scm, &rsta, &rmon);
	scsun = scs / rsta / rsun;
	scmon = scm / rsta / rmon;
	/* ** computation of new h2 and l2 */
	cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
	t2 = (1. - cosphi * 1.5 * cosphi);
	h2 = h20 - t2 * 6e-4;
	l2 = l20 + t2 * 2e-4;
	/* ** p2-term */
	t2 = (h2 / 2. - l2) * 3;
	p2sun = t2 * scsun * scsun - h2 / 2.;
	p2mon = t2 * scmon * scmon - h2 / 2.;
	/* ** p3-term */
	t2 = (h3 - l3 * 3.) * 2.5;
	p3sun = t2 * (scsun * scsun * scsun) + (l3 - h3) * 1.5 * scsun;
	p3mon = t2 * (scmon * scmon * scmon) + (l3 - h3) * 1.5 * scmon;
	/* ** term in direction of sun/moon vector */
	x2sun = l2 * 3. * scsun;
	x2mon = l2 * 3. * scmon;
	x3sun = l3 * 1.5 * (scsun * 5 * scsun - 1);
	x3mon = l3 * 1.5 * (scmon * 5 * scmon - 1);
	/* ** factors for sun/moon */
	mass_ratio_sun = 332945.943062;
	mass_ratio_moon = 0.012300034;
	re = 6378136.55;
	re_over_rsun = re / rsun;
	fac2sun = mass_ratio_sun * re * re_over_rsun * re_over_rsun * re_over_rsun;
	re_over_rmon = re / rmon;
	fac2mon = mass_ratio_moon * re * re_over_rmon * re_over_rmon * re_over_rmon;
	fac3sun = fac2sun * re_over_rsun;
	fac3mon = fac2mon * re_over_rmon;
	/* ** total displacement */
	inv_rsta = 1 / rsta;
	for (i = 0; i < 3; i++) {
		t2 = xsta[i] * inv_rsta;
		dxtide[i] = fac2sun * (x2sun * xsun[i] / rsun + p2sun * t2) +
					fac2mon * (x2mon * xmon[i] / rmon + p2mon * t2) +
					fac3sun * (x3sun * xsun[i] / rsun + p3sun * t2) +
					fac3mon * (x3mon * xmon[i] / rmon + p3mon * t2);
	}
	xcorsta[0] = xcorsta[1] = xcorsta[2] = 0;
	/* ** corrections for the out-of-phase part of love numbers */
	/* **     (part h_2^(0)i and l_2^(0)i ) */
	/* ** first, for the diurnal band */
	st1idiu(&xsta[0], &xsun[0], &xmon[0], fac2sun, fac2mon, xcorsta);
	dxtide[0] += xcorsta[0];
	dxtide[1] += xcorsta[1];
	dxtide[2] += xcorsta[2];
	/* ** second, for the semi-diurnal band */
	st1isem(&xsta[0], &xsun[0], &xmon[0], fac2sun, fac2mon, xcorsta);
	dxtide[0] += xcorsta[0];
	dxtide[1] += xcorsta[1];
	dxtide[2] += xcorsta[2];
	/* ** corrections for the latitude dependence of love numbers (part l^(1) ) */
	st1l1(&xsta[0], &xsun[0], &xmon[0], fac2sun, fac2mon, xcorsta);
	dxtide[0] += xcorsta[0];
	dxtide[1] += xcorsta[1];
	dxtide[2] += xcorsta[2];
	/* ** consider corrections for step 2 */
	/* ** corrections for the diurnal band: */
	/* **  first, we need to know the date converted in julian centuries */
	/* **  this is now handled at top of code   (also convert to TT time system) */
	/* **** t=(dmjd-51545.)/36525. */
	/* **** fhr=dmjd-int(dmjd)             !*** this is/was a buggy line (day vs. hr) */
	/* **  second, the diurnal band corrections, */
	/* **   (in-phase and out-of-phase frequency dependence): */
	step2diu(&xsta[0], fhr, t, xcorsta);
	dxtide[0] += xcorsta[0];
	dxtide[1] += xcorsta[1];
	dxtide[2] += xcorsta[2];
	/* **  corrections for the long-period band, */
	/* **   (in-phase and out-of-phase frequency dependence): */
	step2lon(&xsta[0], t, xcorsta);
	dxtide[0] += xcorsta[0];
	dxtide[1] += xcorsta[1];
	dxtide[2] += xcorsta[2];
	/* ** consider corrections for step 3 */
	/* -----------------------------------------------------------------------
	 * The code below is commented to prevent restoring deformation
	 * due to permanent tide.  All the code above removes
	 * total tidal deformation with conventional Love numbers.
	 * The code above realizes a conventional tide free crust (i.e. ITRF).
	 * This does NOT conform to Resolution 16 of the 18th General Assembly
	 * of the IAG (1983).  This resolution has not been implemented by
	 * the space geodesy community in general (c.f. IERS Conventions 2003).
	 * -----------------------------------------------------------------------
	 * ** uncorrect for the permanent tide  (only if you want mean tide system)
	 * **   pi=3.141592654
	 * **   sinphi=xsta(3)/rsta
	 * **   cosphi=dsqrt(xsta(1)**2+xsta(2)**2)/rsta
	 * **   cosla=xsta(1)/cosphi/rsta
	 * **   sinla=xsta(2)/cosphi/rsta
	 * **   dr=-dsqrt(5./4./pi)*h2*0.31460*(3./2.*sinphi**2-0.5)
	 * **   dn=-dsqrt(5./4./pi)*l2*0.31460*3.*cosphi*sinphi
	 * **   dxtide(1)=dxtide(1)-dr*cosla*cosphi+dn*cosla*sinphi
	 * **   dxtide(2)=dxtide(2)-dr*sinla*cosphi+dn*sinla*sinphi
	 * **   dxtide(3)=dxtide(3)-dr*sinphi      -dn*cosphi
	 */
}

/* ******************************************************************************* */
static void getghar(int mjd, double fmjd, double *ghar) {
	/* convert mjd/fmjd in UTC time to Greenwich hour angle (in radians)
	 * "satellite orbits: models, methods, applications" montenbruck & gill(2000)
	 * section 2.3.1, pg. 33
	 * need UTC to get sidereal time  ("astronomy on the personal computer", 4th ed)
	 *                               (pg.43, montenbruck & pfleger, springer, 2005)
	 */
	int i;
	double ghad, fmjdutc, tsecutc;

	tsecutc = fmjd * 86400.;			/* UTC time (sec of day) */
	fmjdutc = tsecutc / 86400.;			/* UTC time (fract. day) */
	/*  d = MJD - 51544.5d0                               !*** footnote
	 *  greenwich hour angle for J2000  (12:00:00 on 1 Jan 2000)
	 *  ghad = 100.46061837504d0 + 360.9856473662862d0*d  !*** eq. 2.85 (+digits)
	 */
	ghad = (mjd - 51544 + (fmjdutc - 0.5)) * 360.9856473662862 + 280.46061837504;	/* days since J2000 */

	/* normalize to 0-360 and convert to radians */
	i = (int)(ghad / 360.);
	*ghar = (ghad - i * 360.) * D2R;

	while (*ghar > TWO_PI)
		*ghar -= TWO_PI;

	while (*ghar < 0.)
		*ghar += TWO_PI;
}

/* ----------------------------------------------------------------------- */
static void rge(double lat, double lon, double *u, double *v, double *w, double x, double y, double z) {
	/* given a rectangular cartesian system (x,y,z) */
	/* compute a geodetic h cartesian sys   (u,v,w) */
	static double cb, cl, sb, sl;

	sb = sin(lat);
	cb = cos(lat);
	sl = sin(lon);
	cl = cos(lon);
	*u = -sb * cl * x - sb * sl * y + cb * z;
	*v = -sl * x + cl * y;
	*w = cb * cl * x + cb * sl * y + sb * z;
}

/* ----------------------------------------------------------------------- */
static void rot1(double theta, double x, double y, double z, double *u, double *v, double *w) {
	/* ** rotate coordinate axes about 1 axis by angle of theta radians */
	/* ** x,y,z transformed into u,v,w */
	static double c, s;

	s = sin(theta);
	c = cos(theta);
	*u = x;
	*v = c * y + s * z;
	*w = c * z - s * y;
}

/* ----------------------------------------------------------------------- */
static void rot3(double theta, double x, double y, double z, double *u, double *v, double *w) {
	/* rotate coordinate axes about 3 axis by angle of theta radians */
	/* x,y,z transformed into u,v,w */
	static double c, s;

	s = sin(theta);
	c = cos(theta);
	*u = c * x + s * y;
	*v = c * y - s * x;
	*w = z;
}

/* ------------------------------------------------------------------------------- */
static void moonxyz(int mjd, double fmjd, double *rm, bool *leapflag) {
	/* get low-precision, geocentric coordinates for moon (ECEF)
	 * UTC Version
	 * input:  mjd/fmjd, is Modified Julian Date (and fractional) in UTC time
	 * output: rm, is geocentric lunar position vector [m] in ECEF
	 *		   lflag  -- leap second table limit flag,  false:flag not raised
	 * 1."satellite orbits: models, methods, applications" montenbruck & gill(2000)
	 * section 3.3.2, pg. 72-73
	 * 2."astronomy on the personal computer, 4th ed." montenbruck & pfleger (2005)
	 * section 3.2, pg. 38-39  routine MiniMoon
	 */
	double d2, d__, f, q, t, t1, t2, t3, el, el0;
	double rm1, rm2, rm3, elp, rse, ghar, oblir, tjdtt;
	double cselat, selatd,  selond, fmjdtt, tsectt, tsecutc;

	/* ** use TT for lunar ephemerides */
	tsecutc = fmjd * 86400;	        		/* UTC time (sec of day) */
	tsectt  = utc2ttt(tsecutc, leapflag);	/* TT  time (sec ofday)  */
	fmjdtt  = tsectt / 86400.;				/* TT  time (fract. day) */

	/* julian centuries since 1.5 january 2000 (J2000) */
	/*   (note: also low precision use of mjd --> tjd) */
	tjdtt = mjd + fmjdtt + 2400000.5;		/* Julian Date, TT */
	t = (tjdtt - 2451545.) / 36525.;		/*  julian centuries, TT */

	/* el0 -- mean longitude of Moon (deg) */
	/* el  -- mean anomaly of Moon (deg) */
	/* elp -- mean anomaly of Sun  (deg) */
	/* f   -- mean angular distance of Moon from ascending node (deg) */
	/* d   -- difference between mean longitudes of Sun and Moon (deg) */
	/* equations 3.47, p.72 */
	el0 = t * 481267.88088 + 218.31617 - t * 1.3972f;
	el  = t * 477198.86753 + 134.96292;
	elp = t * 35999.04944  + 357.52543;
	f   = t * 483202.01873 + 93.27283;
	d__ = t * 445267.11135 + 297.85027;
	d2  = 2 * d__;
	/* ** longitude w.r.t. equinox and ecliptic of year 2000 */
	selond = el0 + sin(el * D2R) * 22640.0/3600. + sin((el + el) * D2R) * 769./3600. -
			 sin((el - d2) * D2R) * 4586.0/3600. + sin(d2 * D2R) * 2370.0/3600. -
			 sin(elp * D2R) * 668.0/3600. - sin((f + f) * D2R) * 412.0/3600. -
			 sin((el + el - d2) * D2R) * 212.0/3600. - sin((el + elp - d2) * D2R) * 206.0/3600. +
			 sin((el + d2) * D2R) * 192.0/3600. - sin((elp - d2) * D2R) * 165.0/3600. +
			 sin((el - elp) * D2R) * 148.0/3600. - sin(d__ * D2R) * 125.0/3600. -
			 sin((el + elp) * D2R) * 110.0/3600. - sin((f + f - d2) * D2R) * 55.0/3600.;
	/* latitude w.r.t. equinox and ecliptic of year 2000 */
	/*  eq 3.48, p.72 */
	q = sin((f + f) * D2R) * 412.0/3600. + sin(elp * D2R) * 541.0/3600.;
	/*  temporary ter */
	selatd = sin((f + selond - el0 + q) * D2R) * 18520.0/3600. -
			 sin((f - d2) * D2R) * 526.0/3600. + sin((el + f - d2) * D2R) * 44.0/3600. -
			 sin((-el + f - d2) * D2R) * 31.0/3600. - sin((-el - el + f) * D2R) * 25.0/3600. -
			 sin((elp + f - d2) * D2R) * 23.0/3600. + sin((-el + f) * D2R) * 21.0/3600. +
			 sin((-elp + f - d2) * D2R) * 11.0/3600.;
	/* distance from Earth center to Moon (m) */
	/*  eq 3.49, p.72 */
	rse = 3.85e8 - cos(el * D2R) * 2.0905e7 - cos((d2 - el) * D2R) * 3.699e6 - cos(d2 * D2R) * 2.956e6
		- cos((el + el) * D2R) * 5.7e5 + cos((el + el - d2) * D2R) * 2.46e5 - cos((elp - d2) * D2R) *
		2.05e5 - cos((el + d2) * D2R) * 1.71e5 - cos((el + elp - d2) * D2R) * 1.52e5;
	/* convert spherical ecliptic coordinates to equatorial cartesian */
	/* precession of equinox wrt. J2000   (p.71) */
	/*  eq 3.50, p.72 */
	selond += t * 1.3972;
	/* position vector of moon (mean equinox & ecliptic of J2000) (EME2000, ICRF) */
	/*                         (plus long. advance due to precession -- eq. above) */
	selatd *= D2R;
	selond *= D2R;
	oblir = 23.43929111 * D2R;			/* obliquity of the J2000 eclipti */
	cselat = cos(selatd);
	t1 = rse * cos(selond) * cselat;		/* meters          !*** eq. 3.51, */
	t2 = rse * sin(selond) * cselat;		/* meters          !*** eq. 3.51, */
	t3 = rse * sin(selatd);					/* meters          !*** eq. 3.51, */
	rot1(-oblir, t1, t2, t3, &rm1, &rm2, &rm3);
	/* convert position vector of moon to ECEF  (ignore polar motion/LOD) */
	/*  eq. 3.51, */
	getghar(mjd, fmjd, &ghar);						/* sec 2.3.1, */
	rot3(ghar, rm1, rm2, rm3, &rm[0], &rm[1], &rm[2]); /* eq. 2.89, */
}

/* ******************************************************************************* */
static void sunxyz(int mjd, double fmjd, double *rs, bool *leapflag) {
	/* get low-precision, geocentric coordinates for sun (ECEF)
	 * input, mjd/fmjd, is Modified Julian Date (and fractional) in UTC time
	 * output, rs, is geocentric solar position vector [m] in ECEF
	 *      	  lflag  -- leap second table limit flag,  false:flag not raised 
	 * 1."satellite orbits: models, methods, applications" montenbruck & gill(2000)
	 * section 3.3.2, pg. 70-71
	 * 2."astronomy on the personal computer, 4th ed." montenbruck & pfleger (2005)
	 * section 3.2, pg. 39  routine MiniSun
	 */
	double r__, t, em, em2, rs1, rs2, rs3, obe;
	double ghar, opod, slon, emdeg, slond, tjdtt, sslon;
	double fmjdtt, tsectt, tsecutc;

	/* ** mean elements for year 2000, sun ecliptic orbit wrt. Earth */
	obe = 23.43929111 * D2R;			/* obliquity of the J2000 ecliptic */
	opod = 282.94;
	/*  use TT for solar ephemerides */
	/*  RAAN + arg.peri.  (deg.) */
	tsecutc = fmjd * 86400.;			/* UTC time (sec of */
	tsectt = utc2ttt(tsecutc, leapflag);/* TT  time (sec of */
	fmjdtt = tsectt / 86400.;
	/* julian centuries since 1.5 january 2000 (J2000) */
	/*   (note: also low precision use of mjd --> tjd) */
	/*  TT  time (fract. */
	tjdtt = mjd + fmjdtt + 2400000.5;	/* Julian Date, TT */
	t = (tjdtt - 2451545.) / 36525.;	/* julian centuries, */
	emdeg = t * 35999.049 + 357.5256;	/* degrees */
	em = emdeg * D2R;
	em2 = em + em;
	/* ** series expansions in mean anomaly, em   (eq. 3.43, p.71) */
	/* *** radians */
	r__ = (149.619 - cos(em) * 2.499 - cos(em2) * 0.021) * 1e9; /* *** m. */
	slond = opod + emdeg + (sin(em) * 6892 + sin(em2) * 72) / 3600.;	/* precession of equinox wrt. J2000   (p.71) */
	slond += t * 1.3972;
	/* position vector of sun (mean equinox & ecliptic of J2000) (EME2000, ICRF) */
	/*                        (plus long. advance due to precession -- eq. above) */
	slon = slond * D2R;
	sslon = sin(slon);
	rs1 = r__ * cos(slon);				/* meters  !*** eq. 3.46, */
	rs2 = r__ * sslon * cos(obe);		/* meters  !*** eq. 3.46, */
	rs3 = r__ * sslon * sin(obe);
	/* ** convert position vector of sun to ECEF  (ignore polar motion/LOD) */
	/* meters             !*** eq. 3.46, */
	getghar(mjd, fmjd, &ghar);			/* sec 2.3.1, */
	rot3(ghar, rs1, rs2, rs3, &rs[0], &rs[1], &rs[2]);		/* eq. 2.89, */
}

/* ----------------------------------------------------------------------- */
static void geoxyz(double lat, double lon, double eht, double *x, double *y, double *z) {
	/* convert geodetic lat, long, ellip ht. to x,y,z */
	double w, w2, en, cla, sla, t;

	sla = sin(lat);
	cla = cos(lat);
	w2 = 1. - ECC2 * sla * sla;
	w = sqrt(w2);
	en = EARTH_RAD / w;
	t = (en + eht) * cla;
	*x = t * cos(lon);
	*y = t * sin(lon);
	*z = (en * (1. - ECC2) + eht) * sla;
}


/* *********************************************************************** */
/* ** time conversion **************************************************** */
/* *********************************************************************** */
static void setjd0(int iyr, int imo, int idy) {
	/* set the integer part of a modified julian date as epoch, mjd0
	   the modified julian day is derived from civil time as in civmjd()
	   allows single number expression of time in seconds w.r.t. mjd0 */
	static int m, y, it1, it2, mjd;

	if (imo <= 2) {
		y = iyr - 1;
		m = imo + 12;
	} else {
		y = iyr;
		m = imo;
	}
	it1 = (int) (y * 365.25);
	it2 = (int) ((m + 1) * 30.6001);
	mjd = it1 + it2 + idy - 679019;
	/* ** now set the epoch for future time computations */
	mjdoff_1.mjd0 = mjd;
}

/* *********************************************************************** */
static void civmjd(int iyr, int imo, int idy, int ihr, int imn, double sec, int *mjd, double *fmjd) {
	/* convert civil date to modified julian date */
	/* imo in range 1-12, idy in range 1-31 */
	/* only valid in range mar-1900 thru feb-2100     (leap year protocols) */
	/* ref: hofmann-wellenhof, 2nd ed., pg 34-35 */
	/* operation confirmed against table 3.3 values on pg.34 */
	int m, y, it1, it2;

	if (imo <= 2) {
		y = iyr - 1;
		m = imo + 12;
	}
	else {
		y = iyr;
		m = imo;
	}
	it1 = (int) (y * 365.25);
	it2 = (int) ((m + 1) * 30.6001);
	*mjd = it1 + it2 + idy - 679019;
	*fmjd = (ihr * 3600 + imn * 60 + sec) / 86400.;
}

static void mjdciv(int mjd, double fmjd, int *iyr, int *imo, int *idy, int *ihr, int *imn, double *sec) {
	/* convert modified julian date to civil date */
	/* imo in range 1-12, idy in range 1-31 */
	/* only valid in range mar-1900 thru feb-2100 */
	/* ref: hofmann-wellenhof, 2nd ed., pg 34-35 */
	/* operation confirmed for leap years (incl. year 2000) */
	static int ia, ib, ic, id, ie, it1, it2, it3;
	static double rjd, tmp;

	rjd = mjd + fmjd + 2400000.5;
	ia = (int)(rjd + .5);
	ib = ia + 1537;
	ic = (int)((ib - 122.1) / 365.25);
	id = (int)(ic * 365.25);
	ie = (int)((ib - id) / 30.6001);
	/* the fractional part of a julian day is fractional mjd + 0.5
	   therefore, fractional part of julian day + 0.5 is fractional mjd */
	it1  = (int)(ie * 30.6001);
	*idy = (int)(ib - id - it1 +fmjd);
	it2  = (int)(ie / 14.);
	*imo = ie - 1 - it2 * 12;
	it3  = (int)((*imo + 7) / 10.);
	*iyr = ic - 4715 - it3;
	tmp  = fmjd * 24.;
	*ihr = (int)tmp;
	tmp  = (tmp - *ihr) * 60.;
	*imn = (int)tmp;
	*sec = (tmp - *imn) * 60.;
}
