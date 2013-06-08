/*--------------------------------------------------------------------
 *	$Id:$
 *
 * \file geoglib_mex.cpp
 * \brief Matlab mex file for evaluating gravity fields
 *
 *	Copyright (c) 2004-2012 by J. Luis
 *
 * 	This program is part of Mirone and is free software and licensed under
 *  the MIT/X11 License. 
 * 
 * 	This program is distributed in the hope that it will be useful,
 * 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 * 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *	Contact info: w3.ualg.pt/~jluis/mirone
 *--------------------------------------------------------------------*/

#include <string>
#include <algorithm>
#include <GeographicLib/Geoid.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/GravityCircle.hpp>
#include <GeographicLib/DMS.hpp>
#include <GeographicLib/Utility.hpp>
#include <mex.h>

#if HAVE_OPENMP
#include <omp.h>
#endif

#define	FALSE	0
#define	TRUE	1
#define CNULL	((char *)NULL)
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#ifndef rint
#define rint(x) (floor((x)+0.5))
#endif
#ifndef irint
#define irint(x) ((int)rint(x))
#endif

using namespace std;
using namespace GeographicLib;

int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (double w, double e, double s, double n);
double ddmmss_to_degree (char *text);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	int		argc = 0, i, j, n, n_row, n_col, n_pts, n_fields = 0, n_arg_no_char = 0;
	int		got_R = FALSE, error = FALSE, circle = FALSE;
	char	**argv, *i_1;
	float	*out4;
	double	*in, *out, *hdr, lon, lat = 0, h = 0, delta_x = 0, delta_y = 0;
	double	west, east, south, north;
	std::string dir;
	std::string model = GravityModel::DefaultGravityName();
	enum {
		GRAVITY = 0,
		DISTURBANCE = 1,
		ANOMALY = 2,
		UNDULATION = 3,
	};
	unsigned mode = UNDULATION;	

#ifdef MIR_TIMEIT
	//#if HAVE_OPENMP
		//double tic = omp_get_wtime( );
	//#else
		clock_t tic = clock();
	//#endif
#endif

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
	argv[0] = "gravity_glib";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			
				case 'R':
					if (decode_R (argv[i], &west, &east, &south, &north)) {
						mexPrintf ("GRAVITY_GLIB: Error in option -R");
						error++;
					}
					got_R = TRUE;
					break;
				case 'A':
					mode = ANOMALY;
					break;
				case 'D':
					mode = DISTURBANCE;
					break;
				case 'G':
					mode = GRAVITY;
					break;
				case 'H':
					mode = UNDULATION;
					break;
				case 'I':
					j = sscanf(&argv[i][2], "%lf/%lf", &delta_x, &delta_y);
					if (j == 1) 
						delta_y = delta_x;
					else if (j != 2) {
						mexPrintf ("GRAVITY_GLIB: Error using option -I");
						error++;
					}
					break;
				case 'c':
					j = sscanf(&argv[i][2], "%lf/%lf", &lat, &h);
					if (j < 1) {
						mexPrintf ("GRAVITY_GLIB: Error using option -c. Need to set at least 'lat'");
						error++;
					}
					circle = TRUE;
					break;
				case 'd':
					dir = &argv[i][2];
					break;
				case 'n':
					model = &argv[i][2];
					break;
				default:
					error = TRUE;
					break;
			}
		}
	}

	if (nrhs == 0 || error)
		mexErrMsgTxt("One input argument required.");
	else if (nlhs > 2)
		mexErrMsgTxt("GRAVITY_GLIB: more than two output arguments specified.");

	if (got_R) {		// A grid request. We need to compute n_row, n_columns and make adjustments
		if (delta_x <= 0 || delta_y <= 0)
			mexErrMsgTxt("GRAVITY_GLIB: Error, option -I not provided or non-sense values transmitted.");
		n_row = irint((north - south) / delta_y) + 1;
		n_col = irint((east - west) / delta_x) + 1;
	}

	if (!got_R) {

		if (!( mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0]) ))
			mexErrMsgTxt("GRAVITY_GLIB: longlat coordinates are not of type double.");

		/* Check that thirth argument contains at least a mx2 table */
		n_pts = mxGetM (prhs[0]);
		n_fields = mxGetN(prhs[0]);
		if (!mxIsNumeric(prhs[0]) || (!circle && mxGetN(prhs[0]) < 2))
			mexErrMsgTxt("GRAVITY_GLIB: first argument must contain the x,y positions.\n");

		/* Read the interpolation points and convert them to double */
		if (mxIsDouble(prhs[0]))
			in = mxGetPr(prhs[0]);
		else if (mxIsSingle(prhs[0]))
			in = (double *)mxGetData(prhs[0]);

		plhs[0] = mxCreateDoubleMatrix(n_pts, mode == UNDULATION ? 1 : 3, mxREAL);
		out = mxGetPr(plhs[0]);
	}
	else {
		plhs[0] = mxCreateNumericMatrix (n_row, n_col, mxSINGLE_CLASS, mxREAL);
		out4 = (float *)mxGetData(plhs[0]);
		if (nlhs == 2) {
			plhs[1] = mxCreateDoubleMatrix(1, 9, mxREAL);
			hdr = mxGetPr(plhs[1]);
			hdr[0] = west;		hdr[1] = west + (n_col - 1) * delta_x;
			hdr[2] = south;		hdr[3] = south + (n_row - 1) * delta_y;
			hdr[6] = 0;			// Grid registration
			hdr[7] = delta_x;	hdr[8] = delta_y;
		}
	}

	const GravityModel g(model, dir);

	unsigned mask = (mode == GRAVITY ? GravityModel::GRAVITY :
						(mode == DISTURBANCE ? GravityModel::DISTURBANCE :
							(mode == ANOMALY ? GravityModel::SPHERICAL_ANOMALY :
								GravityModel::GEOID_HEIGHT))); // mode == UNDULATION
      
	const GravityCircle c(circle ? g.Circle(lat, h, mask) : GravityCircle());

	switch (mode) {
		case GRAVITY:
			for (i = 0; i < n_pts; i++) {
				if (circle)
					c.Gravity(in[i], out[i], out[i+n_pts], out[i+2*n_pts]);		// out = [gx gy gz]
				else {
					if (n_fields == 3) h = in[i+2*n_pts];
					g.Gravity(in[i+n_pts], in[i], h, out[i], out[i+n_pts], out[i+2*n_pts]);
				}
			}
			break;
		case DISTURBANCE:
			double deltax, deltay, deltaz;
			for (i = 0; i < n_pts; i++) {
				if (circle)
					c.Disturbance(in[i], deltax, deltay, deltaz);
				else {
					if (n_fields == 3) h = in[i+2*n_pts];
					g.Disturbance(in[i+n_pts], in[i], h, deltax, deltay, deltaz);
				}
				out[i]         = deltax * 1e5;		// Convert to mGals
				out[i+n_pts]   = deltay * 1e5;
				out[i+2*n_pts] = deltaz * 1e5;
			}
			break;
		case ANOMALY:
			double Dg01, xi, eta;
			for (i = 0; i < n_pts; i++) {
				if (circle)
					c.SphericalAnomaly(in[i], Dg01, xi, eta);
				else {
					if (n_fields == 3) h = in[i+2*n_pts];
					g.SphericalAnomaly(in[i+n_pts], in[i], h, Dg01, xi, eta);
				}
				out[i]         = Dg01 * 1e5;      // Convert to mGals
				out[i+n_pts]   = xi * 3600;       // Convert to arcsecs
				out[i+2*n_pts] = eta * 3600;
			}
			break;
		case UNDULATION:
		default:
			if (!got_R) {
				for (i = 0; i < n_pts; i++)
					out[i] = circle ? c.GeoidHeight(in[i]) : g.GeoidHeight(in[i+n_pts], in[i]);
			}
			else {
#if HAVE_OPENMP
#pragma omp parallel for private(j)
#endif
				for (i = 0; i < n_row; i++) { // Loop over latitudes
					lat = south + i * delta_y;
					GravityCircle c(g.Circle(lat, 0, GravityModel::GEOID_HEIGHT));
					for (j = 0; j < n_col; j++) { // Loop over longitudes
						lon = west + j * delta_x;
						out4[i + j*n_row] = float(c.GeoidHeight(lon));
					}
				}				
			}
			break;
	}

#ifdef MIR_TIMEIT
	//#if HAVE_OPENMP
		//double toc = omp_get_wtime( );
	//#else
		clock_t toc = clock();
	//#endif

	mexPrintf("GRAVITY_GEOGLIB: CPU secs/ticks = %.3f\n", (double)(toc - tic));
#endif

}

/* -------------------------------------------------------------------- */
int decode_R (char *item, double *w, double *e, double *s, double *n) {
	char *text, string[BUFSIZ];
	
	/* Minimalist code to decode option -R extracted from GMT_get_common_args */
	
	int i, error = 0;
	double *p[4];
	
	p[0] = w;	p[1] = e;	p[2] = s;	p[3] = n;
			
	i = 0;
	strcpy (string, &item[2]);
	text = strtok (string, "/");
	while (text) {
		*p[i] = ddmmss_to_degree (text);
		i++;
		text = strtok (CNULL, "/");
	}
	if (item[strlen(item)-1] == 'r')	/* Rectangular box given, but valid here */
		error++;
	if (i != 4 || check_region (*p[0], *p[1], *p[2], *p[3]))
		error++;
	w = p[0];	e = p[1];
	s = p[2];	n = p[3];
	return (error);
}

/* -------------------------------------------------------------------- */
int check_region (double w, double e, double s, double n) {
	/* If region is given then we must have w < e and s < n */
	return ((w >= e || s >= n));
}

/* -------------------------------------------------------------------- */
double ddmmss_to_degree (char *text) {
	int i, colons = 0, suffix;
	double degree, minute, degfrac, second;

	for (i = 0; text[i]; i++) if (text[i] == ':') colons++;
	suffix = (int)text[i-1];	/* Last character in string */
	if (colons == 2) {	/* dd:mm:ss format */
		sscanf (text, "%lf:%lf:%lf", &degree, &minute, &second);
		degfrac = degree + Loc_copysign (minute / 60.0 + second / 3600.0, degree);
	}
	else if (colons == 1) {	/* dd:mm format */
		sscanf (text, "%lf:%lf", &degree, &minute);
		degfrac = degree + Loc_copysign (minute / 60.0, degree);
	}
	else
		degfrac = atof (text);
	if (suffix == 'W' || suffix == 'w' || suffix == 'S' || suffix == 's') degfrac = -degfrac;	/* Sign was given implicitly */
	return (degfrac);
}
