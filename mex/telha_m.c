/*--------------------------------------------------------------------
 *	$Id:$
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
 * 
 * read a file of lon, lat, z and ...
 * 
 * Author:	J. Luis
 * Date:	24 Janeiro, 2000
 * Revision:	30/11/00 (some options slightly changed)
 *          	 2/12/00 Added -C optin which allows to select between 
 *          	         different isochron time scales
 *          	 6/12/00 (bug in -C) 
 *          	 26/12/04
 *          	 18/09/05 Realy implemented -B option. It must be called like that:
 *			  [out_x,out_y] = telha_m(linha, opt_E, opt_I, '-B', opt_T, opt_N);
 *				or
 *			  [out_x,out_y,first_mag] = telha_m(linha, opt_E, opt_I, '-B', opt_T, opt_N);
 *
 *			  Where out_x & out_y are [npts x n_flow] matrices
 *          	 22/11/05 Take also into account the brick's portions between t_stat & upper_age
 *			  and their's closest neighbors in the chronological table.
 *
 *			  Changed the default time scale for that of Cand & Kent 95
 *			  The old Cand scale was renamed isoc_o
 *          	 01/07/08 Killed all the globals
 *
 *          	 22/07/10 Convert latitudes to geocentric -> rotations -> back to geodetics
 */

#include "telha.h"
#include "mex.h"
#include <string.h>

int read_time_scale (FILE *fp_ts, struct ISOC_SCALE *isoc_scale);
int int_perf (double *c, int i_min, int n_flow, int k, int age_flow, int multi_flow,
	      double ddeg, double filter_width, struct ISOC_SCALE *isoc_scale);
int intp_lin (double *x, double *y, int n, int m, double *u, double *v, int mode);
int spotter_init (char *file, struct EULER **p, int flowline, double t_max);
int spotter_backtrack (double xp[], double yp[], double tp[], int np, struct EULER p[], int ns, double t_zero);
int do_the_filter (double *fm, int n_rows, int n_cols, double half_width, double *mag_intp, int n_f_wts,
		   int half_n_f_wts, double dt, double t_start, double t_stop, double *f_wt, double *f_in_m);
int flow_in_meters (double *x, double *y, int n_pts, double *f_in_m);
void no_sys_mem (char *where, int n);
int intp_along_seg(int n_seg, double dl);
double	gaussian_weight (double radius, double half_width);
double	set_up_filter(double t_i, double t_f, int n_pts, double filter_width,
		      double *dt, double *t_start, double *t_stop, int *half_n_f_wts, 
		      int *n_f_wts, double *f_wt);


/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	struct	EULER *p;		/* Pointer to array of stage poles */
	struct	ISOC_SCALE *isoc_scale;	/* Pointer to array of time scale */
	int	error = FALSE, brick = FALSE, age_flow = FALSE;
	int	multi_flow = FALSE, first_seg = TRUE, invert_rot = FALSE, linear_age = FALSE;
	int	np_in_seg, switch_xy = FALSE, first = TRUE, n_flow, isoc_c = TRUE;
	int	n_seg;			/* Number of segments in input data */
	int	m, i_min = 0, int_seg = FALSE, ts_file = FALSE, isoc_p = FALSE, isoc_o = FALSE;
	int	data_in_input = FALSE, ages_given = FALSE, fake_origin = FALSE;
	int     i, ii, j, jj, k, l, ll, ndata, n_alloc, ix = 0, iy = 1, ns;
	int	n_stages;		/* Number of stage poles */
	int	n_isoc;			/* Number of isochrons in time scale */
	int 	n_chunk;		/* Total length or array returned by libeuler functions */
	int	argc, n_arg_no_char = 0, np_in, np_out;
	double	in[2], lon, lat, *tmp, *out, *out2, *out_n_data, *out_n_seg, *out_n_flow;
	double	upper_age = 0.0;	/* Extend oldest age back to this time, in Ma */
	double	t_zero = 0.0;		/* Current age in Ma */
	double	age = 0.0;		/* Age of isochron, in Ma */
	double	age_shift = 0.0;	/* Fake the chron table by shifting it by this amount, in Ma */
	double	*c;			/* Array com os pts de todos os tijolos de cada seg */
	double	filter_width = 0.001;	/* Full width of filter in user's units */
	double	dl = 0.05;		/* Delta resolution in degrees for along segment interpolation */
	double	ddeg = 0.005;		/* Delta resolution in degrees for flow interpolation */
	double	age_lin_start, age_lin_stop, age_lin_inc, *first_mag;
	double	ecc2_1 = (1 - 0.0818191908426215 * 0.0818191908426215);	/* Parameter for convertion to geocentrics in WGS84 */
	char	line[MAXCHAR], *newseg = ">";
	char	*euler_file = NULL;	/* Name pointer for file with stage poles */
	char	*t_scale = NULL;	/* Name pointer for file with time scale */
	char	**argv;
	FILE    *fp = NULL, *fp_ts = NULL;
	
	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(!mxIsChar(prhs[i])) {
			argc--;
			n_arg_no_char++;	/* Number of arguments that have a type other than char */
		}
	}
	argc++;					/* to account for the program's name to be inserted in argv[0] */

	/* get the length of the input string */
	argv=(char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "telha";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

        for (i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        switch (argv[i][1]) {
 		
				/* Common parameters */
			
				case '\0':
					error = TRUE;
					break;
				case ':':
					switch_xy = TRUE;
					break;
				case 'A':
					age_flow = TRUE;
					break;
				case 'B':
					brick   = TRUE;
					break;
				case 'C':
					if (strlen (&argv[i][2]) == 1) { /* use one of default time scales */
						if (strcmp (&argv[i][2], "P") == 0) /* Patriat */
							isoc_p = TRUE;
						else if (strcmp (&argv[i][2], "C") == 0) /* Cande & Kent 95 */
							isoc_c = TRUE;
						else if (strcmp (&argv[i][2], "O") == 0) /* Cande (old scale given by Patriat)*/
							isoc_o = TRUE;
						else
							mexPrintf ("SYNTAX ERROR -C option: Must choose only between -CP (Patriat) or -CC (Cande)\n");
					}
					else {
						t_scale = &argv[i][2];
						ts_file = TRUE;
					}
					break;
				case 'D':
					ddeg = atof (&argv[i][2]);
					break;
				case 'E':	/* File with stage poles */
					euler_file  = &argv[i][2];
					break;
				case 'F':
					filter_width = atof (&argv[i][2]);
					break;
				case 'G':
                                        multi_flow = TRUE;
					break;
				case 'I':
                                        invert_rot = TRUE;
					break;
				case 'L':
					sscanf (&argv[i][2], "%lf/%lf/%lf", &age_lin_start, &age_lin_stop, &age_lin_inc);
                                        linear_age = TRUE;
					break;
				case 'N':	/* Extend oldest stage back to this time [no extension] */
					upper_age = atof (&argv[i][2]);
					break;
				case 'O':
					age_shift = fabs(atof (&argv[i][2]));
                                        fake_origin = TRUE;
                                        break;
				case 'P':
                                        ages_given = TRUE;
                                        break;
				case 'S':
					dl = atof (&argv[i][2]);
                                        int_seg = TRUE;
					break;
				case 'T':	/* Current age [0 Ma] */
					t_zero = atof (&argv[i][2]);
					break;
                                default:
                                        error = TRUE;
                                        break;
                        }
                }
                else {
			if (!data_in_input && argv[i][0] != ' ') {
                        	if ((fp = fopen(argv[i], "r")) == NULL) {
                                	mexPrintf("%s:  Cannot open %s\n", "telha", argv[i]);
                                	error = TRUE;
                        	}
                        }
                }
        }

	if (argc == 1) {	/* Display usage */
 		mexPrintf ("Telha - 3D magnetic modeling program\n");
		mexPrintf("usage:  telha infile -E<euler> -N<upper_age> [-A] [-T<t_zero>] [-F<filterwidth>] [-G] [-I] [-B] [-D<inc_deg>] [-S<inc_deg>] [-C[P|C][<time_scale>]] [-L<start/stop/inc>] [-P] [-:]\n\n");
		mexPrintf ("	infile (ASCII) has 2 columns defining segment limits.\n");
 		
		mexPrintf ("\t-E specifies the stage pole file to be used\n");
		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("\t-A create age flow lines (instead of mag) starting at points as \n\t  determined by -S and -G options\n");
		mexPrintf ("\t-B output a file with direct/reverse polarity blocks\n");
		mexPrintf ("\t-C choose between different isochron time scales. Append C (-CC) to\n");
		mexPrintf ("\t   select Cande & Kent 95 time scale (no need, is the Default), O (-CO)\n");
		mexPrintf ("\t   to select the old Cand time scale (used to be the default), P (-CP)\n");
		mexPrintf ("\t   to select Patriat time scale. Altervatively give a file name\n");
		mexPrintf ("\t   (-Ctime_scale) with a user provided time scale.\n");
		mexPrintf ("\t-D<inc_deg> distance in degrees used to reinterpolate the flow line \n\t  equidistantly [default 0.005]\n");
		mexPrintf ("\t-F<filterwidth> apply a gaussian filter (width in meters)\n");
		mexPrintf ("\t-G create multiple flow lines starting at points defining the \n\t  segment OR at interpolated positions controled by -S option\n");
		mexPrintf ("\t-I invert sense of rotation as defined in the stage pole file\n");
		mexPrintf ("\t-L compute rotated points along a linear 'age-line'.\n");
		mexPrintf ("\t-P It means that age points are transmited in input (They MUST be).\n");
		mexPrintf ("\t-N<age> extends earliest stage pole back to <upper_age> [no extension]\n");
		mexPrintf ("\t-S<inc_deg> reinterpolate the segment at inc_deg spacing\n");
		mexPrintf ("\t-T<age> sets the starting age in Ma [0]\n");
		mexPrintf("\t-: Expect lat/lon input rather than lon/lat [OFF].\n");
		mexErrMsgTxt("\n");
        }

	if (age_flow && filter_width > 0.001) {
		mexPrintf ("%s: Warning -A option overuns -F<filterwidth>\n", "telha");
		filter_width = 0.001;
	}
	if (multi_flow && brick) {
		mexPrintf ("%s: Warning -B option overuns -G\n", "telha");
	}
	if (ts_file) {
		if ((fp_ts = fopen (t_scale, "r")) == NULL) {
			mexPrintf ("%s: Cannot open file %s\n", argv[0], t_scale);
			mexErrMsgTxt("\n");
		}
	}
	if (n_arg_no_char == 0 && ages_given) {		/* This will fail if -P and data_in_input */
		mexPrintf ("%s: ERROR: -P option requires ages to be transmited in input\n", "telha");
		mexErrMsgTxt("\n");
	}

	if (error) mexErrMsgTxt("\n");

	if (switch_xy) {iy = 0; ix = 1;}

	if (n_arg_no_char == 1 && !ages_given) {		/* Not all combinations are tested */
		tmp = mxGetPr(prhs[0]);
		np_in = mxGetM(prhs[0]);
		data_in_input = TRUE;
	}
	else if (n_arg_no_char == 1 && ages_given) {		/* prhs[0] MUST contain the "rotation ages" */
		tmp = mxGetPr(prhs[0]);
		n_isoc = MAX(mxGetM(prhs[0]),mxGetN(prhs[0]));
		isoc_scale = (struct ISOC_SCALE *) mxCalloc ((size_t) n_isoc, sizeof(struct ISOC_SCALE));
		for (i = 0; i < n_isoc; i++)
			isoc_scale[i].age = tmp[i];
		linear_age = ts_file = isoc_p = isoc_o = isoc_c = FALSE;
	}
	else if (n_arg_no_char == 2) {		/* Here, prhs[0] MUST contain the object to rotate */
						/* and prhs[1] the "rotation ages" */ 
		/* The "rotation ages" */
		tmp = mxGetPr(prhs[1]);
		n_isoc = MAX(mxGetM(prhs[1]),mxGetN(prhs[1]));
		isoc_scale = (struct ISOC_SCALE *) mxCalloc ((size_t) n_isoc, sizeof(struct ISOC_SCALE));
		for (i = 0; i < n_isoc; i++)
			isoc_scale[i].age = tmp[i];
		linear_age = ts_file = isoc_p = isoc_o = isoc_c = FALSE;
		/* And now the the object to rotate */
		tmp = mxGetPr(prhs[0]);
		np_in = mxGetM(prhs[0]);
		data_in_input = TRUE;
	}

	if (linear_age) {
		n_isoc = (int)floor(fabs(age_lin_stop - age_lin_start) / age_lin_inc + 0.5) + 1;
		isoc_scale = (struct ISOC_SCALE *) mxCalloc ((size_t) n_isoc, sizeof(struct ISOC_SCALE));
		isoc_scale[0].age = age_lin_start;
		for (i = 1; i < n_isoc; i++)
			isoc_scale[i].age = isoc_scale[i-1].age + age_lin_inc;
	}
	else if (ts_file) /* read time scale provided by user */
		n_isoc = read_time_scale (fp_ts, isoc_scale);
	else if (isoc_p) {	/* use Patriat time scale */ 
		n_isoc = sizeof isoc_scale_p / sizeof (struct ISOC_SCALE);
		isoc_scale = (struct ISOC_SCALE *) mxCalloc ((size_t)n_isoc, sizeof(struct ISOC_SCALE));
		isoc_scale = isoc_scale_p;
	}
	else if (isoc_c) {     /* use Cande & Kent 95 time scale (DEFAULT) */
		n_isoc = sizeof isoc_scale_c / sizeof (struct ISOC_SCALE);
		isoc_scale = (struct ISOC_SCALE *) mxCalloc ((size_t)n_isoc, sizeof(struct ISOC_SCALE));
		isoc_scale = isoc_scale_c;
	}
	else if (isoc_o) {     /* use old Cande time scale (by Patriat) */
		n_isoc = sizeof(isoc_scale_o) / sizeof (struct ISOC_SCALE);
		isoc_scale = (struct ISOC_SCALE *) mxCalloc ((size_t)n_isoc, sizeof(struct ISOC_SCALE));
		isoc_scale = isoc_scale_o;
	}
	else if (ages_given);
	else {
		mexPrintf ("%s: ERROR, should not pass here %s\n", argv[0]);
		mexErrMsgTxt("\n");
	}


        n_alloc = CHUNK;
	ndata = 0;	n_seg = 0;	np_in_seg = 0;

	if (data_in_input) {	/* Data was transmited as arguments. */
		data = (struct DATA *) mxCalloc ((size_t) np_in, sizeof(struct DATA));
		for (i = 0; i < np_in; i++) {
			if (!mxIsNaN (tmp[i])) {
				data[ndata].x = tmp[i];
				/* Convert to geocentrics */
				data[ndata].y = atan2( ecc2_1 * sin(tmp[i+np_in]*D2R), cos(tmp[i+np_in]*D2R) ) / D2R;
				ndata++;
				np_in_seg++;
			}
			else if (first) {
				first = FALSE;
				n_seg++;
			}
			else {
				data[n_seg-1].np = np_in_seg; 
				n_seg++;
				np_in_seg = 0;
			}
		}
		if (first) n_seg++; 		/* Single segment file without any NaNs */
		data[n_seg-1].np = np_in_seg;	/* Save npoints of the last segment */
		np_in = ndata + n_seg - 1;
	}
	else {			/* Data was not transmited as arguments. So we must read it from file */
		data = (struct DATA *) mxCalloc ((size_t) n_alloc, sizeof(struct DATA));
		while (fgets (line, MAXCHAR, fp)) { /* Reading file */
			if (ns = strncmp(newseg, line, 1)) { /* Multisegment file */
				if (sscanf (line, "%lg %lg", &in[0], &in[1]) < 2)
					mexPrintf ("ERROR deciphering line %d\n",ndata+1);
                		if (ndata == n_alloc) {
                        		n_alloc += CHUNK;
					data = (struct DATA *) mxRealloc ((void *)data, (size_t)(n_alloc * sizeof(struct DATA)));
                		}
				data[ndata].x = in[ix];
				data[ndata].y = in[iy];
				ndata++;
				np_in_seg++;
			}
			else if (first) {
				first = FALSE;
				n_seg++;
			}
			else {
				data[n_seg-1].np = np_in_seg;
				n_seg++;
				np_in_seg = 0;
			}
		}
		if (first) n_seg++; 		/* Single segment file without any ">" */
		data[n_seg-1].np = np_in_seg;	/* Save npoints of the last segment */
		fclose(fp);
		np_in = ndata + n_seg - 1;
	}
	
	if (int_seg) intp_along_seg(n_seg, dl);

	/* I'm not sure that the user realy wants this all the times */
	if (upper_age == 0.0 && linear_age) upper_age = age_lin_stop;

	/* Load in the stage poles */

	n_stages = spotter_init (euler_file, &p, 0, upper_age);
	if (upper_age == 0.0) upper_age = p[0].t_start;

	if (invert_rot) /* Want to revert sense of rotation */
		for (i = 0; i < n_stages; i++) {
			p[i].omega *= -1;	p[i].omega_r *= -1;
		}

	if (p[n_stages - 1].t_stop > 0) {
		mexPrintf ("ERROR: A idade mais nova do polo mais recente \ntem que ser 0 nao %lg como e o caso\n", p[n_stages - 1].t_stop);
		return;
	}
	if (t_zero > 0)
		while (isoc_scale[i_min].age < t_zero) i_min++;


	if (linear_age || ages_given) {
		n_flow = 0;	l = 0;
		/* Count number of rotations for each point in segment */
		for (jj = 0; jj < n_isoc && isoc_scale[jj].age <= upper_age; jj++) n_flow++;
		n_flow -= i_min;
		np_out = np_in * n_flow + n_flow;
		plhs[0] = mxCreateDoubleMatrix (np_out, 2, mxREAL); 
		out = mxGetPr(plhs[0]); 
		for (j = i_min; j < n_flow + i_min; j++) {	/* Loop over number of rotations*/
			age = isoc_scale[j].age;
			m = 0;
			for (k = 0; k < n_seg; k++) {	/* Loop over segments */
				ii = 0;
				c = mxCalloc (data[k].np * 2, sizeof (double));
				for (i = 0; i < data[k].np; i++) { /* Loop over np in segment */
					lon =  data[m].x * D2R;	
					lat = data[m].y * D2R;
					n_chunk = spotter_backtrack (&lon, &lat, &age, 1, p, n_stages, 0);
					lon *= R2D;	lat *= R2D;
					if (fabs(lon / 180.) > 1) lon -= 360.;
					c[ii++] = lon;	c[ii++] = lat;
					m++;
				}
				if (n_flow > 1 || n_seg > 1) {
					out[l] = mxGetNaN();	out[l+np_out] = mxGetNaN();
					l++;
				}
				for (ll = 0; ll < data[k].np; ll++) {
					out[l] = c[ll*2];	out[l+np_out] = c[ll*2+1];
					l++;
				}
				mxFree(c);
			}
		}
		if (!data_in_input) mxFree ((void *)data);
	}
	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix (1, 1, mxREAL); 
		out_n_data = mxGetPr(plhs[1]);	out_n_data[0] = ndata; 
	}
	if (nlhs > 2) {
		plhs[2] = mxCreateDoubleMatrix (1, 1, mxREAL); 
		out_n_seg = mxGetPr(plhs[2]);	out_n_seg[0] = n_seg; 
	}
	if (nlhs > 3) {
		plhs[3] = mxCreateDoubleMatrix (1, 1, mxREAL); 
		out_n_flow = mxGetPr(plhs[3]);	out_n_flow[0] = n_flow; 
	}

	if (!brick)		/* Not very clean code, but we are done here */
		return;


	n_flow = 0;
	/* Count number of rotations for each point in segment */
	for (jj = 0; jj < n_isoc && isoc_scale[jj].age <= upper_age; jj++) n_flow++;
	n_flow -= i_min;

	/* NOTE. Since from the above loop jj is allways < n_isoc, but it also is incremented one
	   too much before the condition fails, we can safely do the following. The intention is that
	   we can have also the portion of the last brick between isoc_scale[jj-1].age and upper_age */
	if (isoc_scale[jj].age != upper_age) {		/* Otherwise we would be repeating last point */
		isoc_scale[jj].age = upper_age;
		n_flow++;
	}

	/* For the same reason, that is to account for the portion of the brick from t_start
	   to the next closest entry in the chronological table, we do this trick. */
	if (t_zero > 0 && i_min > 0) {
		i_min--;
		isoc_scale[i_min].age = t_zero;
		n_flow++;
	}

	for (k = m = 0; k < n_seg; k++) {	/* Loop over segments */
		c = mxCalloc ((size_t)(n_flow * data[k].np * 2), sizeof(double));
		ii = 0;		m = 0;
		for (i = 0; i < data[k].np; i++) { /* Loop over np in segment */
			for (j = i_min; j < n_flow + i_min; j++) { /* Loop over number of rotations*/ 
				lon = data[i].x * D2R;	
				lat = data[i].y * D2R;
				age = isoc_scale[j].age;
				n_chunk = spotter_backtrack (&lon, &lat, &age, 1, p, n_stages, 0);
				lon *= R2D;
				lat = atan2( sin(lat), ecc2_1 * cos(lat) ) * R2D;	/* Convert back to geodetics */
				if (fabs(lon / 180.) > 1) lon -= 360.;
				c[ii++] = lon;	c[ii++] = lat;
			}
			m++;
		}
		if (brick) {
			i = 0;
			if (first_seg) {
				plhs[0] = mxCreateDoubleMatrix (np_in * 2, n_flow, mxREAL); 
				plhs[1] = mxCreateDoubleMatrix (np_in * 2, n_flow, mxREAL); 
				out = mxGetPr(plhs[0]); 
				out2 = mxGetPr(plhs[1]);
				if (nlhs > 2) {
					plhs[2] = mxCreateDoubleMatrix (1, 1, mxREAL); 
					first_mag = mxGetPr(plhs[2]);
					first_mag[0] = isoc_scale[i_min].mag;
				}
			}

			for (l = 0; l < n_flow-1; l++) {
				for (ll = 0; ll < data[k].np; ll++) {
					out[i] = c[(l+ll*n_flow)*2];
					out2[i] = c[(l+ll*n_flow)*2+1];
					i++;
				}
				for (ll = data[k].np; ll > 0; ll--) {
					out[i] = c[((ll-1)*n_flow+l+1)*2];
					out2[i] = c[((ll-1)*n_flow+l+1)*2+1];
					i++;
				}
			}
			first_seg = FALSE;
		}
		else
			int_perf (c, i_min, n_flow, k, age_flow, multi_flow, ddeg, filter_width, isoc_scale);
		mxFree ((void *)c);
	}
	mxFree ((void *)data);
	mxFree ((void *)p);
	mxFree ((void *)argv);
}

int read_time_scale (FILE *fp_ts, struct ISOC_SCALE *isoc_scale) {
	/* Read file with a isochron time scale */
	int n_alloc, ndata = 0;
	double in[2];
	char line[512];

       	n_alloc = SMALL_CHUNK;

	if ((isoc_scale = (struct ISOC_SCALE *) mxCalloc ((size_t) n_alloc, sizeof(struct ISOC_SCALE)) ) == NULL) 
			{no_sys_mem("telha --> (read_isoc_scale)", n_alloc);}

	while (fgets (line, 512, fp_ts)) { 
		if(sscanf (line, "%lg %lg", &in[0], &in[1]) !=2)
			mexPrintf ("ERROR deciphering line %d of isoc_scale file\n", ndata+1);
               	if (ndata == n_alloc) {
                       	n_alloc += SMALL_CHUNK;
			if ((isoc_scale = (struct ISOC_SCALE *) mxRealloc ((void *)isoc_scale, (size_t)(n_alloc * sizeof(struct ISOC_SCALE)))) == NULL) 
				{no_sys_mem("telha --> (read_isoc_scale)", n_alloc);}
               	}
		isoc_scale[ndata].age = in[0];
		isoc_scale[ndata].mag = in[1];
		ndata++;
	}
	fclose(fp_ts);
	return (ndata);
}

int int_perf (double *c, int i_min, int n_flow, int k, int age_flow, int multi_flow,
	      double ddeg, double filter_width, struct ISOC_SCALE *isoc_scale) {
	int l, ll, i, j, n_pts, ii, go_x, sign;
	double x1, y1, x2, y2, *xx, *yy, *m_a, *u, *v, *vm, dx, dy, dr;
	int	n_f_wts;		/* Number of filter weights  */
	int	half_n_f_wts;		/* Half the number of filter weights  */
	double	*mag_intp = NULL;	/* Pointer for array of magnetizations */
	double  t_start;		/* x-value of first output point */
	double  t_stop;			/* x-value of last output point */
	double	half_width, dt;	
	double	*f_wt = NULL;		/* Pointer for array of filter coefficients  */
	double	*f_in_m = NULL;		/* Pointer for array of flow line in meters */
	/* Vai ser preciso mais tarde decidir quem sao os xx e yy transmitidos para intp_lin	*/

	for (ll = 0; ll < data[k].np; ll++) {
		if ((xx = (double *) mxCalloc ((size_t)(2*n_flow), sizeof(double)) ) == NULL) 
			{no_sys_mem("telha --> int_perf (xx)", 2*n_flow);} 
		if ((yy = (double *) mxCalloc ((size_t)(2*n_flow), sizeof(double)) ) == NULL) 
			{no_sys_mem("telha --> int_perf (yy)", 2*n_flow);} 
		if ((m_a = (double *) mxCalloc ((size_t)(2*n_flow), sizeof(double)) ) == NULL) 
			{no_sys_mem("telha --> int_perf (m_a)", 2*n_flow);} 
		i = i_min;
		for (l = 0; l < n_flow; l++, i++) {
			x1 = c[(l+ll*n_flow)*2]; y1 = c[(l+ll*n_flow)*2+1];
			x2 = c[(l+1+ll*n_flow)*2]; y2 = c[(l+1+ll*n_flow)*2+1];
			xx[2*l] = x1;	xx[2*l+1] = x2;	
			yy[2*l] = y1;	yy[2*l+1] = y2;	
			if (!age_flow) {
				m_a[2*l] = isoc_scale[i].mag;	
				m_a[2*l+1] = isoc_scale[i].mag;	
			}
			else {
				m_a[2*l] = isoc_scale[i].age;	
				m_a[2*l+1] = isoc_scale[i+1].age; /*Pode dar merda */
			}
		}
		xx[2*n_flow-1] = xx[2*n_flow-2];	
		yy[2*n_flow-1] = yy[2*n_flow-2];	/* Last point is not defined */
		dx = xx[2*n_flow-1] - xx[0];	dy = yy[2*n_flow-1] - yy[0];
		dr = (fabs(dx) > fabs(dy)) ? dx : dy;
		go_x = (fabs(dx) > fabs(dy)) ? 1 : 0;
		n_pts = (int)fabs(dr / ddeg) + 1;
		if ((u = (double *) mxCalloc ((size_t)(n_pts+1), sizeof(double)) ) == NULL) 
			{no_sys_mem("telha --> int_perf (u)", n_pts);}
		if ((v= (double *) mxCalloc ((size_t)(n_pts+1), sizeof(double)) ) == NULL) 
			{no_sys_mem("telha --> int_perf (vy)", n_pts);} 
		if ((vm= (double *) mxCalloc ((size_t)(n_pts+1), sizeof(double)) ) == NULL) 
			{no_sys_mem("telha --> int_perf (vm)", n_pts);} 

		if (go_x) {
			sign = (dx > 0.) ? 1 : -1;
			u[0] = xx[0];
			for (j = 1; j < n_pts; j++) u[j] = u[j-1] + ddeg * sign;
			intp_lin (xx, yy, 2*n_flow, n_pts, u, v, 0);
			intp_lin (xx, m_a, 2*n_flow, n_pts, u, vm, 0);
		}
		else {
			sign = (dy > 0.) ? 1 : -1;
			u[0] = yy[0];
			for (j = 1; j < n_pts; j++) u[j] = u[j-1] + ddeg * sign;
			intp_lin (yy, xx, 2*n_flow, n_pts, u, v, 0);
			intp_lin (yy, m_a, 2*n_flow, n_pts, u, vm, 0);
		}

		mxFree ((void *) xx);
		mxFree ((void *) yy);
		mxFree ((void *) m_a);

		if (age_flow) {
			if (multi_flow) mexPrintf ("%>\n");
			for (j = 0; j < n_pts; j++)
				if (go_x)
					mexPrintf ("%lf %lf %.3f\n", u[j], v[j], vm[j]);
				else
					mexPrintf ("%lf %lf %.3f\n", v[j], u[j], vm[j]);
		}
		else {

			if ((f_in_m = (double *) mxCalloc ((size_t)(n_pts+1), sizeof(double))) == NULL) 
			{no_sys_mem("telha --> int_perf (f_in_m)", n_pts+1);} 
			if ((mag_intp = (double *) mxCalloc ((size_t)(n_pts+1), sizeof(double))) == NULL) 
			{no_sys_mem("telha --> int_perf (mag_intp)", n_pts+1);} 
			/* Now covert to distance in m from begining of flow line */
			flow_in_meters (u, v, n_pts, f_in_m);	/* output in pointer f_in_m */
			/*half_width = set_up_filter(f_in_m[0], f_in_m[n_pts-1], n_pts, filter_width);*/
			half_width = set_up_filter(f_in_m[0], f_in_m[n_pts-1], n_pts, filter_width,
						&dt, &t_start, &t_stop, &half_n_f_wts, &n_f_wts, f_wt);

			do_the_filter (vm, n_pts, 2, half_width, mag_intp, n_f_wts, half_n_f_wts, dt, t_start, t_stop, f_wt, f_in_m);

			if (multi_flow) mexPrintf ("%>\n");
			for (j = 0; j < n_pts; j++)
				if (go_x)
					mexPrintf ("%lf %lf %.3f\n", u[j], v[j], mag_intp[j]);
				else
					mexPrintf ("%lf %lf %.3f\n", v[j], u[j], mag_intp[j]);

			mxFree ((void *)f_wt);
			mxFree ((void *)f_in_m);
			mxFree ((void *)mag_intp);
		}
		mxFree ((void *) u);
		mxFree ((void *) v);
		mxFree ((void *) vm);
	}
	return (0);
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
	double dx, *c;
	
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

double	set_up_filter(double t_i, double t_f, int n_pts, double filter_width,
		      double *dt, double *t_start, double *t_stop, int *half_n_f_wts, 
		      int *n_f_wts, double *f_wt) {
	int	i, i1, i2;
	double	time, w_sum, half_width;
	
	*dt = (t_f - t_i) / (n_pts - 1);

	half_width = 0.5 * filter_width;
	*half_n_f_wts = (int)floor (half_width / (*dt));
	*n_f_wts = 2 * (*half_n_f_wts) + 1;
		
	f_wt = (double *) mxCalloc ((size_t)(*n_f_wts), sizeof(double));
	for (i = 0; i <= *half_n_f_wts; i++) {
		time = i * (*dt);
		i1 = *half_n_f_wts - i;
		i2 = *half_n_f_wts + i;
		f_wt[i1] = f_wt[i2] = gaussian_weight(time, half_width);
	}

	*t_start = t_i;		*t_stop = t_f;	/* Initialize start/stop time */

	return(half_width);
}


double	gaussian_weight (double radius, double half_width) {
	double weight, gauss_constant;
	gauss_constant = -4.5 / (half_width * half_width);
	if (radius > half_width)
		weight = 0.0;
	else
		weight = exp (radius * radius * gauss_constant);
		
	return (weight);
}

int spotter_init (char *file, struct EULER **p, int flowline, double t_max) {
	/* file;	Name of file with backward stage poles */
	/* p;		Pointer to stage pole array */
	/* flowline;	TRUE if flowlines rather than hotspot-tracks are needed */
	/* t_max;	Extend earliest stage pole back to this age */
	FILE *fp;
	struct EULER *e;
	char  buffer[MAXCHAR];
	int n, i = 0, n_alloc = SMALL_CHUNK;
	double x, y, last_t = DBL_MAX;

	e = (struct EULER *) mxCalloc ((size_t) n_alloc, sizeof(struct EULER));

	if ((fp = fopen (file, "r")) == NULL) {
		mexPrintf ("EULER: Cannot open stage pole file: %s\n", file);
		mexErrMsgTxt("\n");
	}

	while (fgets (buffer, 512, fp) != NULL) { /* Expects lon lat t0 t1 ccw-angle */
		if (buffer[0] == '#' || buffer[0] == '\n') continue;

		sscanf (buffer, "%lf %lf %lf %lf %lf", &e[i].lon, &e[i].lat, &e[i].t_start, &e[i].t_stop, &e[i].omega);

		if (e[i].t_stop >= e[i].t_start) {
			mexPrintf ("ERROR: Stage rotation %d has start time younger than stop time\n", i);
			mexErrMsgTxt("\n");
		}
		if (e[i].t_stop > last_t) {
			mexPrintf ("ERROR: Stage rotations must go from oldest to youngest\n");
			mexErrMsgTxt("\n");
		}
		last_t = e[i].t_stop;

		e[i].duration = e[i].t_start - e[i].t_stop;
		e[i].omega /= e[i].duration;
		e[i].omega_r = e[i].omega * D2R;
		e[i].sin_lat = sin (e[i].lat * D2R);
		e[i].cos_lat = cos (e[i].lat * D2R);
		x = e[i].lon * D2R;	y = e[i].lat * D2R;
		e[i].lon_r = x;		e[i].lat_r = y;
		i++;
		if (i == n_alloc) {
			n_alloc += SMALL_CHUNK;
			e = (struct EULER *) mxRealloc ((void *)e, (size_t)(n_alloc * sizeof(struct EULER)));

		}
	}
	n = i;
	fclose(fp);

	/* Extend oldest stage pole back to t_max Ma */

	if (t_max > 0.0 && e[0].t_start < t_max) {
		/*mexPrintf ("Extending oldest stage pole back to %lg Ma\n", t_max);*/
		e[0].t_start = t_max;
		e[0].duration = e[0].t_start - e[0].t_stop;
	}

	e = (struct EULER *) mxRealloc ((void *)e, (size_t)(n * sizeof(struct EULER)));
	*p = e;

	return (n);
}

/* spotter_backtrack: Given a seamount location and age, trace the
 *	hotspot-track between this seamount and a seamount of 
 *	age t_zero.  For t_zero = 0 this means the hotspot
 */

int spotter_backtrack (double xp[], double yp[], double tp[], int np, struct EULER p[], int ns, double t_zero) {
/* xp, yp;	Points, in RADIANS */
/* tp;		Age of feature in m.y. */
/* np;		# of points */
/* p;		Stage poles */
/* ns;		# of stage poles */
/* t_zero;	Backtrack up to this age */
 
	int j, k, kk = 0, start_k, nd, n_alloc = 2 * CHUNK;
	double t, tt, dt, d_lon, tlon, dd, i_km, xnew, xx, yy, *track;
	double s_lat, c_lat, s_lon, c_lon, cc, ss, cs;

		t = tp[0];
		while (t > 0) {	/* As long as we're not back at zero age */
			j = 0;
			while (j < ns && t <= p[j].t_stop){ j++;	/* Find first applicable stage pole */
			}
			dt = MIN (p[j].duration, t - MAX(p[j].t_stop, t_zero));
			d_lon = p[j].omega_r * dt;
			xnew = xp[0] - p[j].lon_r;
			s_lat = sin (yp[0]);
			c_lat = cos (yp[0]);
			s_lon = sin (xnew);
			c_lon = cos (xnew);
			cc = c_lat * c_lon;
			tlon = d_atan2 (c_lat * s_lon, p[j].sin_lat * cc - p[j].cos_lat * s_lat);
			s_lat = p[j].sin_lat * s_lat + p[j].cos_lat * cc;
			c_lat = sqrt (1.0 - s_lat * s_lat);
			ss = p[j].sin_lat * s_lat;
			cs = p[j].cos_lat * s_lat;
			xnew = tlon + d_lon;
			s_lon = sin (xnew);
			c_lon = cos (xnew);
			cc = c_lat * c_lon;
			yp[0] = d_asin (ss - p[j].cos_lat * cc);
			xp[0] = p[j].lon_r + d_atan2 (c_lat * s_lon, p[j].sin_lat * cc + cs);

			if (xp[0] < 0.0) xp[0] += TWO_PI;
			if (xp[0] >= TWO_PI) xp[0] -= TWO_PI;
			t -= dt;
		}
	return (np);
}

int	do_the_filter (double *fm, int n_rows, int n_cols, double half_width, double *mag_intp, 
		       int n_f_wts, int half_n_f_wts, double dt, double t_start, double t_stop, 
		       double *f_wt, double *f_in_m) {
/*	n_rows	Number of array points	*/
/*	n_cols	Should allways be = 2 in this program 	*/

	int	i_row, i_col, left, right, i_f_wt, j = 0;
	int	i_t_output, n_in_filter, n_good_ones;
	int	*n_this_col;		/* Pointer to array of counters [one per column]  */
	double	time, delta_time, *outval, wt, val;
	double	*wt_sum;		/* Pointer for array of weight sums [each column]  */
	double	*data_sum;		/* Pointer for array of data * weight sums [columns]  */
	int	*good_one;		/* Pointer to array of logicals [one per column]  */

	outval = (double *)mxCalloc ((size_t)n_cols, sizeof (double));
	n_this_col = (int *)mxCalloc ((size_t)n_cols, sizeof(int));
	wt_sum = (double *)mxCalloc ((size_t)n_cols, sizeof(double));
	data_sum = (double *)mxCalloc ((size_t)n_cols, sizeof(double));
	good_one = (int *)mxCalloc ((size_t)n_cols, sizeof(int));
	
	/* Position i_t_output at first output time  */
	for (i_t_output = 0; f_in_m[i_t_output] < t_start; i_t_output++);

	time = t_start;
	left = right = 0;		/* Left/right end of filter window */

	while (time <= t_stop) {
		while ((time - f_in_m[left]) > half_width) left++;
		while (right < n_rows && (f_in_m[right] - time) <= half_width) right++;
		n_in_filter = right - left;
		if (!(n_in_filter)) {
			i_t_output++;
			time = (i_t_output < n_rows) ? f_in_m[i_t_output] : t_stop + 1.0;
			continue;
		}
			
		for (i_col = 0; i_col < n_cols; i_col++) {
			n_this_col[i_col] = 0;
			wt_sum[i_col] = 0.0;
			data_sum[i_col] = 0.0;
			if (i_col == 0) {
				good_one[i_col] = FALSE;
			}
			else {
				good_one[i_col] = TRUE;
			}
		}
			
		for (i_row = left; i_row < right; i_row++) {
			delta_time = time - f_in_m[i_row];
			i_f_wt = half_n_f_wts + (int)floor(0.5 + delta_time / dt);
			if ( (i_f_wt < 0) || (i_f_wt >= n_f_wts) ) continue;
				
				i_col = 1;
				if ((good_one[i_col])) {;
					wt = f_wt[i_f_wt];
					val = fm[i_row];
					wt_sum[i_col] += wt;
					data_sum[i_col] += (wt * val);
					n_this_col[i_col]++;
				}

		}
		for (i_col = 0; i_col < n_cols; i_col++) {
			if (i_col == 0)
				outval[i_col] = time;
			else if (good_one[i_col])
				outval[i_col] = FALSE ? data_sum[i_col] : data_sum[i_col] / wt_sum[i_col];
			else
				/*outval[i_col] = GMT_d_NaN;*/
				mexPrintf ("Nunca devia ter passado aqui\n");
		}
		mag_intp[j++] = -outval[1];

		/* Go to next output time */

		i_t_output++;
		time = (i_t_output < n_rows) ? f_in_m[i_t_output] : t_stop + 1.0;
	}
	
	mxFree ((void *)outval);
	mxFree ((void *)n_this_col);
	mxFree ((void *)wt_sum);
	mxFree ((void *)data_sum);
	mxFree ((void *)good_one);
	
	return(0);
}

int	flow_in_meters (double *x, double *y, int n_pts, double *f_in_m) {

	int i;
	double dx, dy;

	for (i = 1; i < n_pts; i++) {
		dx = (x[i] - x[i-1]);
		dy = (y[i] - y[i-1]);
		dy *= M_PR_DEG;
		dx *= (M_PR_DEG * cos (0.5 * (x[i] + x[i-1]) * D2R));
		f_in_m[i] = f_in_m[i-1] + hypot(dx,dy);
	}
	return (0);
}

void	no_sys_mem (char *where, int n) {	
		mexPrintf ("Fatal Error: %s could not allocate memory, n = %d\n", where, n);
}

int	intp_along_seg(int n_seg, double dl) {
	int i, m, j, n_ptsx, n_ptsy, n_pts, go_x, ndata = 0, n_alloc = CHUNK;
	int sign, size_data;
	double *xx, *yy, *u, *v, dx, dy, dr;
	struct DATA *p;

	if ((p = (struct DATA *)mxCalloc ((size_t)n_alloc, sizeof(struct DATA))) == NULL) 
			{no_sys_mem("telha --> (intp_along_seg)", n_alloc);}

	for (i = m = 0; i < n_seg; i++) {	/* Loop over segments */
		xx = (double *) mxCalloc ((size_t)(data[i].np), sizeof(double));
		yy = (double *) mxCalloc ((size_t)(data[i].np), sizeof(double));
		for (j = 0; j < data[i].np; j++, m++) {
			xx[j] = data[m].x;	yy[j] = data[m].y;
		}
		dx = xx[j-1] - xx[0];	dy = yy[j-1] - yy[0];
		dr = (fabs(dx) > fabs(dy)) ? dx : dy;
		go_x = (fabs(dx) > fabs(dy)) ? 1 : 0;
		n_pts = (int)fabs(dr / dl) + 1;
		u = (double *) mxCalloc ((size_t)(n_pts+1), sizeof(double));
		v = (double *) mxCalloc ((size_t)(n_pts+1), sizeof(double));

		if (go_x) {
			sign = (dx > 0.) ? 1 : -1;
			u[0] = xx[0];
			for (j = 1; j < n_pts; j++) u[j] = u[j-1] + dl * sign;
			if (u[n_pts-1] != xx[data[i].np-1]) {
				u[n_pts] = xx[data[i].np-1];
				n_pts++;
			}
			intp_lin (xx, yy, data[i].np, n_pts, u, v, 0);
		}
		else {
			sign = (dy > 0.) ? 1 : -1;
			u[0] = yy[0];
			for (j = 1; j < n_pts; j++) u[j] = u[j-1] + dl * sign;
			if (u[n_pts-1] != yy[data[i].np-1]) {
				u[n_pts] = yy[data[i].np-1];
				n_pts++;
			}
			intp_lin (yy, xx, data[i].np, n_pts, u, v, 0);
		}
		for (j = 0; j < n_pts; j++) {
			if (go_x) {
				p[ndata].x = u[j];
				p[ndata].y = v[j];
				ndata++;
               			if (ndata == n_alloc) {
                       			n_alloc += CHUNK;
					if ((p = (struct DATA *) mxRealloc ((void *)p, (size_t)(n_alloc * sizeof(struct DATA)))) == NULL) {no_sys_mem("telha --> (intp_along_seg)", n_alloc);}
               			}
			}
			else {
				p[ndata].x = v[j];
				p[ndata].y = u[j];
				ndata++;
               			if (ndata == n_alloc) {
                       			n_alloc += CHUNK;
					if ((p = (struct DATA *) mxRealloc ((void *)p, (size_t)(n_alloc * sizeof(struct DATA)))) == NULL) {no_sys_mem("telha --> (intp_along_seg)", n_alloc);} 
               			}
			}
		}
		p[i].np = n_pts;
		mxFree ((void *) xx);
		mxFree ((void *) yy);
		mxFree ((void *) u);
		mxFree ((void *) v);
	}
	size_data = sizeof data / sizeof (struct DATA);
	if (size_data < n_alloc) { /* sizeof p could be > sizeof data */
		if ((data = (struct DATA *) mxRealloc ((void *)data, (size_t)(n_alloc * sizeof(struct DATA)))) == NULL) 
			{no_sys_mem("telha --> (intp_along_seg)", n_alloc);}
	}
	data = p;
	return (0);
}
