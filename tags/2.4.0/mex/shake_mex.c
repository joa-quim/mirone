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

/*--------------------------------------------------------------------
 * Mex function to compute Peak Ground Acceleration/Velocity and Intensity
 * Author:	J. Luis, after the FORTRAN version from J M Miranda
 * Date: 	06-Apr-2010
 *
 * Usage
 * out[,out2[,out3]] = shake_mex(vs30, head, line, [-F<mecatype>], ['-M<mag>'], ['-O<avi>']);\n");
 *
 */ 

#include "gmt.h"
#include "mex.h"

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	GMT_LONG	greenwich = 0;
	GMT_LONG	proj_type = 0;		/* Geographic */
	int	i, j, k, nx, ny, argc = 0, n_arg_no_char = 0, imeca = 0;
	int	no[3], do_PGA = FALSE, do_PGV = FALSE, do_INT = FALSE, error = FALSE;
	
	char	*line_file = "lixo_line.xy";
	char	**argv;
	float	*vs30, *pga, *pgv, *pint;
	double	west = 0.0, east = 0.0, south = 0.0, north = 0.0, dist, xnear, ynear, d_scale;
	double	rmh_pga, rmh_pgv, r_pga, r_pgv, rfd_pga, rfd_pgv, rfd_pga4nl, k1, k2, k3;
	double	flin_pga, flin_pgv, pga4nl, bnl_pga, bnl_pgv, fnl_pga, fnl_pgv, fs_pga, fs_pgv;
	double	rfm_pgv, rfm_pga, c_pga, d_pga, c_pgv, d_pgv, tmp, tmp1, tmp2, tmp3, rlon, rlat;
	double	rmw = 7, u = 1, ss = 0, rns = 0, rs = 0;
	double	*fault, *head;

	struct	GMT_TABLE *xyline;
	FILE	*fp;
	PFD	distance_func;
	PFI	near_a_line;
	
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
	argv[0] = "shake_mex";
	for (i = 1; i < argc; i++)
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);

	argc = GMT_begin (argc, argv);

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
		
				/* Supplemental parameters */

				case 'F':
					imeca = atoi(&argv[i][2]);
					if (imeca < 1 || imeca > 4)
						mexErrMsgTxt("SHAKE_MEX: -F<type> error. 'type' must be in [1 4].");
					break;
				case 'M':
					rmw = atof(&argv[i][2]);
					break;
				case 'O':
					for (j = 2; argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'a':		/* PGA is requested */
								do_PGA = TRUE;
								no[0] = j-2;
								break;
							case 'v':		/* PGV is requested */
								do_PGV = TRUE;
								no[1] = j-2;
								break;
							case 'i':		/* Intensity is requested */
								do_INT = TRUE;
								no[2] = j-2;
								break;
						}
					}
					break;
				default:
					mexPrintf("%s: GMT SYNTAX ERROR.\n", "mapproject");
					error++;
					break;
			}
		}
	}
	
	if ((argc == 1 && n_arg_no_char == 0) || error) {
		mexPrintf("shake_mex - Compute Peak Ground Acceleration/Velocity and Intensity.\n\n");
		mexPrintf("usage: out = shake_mex(vs30, head, line, [-F<mecatype>], ['-M<mag>'], ['-O<avi>']);\n");

		mexPrintf("\tout  One to three arrays (see -O) with results (e.g. [int,pga = ...].\n");
		mexPrintf("\tvs30 Array VS30 velocities.\n");
		mexPrintf("\thead Row vector with the (9) parameter defining the grid's header.\n");
		mexPrintf("\tline A Mx2 array with the fault's trace vertices.\n");
		
		mexPrintf("\n\tOPTIONS:\n");

		mexPrintf("\t-F Select focal mechanism type (e.g. -F1 or -F2 ...).\n");
		mexPrintf("\t   1  unknown [Default].\n");
		mexPrintf("\t   2  strike-slip.\n");
		mexPrintf("\t   3  normal.\n");
		mexPrintf("\t   4  thrust.\n");
		mexPrintf("\t-M Select magnitude [Default is 7].\n");
		mexPrintf("\t-O avi is a string made up of 1 or more of these characters:\n");
		mexPrintf("\t   a  Peak Ground Acceleration [Default].\n");
		mexPrintf("\t   v  Peak Ground Velocity.\n");
		mexPrintf("\t   i  Intensity (computed from PGV) [Default].\n");
		mexPrintf("\t	The data is written out in the order specified in <avi>\n");
		return;
	}

	/* Find out in which data type was given the input array */
	if (mxIsSingle(prhs[0]))
		vs30 = (float *)mxGetData(prhs[0]);
	else {
		mexPrintf("SHAKE_MEX ERROR: Invalid input data type.\n");
		mexErrMsgTxt("Only valid data type is single.");
	}

	nx = mxGetN (prhs[0]);
	ny = mxGetM (prhs[0]);
	if (!mxIsNumeric(prhs[0]) || ny < 10 || nx < 10)
		mexErrMsgTxt("SHAKE_MEX ERROR: First arg must contain a decent array");

	if (!mxIsDouble(prhs[1]) || (mxGetN(prhs[1]) != 9))
		mexErrMsgTxt("SHAKE_MEX ERROR: Second arg must contain a header row vector with 9 doubles");

	head  = mxGetPr(prhs[1]);		/* Get header info */


	/* Check that third argument contains the line location */
	if ( !mxIsNumeric(prhs[2]) || (mxGetM(prhs[2]) < 2) || (mxGetN(prhs[2]) != 2) )
		mexErrMsgTxt("SHAKE_MEX ERROR: third arg must contain a Mx2 double array with the fault location.");

	fault = mxGetPr(prhs[2]);

	if ((do_PGA + do_PGV + do_INT) == 0) {do_INT = TRUE;	no[2] = 0;}	/* Default */

	if ( nlhs == 0 )
		mexErrMsgTxt("SHAKE_MEX ERROR: Must provide at least one output.");
	if ( nlhs > (do_PGA + do_PGV + do_INT) )
		mexErrMsgTxt("SHAKE_MEX ERROR: number of requested ouputs exceeds -O option selection.");

	/* -------------------------------------------------------------------------------
	 Since it's actually so complicated to fill in a GMT_TABLE struct, easier to write
	 the line vertices in a tmp file, and let GMT_import_table() do the job.
	 ------------------------------------------------------------------------------ */
	if ((fp = fopen (line_file, "w")) == NULL)
		mexErrMsgTxt ("SHAKE_MEX: Unable to open temporary file for writing\n");

	for (k = 0; k < mxGetM(prhs[2]); k++) {
		fprintf (fp, "%f\t%f\n", fault[k], fault[k+mxGetM(prhs[2])]);
	}
	fclose ((void *) fp);

	GMT_import_table ((void *)line_file, GMT_IS_FILE, &xyline, 0.0, greenwich, FALSE, FALSE);
	remove (line_file);
	/* ------------------------------------------------------------------------------------ */

	GMT_get_dist_scale ('k', &d_scale, &proj_type, &distance_func);

	if (project_info.projection == GMT_NO_PROJ) {	/* Supply dummy linear proj */
		GMT_parse_J_option ("x1d");
		if (!project_info.region_supplied) {
			west = 0.0;	east = 360.0;
			south = -90.0;	north = 90.0;
		}
	}

	if (west == east && south == north)
		mexErrMsgTxt ("GMT Fatal Error: No region selected - Aborts!\n");

	GMT_map_setup (west, east, south, north);
	
	if (proj_type == 0) {	/* Geographic data */
		GMT_distance_func = (PFD) GMT_great_circle_dist;
		near_a_line = (PFI) GMT_near_a_line_spherical;
		greenwich = (west < 0.0 && east > 0.0);
	}
	else {
		GMT_distance_func = (PFD) GMT_cartesian_dist;
		near_a_line = (PFI) GMT_near_a_line_cartesian;
	}

	/*  -------------------------------------------------------------
	! Allocate arrays for output ... respecting -O (or default) order
	! ------------------------------------------------------------- */
	if (do_PGA) {
		plhs[no[0]] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		pga = (float *)mxGetData(plhs[no[0]]);
	}
	if (do_PGV) {
		plhs[no[1]] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		pgv = (float *)mxGetData(plhs[no[1]]);
	}
	if (do_INT) {
		plhs[no[2]] = mxCreateNumericMatrix (ny,nx,mxSINGLE_CLASS,mxREAL);
		pint = (float *)mxGetData(plhs[no[2]]);
		if (!do_PGV) {		/* PGV not selected but we need it anyway for INT */
			pgv = mxCalloc(nx * ny, sizeof (float));
		}
	}

	/* ------------------------------------------------------------- */

	if (imeca == 2) {
		u = rns = rs = 0;	ss = 1;
	}
	else if (imeca == 3) {
		u = ss = rs = 0;	rns = 1;
	}
	else if (imeca == 4) {
		u = ss = rns = 0;	rs = 1;
	}

	/* - ------------------------------------------------------------
	! calcula o termo de magnitude scaling para o pga e pgv
	! nao depende da distancia mas apenas da magnitude
	! ------------------------------------------------------------- */
	rmh_pga = 6.75;
	rmh_pgv = 8.50;
	rfm_pga = -0.53804*u - 0.50350*ss - 0.75472*rns - 0.50970*rs;
	if (rmw <= rmh_pga)
		rfm_pga += (0.28805*(rmw-rmh_pga) - 0.10164*(rmw-rmh_pga)*(rmw-rmh_pga));

	rfm_pgv = 5.00121*u + 5.04727*ss + 4.63188*rns + 5.08210*rs;
	if(rmw <= rmh_pgv)
		rfm_pgv += (0.18322*(rmw-rmh_pgv) - 0.12736*(rmw-rmh_pgv)*(rmw-rmh_pgv));

	/* - ------------------------------------------------------------
	! itera na matriz de posiçoes pga (i,j), pgv (i,j) e pint (i,j)
	! para calcular o termo de distance scaling
	! ------------------------------------------------------------- */
	/* The k2 & k3 coeff bellow is to simplify the computations of origianl equations such:
	   c_pga = +bnl_pga * (3.*0.40546510 - 1.0986123)/(1.0986123 * 1.0986123); */
	k2 = (3.*0.40546510 - 1.0986123) / (1.0986123 * 1.0986123);
	k3 = (2.*0.40546510 - 1.0986123) / (1.0986123 * 1.0986123 * 1.0986123);
	k1 = log(0.6);
	for (i = k = 0; i < nx; i++) {
		rlon = head[0] + i * head[7];
		for (j = 0; j < ny; j++, k++) {
			rlat = head[2] + j * head[8];

			(void) near_a_line (rlon, rlat, xyline, 2, &dist, &xnear, &ynear);
			dist *= d_scale;

			r_pga = sqrt (dist*dist + 1.35*1.35);
			r_pgv = sqrt (dist*dist + 2.54*2.54);
			rfd_pga = (-0.66050+0.11970*(rmw-4.5))*log(r_pga/1.0) - 0.01151*(r_pga-1.0);
			rfd_pgv = (-0.87370+0.10060*(rmw-4.5))*log(r_pgv/1.0) - 0.00334*(r_pgv-1.0);
			rfd_pga4nl = (-0.66050+0.11970*(rmw-4.5))*log(r_pga/5.0) - 0.01151*(r_pga-5.0);

			/* -------------------------------------------------------------
			! calcula o efeito de sitio
			! ------------------------------------------------------------- */
			tmp = log(vs30[k] / 760);
			flin_pga = -0.360 * tmp;
			flin_pgv = -0.600 * tmp;

			pga4nl = exp(rfm_pga + rfd_pga4nl);

			if (vs30[k] <= 180.) {
				bnl_pga = -0.640;
				bnl_pgv = -0.600;
			}
			else if (vs30[k] <= 300 && vs30[k] > 180.) {
				tmp = log(vs30[k]/300.) / log(180./300.);
				bnl_pga = (-0.64+0.14) * tmp - 0.14;
				bnl_pgv = (-0.60+0.06) * tmp - 0.14;
			}
			else if (vs30[k] <= 760 && vs30[k] > 300.) {
				tmp = log(vs30[k]/760.) / log(300./760.);
				bnl_pga = -0.14 * tmp;
				bnl_pgv = -0.06 * tmp;
			}
			else {
				bnl_pga = bnl_pgv = 0.;
			}

			/* ---- verificar se as condicoes são sempre em pga !!!!!!!!!!!!!!!!!!!!! */
			c_pga = +bnl_pga * k2;		d_pga = -bnl_pga * k3;
			c_pgv = +bnl_pgv * k2;		d_pgv = -bnl_pgv * k3;

			if (pga4nl <= 0.03) {
				fnl_pga = bnl_pga * k1;
				fnl_pgv = bnl_pgv * k1;
			}
			else if (pga4nl > 0.03 && pga4nl <= 0.09) {
				tmp3 = log(pga4nl/0.03);	tmp1 = tmp3 * tmp3;	tmp2 = tmp1 * tmp3;
				fnl_pga = bnl_pga * k1 + c_pga * tmp1 + d_pga * tmp2;
				fnl_pgv = bnl_pgv * k1 + c_pgv * tmp1 + d_pgv * tmp2;
			}
			else if (pga4nl > 0.09) {
				tmp = log (pga4nl/0.1);
				fnl_pga = bnl_pga * tmp;
				fnl_pgv = bnl_pgv * tmp;
			}

			fs_pga = flin_pga + fnl_pga;
			fs_pgv = flin_pgv + fnl_pgv;

			if (do_PGA)
				pga[k] = (float)(exp (rfm_pga + rfd_pga + fs_pga) * 980.);

			if (do_PGV || do_INT)
				pgv[k] = (float)exp (rfm_pgv + rfd_pgv + fs_pgv);
		}
	}

	if (do_INT) {
		for (k = 0; k < nx * ny; k++) {
			tmp = log10(pgv[k]);
            		pint[k] = (float)(4.398 + 1.916 * tmp + 0.280 * tmp*tmp);
		}
	}

	GMT_end (argc, argv);

	mxFree(argv);
	if (do_INT && !do_PGV) mxFree(pgv);
}
