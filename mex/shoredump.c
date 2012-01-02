/*--------------------------------------------------------------------
 *	$Id:$
 *
 *	Coffeeright (c) 2004-2012 by J. Luis
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
 *
 *    shoredump extracts segments from the shoreline, rivers, and
 *    border databases.
 *
 *--------------------------------------------------------------------*/
/*
 * Original Author:	Paul Wessel
 * Date:	21-JUN-1995
 * Version:	3.0
 *
 * Mexified by	Joaquim Luis
 * Date:	16-APR-2004
 * Modified:	15-JUIN-2005		Now uses the -R by GMT_get_commom_args 
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		17/06/09 J Luis, Updated to compile with version 4.5.0
 */

#include "gmt.h"
#include "mex.h"
#include "mwsize.h"

#define LAKE	0
#define RIVER	1

mwSize prep_polygons(struct GMT_GSHHS_POL **p_old, mwSize np, int greenwich, int sample, double step, mwSize anti_bin);
int getpathname (char *name);

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

/* Matlab Gateway routine */

void mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[]) {

	GMT_LONG	i, np, ind, bin, base = 3, max_level = GMT_MAX_GSHHS_LEVEL, direction = 1, min_level = 0;
	GMT_LONG	blevels[GMT_N_BLEVELS], n_blevels = 0, rlevels[GMT_N_RLEVELS], n_rlevels = 0;
	mwSize is = 0, ir = 0, ib = 0, p_alloc = 0, bytes_to_copy, k, j;
	mwSize argc, dims[] = {0,0}, n;
	
	int	error = FALSE, get_river = FALSE, shift = FALSE, first_shore = TRUE, first_river = TRUE;
	int	greenwich = FALSE, get_shore = FALSE, get_border = FALSE, first_border = TRUE, test = FALSE;
	
	double	west = 0.0, east = 0.0, south = 0.0, north = 0.0, edge = 720.0, left, right, bsize;
	double min_area = 0.0, *p_shore, *p_river, *p_border;
	double west_border, east_border, nan, *pdata, step = 0.01;
	
	/*struct POL *p;*/
	struct GMT_GSHHS_POL *p;
	struct GMT_SHORE c;
	struct GMT_BR b, r;
	char **argv, res = 'l';
	char *shore_resolution[5] = {"full", "high", "intermediate", "low", "crude"};
#ifdef GMT_MINOR_VERSION
	struct GMT_SHORE_SELECT Ainfo;
#endif

	if (nrhs < 1 || nrhs > 6) {
		mexPrintf ("shoredump - Extract shorelines, rivers, or borders\n\n");
		mexPrintf ("usage: shoredump -R<west>/<east>/<south>/<north> [-A<min_area>[/<min_level>/<max_level>]]\n");
		mexPrintf ("	 [-D<resolution>] [-I<feature>] [-N<feature>] [-S] \n");

		mexPrintf ("\n\tOPTIONS:\n");
		mexPrintf ("	-A features smaller than <min_area> (in km^2) or of levels (0-4) outside min-max levels\n");
		mexPrintf ("	will be skipped [0/4 (4 means lake inside island inside lake)]\n");
		mexPrintf ("	-D Choose one of the following resolutions:\n");
		mexPrintf ("	   f - full resolution (may be very slow for large regions)\n");
		mexPrintf ("	   h - high resolution (may be slow for large regions)\n");
		mexPrintf ("	   i - intermediate resolution\n");
		mexPrintf ("	   l - low resolution [Default]\n");
		mexPrintf ("	   c - crude resolution\n");
		mexPrintf ("	-I extract rIvers.  Choose features below, repeat -I as many times as needed\n");
		mexPrintf ("	       1 = Permanent major rivers\n");
		mexPrintf ("	       2 = Additional major rivers\n");
		mexPrintf ("	       3 = Additional rivers\n");
		mexPrintf ("	       4 = Minor rivers\n");
		mexPrintf ("	       5 = Intermittent rivers - major\n");
		mexPrintf ("	       6 = Intermittent rivers - additional\n");
		mexPrintf ("	       7 = Intermittent rivers - minor\n");
		mexPrintf ("	       8 = Major canals\n");
		mexPrintf ("	       9 = Minor canals\n");
		mexPrintf ("	      10 = Irrigtion canals\n");
		mexPrintf ("	       a = All rivers and canals (1-10)\n");
		mexPrintf ("	       r = All permanent rivers (1-4)\n");
		mexPrintf ("	       i = All intermittent rivers (6-8)\n");
		mexPrintf ("	       c = All canals (9-10)\n");
		mexPrintf ("	-N extract boundaries.  Choose features below, repeat -N as many times as needed\n");
		mexPrintf ("	       1 = National boundaries\n");
		mexPrintf ("	       2 = State boundaries within the Americas\n");
		mexPrintf ("	       3 = Marine boundaries\n");
		mexPrintf ("	       a = All boundaries (1-3)\n");
		mexPrintf ("	-S extract shorelines [Default].\n");
		return;
	}

	if (nlhs < 1)
		mexErrMsgTxt ("ERROR: Need to specify at least one output;\n");

	if (nlhs > 3)
		mexErrMsgTxt ("ERROR: Pssible outputs are: shore, borders & rivers. What else do you want more?\n");

	nan = mxGetNaN();
	memset((char *)rlevels, 0, GMT_N_RLEVELS * sizeof(int));
	memset((char *)blevels, 0, GMT_N_BLEVELS * sizeof(int));

	argc = nrhs;
	for (i = 0; i < argc; i++) {		/* Check input to be sure it is of type char. */
		if(!mxIsChar(prhs[i]))
			mexErrMsgTxt("Input must be of type char.");
	}
	argc++;			/* to account for the program's name to be inserted in argv[0] */

	/* get the length of the input string */
	argv=(char **)mxCalloc(argc, sizeof(char *));
	argv[0] = "shoredump";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i-1]);
	}

	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);*/
	argc = GMT_begin (argc, argv);

	/* Check and interpret the command line arguments */
	
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
		
				/* Common parameters */
				
				case 'R':
					error += GMT_parse_common_options (argv[i], &west, &east, &south, &north);
					break;
				
				/* Supplemental parameters */

				case 'A':
#ifdef GMT_MINOR_VERSION
					Ainfo.fraction = Ainfo.flag = Ainfo.low = 0;
					Ainfo.high = GMT_MAX_GSHHS_LEVEL;
					Ainfo.area = 0.;
					GMT_set_levels (&argv[i][2], &Ainfo);
#else
					j = sscanf (&argv[i][2], "%lf/%d/%d", &min_area, &min_level, &max_level);
					if (j == 1) min_level = 0, max_level = GMT_MAX_GSHHS_LEVEL;
#endif
					break;
				case 'D':
					res = argv[i][2];
					switch (res) {
						case 'f':
							base = 0;
							break;
						case 'h':
							base = 1;
							break;
						case 'i':
							base = 2;
							break;
						case 'l':
							base = 3;
							break;
						case 'c':
							base = 4;
							break;
						default:
							mexPrintf("SHOREDUMP SYNTAX ERROR -D option:  Unknown modifier %c\n", argv[i][2]);
							error++;
							break;
					}
					break;
				case 'N':
					if (!argv[i][2]) {
						mexPrintf ("SHOREDUMP SYNTAX ERROR:  -N option takes one argument\n");
						error++;
						continue;
					}
					get_border = TRUE;
					switch (argv[i][2]) {
						case 'a':
							for (k = 0; k < GMT_N_BLEVELS; k++) blevels[k] = TRUE;
							break;
						default:
							k = atoi (&argv[i][2]) - 1;
							if (k < 0 || k >= GMT_N_BLEVELS) {
								mexPrintf ("SHOREDUMP SYNTAX ERROR -N option: Feature not in list!\n");
								error++;
							}
							else
								blevels[k] = TRUE;
							break;
					}
					break;
					
				case 'I':
					if (!argv[i][2]) {
						mexPrintf ("SHOREDUMP SYNTAX ERROR:  -I option takes one argument\n");
						error++;
						continue;
					}
					get_river = TRUE;
						
					switch (argv[i][2]) {
						case 'a':
							for (k = 0; k < GMT_N_RLEVELS; k++) rlevels[k] = TRUE;
							break;
						case 'r':
							for (k = 0; k < 4; k++) rlevels[k] = TRUE;
							break;
						case 'i':
							for (k = 5; k < 8; k++) rlevels[k] = TRUE;
							break;
						case 'c':
							rlevels[8] = rlevels[9] = TRUE;
							break;
						default:
							k = atoi (&argv[i][2]) - 1;
							if (k < 0 || k >= GMT_N_RLEVELS) {
								mexPrintf ("SHOREDUMP SYNTAX ERROR -I option: Feature not in list!\n");
								error++;
							}
							else
								rlevels[k] = TRUE;
							break;
					}
					break;
					
				case 'S':
					get_shore = TRUE;
					break;
				default:		/* Options not recognized */
					error = TRUE;
					break;
			}
		}
	}
	if (error) return;

	/* Check that the options selected are mutually consistant */
	if ((get_shore + get_border + get_river) == 0) get_shore = TRUE;	/* Default */

	if (!project_info.region_supplied)
		mexErrMsgTxt ("SHOREDUMP SYNTAX ERROR:  Must specify -R option\n");

	if (get_river) {
		for (k = 0; k < GMT_N_RLEVELS; k++) {
			if (!rlevels[k]) continue;
			rlevels[n_rlevels] = k + 1;
			n_rlevels++;
		}
	}
	
	if (get_border) {
		for (k = 0; k < GMT_N_BLEVELS; k++) {
			if (!blevels[k]) continue;
			blevels[n_blevels] = k + 1;
			n_blevels++;
		}
	}

	if (east > 360.0) {
		west -= 360.0;	east -= 360.0;
	}
	
#ifdef GMT_MINOR_VERSION
	if (get_shore && GMT_init_shore(res, &c, west, east, south, north, &Ainfo)) {
#else
	if (get_shore && GMT_init_shore(res, &c, west, east, south, north))  {
#endif
		mexPrintf ("SHOREDUMP: %s resolution shoreline data base not installed\n", shore_resolution[base]);
		mexErrMsgTxt("");
	}
	
	if (get_border && GMT_init_br('b', res, &b, west, east, south, north)) {
		mexPrintf ("SHOREDUMP: %s resolution political boundary data base not installed\n", shore_resolution[base]);
		mexErrMsgTxt("");
	}
	
	if (get_river && GMT_init_br('r', res, &r, west, east, south, north)) {
		mexPrintf ("SHOREDUMP: %s resolution river data base not installed\n", shore_resolution[base]);
		mexErrMsgTxt("");
	}

	if (west < -180 && east > 180) {	/* Patch for a probable GMT bug */
		west = -180;
		project_info.w = west;
		east = 180;
		project_info.e = east;
	}

	if (get_shore)			/* Get the bin size from the first one available */
		bsize = c.bsize;
	else if (get_border)
		bsize = b.bsize;
	else
		bsize = r.bsize;

	if (west < 0.0 && east > 0.0) greenwich = TRUE;
	if ((360.0 - fabs (project_info.e - project_info.w) ) < bsize)
		GMT_world_map = TRUE;

	if (GMT_world_map && greenwich) {
		edge = 180;
		shift = TRUE;					/* I WANT LONGS <-180;+180>*/
	}
	else if (!GMT_world_map && greenwich) {
		shift = TRUE;
		edge = west;	if (edge < 0.0) edge += 360.0;
		if (edge > 180.0) edge = 180.0;
	}
	else if (!GMT_world_map && (west < 0.0 && east <= 0.0)) {	/* I WANT LONGS <-180;+180>*/
		shift = TRUE;
		edge = west;	if (edge < 0.0) edge += 360.0;
		if (edge > 180.0) edge = 180.0;
	}

	west_border = floor (project_info.w / bsize) * bsize;
	east_border = ceil (project_info.e / bsize) * bsize;

	if (get_shore) {	/* Read shoreline file and extract lines */
		n = 0;
		for (ind = 0; ind < c.nb; ind++) {	/* Loop over necessary bins only */
			bin = c.bins[ind];
#ifdef GMT_MINOR_VERSION
			GMT_get_shore_bin (ind, &c);
#else
			GMT_get_shore_bin (ind, &c, min_area, min_level, max_level);
#endif
			if (c.ns == 0) continue;
			
#ifdef GMT_MINOR_VERSION
			if ((np = GMT_assemble_shore (&c, direction, FALSE, greenwich, west_border, east_border, &p)) == 0)
#else
			if ((np = GMT_assemble_shore (&c, direction, min_level, FALSE, greenwich, west_border, east_border, &p)) == 0)
#endif
				continue;

			if (first_shore) {
				for (i = 0; i < np; i++) for (k = 0; k < p[i].n; k++) n++;
				p_alloc = 2*(n + np);
				p_shore = mxCalloc (p_alloc, sizeof (double));
				n = 0;
				first_shore = FALSE;
			}
			else {
				for (i = 0; i < np; i++) for (k = 0; k < p[i].n; k++) n++;
				p_alloc += 2*(n + np);
				n = 0;
				if ((p_shore = mxRealloc(p_shore, p_alloc * sizeof (double))) == 0) {
					mexPrintf("SHOREDUMP ERROR: Could not reallocate memory\n");
					mxFree(p_shore);
					GMT_free_polygons(p, np);	GMT_free((void *)p);
					GMT_free_shore(&c);		GMT_shore_cleanup(&c);
					return;
				}
			}
		
			for (i = 0; i < np; i++) {
				for (k = 0; k < p[i].n; k++) {
					p_shore[is++] = p[i].lon[k];
					p_shore[is++] = p[i].lat[k];
				}
				p_shore[is++] = nan;		/* Separate blocks */
				p_shore[is++] = nan;
			}

			/* Free up memory */
			GMT_free_polygons(p, np);
			GMT_free((void *)p);
			GMT_free_shore(&c);
		}
		GMT_shore_cleanup(&c);
	}

	if (get_border) {	/* Read borders file */
		n = 0;
		for (ind = 0; ind < b.nb; ind++) {	/* Loop over necessary bins only */
			bin = b.bins[ind];
			GMT_get_br_bin(ind, &b, blevels, n_blevels);
			
			if (b.ns == 0) continue;
			
			if (FALSE && GMT_world_map && greenwich) {	/* Don't understand what is this for. So short-circuit it */
				left = b.bsize * (bin % b.bin_nx);	right = left + b.bsize;
				shift = ((left - edge) <= 180.0 && (right - edge) > 180.0);
			}

			if ((np = GMT_assemble_br (&b, shift, edge, &p)) == 0) continue;

			if (first_border) {
				for (i = 0; i < np; i++) for (k = 0; k < p[i].n; k++) n++;
				p_alloc = 2*(n + np);
				p_border = mxCalloc (p_alloc, sizeof (double));
				n = 0;
				first_border = FALSE;
			}
			else {
				for (i = 0; i < np; i++) for (k = 0; k < p[i].n; k++) n++;
				p_alloc += 2*(n + np);
				n = 0;
				if ((p_border = mxRealloc(p_border, p_alloc * sizeof (double))) == 0) {
					mexPrintf("SHOREDUMP ERROR: Could not reallocate memory\n");
					mxFree(p_border);
					GMT_free_polygons(p, np);
					GMT_free((void *)p);
					GMT_free_br(&b);
					GMT_br_cleanup(&b);
					return;
				}
			}

			for (i = 0; i < np; i++) {
				for (k = 0; k < p[i].n; k++) {
					p_border[ib] = p[i].lon[k];
					p_border[ib+1] = p[i].lat[k];
					ib += 2;
				}
				p_border[ib++] = nan;		/* Separate blocks */
				p_border[ib++] = nan;
			}
			
			/* Free up memory */
			GMT_free_polygons(p, np);
			GMT_free((void *)p);
			GMT_free_br(&b);
		}
		GMT_br_cleanup(&b);
	}

	if (get_river) {	/* Read rivers file */
		n = 0;
		for (ind = 0; ind < r.nb; ind++) {	/* Loop over necessary bins only */
			bin = r.bins[ind];
			GMT_get_br_bin (ind, &r, rlevels, n_rlevels);
			
			if (r.ns == 0) continue;
			
			if (GMT_world_map && greenwich) {
				left = r.bsize * (bin % r.bin_nx);	right = left + r.bsize;
				shift = ((left - edge) <= 180.0 && (right - edge) > 180.0);
			}
			
			if ((np = GMT_assemble_br (&r, shift, edge, &p)) == 0) continue;

			if (first_river) {
				for (i = 0; i < np; i++) for (k = 0; k < p[i].n; k++) n++;
				p_alloc = 2*(n + np);
				p_river = mxCalloc (p_alloc, sizeof (double));
				n = 0;
				first_river = FALSE;
			}
			else {
				for (i = 0; i < np; i++) for (k = 0; k < p[i].n; k++) n++;
				p_alloc += 2*(n + np);
				n = 0;
				if ((p_river = mxRealloc(p_river, p_alloc * sizeof(double))) == 0) {
					mexPrintf("SHOREDUMP ERROR: Could not reallocate memory\n");
					mxFree(p_river);
					GMT_free_polygons(p, np);
					GMT_free((void *)p);
					GMT_free_br(&r);
					GMT_br_cleanup(&r);
					return;
				}
			}
			
			for (i = 0; i < np; i++) {
				for (k = 0; k < p[i].n; k++) {
					p_river[ir] = p[i].lon[k];
					p_river[ir+1] = p[i].lat[k];
					ir += 2;
				}
				p_river[ir++] = nan;		/* Separate blocks */
				p_river[ir++] = nan;
			}

			/* Free up memory */
			GMT_free_br(&r);
			GMT_free_polygons(p, np);
			GMT_free((void *)p);
		}
		GMT_br_cleanup(&r);
	}

	/* Create and populate matrices for the return array(s) */
	if (nlhs == 1) {
		if (get_shore) {
			plhs[0] = mxCreateDoubleMatrix (2,(int)(is/2), mxREAL);
			pdata = mxGetPr(plhs[0]);
			bytes_to_copy = is * mxGetElementSize(plhs[0]);
			memcpy(pdata, p_shore, bytes_to_copy);
			mxFree(p_shore);
		}
		else if (get_river) {
			plhs[0] = mxCreateDoubleMatrix (2,(int)(ir/2), mxREAL);
			pdata = mxGetPr(plhs[0]);
			bytes_to_copy = ir * mxGetElementSize(plhs[0]);
			memcpy(pdata, p_river, bytes_to_copy);
			mxFree(p_river);
		}
		else {	/* Has to be a border */
			plhs[0] = mxCreateDoubleMatrix (2,(int)(ib/2), mxREAL);
			pdata = mxGetPr(plhs[0]);
			bytes_to_copy = ib * mxGetElementSize(plhs[0]);
			memcpy(pdata, p_border, bytes_to_copy);
			mxFree(p_border);
		}
	}
	else if (nlhs == 2) {
		if (get_shore && get_border) {
			plhs[0] = mxCreateDoubleMatrix (2,(int)(is/2), mxREAL);
			pdata = mxGetPr(plhs[0]);
			bytes_to_copy = is * mxGetElementSize(plhs[0]);
			memcpy(pdata, p_shore, bytes_to_copy);
			mxFree(p_shore);

			plhs[1] = mxCreateDoubleMatrix (2,(int)(ib/2), mxREAL);
			pdata = mxGetPr(plhs[1]);
			bytes_to_copy = ib * mxGetElementSize(plhs[1]);
			memcpy(pdata, p_border, bytes_to_copy);
			mxFree(p_border);
		}
		else if (get_shore && get_river) {
			plhs[0] = mxCreateDoubleMatrix (2,(int)(is/2), mxREAL);
			pdata = mxGetPr(plhs[0]);
			bytes_to_copy = is * mxGetElementSize(plhs[0]);
			memcpy(pdata, p_shore, bytes_to_copy);
			mxFree(p_shore);

			plhs[1] = mxCreateDoubleMatrix (2,(int)(ir/2), mxREAL);
			pdata = mxGetPr(plhs[1]);
			bytes_to_copy = ir * mxGetElementSize(plhs[1]);
			memcpy(pdata, p_river, bytes_to_copy);
			mxFree(p_river);
		}
		else if (get_border && get_river) {
			plhs[0] = mxCreateDoubleMatrix (2,(int)(ib/2), mxREAL);
			pdata = mxGetPr(plhs[0]);
			bytes_to_copy = ib * mxGetElementSize(plhs[0]);
			memcpy(pdata, p_border, bytes_to_copy);
			mxFree(p_border);

			plhs[1] = mxCreateDoubleMatrix (2,(int)(ir/2), mxREAL);
			pdata = mxGetPr(plhs[1]);
			bytes_to_copy = ir * mxGetElementSize(plhs[1]);
			memcpy(pdata, p_river, bytes_to_copy);
			mxFree(p_river);
		}
	}
	else if (nlhs == 3) {		/* All three were demanded */
		plhs[0] = mxCreateDoubleMatrix (2,(int)(is/2), mxREAL);
		pdata = mxGetPr(plhs[0]);
		bytes_to_copy = is * mxGetElementSize(plhs[0]);
		memcpy(pdata, p_shore, bytes_to_copy);
		mxFree(p_shore);

		plhs[1] = mxCreateDoubleMatrix (2,(int)(ib/2), mxREAL);
		pdata = mxGetPr(plhs[1]);
		bytes_to_copy = ib * mxGetElementSize(plhs[1]);
		memcpy(pdata, p_border, bytes_to_copy);
		mxFree(p_border);

		plhs[2] = mxCreateDoubleMatrix (2,(int)(ir/2), mxREAL);
		pdata = mxGetPr(plhs[2]);
		bytes_to_copy = ir * mxGetElementSize(plhs[2]);
		memcpy(pdata, p_river, bytes_to_copy);
		mxFree(p_river);
	}

	GMT_end(argc, argv);
	mxFree(argv);
}

mwSize prep_polygons(struct GMT_GSHHS_POL **p_old, mwSize np, int greenwich, int sample, double step, mwSize anti_bin) {
	/* This function will go through each of the polygons and determine
	 * if the polygon is clipped by the map boundary, and if so if it
	 * wraps around to the other side due to 360 degree periodicities
	 * A wrapped polygon will be returned as two new polygons so that
	 * this function may return more polygons that it receives.
	 * Upon return the polygons are in x,y inches, not degrees.
	 *
	 * *p is the array of np polygons
	 * greenwich is TRUE if area crosses Greenwich
	 * sample is TRUE if we need to resample the polygons to reduce point spacing
	 * step is the new maximum point separation in degrees
	 * anti_bin, if >= 0, indicates a possible problem bin at the antipole using -JE only
	 */

	mwSize k, np_new, n_use, n, start;
	double *xtmp, *ytmp;
	struct GMT_GSHHS_POL *p;

	p = *p_old;
	np_new = np;
		
	xtmp = (double *)0;	ytmp = (double *)0;

	for (k = 0; k < np; k++) {
			
		if (sample) p[k].n = GMT_fix_up_path (&p[k].lon, &p[k].lat, p[k].n, greenwich, step);
				
		/* Clip polygon against map boundary if necessary and return plot x,y in inches */
				
		if ((n = GMT_clip_to_map (p[k].lon, p[k].lat, p[k].n, &xtmp, &ytmp)) == 0) {	/* Completely outside */
			p[k].n = 0;
			continue;
		}
			
		/* Must check if polygon must be split and partially plotted at both edges of map */
				
		if ((*GMT_will_it_wrap) (xtmp, ytmp, n, &start)) {	/* Polygon does indeed wrap */
				
			/* First truncate agains left border */
						
			GMT_n_plot = (*GMT_truncate) (xtmp, ytmp, n, start, -1);
			n_use = GMT_compact_line (GMT_x_plot, GMT_y_plot, GMT_n_plot, FALSE, 0);
			if (project_info.three_D) GMT_2D_to_3D (GMT_x_plot, GMT_y_plot, project_info.z_level, GMT_n_plot);
			p[k].lon = (double *) GMT_memory ((void *)p[k].lon, (size_t)n_use, sizeof (double), GMT_program);
			p[k].lat = (double *) GMT_memory ((void *)p[k].lat, (size_t)n_use, sizeof (double), GMT_program);
			memcpy ((void *)p[k].lon, (void *)GMT_x_plot, (size_t)(n_use * sizeof (double)));
			memcpy ((void *)p[k].lat, (void *)GMT_y_plot, (size_t)(n_use * sizeof (double)));
			p[k].n = n_use;
								
			/* Then truncate agains right border */
						
			GMT_n_plot = (*GMT_truncate) (xtmp, ytmp, n, start, +1);
			n_use = GMT_compact_line (GMT_x_plot, GMT_y_plot, GMT_n_plot, FALSE, 0);
			if (project_info.three_D) GMT_2D_to_3D (GMT_x_plot, GMT_y_plot, project_info.z_level, GMT_n_plot);
			p = (struct GMT_GSHHS_POL *) GMT_memory ((void *)p, (size_t)(np_new + 1), sizeof (struct GMT_GSHHS_POL), GMT_program);
			p[np_new].lon = (double *) GMT_memory (VNULL, (size_t)n_use, sizeof (double), GMT_program);
			p[np_new].lat = (double *) GMT_memory (VNULL, (size_t)n_use, sizeof (double), GMT_program);
			memcpy ((void *)p[np_new].lon, (void *)GMT_x_plot, (size_t)(n_use * sizeof (double)));
			memcpy ((void *)p[np_new].lat, (void *)GMT_y_plot, (size_t)(n_use * sizeof (double)));
			p[np_new].n = n_use;
			p[np_new].interior = p[k].interior;
			p[np_new].level = p[k].level;
			np_new++;
		}
		else {
			n_use = GMT_compact_line (xtmp, ytmp, n, FALSE, 0);
			if (project_info.three_D) GMT_2D_to_3D (xtmp, ytmp, project_info.z_level, n_use);
			if (anti_bin > 0 && step == 0.0) {	/* Must warn for donut effect */
				if (gmtdefs.verbose) mexPrintf ("%s: GMT Warning: Antipodal bin # %d not filled!\n", GMT_program, anti_bin);
				GMT_free ((void *)xtmp);
				GMT_free ((void *)ytmp);
				continue;
			}
					
			else {
				p[k].lon = (double *) GMT_memory ((void *)p[k].lon, (size_t)n_use, sizeof (double), GMT_program);
				p[k].lat = (double *) GMT_memory ((void *)p[k].lat, (size_t)n_use, sizeof (double), GMT_program);
				memcpy ((void *)p[k].lon, (void *)xtmp, (size_t)(n_use * sizeof (double)));
				memcpy ((void *)p[k].lat, (void *)ytmp, (size_t)(n_use * sizeof (double)));
				p[k].n = n_use;
			}
		}
				
		GMT_free ((void *)xtmp);
		GMT_free ((void *)ytmp);
	}

	*p_old = p;

	return (np_new);
}
