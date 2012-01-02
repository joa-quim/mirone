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
 * Extracts country polygons as ASCII multi-segment from binary files
 * countries.bin or its DP reduced version countries_dp5.bin
 *
 *	Author Joaquim Luis
 *	Version: ??
 *	Date:	 25-06-2005
 *
 *	NOTE: this is the mex version. A similar program, country_extract, is
 *	also provided and may be run alone. Its purpuse is to be used from scripts.
 */

#define	TRUE	1
#define	FALSE	0

/* In non-Windows this is may not be necessary (or guive conflicts) */
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#include "countries.h"
#include "mex.h"

int get_country_id(char *country);
int get_countries_by_continent(char *continent, short int *id);
int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (float w, float e, float s, float n);
int check_in_region (float w, float e, float s, float n, float tol, float x, float y);
int check_polygon_in_region (float w, float e, float s, float n, float x0, float x1, float y0, float y1);
double ddmmss_to_degree (char *text);


/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	FILE	*fp, *fpc = NULL;
	struct	COUNTRY_HEADER h;
	struct	POINT p;
	short	id[N_COUNTRIES], id_prev, id_curr, id_cont = -1, id_c;
	int	i, j, n_countries = N_COUNTRIES, n_read, status, totaly_inside, partialy_inside = -1;
	int	error = FALSE, name_given = FALSE, list_countries = FALSE, got_one = FALSE;
	int	id_given = FALSE, list_given = FALSE, do_continent = FALSE;
	int	list_continents = FALSE, list_continent_countries = FALSE, got_partialy = FALSE;
	char	line[64], *country = NULL, *continent = NULL, dumbc[32], first_alloc = TRUE;
	char	*infile = NULL;
	float	max_east = 180, area = 0, tol = 5.0;
	float	west = -180.0, east = 180.0, south = -90.0, north = 90.0;
	double	west_d, east_d, south_d, north_d, *p_country, nan, *pdata, tmp_cent[2], lims[4];
	int	argc, p_alloc, is, id_diff, n_arg_no_char = 0;
	char	**argv;

	mxArray *mxCountry;
	mxArray *mxCountryName;
	mxArray *mxCountryCenter;
	mxArray *mxPolygonLims;
	mxArray *countries_struct;
	mxArray	*ptr, *ptr2;
	char	*fieldnames[4];		/* this array contains the names of the fields of the countries structure. */
	int	num_struct_fields;
	int	field_num;
	float	max_area = 0;

	id[0] = -1;	/* Default for error checking */

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
	argv[0] = "country_select";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'R':
					error += decode_R (argv[i], &west_d, &east_d, &south_d, &north_d);
					if (error) mexPrintf("%s ERROR in -R option\n", argv[0]);
					west = (float)west_d;	east = (float)east_d;
					south = (float)south_d;	north = (float)north_d;
					break;
				case 'A':
					area = (float)atof(&argv[i][2]);
					break;
					/*j = sscanf (&argv[i][2], "%f/%f", &min_area, &max_area);*/
					/*if (j == 1) max_level = MAX_LEVEL;*/
				case 'I':
					id[0] = (short int)atoi (&argv[i][2]);
					id_given = TRUE;
					break;
				case 'L':
					if (argv[i][2] == 'c')
						list_continents = TRUE;
					else if ((id_cont = atoi(&argv[i][2]))) {
						if (id_cont < 1 || id_cont > N_CONTINENTS) {
							mexPrintf("%s Wrong continent ID. Valid numbers are in range [1 - %d] %d\n", argv[0], N_CONTINENTS, id_cont);
							error = TRUE;
						}
						list_continent_countries = TRUE;
					}
					else
						list_countries = TRUE;
					break;
				case 'P':	/* Country, or country file list given */
					if ((fpc = fopen (&argv[i][2], "r")) == NULL) {
	 					country = &argv[i][2];
						name_given = TRUE;
					}
					else
						list_given = TRUE;
					break;
				case 'T':
					continent = &argv[i][2]; 
					do_continent = TRUE;
					break;
				default:
                                        error = TRUE;
					break;
			}
		}
                else if (argv[i][0] != ' ')
                        infile = argv[i];
	}

	if (argc == 1 || error) {
		mexPrintf("usage:  country_extract binary_file [-A<area>] [-I<id>] [-L] [-P<name>] [-V] > ...\n");
		mexPrintf("\t -A do not extract polygons with area inferior to <area>\n");
		mexPrintf("\t -I extract country whose ID is <id>\n");
		mexPrintf("\t -L[c|id] list all country names and IDs in binary_file\n");
		mexPrintf("\t -P extract country <name>, or if name is a file extract all countries listed there.\n");
		return;
	}

	if (list_countries) {		/* Just list country names and ID and exit */
		for (i = 0; i < N_COUNTRIES; i++)
			mexPrintf("Country name: %18s\tID = %d\n", country_list[i], i);
		return;
	}

	if (list_continents) {		/* Just list continent names and ID and exit */
		for (i = 0; i < N_CONTINENTS; i++)
			mexPrintf("Continent name: %18s\tID = %d\n", continent_list[i], i);
		return;
	}

	if (list_continent_countries) {	/* Just list coutries by continent and ID and exit */
		if (id_cont == 0)
			for (i = 0; i < N_AFR; i++)
				mexPrintf("Country name: %18s\tID = %d\n", afric_country_list[i], i);
		else if (id_cont == 1)
			for (i = 0; i < N_ANTARCTICA; i++)
				mexPrintf("Country name: %18s\tID = %d\n", antarctica_country_list[i], i);
		else if (id_cont == 2)
			for (i = 0; i < N_ASIA; i++)
				mexPrintf("Country name: %18s\tID = %d\n", asia_country_list[i], i);
		else if (id_cont == 3)
			for (i = 0; i < N_AUSTRALASIA; i++)
				mexPrintf("Country name: %18s\tID = %d\n", australasia_country_list[i], i);
		else if (id_cont == 4)
			for (i = 0; i < N_CENTRAL_AM; i++)
				mexPrintf("Country name: %18s\tID = %d\n", central_america_country_list[i], i);
		else if (id_cont == 5)
			for (i = 0; i < N_SOUTH_AM; i++)
				mexPrintf("Country name: %18s\tID = %d\n", south_america_country_list[i], i);
		else if (id_cont == 6)
			for (i = 0; i < N_EURO; i++)
				mexPrintf("Country name: %18s\tID = %d\n", euro_country_list[i], i);
		else if (id_cont == 7)
			for (i = 0; i < N_NORTH_AM; i++)
				mexPrintf("Country name: %18s\tID = %d\n", north_america_country_list[i], i);
		return;
	}

	/* ------------------------------------------------------------------------------------ */

	if ((fp = fopen (infile, "rb")) == NULL ) {
		mexPrintf ("country_extract:  Could not find file %s.\n", infile);
		return;
	}

	if (do_continent) {	/* Extract all countries in continent */
		status = get_countries_by_continent(continent,id);
		if (status < 0) {
			mexPrintf("%s: continent %s not found in data base\n", argv[0], continent);
			return;
		}
	}

	if (name_given) {
		id[0] = (short int)get_country_id(country);
		if (id[0] == -1) {
			mexPrintf("%s: country %s not found in data base\n", argv[0], country);
			return;
		}
	}

	if (list_given) {	/* Read country names from file */
		n_countries = 0;
		while (fgets (line, BUFSIZ, fpc)) {
			sscanf (line, "%s", dumbc);
			id[n_countries] = (short int)get_country_id(dumbc);
			n_countries++;
		}
		fclose(fpc);
	}

	if (name_given + list_countries + id_given + list_given + do_continent == 0) {	/* extract all countries in file */
		for (i = 0; i < N_COUNTRIES; i++) id[i] = (short int)i;
	}
	else if (name_given || id_given) {
		for (i = 1; i < N_COUNTRIES; i++) id[i] = -1;
	}

	/* ---------------------------------------------------------------------------------------- */
	/* We are ready to go */

	/* Create the metadata structure. Just one element, with XXX fields. */
	num_struct_fields = 4;
	fieldnames[0] = strdup ("Country");
	fieldnames[1] = strdup("Tag");
	fieldnames[2] = strdup("Centroide");
	fieldnames[3] = strdup("Limits");
	countries_struct = mxCreateStructMatrix ( N_COUNTRIES, 1, num_struct_fields, (const char **)fieldnames );

	id_c = is = 0;
	n_read = fread ((void *)&h, (size_t)sizeof (struct COUNTRY_HEADER), (size_t)1, fp);
	tmp_cent[0] = (double) h.centroid.x;		/* Initialize these */
	tmp_cent[1] = (double) h.centroid.y;
	nan = mxGetNaN();

	while (n_read == 1) {
		if (h.country_id == id[id_c]) {
			id_curr = id[id_c];
			if (!strcmp(country_list[id_c],"russia"))	tol = 20;
			else						tol = 2;
			totaly_inside = check_polygon_in_region(west, east, south, north, h.west, h.east, h.south, h.north);
			if (totaly_inside == 0) {
				partialy_inside = check_polygon_in_region (west, east, south, north, h.west, h.east, h.south, h.north);
			}
			if ((h.area >= area) && totaly_inside == 1) {	/* Easy case */
				got_one = TRUE;
				if (first_alloc) {
					p_alloc = (2 * (h.n + 1));
					p_country = mxCalloc (p_alloc, sizeof (double));
					first_alloc = FALSE;
				}
				else {
					p_alloc += (2 * (h.n + 1));
					p_country = mxRealloc(p_country, p_alloc * sizeof (double));
				}
				for (j = 0; j < h.n; j++) {
					if (fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp) != 1) {
						mexPrintf ("Error reading file %s for country %d, point %d.\n", argv[1], h.country_id, j);
						return;
					}
					p_country[is++] = p.x;	p_country[is++] = p.y;
				}
				p_country[is++] = nan;	p_country[is++] = nan; 	/* Separate polygons */
				if (max_area < h.area) {	/* Keep track of maximum polygon area to report the centroide of it */
					max_area = h.area;
					tmp_cent[0] = (double) h.centroid.x;
					tmp_cent[1] = (double) h.centroid.y;
					lims[0] = (double) h.west;
					lims[1] = (double) h.east;
					lims[2] = (double) h.south;
					lims[3] = (double) h.north;
				}
			}
			else if ((h.area >= area) && partialy_inside == 0) {	/* We have a partial contained polygon inside region */
				got_one = TRUE;
				fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp);
				if (first_alloc) {
					p_alloc = 2 * (h.n + 1);
					p_country = mxCalloc (p_alloc, sizeof (double));
					first_alloc = FALSE;
				}
				else {
					p_alloc += 2 * (h.n + 1);
					p_country = mxRealloc(p_country, p_alloc * sizeof (double));
				}
				if (check_in_region(west,east,south,north,tol,p.x,p.y)) {
					p_country[is++] = p.x;	p_country[is++] = p.y;
					got_partialy = TRUE;
				}
				for (j = 0; j < h.n-1; j++) {	/* n-1 because we already read one record above */
					if (fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp) != 1) {
						mexPrintf ("Error reading file %s for country %d, point %d.\n", argv[1], h.country_id, j);
						return;
					}
					if (check_in_region(west,east,south,north,tol,p.x,p.y)) {
						p_country[is++] = p.x;	p_country[is++] = p.y;
						got_partialy = TRUE;
					}
				}
				if (got_partialy)	/* We got an entry in the middle of polygon. So we must end it with NaNs */
					p_country[is++] = nan;	p_country[is++] = nan; 	/* Separate polygons */
				if (max_area < h.area) {	/* Keep track of maximum polygon area to report the centroide of it */
					max_area = h.area;
					tmp_cent[0] = (double) h.centroid.x;
					tmp_cent[1] = (double) h.centroid.y;
					lims[0] = (double) h.west;
					lims[1] = (double) h.east;
					lims[2] = (double) h.south;
					lims[3] = (double) h.north;
				}
			}
			else	/* Skip non-wanted polygon */
				fseek (fp, (long)(h.n * sizeof(struct POINT)), SEEK_CUR);
		}
		else {		/* Skip data of non-wanted country */
			fseek (fp, (long)(h.n * sizeof(struct POINT)), SEEK_CUR);
			is = 0;
		}
		n_read = fread ((void *)&h, (size_t)sizeof (struct COUNTRY_HEADER), (size_t)1, fp);
		id_prev = id_curr;
		id_curr = h.country_id;
		if (id_curr != id_prev || id_prev == N_COUNTRIES-1) {	/* Got next country */
			if (got_one) {		/* Got one output. Fill in the structure */

				mxCountry = mxCreateDoubleMatrix (2,(int)(is/2), mxREAL);
				pdata = mxGetPr(mxCountry);
				memcpy(pdata, p_country, is * sizeof(double));
				mxSetField(countries_struct, id_prev, "Country" , mxCountry);

				mxCountryName = mxCreateString(country_list[id[id_c]]);
				mxSetField(countries_struct, id_prev, "Tag", mxCountryName);

				mxCountryCenter = mxCreateDoubleMatrix (1, 2, mxREAL);
				pdata = mxGetPr(mxCountryCenter);
				memcpy(pdata, tmp_cent, 2 * sizeof(double));
				mxSetField(countries_struct, id_prev, "Centroide" , mxCountryCenter);

				mxPolygonLims = mxCreateDoubleMatrix (1, 4, mxREAL);
				pdata = mxGetPr(mxPolygonLims);
				memcpy(pdata, lims, 4 * sizeof(double));
				mxSetField(countries_struct, id_prev, "Limits" , mxPolygonLims);

				mxFree(p_country);
				got_one = FALSE;
			}
			/* Re-Initialize those */
			is = 0;
			partialy_inside = -1;
			max_area = 0;
			tmp_cent[0] = (double) h.centroid.x;
			tmp_cent[1] = (double) h.centroid.y;
			got_partialy = FALSE;
			first_alloc = TRUE;
			/*if (h.country_id >= 112 && h.country_id < 127) mexPrintf("%s\tcurr=%d prev=%d id_c=%d %s", country_list[h.country_id], id_curr, id_prev, id_c, country_list[id[id_c]]);*/
			if (h.country_id > id[id_c]) {		/* This is needed with decimated versions of the database */
				/*id_diff = id_curr - id_prev;*/
				id_diff = id_curr - id[id_c];
				if (id_diff > 1)
					id_c += id_diff;
				else
					id_c++;
				if (id[id_c] < 0)
					break;	/* No more countries to extract. We are done */
			}
			/*if (h.country_id >= 112 && h.country_id < 127) mexPrintf("\tid_c=%d %s\n", id_c, country_list[id[id_c]]); */
		}
	}
	fclose (fp);

	/*for (i = 0; i < N_COUNTRIES; i++) {
		field_num = mxGetFieldNumber(countries_struct, country_list[i]); 
		mexPrintf("MMM %s %d\n", country_list[i], field_num);
		ptr = mxGetFieldByNumber(countries_struct, i, 0);
		ptr2 = mxGetFieldByNumber(countries_struct, i, 1);
		if (ptr == NULL) {
			mxFree(ptr);
			mxFree(ptr2);
			mexErrMsgTxt("MMMMMMMaaaaaaa \n");
		}
	}*/

	plhs[0] = countries_struct;
}

int get_country_id(char *country) {
	int i, id = -1;
	for (i = 0; i < N_COUNTRIES; i++) {
		if (!strcmp(country,country_list[i])) {
			id = i;
			break;
		}
	}
	return ((short)id);
} 

int get_countries_by_continent(char *continent, short int *id) {
	int i, k = -1;

	for (i = 0; i < N_CONTINENTS; i++) {	/* Find the continent first */
		if (!strcmp(continent,continent_list[i])) {
			k = i;
			break;
		}
	}
	if (k < 0) return (k);		/* continent not found */
	if (!strcmp(continent_list[k],"africa"))
		for (i = 0; i < N_AFR; i++)
			id[i] = (short int)get_country_id(afric_country_list[i]);
	else if (!strcmp(continent_list[k],"antarctica"))
		for (i = 0; i < N_ANTARCTICA; i++)
			id[i] = (short int)get_country_id(antarctica_country_list[i]);
	else if (!strcmp(continent_list[k],"asia"))
		for (i = 0; i < N_ASIA; i++)
			id[i] = (short int)get_country_id(asia_country_list[i]);
	else if (!strcmp(continent_list[k],"australasia"))
		for (i = 0; i < N_AUSTRALASIA; i++)
			id[i] = (short int)get_country_id(australasia_country_list[i]);
	else if (!strcmp(continent_list[k],"central_america"))
		for (i = 0; i < N_CENTRAL_AM; i++)
			id[i] = (short int)get_country_id(central_america_country_list[i]);
	else if (!strcmp(continent_list[k],"south_america"))
		for (i = 0; i < N_SOUTH_AM; i++)
			id[i] = (short int)get_country_id(south_america_country_list[i]);
	else if (!strcmp(continent_list[k],"europe"))
		for (i = 0; i < N_EURO; i++)
			id[i] = (short int)get_country_id(euro_country_list[i]);
	else if (!strcmp(continent_list[k],"north_america"))
		for (i = 0; i < N_NORTH_AM; i++)
			id[i] = (short int)get_country_id(north_america_country_list[i]);
	return (0);
}

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

int check_region (float w, float e, float s, float n) {
	/* If region is given then we must have w < e and s < n */
	return ((w >= e || s >= n));
}

int check_in_region (float w, float e, float s, float n, float tol, float x, float y) {
	/* See if the point (x,y) is inside the rectangle widened by tol */
	return ( (x >= (w - tol)) && (x <= (e + tol)) && (y >= (s - tol)) && (y <= (n + tol)) );
}

int check_polygon_in_region (float w, float e, float s, float n, float x0, float x1, float y0, float y1) {
	/* See if the rectangle deffined by -Rx0/x1/y0/y1 is totaly or partialy inside the region -Rw/e/s/n
	   If -Rx0/x1/y0/y1 is totaly inside return 1
	   If -Rx0/x1/y0/y1 is partially inside return 0
	   Else return -1 */

	if (x0 >= w && x1 <= e && y0 >= s && y1 <= n)
		return (1);
	else if ( ((x0 >= w && x0 <= e) || (x1 >= w && x1 <= e)) && ((y0 >= s && y0 <= n) || (y1 >= s && y1 <= n)) )
		return (0);
	return (-1);
}
