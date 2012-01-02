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
 *	NOTE: this is the stand-alone version of the mex country_select
 */

#define	TRUE	1
#define	FALSE	0

/* In non-Windows this is may not be necessary (or guive conflicts) */
#define copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#include "countries.h"

int get_country_id(char *country);
void get_rand_color(int seed, char *opt_G);
int get_countries_by_continent(char *continent, short int *id);
int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (float w, float e, float s, float n);
int check_in_region (float w, float e, float s, float n, float tol, float x, float y);
int check_polygon_in_region (float w, float e, float s, float n, float x0, float x1, float y0, float y1);
int time();
double ddmmss_to_degree (char *text);

int main (int argc, char **argv) {
	FILE	*fp, *fpc = NULL;
	struct	COUNTRY_HEADER h;
	struct	POINT p;
	short	id[N_COUNTRIES], id_prev, id_curr, id_cont = -1;
	int	i, j, n_countries = N_COUNTRIES, n_read, status, flaged, totaly_inside, partialy_inside;
	int	error = FALSE, name_given = FALSE, list_countries = FALSE, do_random_color = FALSE;
	int	verbose = FALSE, id_given = FALSE, list_given = FALSE, do_continent = FALSE;
	int	list_continents = FALSE, list_continent_countries = FALSE;
	char	line[64], *country = NULL, *continent = NULL, dumbc[32];
	char	*infile = NULL, opt_G[16];
	float	lon, max_east = 180, area = 0, tol = 5.0;
	float	west = -180.0, east = 180.0, south = -90.0, north = 90.0;
	double	west_d, east_d, south_d, north_d;

	id[0] = -1;	/* Default for error checking */

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'R':
					error += decode_R (argv[i], &west_d, &east_d, &south_d, &north_d);
					if (error) fprintf(stderr, "%s ERROR in -R option\n", argv[0]);
					west = (float)west_d;	east = (float)east_d;
					south = (float)south_d;	north = (float)north_d;
					break;
				case 'A':
					area = (float)atof(&argv[i][2]);
					break;
				case 'C':
					do_random_color = TRUE;
					break;
				case 'I':
					id[0] = (short int)atoi (&argv[i][2]);
					id_given = TRUE;
					break;
				case 'L':
					if (argv[i][2] == 'c')
						list_continents = TRUE;
					else if ((id_cont = atoi(&argv[i][2]))) {
						if (id_cont < 1 || id_cont > N_CONTINENTS) {
							fprintf(stderr, "%s Wrong continent ID. Valid numbers are in range [1 - %d] %d\n", argv[0], N_CONTINENTS, id_cont);
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
				case 'V':
					verbose = TRUE;
					break;
				default:
                                        error = TRUE;
					break;
			}
		}
                else
                        infile = argv[i];
	}

	if (argc == 1 || error) {
		fprintf(stderr,"usage:  country_extract binary_file [-A<area>] [-C] [-I<id>] [-L] [-P<name>] [-V] > ...\n");
		fprintf(stderr,"\t -A do not extract polygons with area inferior to <area>\n");
		fprintf(stderr,"\t -C use a random color schema\n");
		fprintf(stderr,"\t -I extract country whose ID is <id>\n");
		fprintf(stderr,"\t -L[c|id] list all country names and IDs in binary_file. Append c to list only continents\n");
		fprintf(stderr,"\t -P extract country <name>, or if name is a file extract all countries listed there.\n");
		fprintf(stderr,"\t -V verbose\n");
		exit(-1);
	}

	if (list_countries) {		/* Just list country names and ID and exit */
		for (i = 0; i < N_COUNTRIES; i++)
			fprintf(stderr, "Country name: %18s\tID = %d\n", country_list[i], i);
		return (0);
	}

	if (list_continents) {		/* Just list continent names and ID and exit */
		for (i = 0; i < N_CONTINENTS; i++)
			fprintf(stderr, "Continent name: %18s\tID = %d\n", continent_list[i], i);
		return (0);
	}

	if (list_continent_countries) {	/* Just list coutries by continent and ID and exit */
		if (id_cont == 0)
			for (i = 0; i < N_AFR; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", afric_country_list[i], i);
		else if (id_cont == 1)
			for (i = 0; i < N_ANTARCTICA; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", antarctica_country_list[i], i);
		else if (id_cont == 2)
			for (i = 0; i < N_ASIA; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", asia_country_list[i], i);
		else if (id_cont == 3)
			for (i = 0; i < N_AUSTRALASIA; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", australasia_country_list[i], i);
		else if (id_cont == 4)
			for (i = 0; i < N_CENTRAL_AM; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", central_america_country_list[i], i);
		else if (id_cont == 5)
			for (i = 0; i < N_SOUTH_AM; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", south_america_country_list[i], i);
		else if (id_cont == 6)
			for (i = 0; i < N_EURO; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", euro_country_list[i], i);
		else if (id_cont == 7)
			for (i = 0; i < N_NORTH_AM; i++)
				fprintf(stderr, "Country name: %18s\tID = %d\n", north_america_country_list[i], i);
		return (0);
	}

	/* ------------------------------------------------------------------------------------ */

	if ((fp = fopen (infile, "rb")) == NULL ) {
		fprintf (stderr, "country_extract:  Could not find file %s.\n", infile);
		exit (EXIT_FAILURE);
	}

	if (do_continent) {	/* Extract all countries in continent */
		status = get_countries_by_continent(continent,id);
		if (status < 0) {
			fprintf(stderr,"%s: continent %s not found in data base\n", argv[0], continent);
			exit (EXIT_FAILURE);
		}
		if (verbose) fprintf(stderr,"Extracting all countries in %s\n", continent);
	}

	if (name_given) {
		id[0] = (short int)get_country_id(country);
		if (id[0] == -1) {
			fprintf(stderr,"%s: country %s not found in data base\n", argv[0], country);
			exit(-1);
		}
	}

	if (list_given) {	/* Read country names from file */
		if (verbose) fprintf(stderr,"Extracting from a country list\n");
		n_countries = 0;
		while (fgets (line, BUFSIZ, fpc)) {
			sscanf (line, "%s", dumbc);
			id[n_countries] = (short int)get_country_id(dumbc);
			n_countries++;
		}
		fclose(fpc);
	}

	if (name_given + list_countries + id_given + list_given + do_continent == 0) {	/* extract all countries in file */
		for (i = 0; i < N_COUNTRIES; i++) id[i] = (short int)(i+1);
		if (verbose) fprintf(stderr,"Extracting all countries\n");
	}
	else if (name_given || id_given) {
		for (i = 1; i < N_COUNTRIES; i++) id[i] = -1;
		if (verbose) fprintf(stderr,"Extracting country:\t%s\n", country_list[id[0]-1]);
	}

	if (do_random_color)
		get_rand_color(500, opt_G);

	/* ---------------------------------------------------------------------------------------- */
	/* We are ready to go */

	i = 0;
	partialy_inside = -1;
	n_read = fread ((void *)&h, (size_t)sizeof (struct COUNTRY_HEADER), (size_t)1, fp);
	while (n_read == 1) {
		if (h.country_id == id[i]) {
			if (!strcmp(country_list[i],"russia"))
				tol = 40;
			else
				tol = 5;
			totaly_inside = check_polygon_in_region(west, east, south, north, h.west, h.east, h.south, h.north);
			if (totaly_inside == 0)
				partialy_inside = check_polygon_in_region (west, east, south, north, h.west, h.east, h.south, h.north);
			if ((h.area >= area) && totaly_inside == 1) {	/* Easy case */
				(do_random_color) ? printf ("> %s\t -W0.5p\n", opt_G) : printf (">\n");
				for (j = 0; j < h.n; j++) {
					if (fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp) != 1) {
						fprintf (stderr, "Error reading file %s for country %d, point %d.\n", argv[1], h.country_id, j);
						exit (EXIT_FAILURE);
					}
					lon = (h.greenwich && p.x > max_east) ? p.x - 360.0 : p.x;
					printf ("%.5f\t%.5f\n", lon, p.y);
				}
			}
			else if ((h.area >= area) && partialy_inside == 0) {	/* We have a partial contained polygon inside region */
				flaged = FALSE;
				fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp);
				if (check_in_region(west,east,south,north,tol,p.x,p.y)) {
					(do_random_color) ? printf ("> %s\t -W0.5p\n", opt_G) : printf (">\n");
					printf ("%.5f\t%.5f\n", p.x, p.y);
					flaged = TRUE;
				}
				for (j = 0; j < h.n-1; j++) {	/* n-1 because we already read one record above */
					if (fread ((void *)&p, (size_t)sizeof(struct POINT), (size_t)1, fp) != 1) {
						fprintf (stderr, "Error reading file %s for country %d, point %d.\n", argv[1], h.country_id, j);
						exit (EXIT_FAILURE);
					}
					lon = (h.greenwich && p.x > max_east) ? p.x - 360.0 : p.x;
					if (check_in_region(west,east,south,north,tol,lon,p.y)) {
						if (!flaged) {	/* Just entering padded region. Time to throw the multisegment flag */
							(do_random_color) ? printf ("> %s\t -W0.5p\n", opt_G) : printf (">\n");
							printf ("%.5f\t%.5f\n", lon, p.y);
							flaged = TRUE;
						}
						else
							printf ("%.5f\t%.5f\n", lon, p.y);
					}
				}

			}
			else	/* Skip non-wanted polygon */
				fseek (fp, (long)(h.n * sizeof(struct POINT)), SEEK_CUR);
		}
		else {		/* Skip data of non-wanted countries */
			fseek (fp, (long)(h.n * sizeof(struct POINT)), SEEK_CUR);
		}
		n_read = fread ((void *)&h, (size_t)sizeof (struct COUNTRY_HEADER), (size_t)1, fp);
		id_prev = id_curr;
		id_curr = h.country_id;
		if (id_curr != id_prev) {	/* got next country */
			if (h.country_id > id[i]) {
				i++;
				if (id[i] < 0) break;	/* No more countries to extract. We are done */
				if (do_random_color) get_rand_color(i, opt_G);
			}
		}
		partialy_inside = -1;
	}
	fclose (fp);

	return(0);
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

void get_rand_color(int seed, char *opt_G) {
	int c1, c2, c3, s;
	s = (unsigned int)time(NULL) % 100000;
	srand(s);		c1 = (int) ((float)rand() / (float)RAND_MAX * 255); 
	srand(s*seed*3);	c2 = (int) ((float)rand() / (float)RAND_MAX * 255); 
	srand(s*seed*7);	c3 = (int) ((float)rand() / (float)RAND_MAX * 255);
	sprintf(opt_G,"-G%03d/%03d/%03d",c1,c2,c3);
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
		degfrac = degree + copysign (minute / 60.0 + second / 3600.0, degree);
	}
	else if (colons == 1) {	/* dd:mm format */
		sscanf (text, "%lf:%lf", &degree, &minute);
		degfrac = degree + copysign (minute / 60.0, degree);
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
