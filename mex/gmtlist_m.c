/*--------------------------------------------------------------------
 *	$Id$
 *
 *    Copyright (c) 1991-2001 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------*/
/*
 * gmtlist produces ASCII listings of <legid>.gmt files. The *.gmt files
 * contains time(s), latitude(y), longitude(x), gravity(g), magnetics(m),
 * and bathymetry(t), and the user may extract any combination of these 6
 * parameters + distance (in km), heading, velocity, and weight by using
 * the option -Fsxygmtdhvw.  The sequence in which the flag characters
 * appear determines the sequence in which the parameters will be printed
 * out. If no options is specified, the default is -Fsxygmtdhvw.  E.g. to
 * create an input file for surface, use -Fxyg (for gravity).
 * If upper case letters are used for gmt (GMT), then only records that have
 * that particular data are written out.  E.g -Fxyg gives
 * lon/lat/grav, whereas -FxyG gives lon/lat/grav where there is gravity data.
 * To select a section of the track, specify the start/endpoints by:
 *	1) Start-time (mm/dd/yyyy/hh:mm) OR start-distance (km)
 *	2) Stop-time (mm/dd/yyyy/hh:mm) OR stop-distance (km)
 * To select data inside an area, use the -R option.
 * To start output with a header string, use -H.
 * Several formats for time is available. The default when using -Fs is
 * seconds from Jan 1 the year the cruise started. -Fsc (calender) gives
 * yyyymmddhhmms output, and -Fsj (julian) gives yyyyjjhhmmss output.
 *
 * Author:	Paul Wessel
 * Date:	19-APR-1988
 * Version:	2.1 1-JUL-1992
 *		3.2 10-MAR-1999
 *		3.3.1 25-JUN-1999
 *
 * -------------------------------------------------------------------------
 *
 * Mexified version of gmtlist, so some of what is stated above isn't true anymore
 * Also added an extra field (info) that contains cruise information like the
 * one provided by gmtinfo.
 * Mexifier:	Joaquim Luis
 * Date:	30-APR-2005
 *	 
 *		28/07/15 J Luis, Make it stand-alone
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 */
 
#include "mex.h"
#include <string.h>
#include <math.h>

#define KMPRDEG 111.1949e-6
#define MAXLEGS 5000
#define S_PR_DAY 86400
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define D2R (M_PI / 180.0)
#define R2D (180.0 / M_PI)
#define CNULL	((char *)NULL)
#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))

#ifndef rint
#define rint(x) (floor((x)+0.5))
#endif
#ifndef irint
#define irint(x) ((int)rint(x))
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))	/* min and max value macros */
#endif
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define GMTMGG_NODATA (-32000)	/* .gmt file NaN proxy */
#define MDEG2DEG	0.000001	/* Convert millidegrees to degrees */
#define NGDC_OLDEST_YY	39	/* Oldest NGDC cruise is from 1939 */

#define GMTMGG_TIME_MAXMONTH	61	/* 5 years is a long time for one cruise */
#define REC_SIZE 40	/* Rec size for xx_base.b file xover-records and struct XOVERS */
#define NODATA (-32000)

#ifdef WIN32	/* Start of Windows setup */
/* fileno and setmode have leading _ under WIN32 */
#include <io.h>
#define R_OK 04

#define fileno(stream) _fileno(stream)
#define setmode(fd,mode) _setmode(fd,mode)
#endif
#define GMT_swab2(data) ((((data) & 0xff) << 8) | ((unsigned short) (data) >> 8))
#define GMT_swab4(data) \
	(((data) << 24) | (((data) << 8) & 0x00ff0000) | \
	(((data) >> 8) & 0x0000ff00) | ((unsigned int)(data) >> 24))

struct GMTMGG_TIME {
  int daymon[GMTMGG_TIME_MAXMONTH];	/* Cumulative number of days up to last month */
  int first_year;			/* The year the cruise started */
};

struct GMTMGG_REC {	/* Format of *.gmt file records */
	int time;
	int lat;
	int lon;
	short int gmt[3];
};

struct LEG {	/* Structure with info about one leg */
	char name[10];			/* Name of leg */
	char agency[10];		/* Collecting agency */
	int year;			/* Year the leg started */
	int n_x_int;			/* Total number of internal cross-over points */
	int n_x_ext;			/* Total number of external cross-over points */
	int n_gmtint[3];		/* Number of internal gravity/magnetics/topography crossovers */
	int n_gmtext[3];		/* Number of external gravity/magnetics/topography crossovers */
	double mean_gmtint[3];		/* Mean gravity/magnetics/topography internal xover value */
	double mean_gmtext[3];		/* Mean gravity/magnetics/topography external xover value */
	double st_dev_gmtint[3];	/* St. Dev. of the internal gravity/magnetics/topography crossovers */
	double st_dev_gmtext[3];	/* Same for external xovers */
	double dc_shift_gmt[3];		/* Best fitting d.c.-shift for gravity/magnetics/topography */
	double drift_rate_gmt[3];	/* Best fitting drift rate for gravity/magnetics/topography */
	struct LEG *next_leg;		/* Pointer to next leg in list */
};

struct CORR {	/* Structure with the corrections for each leg */
	char name[10];			/* Name of leg */
	short int year;			/* Year the leg started */
	float dc_shift_gmt[3];		/* Best fitting d.c.-shift for gravity, magnetics, and topo */
	float drift_rate_gmt[3];	/* Best fitting drift-rate for gravity, magnetics, and topo */
};
struct CORR **bin;

int binsize = sizeof(struct CORR);
int nlegs = 0;

int get_id (char *name);
int gmtmgg_date (int time, int *year, int *month, int *day, int *hour, int *minute, int *second, struct GMTMGG_TIME *gmt_struct);
int gmtmgg_time (int *time, int year, int month, int day, int hour, int minute, int second, struct GMTMGG_TIME *gmt_struct);
struct GMTMGG_TIME *gmtmgg_init (int year1);
int gmtmggpath_func (char *leg_path, char *leg);
int decode_R (char *item, double *w, double *e, double *s, double *n);
int check_region (double w, double e, double s, double n);
double ddmmss_to_degree (char *text);

/* --------------------------------------------------------------------------- */
/* Matlab Gateway routine */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int leg_year, rec, i, j, n_records, argno, id, no[10], nval = 0;
	int begin_day, begin_mo, begin_yr, end_day, end_mo, end_yr;
	int hour, minute, second, n_grv, n_mag, n_top;
	int start_time = 0, stop_time = 2000000000, dlon, last_lon = 0;
	int mon1, day1 = 0, year1, hour1, min1, time, last_lat = 0, mon2, day2 = 0, year2, hour2, min2;
	int dt, last_time = 0, n_cruises = 0;
	
	int error = FALSE, wantgmt, want_all = FALSE, geodetic = TRUE, swapa = FALSE;
	int correct = FALSE, tsec = FALSE, calender = FALSE, do_heading = FALSE;
	int no_g, no_m, no_t, greenwich = FALSE, do_speed = FALSE, binary = FALSE;
	
	char gmtfile[BUFSIZ], agency[10], g, m, t, greenwich_in;
	
	double lat, lon, grv, mag, top, dist, start_dist = 0., stop_dist = 1.0E100;
	double ds, dx, dy, west = 0.0, east = 360.0, south = -90.0, north = 90.0;
	double heading, speed, weight = 1.0;
	double xmin1, xmax1, xmin2, xmax2, ymin, ymax, xmin, xmax, old_x = 0.0;

	int	argc, n_arg_no_char = 0, num_struct_fields;
	double	*pdata, *p_out[10], nan, tmp, info[14];
	char	**argv;
	char	*fieldnames[13];	/* this array contains the names of the fields of the gmtlist structure. */
	mxArray *mxOut, *mxStr;
	mxArray *gmtlist_struct;
	
	struct GMTMGG_TIME *gmt;
	
	FILE *fp;
	
	struct GMTMGG_REC record;
	
	g = m = t = no_g = no_m = no_t = FALSE;
	
	/* Check and interpret the command line arguments */

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
	argv[0] = "gmtlist_m";
	for (i = 1; i < argc; i++) {
		argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	for (i = 1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
			
				case 'R':
				case '\0':
					error += decode_R (argv[i], &west, &east, &south, &north);
					break;
					
				case 'D':		/* Assign start/stop times for sub-section */
					if (argv[i][2] == 'a') {	/* Start date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d", &mon1, &day1, &year1, &hour1, &min1);
					}
					else if (argv[i][2] == 'b')	 {	/* Stop date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d", &mon2, &day2, &year2, &hour2, &min2);
					}
					else
						error = TRUE;
					break;

				case 'F':	/* Selected output fields */
					for (j = 2; argv[i][j]; j++) {
						switch (argv[i][j]) {
							case 'G':	/* Records with gravity != GMTMGG_NODATA requested */
								no_g = TRUE;
	       						case 'g':		/* Gravity is requested */
								no[nval++] = 3;
								g = 1;
								break;
							case 'M':	/* Records with magnetics != GMTMGG_NODATA requested */
								no_m = TRUE;
							case 'm':		/* Magnetics is requested */
								no[nval++] = 4;
								m = 1;
								break;
							case 'T':	/* Records with topo != GMTMGG_NODATA requested */
								no_t = TRUE;
							case 't':		/* Topography is requested */
								no[nval++] = 5;
								t = 1;
								break;
							case 'x':		/* Longitude is requested */
								no[nval++] = 1;
								break;
							case 'y':		/* Latitude is requested */
								no[nval++] = 2;
								break;
							case 's':		/* Time (in sec) is requested */
								no[nval++] = 0;
								if (argv[i][j+1] == 'c') {
									calender = TRUE;
									j++;
								}
								else if (argv[i][j+1] == 'j')
									j++;
								else
									tsec = TRUE;
								break;
							case 'd':		/* Distance (in km) is requested */
								no[nval++] = 6;
								break;
							case 'h':		/* Heading is requested */
								no[nval++] = 7;
								do_heading = TRUE;
								break;
							case 'v':		/* velocity (in m/s) is requested */
								no[nval++] = 8;
								do_speed = TRUE;
								break;
							case 'w':		/* weights (Set with -W) is requested */
								no[nval++] = 9;
								break;
							default:
								error = TRUE;
								break;
						}
					}
					break;
						
				case 'S':		/* Assign start/stop position for sub-section */
					if (argv[i][2] == 'a')	/* Start position */
						start_dist = atof(&argv[i][3]);
					else if (argv[i][2] == 'b')	/* Stop position */
						stop_dist = atof(&argv[i][3]);
					else
						error = TRUE;
					break;
					
				case 'G':
					geodetic = FALSE;
					break;
				case 'Y':
					swapa = TRUE;
					break;
				case 'W':		/* Assign a weight to these data */
					weight = atof (&argv[i][2]);
					break;
					
				default:		/* Options not recognized */
					error = TRUE;
					break;
			}
		}
		else
			n_cruises++;
	}
	
	/* Check that the options selected are mutually consistent */
	
	if (nval > 10) error = TRUE;
	if (start_dist > stop_dist || start_time > stop_time) error = TRUE;
	if ((day1 > 0 && start_dist > 0.) || (day2 > 0 && stop_dist < 1.0e100)) error = TRUE;
	if (east < west || south > north) error = TRUE;
	if (n_cruises == 0) error = TRUE;
	if (weight <= 0.0) error = TRUE;

	if (error || argc == 1) {	/* Display usage */
		mexPrintf("usage: OUT = gmtlist_m(<cruise(s)>, '[-Da<startdate>]', '[-Db<stopdate>]', '[-F<dataflags>]',\n");
		mexPrintf("	'[-G]', '[-R<west>/<east>/<south>/<north>]', '[-Sa<startdist>]', '[-Sb<stopdist>]', '[-W<Weight>]', '[-Y]')\n\n");
         
		mexPrintf("	OUT is a MxN structure matrix where M = n_cruises and N = n_dataflags + 3 (-F<option> + 3) \n");
		mexPrintf("	the three extra fields contain the cruises''s YEAR, AGENCY info and INFO (lower cases).\n");
		mexPrintf("	The field names of the OUT structure follow the <dataflags> nomenclature\n");
		mexPrintf("	The INFO field is a 1x14 row vector with:\n");
		mexPrintf("	[n_records,n_grv,n_mag,n_top,xmin,xmax,ymin,ymax,...\n");
		mexPrintf("	 begin_day,begin_mo,begin_yr,end_day,end_mo,end_yr]\n\n");

		mexPrintf("	<cruises> is one or more legnames, e.g. c2104 v3206 etc.\n");
		mexPrintf("	OPTIONS:\n\n");
		mexPrintf("	-Da<date> lists from date (given as mm/dd/yr/hh:mm)\n");
		mexPrintf("	-Db<date> lists up to date (given as mm/dd/yr/hh:mm)\n");
		mexPrintf("	-F Dataflags is a string made up of 1 or more of these characters:\n");
		mexPrintf("	  s means list time in seconds, sc gives dates, sj gives Julian day\n");
		mexPrintf("	  x means list longitude (degrees)\n");
		mexPrintf("	  y means list latitude (degrees)\n");
		mexPrintf("	  g means list gravity (mGal)\n");
		mexPrintf("	  m means list magnetics (nTesla)\n");
		mexPrintf("	  t means list topography (m)\n");
		mexPrintf("	  d means list distance (km)\n");
		mexPrintf("	  h means list heading (Degrees east from north)\n");
		mexPrintf("	  v means list velocity (m/s)\n");
		mexPrintf("	  w means list weight (see -W)\n");
		mexPrintf("	  If G, M, or T is used instead of g, m, or t, then only the records\n");
		mexPrintf("	  that have that combination of data will be listed\n");
		mexPrintf("	  The data is written out in the order specified in <dataflags>\n");
		mexPrintf("	  [Default is -Fsxygmtdhvw and all records]\n");
		mexPrintf("	-G force geographical longitudes (-180/+180) [Default is 0-360]\n");
		mexPrintf("	-R only return data inside the specified region\n");
		mexPrintf("	-Sa<dist> lists from dist (in km)\n");
		mexPrintf("	-Sb<dist> lists up to dist (in km)\n");
		mexPrintf("	-W sets weight for these data\n");
		mexPrintf("	-Y swapp bytes. Use this if file was created in a machine with different endianess.\n");
		return;
	}

	if ((west < 0.0 && east > 0.0) || (west < 360.0 && east > 360.0)) greenwich = TRUE;
	if (!geodetic) greenwich = TRUE;
	
	/* Sort the  order in which the parameters appear */
	
	if (nval == 0) {		/* Nothing selected, default used */
		g = m = t = TRUE;	/* No data was specified so all data [default] is output */
		for (i = 0; i < 10; i++) no[i] = i;
		nval = 10;
		want_all = tsec = TRUE;
		do_heading = do_speed = TRUE;
	}

	/* Structure field names. If they aren't used the corresponding field will be empty */
	fieldnames[0] = strdup ("time");
	fieldnames[1] = strdup ("longitude");
	fieldnames[2] = strdup ("latitude");
	fieldnames[3] = strdup ("gravity");
	fieldnames[4] = strdup ("magnetics");
	fieldnames[5] = strdup ("topography");
	fieldnames[6] = strdup ("distance");
	fieldnames[7] = strdup ("heading");
	fieldnames[8] = strdup ("velocity");
	fieldnames[9] = strdup ("weights");
	fieldnames[10] = strdup ("year");
	fieldnames[11] = strdup ("agency");
	fieldnames[12] = strdup ("info");

	nan = mxGetNaN();
	num_struct_fields = 13;
	gmtlist_struct = mxCreateStructMatrix (n_cruises, 1, num_struct_fields, (const char **)fieldnames );
	n_cruises = 0;
	 
	for (argno = 1; argno < argc; argno++) {	/* Loop over all the files */
	
		if (argv[argno][0] == '-') continue;
  		if (gmtmggpath_func (gmtfile, argv[argno])) {
   			mexPrintf ( "gmtlist : Cannot find leg %s\n", argv[argno]);
     			continue;
  		}
		if ((fp = fopen (gmtfile, "rb")) == NULL) {
			mexPrintf ("gmtlist: Could not open %s\n", gmtfile);
			continue;
		}

		/* Read first record of file containing start-year, n_records and info */
		
		if (fread ((void *)(&leg_year), (size_t)4, (size_t)1, fp) != 1) {
			mexPrintf ("gmtlist: Error while reading first year\n");
			return;
		}
		if (fread ((void *)(&n_records), (size_t)4, (size_t)1, fp) != 1) {
			mexPrintf ("gmtlist: Error while reading no of records\n");
			return;
		}
		if (fread ((void *)agency, (size_t)10, (size_t)1, fp) != 1) {
			mexPrintf ("gmtlist: Error while reading info-header\n");
			return;
		}

		if (n_records < 0 && leg_year < 0) {	/* Almost sure this happens only with big endian files */
			swapa = TRUE;
		}
		if (swapa) {
			leg_year  = GMT_swab4(leg_year);
			n_records = GMT_swab4(n_records);
		}
	
		gmt = gmtmgg_init (leg_year);	/* Initialize gmt_structure */
		if (correct)
			id = get_id(argv[argno]);
		else
			id = -1;
	
		/* Decode date to time in sec if needed */
	
		if (day1 > 0) gmtmgg_time (&start_time, year1, mon1, day1, hour1, min1, 0, gmt);
		if (day2 > 0) gmtmgg_time (&stop_time, year2, mon2, day2, hour2, min2, 0, gmt);
	
		dist = 0.0;
		time = 0;
	
		wantgmt = (g || m || t) ? TRUE : FALSE;
		if (want_all) wantgmt = FALSE;
		
		/* Alloc space for selected outputs */
		for (i = 0; i < nval; i++)
			p_out[i] = mxCalloc (n_records, sizeof (double));

		nan = mxGetNaN();

		/* Start reading data from file */
	
		xmin1 = xmin2 = 360.0;
		xmax1 = xmax2 = -360.0;
		ymin = 180.0;
		ymax = -180.0;
		n_grv = n_mag = n_top = 0;
		greenwich_in = FALSE;
		for (rec = 0; rec < n_records && dist < stop_dist && time < stop_time; rec++) {
			if (fread ((void *)(&record), (size_t)18, (size_t)1, fp) != 1) {
				mexPrintf ("gmtlist: Error reading data record no %d\n",rec);
				continue;
			}

			if (swapa) {
				record.lon = GMT_swab4 (record.lon);
				record.lat = GMT_swab4 (record.lat);
				record.time = GMT_swab4 (record.time);
				record.gmt[0] = GMT_swab2 (record.gmt[0]);
				record.gmt[1] = GMT_swab2 (record.gmt[1]);
				record.gmt[2] = GMT_swab2 (record.gmt[2]);
			}
		
			/* Compute accumulated distance along track (Flat Earth) */
		
			if (rec == 0) {
				last_lon = record.lon;
				last_lat = record.lat;
				last_time = record.time;
				ds = 0.0;
				dt = 0;
				heading = speed = GMTMGG_NODATA;
				gmtmgg_date(record.time,&begin_yr,&begin_mo,&begin_day,&hour,&minute,&second,gmt);
			}
			else {
				dlon = record.lon - last_lon;
				if (abs (dlon) > 180000000) dlon = irint(copysign ((double) (360000000 - abs (dlon)), (double)dlon));
				dx = (double) dlon * cos (0.5e-06*D2R*(double)(record.lat+last_lat));
				dy = (double) (record.lat - last_lat);
				ds = KMPRDEG * hypot (dx, dy);
				dt = record.time - last_time;
				if (do_heading) {
					heading = (dx == 0.0 && dy == 0.0) ? GMTMGG_NODATA : 90.0 - R2D * atan2 (dy, dx);
					if (heading < 0.0) heading += 360.0;
				}
				if (do_speed) speed = (dt == 0) ? GMTMGG_NODATA : 1000.0 * ds / dt;
				last_lon = record.lon;
				last_lat = record.lat;
				last_time = record.time;
			}
			dist += ds;
			
			/* Check if record has the required fields */
			if (record.gmt[0] != GMTMGG_NODATA) n_grv++;
			if (record.gmt[1] != GMTMGG_NODATA) n_mag++;
			if (record.gmt[2] != GMTMGG_NODATA) n_top++;
			
			if (no_g && record.gmt[0] == GMTMGG_NODATA) continue;
			if (no_m && record.gmt[1] == GMTMGG_NODATA) continue;
			if (no_t && record.gmt[2] == GMTMGG_NODATA) continue;
			
			/* Check if time or dist falls outside specified range */

			if (dist < start_dist) continue;
			if (dist > stop_dist) continue;
			if (record.time < start_time) continue;
			if (record.time > stop_time) continue;
		
			time = record.time;
			lat = (double) record.lat*0.000001;
			lon = (double) record.lon*0.000001;

			/* Compute limits for info purposes */
			if (lon < 180.0) {
				xmin1 = MIN (lon, xmin1);
				xmax1 = MAX (lon, xmax1);
			}
			else {
				xmin2 = MIN (lon, xmin2);
				xmax2 = MAX (lon, xmax2);
			}
			ymin = MIN (lat, ymin);
			ymax = MAX (lat, ymax);
			if (rec > 0 && (fabs(old_x-lon) > 180.0))
				greenwich_in = TRUE;
			old_x = lon;
		
			/* Check if lat/lon is outside specified area */
		
			if (lat < south || lat > north) continue;
			while (lon > east) lon -= 360.0;
			while (lon < west) lon += 360.0;
			if (lon > east) continue;
			while (lon > 360.0) lon -= 360.0;
		
			grv = (record.gmt[0] != GMTMGG_NODATA) ? (double) record.gmt[0]*0.1 : GMTMGG_NODATA;
			mag = record.gmt[1];
			top = record.gmt[2];
			
			/* This record will now be printed out */
		
			for (i = 0; i < nval; i++) {
				switch(no[i]) {
					case 0:	/* Print out time */
						p_out[i][rec] = record.time;
						break;
					case 1:	/* Print out longitude */
						if (lon > 180.0 && greenwich) lon -= 360.0;
						p_out[i][rec] = lon;
						break;
					case 2:	/* Print out latitude */
						p_out[i][rec] = lat;
						break;
					case 3:	/* Print out gravity */
						if (record.gmt[0] == GMTMGG_NODATA)
							p_out[i][rec] = nan;
						else {
							if (id >= 0) grv -= bin[id]->dc_shift_gmt[0] + bin[id]->drift_rate_gmt[0] * record.time;
							p_out[i][rec] = grv;
						}
		 				break;
					case 4:	/* Print out magnetics */
						if (record.gmt[1] == GMTMGG_NODATA)
							p_out[i][rec] = nan;
						else {
							if (id >= 0) mag -= bin[id]->dc_shift_gmt[1] + bin[id]->drift_rate_gmt[1] * record.time;
							p_out[i][rec] = mag;
						}
						break;
					case 5:	/* Print out bathymetry */
						if (record.gmt[2] == GMTMGG_NODATA)
							p_out[i][rec] = nan;
						else {
							if (id >= 0) top -= bin[id]->dc_shift_gmt[2] + bin[id]->drift_rate_gmt[2] * record.time;
							p_out[i][rec] = top;
						}
						break;
					case 6:	/* Print out distance */
						p_out[i][rec] = dist;
						break;
					case 7:	/* Print out heading */
						if (heading == GMTMGG_NODATA)
							p_out[i][rec] = nan;
						else
							p_out[i][rec] = heading;
						break;
					case 8:	/* Print out velocity */
						if (speed == GMTMGG_NODATA)
							p_out[i][rec] = nan;
						else
							p_out[i][rec] = speed;
						break;
					case 9:	/* Print out weight */
						p_out[i][rec] = weight;
						break;
				}
			}
		}
		fclose (fp);
		gmtmgg_date (record.time,&end_yr,&end_mo,&end_day,&hour,&minute,&second,gmt);
		mxFree(gmt);

		if (greenwich_in) {
			xmin = MAX(xmin1,xmin2);
			xmax = MIN(xmax1,xmax2);
		}
		else {
			xmin = MIN(xmin1,xmin2);
			xmax = MAX(xmax1,xmax2);
		}
		if (xmin > xmax) xmin -= 360.0;
		
		/* Fill struct fields and free allocated space for selected outputs */

		for (i = 0; i < nval; i++) {
			mxOut = mxCreateDoubleMatrix (1,n_records, mxREAL);
			pdata = mxGetPr(mxOut);
			memcpy(pdata, p_out[i], n_records * sizeof(double));
			mxSetField(gmtlist_struct, n_cruises, fieldnames[no[i]], mxOut);
			mxFree(p_out[i]);
		}
		mxOut = mxCreateDoubleMatrix (1,1, mxREAL);
		pdata = mxGetPr(mxOut);
		tmp = (double)leg_year;
		memcpy(pdata, &tmp, sizeof(double));
		mxSetField(gmtlist_struct, n_cruises, fieldnames[10], mxOut);
		/* I'm really pissed off with this so I'll use all the letters.
		   Why the f... this crashs Matlab? (not because leg_year is a int)
			mxSetField(gmtlist_struct, n_cruises, fieldnames[10], &leg_year); */

		mxStr = mxCreateString(agency);
		mxSetField(gmtlist_struct, n_cruises, fieldnames[11], mxStr);

		mxOut = mxCreateDoubleMatrix (1,14, mxREAL);
		pdata = mxGetPr(mxOut);
		info[0] = (double)n_records;	info[1] = (double)n_grv;	info[2] = (double)n_mag;
		info[3] = (double)n_top;	info[4] = xmin;			info[5] = xmax;
		info[6] = ymin;			info[7] = ymax;
		info[8] = (double)begin_day;	info[9] = (double)begin_mo;	info[10] = (double)begin_yr;
		info[11] = (double)end_day;	info[12] = (double)end_mo;	info[13] = (double)end_yr;
		memcpy(pdata, info, 14 * sizeof(double));
		mxSetField(gmtlist_struct, n_cruises, fieldnames[12], mxOut);

		n_cruises++;

	}
	plhs[0] = gmtlist_struct;
	
}

int get_id (char *name) {
	int left, right, mid, cmp;
	
	left = 0;
	right = nlegs-1;
	while (left <= right) {
		mid = (left + right)/2;
		cmp = strcmp(name, bin[mid]->name);
		if (cmp < 0)
			right = mid-1;
		else if (cmp > 0)
			left = mid+1;
		else
			return (mid);
	}
	return (-1);
}


int gmtmggpath_func (char *leg_path, char *leg) {
	int id;
	char geo_path[BUFSIZ];

	sprintf (geo_path, "%s.gmt", leg);
	if (!access(geo_path, R_OK)) {
		strcpy(leg_path, geo_path);
		return (0);
	}
	return(1);
}

/*
 *	GMT subroutine gmtmgg_init sets up the structure GMT  which is
 *	used by other gmt routines (gmtmgg_time,gmtmgg_date) to convert
 *	times.  Daymon[month] contains the cumulative number of
 *	days from Jan 1 in first_year through the months PRIOR to the
 *	value of month.  0 <= month <= 60.  month = 0 only occurs
 *	during initializing in this routine. The user must declare
 *	a pointer to the struct GMTMGG_TIME in the main program and pass it
 *	when calling the gmt_* functions. To define the GMT structure,
 *	include the file gmt.h
 *
 *	Paul Wessel
 *	12-JUL-1987
 *
 */

struct GMTMGG_TIME *gmtmgg_init (int year1) {
	struct GMTMGG_TIME *gmt_struct = NULL;
	int dm[12];	/* No of days in each month */
	int year, this_year, month, m;
	gmt_struct = (struct GMTMGG_TIME *) mxCalloc((size_t)1, sizeof(struct GMTMGG_TIME));
	gmt_struct->first_year = year1;
	/* initialize days of the month etc. */
	dm[0] = 0;
	dm[1] = 31;
	dm[2] = 28;
	dm[3] = 31;
	dm[4] = 30;
	dm[5] = 31;
	dm[6] = 30;
	dm[7] = 31;
	dm[8] = 31;
	dm[9] = 30;
	dm[10] = 31;
	dm[11] = 30;
	gmt_struct->daymon[0] = 0;
	for (year = 0, month = 0; year < 5; year++) {
		this_year = gmt_struct->first_year + year;
		if (this_year%4 == 0 && !(this_year%400 == 0)) dm[2] = 29;
		for (m = 1; m <= 12; m++) {
			month++;
			gmt_struct->daymon[month] = gmt_struct->daymon[month - 1] + dm[m - 1];
		}
   		dm[2] = 28;
   		dm[0] = 31;
   	}

  	return (gmt_struct);
}

/* GMT function gmtmgg_time returns the number of seconds from
 * first_year calculated from (hr/mi/sc/dd/mm/yy). The pointer
 * to the GMT structure is passed along with the arguments.
 */
 /* MODIFIED 10 July, 1987 by W. Smith  --  I killed a bug in month calculation */

int gmtmgg_time (int *time, int year, int month, int day, int hour, int minute, int second, struct GMTMGG_TIME *gmt_struct) {
	int mon, n_days, bad = 0;
	if ((mon = (year - gmt_struct->first_year)) > 4) {
		mexPrintf ("gmtmgg_time:  Year - first_year > 4\n");
		return(-1);
	}
	if (month < 1 || month > 12) mexPrintf ("GMT WARNING: in gmtmgg_time: Month out of range [1-12]: %d\n", month), bad++;
	if (day < 1 || day > 31) mexPrintf("GMT WARNING: in gmtmgg_time: Day out of range [1-31]: %d\n", day), bad++;
	if (hour < 0 || hour > 24) mexPrintf("GMT WARNING: in gmtmgg_time: Hour out of range [0-24]: %d\n", hour), bad++;
	if (minute < 0 || minute > 60) mexPrintf("GMT WARNING: in gmtmgg_time: Minute out of range [0-60]: %d\n", minute), bad++;
	if (second < 0 || second > 60) mexPrintf("GMT WARNING: in gmtmgg_time: Second out of range [0-60]: %d\n", second), bad++;
	if (bad) return (-1);	/* When we got garbage input */
	mon = mon * 12 + month;
	n_days = gmt_struct->daymon[mon] + day - 1;
	*time = n_days * 86400 + hour * 3600 + minute * 60 + second;
	return (*time);
}


/* GMT function gmtmgg_date computes the date (hr/mi/sec/dd/mm/yy) based
 * on the total time in seconds since the beginning of first_year.
 * The pointer to the GMT structure is passed allong with the other
 * arguments. The Julian day is returned. the yymmddhhmmss is passed
 * through the argument list.
 */

int gmtmgg_date (int time, int *year, int *month, int *day, int *hour, int *minute, int *second, struct GMTMGG_TIME *gmt_struct) {
	int day_time, julian_day;
	day_time = time/86400;
	*month = day_time / 31 + 1;	/* Only approximately, may be smaller */

	if ((*month) < 0 || (*month) >= GMTMGG_TIME_MAXMONTH) {
		mexPrintf ("GMT ERROR: in gmtmgg_date: Month outside valid range [0-%d>: %d\n", GMTMGG_TIME_MAXMONTH, *month);
		return (EXIT_FAILURE);
	}
	while (gmt_struct->daymon[*month +1] <= day_time) {
		(*month)++;
		if ((*month) < 0 || (*month) > GMTMGG_TIME_MAXMONTH) {
			mexPrintf ("GMT ERROR: in gmtmgg_date: Month outside valid range [0-%d>: %d\n", GMTMGG_TIME_MAXMONTH, *month);
			return (EXIT_FAILURE);
		}
	}
	*year = (*month  - 1) / 12 + gmt_struct->first_year;
	*day = day_time - gmt_struct->daymon[*month] + 1;
	julian_day = (*month > 12) ?
		gmt_struct->daymon[*month] - gmt_struct->daymon[(*month - (*month)%12)] + *day :
		gmt_struct->daymon[*month] + *day;
	*month  = (*month-1)%12 + 1;
	time %= 86400;
	*hour = time / 3600;
	*minute = (time%3600) / 60;
	*second = time - *hour * 3600 - *minute * 60;
	return (julian_day);
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

