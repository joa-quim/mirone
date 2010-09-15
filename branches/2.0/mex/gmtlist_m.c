/*--------------------------------------------------------------------
 *	$Id: gmtlist.c,v 1.7 2005/03/04 21:00:54 remko Exp $
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
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		14/10/06 J Luis, Now includes the memory leak solving solution
 */
 
#include "mex.h"
#include "gmt.h"
#include "gmt_mgg.h"
#include "x_system.h"

#define KMPRDEG 111.1949e-6
#define MAXLEGS 5000
#define S_PR_DAY 86400

struct CORR **bin;

int binsize = sizeof(struct CORR);
int nlegs = 0;

int get_id (char *name);

/* int GMTisLoaded = FALSE;	/* Used to know wether GMT stuff is already in memory or not */

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
	
	/*if (!GMTisLoaded) {
		argc = GMT_begin (argc, argv);
		GMTisLoaded = TRUE;
	}
	else
		argc = GMT_short_begin (argc, argv);*/
	argc = GMT_begin (argc, argv);

	for (i =1; !error && i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
			
				case 'R':
				case '\0':
					error += GMT_get_common_args (argv[i], &west, &east, &south, &north);
					break;
					
				case 'D':		/* Assign start/stop times for sub-section */
					if (argv[i][2] == 'a') {	/* Start date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d",
							&mon1, &day1, &year1, &hour1, &min1);
					}
					else if (argv[i][2] == 'b')	 {	/* Stop date */
						sscanf(&argv[i][3], "%d/%d/%d/%d:%d",
							&mon2, &day2, &year2, &hour2, &min2);
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

	gmtmggpath_init(GMT_SHAREDIR);

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
		GMT_free ((void *)gmt);

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
	/*GMT_end_for_mex (argc, argv); */
	GMT_end (argc, argv);
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
