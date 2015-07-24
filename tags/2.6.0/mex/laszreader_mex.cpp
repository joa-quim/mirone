/*--------------------------------------------------------------------
 *	$Id$
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
 * Purpose: Read LASZ/LAS LIDAR data into Matlab (Needs to link against LASlib_i)
 *
 * Author:	Joaquim Luis
 * Date:	14-Nov-2012
 *      Contact info: w3.ualg.pt/~jluis
 *--------------------------------------------------------------------*/

#include "mex.h"
#include "mwsize.h"
#include <time.h>
#include <math.h>
#include "lasreader.hpp"
#include "laswaveform13reader.hpp"

void print_header(LASheader *header, int bSkipVLR);
void plain_xyz(mxArray *plhs[], LASreader *reader, unsigned int nPoints);
int  conditional_xyz ( mxArray *plhs[], LASreader *reader, LASheader *header,
		double ang, int classif, int intens, int nRet, int srcID, int got_R,
		double west, double east, double south, double north, double z_min, double z_max);
void get_classification_list ( mxArray *plhs[],  LASreader *reader);
void get_ID_list ( mxArray *plhs[],  LASreader *reader);
double ddmmss_to_degree (char *text);
int decode_R (char *item, double *w, double *e, double *s, double *n, double *z_min, double *z_max);
int check_region (double w, double e, double s, double n);

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	verbose = FALSE, got_R = FALSE, scanC = FALSE, scanD = FALSE, get_BB_only = FALSE;
	int	i, argc = 0, n_arg_no_char = 0, classif = 0, intens = 0, nRet = 0, srcID = 0;
	unsigned int nPoints;
	char	**argv, *fname = NULL, *parse_string = "xyz";
	double	*bbox, west, east, south, north, z_min = -1001, z_max, angle = 0;
 
	if (nrhs == 0) {
		mexPrintf ("usage: [xyz, bbox] = lasreader_mex ('filename', ['-A<ang>'], ['-C<class>'], ['-D<id>']\n");
		mexPrintf ("       [-I<intens>]', ['-N<return>'], ['-R<x_min/x_max/y_min/y_max[/z_min/z_max]>']);\n");
		mexPrintf ("  OR\n");
		mexPrintf ("       [class, bbox] = lasreader_mex ('filename', '-S'),\n\n");
		mexPrintf ("  OR\n");
		mexPrintf ("       [IDs, bbox] = lasreader_mex ('filename', '-D'),\n\n");
		mexPrintf ("  OR\n");
		mexPrintf ("       bbox = lasreader_mex ('filename', '-B'),\n\n");

		mexPrintf ("First case usage:\n");
		mexPrintf ("where xyz is 3xN array with the XYZ point coordinates:\n");
		mexPrintf ("and bbox is a 1x6 vector with [xmin xmax ymin ymax zmin zmax]\n");
		mexPrintf ("  -A<ang> Clip out points with scan angle > ang\n");
		mexPrintf ("  -C<class> Retain only points with classification = class\n");
		mexPrintf ("  -D<id> Retain only points with Source IDs = id (DO NOT CONFUSE WITH -D)\n");
		mexPrintf ("  -D Scan file for a list of Source IDs (see below) (DO NOT CONFUSE WITH -D<id>).\n");
		mexPrintf ("  -I<intens> Clip out points with intensity < intens\n");
		mexPrintf ("  -N<return> Select first return (-N1) or last return (-N10)\n");
		mexPrintf ("  -R<x_min/x_max/y_min/y_max> - Clip to bounding box.\n");
		mexPrintf ("    Optionaly add z_min/z_max to make a 3D bounding box.\n\n");
		mexPrintf ("  -S Scan file for a list of Classifications (see below).\n");
		mexPrintf ("  -V Prints header contents info on ML shell.\n");

		mexPrintf ("Second/Third cases:\n");
		mexPrintf ("class|ID is a 1xN vector with a list of different classificatins|IDs\n");
		mexPrintf ("present in file. No other data (except optional BBox) is return here.\n");

		mexPrintf ("Fourth case:\n");
		mexPrintf ("\tReturns only the bonding box 1x6 vector.\n");
		return;
	}

	if (!mxIsChar(prhs[0])) mexErrMsgTxt ("First arg must contain the filename string.\n");
	if (nlhs != 1 && get_BB_only) mexErrMsgTxt ("-B option implies one output only.\n"); 
	fname = (char *) mxArrayToString (prhs[0]);	/* Load the file name into a char string */
 
	if (nrhs > 1) {
		argc = nrhs;
		for (i = 1; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
			if(!mxIsChar(prhs[i])) {
				argc--;
				n_arg_no_char++;	/* Number of arguments that have a type other than char */
			}
		}
		argc++;			/* to account for the program's name to be inserted in argv[0] */

		/* get the length of the input string */
		argv = (char **)mxCalloc(argc, sizeof(char *));
		argv[0] = "LASreader";
		for (i = 1; i < argc; i++)
			argv[i] = (char *)mxArrayToString(prhs[i+n_arg_no_char-1]);
	}

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {

				case 'A':
					angle = atof(&argv[i][2]);
					break;
				case 'B':
					get_BB_only = TRUE;
					break;
				case 'C':
					classif = atoi(&argv[i][2]);
					break;
				case 'D':
					if (argv[i][2])
						srcID = atoi(&argv[i][2]);
					else
						scanD = TRUE;
					break;
				case 'I':
					intens = atoi(&argv[i][2]);
					break;
				case 'N':
					nRet = atoi(&argv[i][2]);
					break;
				case 'R':
					if (decode_R (argv[i], &west, &east, &south, &north, &z_min, &z_max))
						mexErrMsgTxt("Error decoding -R option!");
					got_R = TRUE;
					break;
				case 'S':
					scanC = TRUE;
					break;
				case 'V':
					verbose = TRUE;
					break;
			}
		}
	}

	LASreadOpener lasreadopener;
	lasreadopener.set_merged(FALSE);
	lasreadopener.set_populate_header(FALSE);
	lasreadopener.set_file_name(fname);

	LASreader *lasreader = lasreadopener.open();
	if (!lasreader) mexErrMsgTxt("LASREADER Error! could not open lasreader!");

	LASheader *header = &(lasreader->header);
	if (!header) mexErrMsgTxt("LASREADER: Unable to fetch header for file");

	if (get_BB_only && (scanC || scanD) )
		mexPrintf("LASREADER WARNING: option -B takes precedence over -C or -D\n");

	else if (scanC && scanD)
		mexPrintf("LASREADER WARNING: option -C takes precedence over -D\n");

	if (get_BB_only) {
		plhs[0] = mxCreateDoubleMatrix(1, 6, mxREAL);
		bbox = mxGetPr(plhs[0]);
		bbox[0] = header->min_x;	bbox[1] = header->max_x;
		bbox[2] = header->min_y;	bbox[3] = header->max_y;
		bbox[4] = header->min_z;	bbox[5] = header->max_z;
		return;
	}

	if (verbose) print_header(header, FALSE);

	if (!(scanC || scanD)) {
		if ((got_R + nRet + intens + classif + angle + srcID) == 0)
			plain_xyz(plhs, lasreader, header->number_of_point_records);
		else
			conditional_xyz(plhs, lasreader, header, angle, classif,
					intens, nRet, srcID, got_R, west, east, south, north, z_min, z_max);
	}
	else if (scanC)		/* Scan file for a list of different classifications */
		get_classification_list ( plhs, lasreader);

	else if (scanD)		/* Scan file for a list of different Source IDs */
		get_ID_list ( plhs, lasreader);

	if (nlhs == 2) {
		plhs[1] = mxCreateDoubleMatrix(1, 6, mxREAL);
		bbox = mxGetPr(plhs[1]);
		bbox[0] = header->min_x;	bbox[1] = header->max_x;
		bbox[2] = header->min_y;	bbox[3] = header->max_y;
		bbox[4] = header->min_z;	bbox[5] = header->max_z;
	}
    
    // close the reader
	lasreader->close();
	delete lasreader;

	return;
}

/* ---------------------------------------------------------------------------------- */
void plain_xyz(mxArray *plhs[], LASreader *reader, unsigned int nPoints) {
	/* Get the XYZ triplets without any condition other than belonging to a valid point */
	unsigned int n = 0;
	double	*ptr;
    
	plhs[0] = mxCreateDoubleMatrix(3, nPoints, mxREAL);
	ptr = mxGetPr(plhs[0]);
    
	while (reader->read_point()) {        
		ptr[n*3]   = reader->point.get_x();
		ptr[n*3+1] = reader->point.get_y();
		ptr[n*3+2] = reader->point.get_z();
		n++;
	}
}

/* ---------------------------------------------------------------------------------- */
int  conditional_xyz ( mxArray *plhs[], LASreader *reader, LASheader *header,
		double ang, int classif, int intens, int nRet, int srcID, int got_R,
		double west, double east, double south, double north, double z_min, double z_max) {
	/* Get the XYZ triplets that escape excluding clausules */

	unsigned int n = 0, n_alloc = 100000;
	int	got_Z = FALSE, first_only, last_only;
	double	*ptr, *tmp, x, y, z, MinZ, MaxZ;
 
	last_only  = (nRet > 9  ) ? TRUE : FALSE;
	first_only = (nRet == 1 ) ? TRUE : FALSE;
	if (z_min > -1000) {
		got_Z = TRUE;
		MinZ = header->min_z;
		MaxZ = header->max_z;
	}

	ptr = (double *)mxMalloc((mwSize)(n_alloc * 3 * sizeof(double)));
 
	while (reader->read_point()) {        
        
		if (classif && reader->point.classification != classif) goto fim;
		if (srcID   && reader->point.point_source_ID != srcID) goto fim;
		if (intens  && reader->point.intensity < intens) goto fim;
		if (last_only  && reader->point.return_number != reader->point.number_of_returns_of_given_pulse) goto fim;
		if (first_only && reader->point.return_number != 1) goto fim;
		if (ang && (reader->point.scan_angle_rank > ang || reader->point.scan_angle_rank < -ang)) goto fim;

		x = reader->point.get_x();		y = reader->point.get_y();
		if (got_R && (x < west || x > east || y < south || y > north)) goto fim; 
		z = reader->point.get_z();
		if (got_Z && (z < MinZ || z > MaxZ)) goto fim; 

		ptr[n*3] = x;
		ptr[n*3+1] = y;
		ptr[n*3+2] = z;
		n++;

		if (n >= n_alloc) {
			n_alloc += 50000;
			ptr = (double *)mxRealloc((void *)ptr, (mwSize)(n_alloc * 3 * sizeof(double)));
		}
fim:
	;
	}

	if (n) {
		n--;
			
		ptr = (double *)mxRealloc((void *)ptr, (mwSize)(n * 3 * sizeof(double)));

		plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL);
		tmp = mxGetPr(plhs[0]);
		mxFree((void *)tmp);
		mxSetPr(plhs[0], ptr);
		mxSetN(plhs[0], (mwSize)n);
	}
	else
		plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);	/* We got no points */

 	return (n);
}

/* ---------------------------------------------------------------------------------- */
void  get_classification_list ( mxArray *plhs[],  LASreader *reader) {
	/* Scan the LAS file for classification and create a vector list of the different ones */
	double	this_classif, last_classif = 0;	/* Use doubles because they are few and makes life easier in ML */
	double	*class_list = NULL;	
	int	nClass = 0, onList, i;

	class_list = (double *)mxMalloc(32 * sizeof(double));
	class_list[0] = 1;

	while (reader->read_point()) {
		this_classif = reader->point.classification;
		if ( this_classif != last_classif) {
			onList = FALSE;
			for (i = 0; i < nClass; i++)	/* Check if we have a new classification index */
				if (this_classif == class_list[i]) onList = TRUE;

			if (!onList) {		/* Got a new classification. Store it */
				class_list[++nClass] = this_classif;
				last_classif = this_classif;
			}
		}
	}

	plhs[0] = mxCreateDoubleMatrix(1, (nClass + 1), mxREAL);
	memcpy(mxGetPr(plhs[0]), class_list, (nClass + 1) * sizeof(double));
	mxFree((void *) class_list);
}

/* ---------------------------------------------------------------------------------- */
void  get_ID_list ( mxArray *plhs[],  LASreader *reader) {
	/* Scan the LAS file for Source ID and create a vector list of the different ones */
	double	this_ID, last_ID = -1;	/* Use doubles because they are few and makes life easier in ML */
	double	*ID_list = NULL;	
	int	nID = 0, onList, i;

	ID_list = (double *)mxMalloc(4096 * sizeof(double));	/* 4096 should be enough no? */
	ID_list[0] = 0;

	while (reader->read_point()) {
		this_ID = reader->point.point_source_ID;
		if ( this_ID != last_ID) { 
			onList = FALSE;
			for (i = 0; i < nID; i++)	/* Check if we have a new classification index */
				if (this_ID == ID_list[i]) onList = TRUE;

			if (!onList) {		/* Got a new classification. Store it */
				ID_list[nID++] = this_ID;
				last_ID = this_ID;
			}
		}
	}

	plhs[0] = mxCreateDoubleMatrix(1, nID, mxREAL);
	memcpy(mxGetPr(plhs[0]), ID_list, nID * sizeof(double));
	mxFree((void *) ID_list);
}

/* ---------------------------------------------------------------------------------- */
int decode_R (char *item, double *w, double *e, double *s, double *n, double *z_min, double *z_max) {
	char *text, string[BUFSIZ];
	
	/* Minimalist code to decode option -R extracted from GMT_get_common_args */
	
	int i, error = 0;
	double *p[6];
	
	p[0] = w;	p[1] = e;	p[2] = s;	p[3] = n;
			
	i = 0;
	strcpy (string, &item[2]);
	text = strtok (string, "/");
	while (text) {
		*p[i] = ddmmss_to_degree (text);
		i++;
		text = strtok (NULL, "/");
	}
	if (item[strlen(item)-1] == 'r')	/* Rectangular box given, but invalid here */
		error++;
	if ((i != 4 && i != 6) || check_region (*p[0], *p[1], *p[2], *p[3]))
		error++;
	w = p[0];	e = p[1];
	s = p[2];	n = p[3];
	if (i == 6) {
		z_min = p[4];		z_max = p[5];
	}
	return (error);
}

/* -------------------------------------------------------------------- */
int check_region (double w, double e, double s, double n) {
	/* If region is given then we must have w < e and s < n */
	return ((w >= e || s >= n));
}

#define Loc_copysign(x,y) ((y) < 0.0 ? -fabs(x) : fabs(x))
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
	return (degfrac);
}

/* ---------------------------------------------------------------------------------- */
void print_header(LASheader *header, int bSkipVLR) {
    
    int i = 0;

    mexPrintf("\n---------------------------------------------------------\n");
    mexPrintf("  Header Summary\n");
    mexPrintf("---------------------------------------------------------\n");

    if (strcmp(header->file_signature, "LASF") != 0) {
        mexPrintf("File signature is not 'LASF'... aborting");
        return;
    }
    mexPrintf("  Version:                   %d.%d\n", 
                 header->version_major, header->version_minor);

    mexPrintf("  Source ID:                 %d\n", header->project_ID_GUID_data_1) ;	// AHHHHHHH
    mexPrintf("  Reserved:                  %d\n", header->global_encoding );
    mexPrintf("  Project ID/GUID:          '%s'\n", header->file_source_id);
    mexPrintf("  System Identifier:        '%s'\n", header->system_identifier);
    mexPrintf("  Generating Software:      '%s'\n", header->generating_software);
    mexPrintf("  File Creation Day/Year:    %d/%d\n", 
                 header->file_creation_day, header->file_creation_year);

    mexPrintf("  Header Size                %d\n", header->header_size);
    mexPrintf("  Offset to Point Data       %d\n", header->offset_to_point_data);
    mexPrintf("  Number Var. Length Records %d\n", header->number_of_variable_length_records);
    mexPrintf("  Point Data Format          %d\n", header->point_data_format);
    mexPrintf("  Point Data Record Length   %d\n", header->point_data_record_length);
    mexPrintf("  Number of Point Records    %d\n", header->number_of_point_records);
    mexPrintf("  Number of Points by Return %u %u %u %u %u\n", 
                 header->number_of_points_by_return[0], 
                 header->number_of_points_by_return[1], 
                 header->number_of_points_by_return[2], 
                 header->number_of_points_by_return[3], 
                 header->number_of_points_by_return[4]);

    mexPrintf("  Scale Factor X Y Z         %.12g %.12g %.12g\n", 
                 header->x_scale_factor, 
                 header->y_scale_factor,
                 header->z_scale_factor);

    mexPrintf("  Offset X Y Z               %.12g %.12g %.12g\n", 
                 header->x_offset, 
                 header->y_offset, 
                 header->z_offset);

    mexPrintf("  Min X Y Z                  %.12g %.12g %.12g\n",
                 header->min_x, 
                 header->min_y, 
                 header->min_z);

    mexPrintf("  Max X Y Z                  %.12g %.12g %.12g\n", 
                 header->max_x, 
                 header->max_y, 
                 header->max_z);    
}
