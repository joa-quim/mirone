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
 * Purpose: Read LAS LIDAR data into Matlab (Needs to link against libLAS)
 *
 * Author:	Joaquim Luis
 * Date:	04-sep-2010
 *      Contact info: w3.ualg.pt/~jluis
 *--------------------------------------------------------------------*/

#include "mex.h"
#include "mwsize.h"
#include "liblas.h"

void print_header(LASHeaderH header, int bSkipVLR);
void plain_xyz(mxArray *plhs[], LASReaderH reader, LASHeaderH header);
void get_classification_list ( mxArray *plhs[],  LASReaderH reader);
void  get_ID_list ( mxArray *plhs[],  LASReaderH reader);
double ddmmss_to_degree (char *text);


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int	verbose = FALSE, got_R = FALSE, scanC = FALSE, scanD = FALSE, get_BB_only = FALSE;
	int	i, argc = 0, n_arg_no_char = 0, classif = 0, intens = 0, nRet = 0, srcID = 0;
	unsigned int nPoints;
	char	**argv, *fname = NULL, *parse_string = "xyz";
	double	*bbox, west, east, south, north, z_min = -1001, z_max, angle = 0;
	LASReaderH reader = NULL;
	LASHeaderH header = NULL;
	LASPointH p = NULL;
 
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


	reader = LASReader_Create(fname);
	if (!reader) mexErrMsgTxt("LASREADER Error! Unable to read file!");

	header = LASReader_GetHeader(reader);
	if (!header) mexErrMsgTxt("LASREADER: Unable to fetch header for file");

	if (get_BB_only && (scanC || scanD) )
		mexPrintf("LASREADER WARNING: option -B takes precedence over -C or -D\n");

	else if (scanC && scanD)
		mexPrintf("LASREADER WARNING: option -C takes precedence over -D\n");

	if (get_BB_only) {
		plhs[0] = mxCreateDoubleMatrix(1, 6, mxREAL);
		bbox = mxGetPr(plhs[0]);
		bbox[0] = LASHeader_GetMinX(header);	bbox[1] = LASHeader_GetMaxX(header);
		bbox[2] = LASHeader_GetMinY(header);	bbox[3] = LASHeader_GetMaxY(header);
		bbox[4] = LASHeader_GetMinZ(header);	bbox[5] = LASHeader_GetMaxZ(header);
		return;
	}

	if (verbose) print_header(header, FALSE);

	if (!(scanC || scanD)) {
		if ((got_R + nRet + intens + classif + angle + srcID) == 0)
			plain_xyz(plhs, reader, header);
		else
			conditional_xyz(plhs, reader, header, angle, classif, intens, nRet,
					srcID, got_R, west, east, south, north, z_min, z_max);
	}
	else if (scanC)		/* Scan file for a list of different classifications */
		get_classification_list ( plhs, reader);

	else if (scanD)		/* Scan file for a list of different Source IDs */
		get_ID_list ( plhs, reader);

	if (nlhs == 2) {
		plhs[1] = mxCreateDoubleMatrix(1, 6, mxREAL);
		bbox = mxGetPr(plhs[1]);
		bbox[0] = LASHeader_GetMinX(header);	bbox[1] = LASHeader_GetMaxX(header);
		bbox[2] = LASHeader_GetMinY(header);	bbox[3] = LASHeader_GetMaxY(header);
		bbox[4] = LASHeader_GetMinZ(header);	bbox[5] = LASHeader_GetMaxZ(header);
	}
    
	LASReader_Destroy(reader);
	LASHeader_Destroy(header);

	return;
}

/* ---------------------------------------------------------------------------------- */
void plain_xyz(mxArray *plhs[], LASReaderH reader, LASHeaderH header) {
	/* Get the XYZ triplets without any condition other than belonging to a valid point */
	unsigned int n = 0, nPoints;
	double	*ptr;
	LASPointH p = NULL;
    
	nPoints = LASHeader_GetPointRecordsCount(header); 
	plhs[0] = mxCreateDoubleMatrix(3, nPoints, mxREAL);
	ptr = mxGetPr(plhs[0]);

	p = LASReader_GetNextPoint(reader);
    
	while (p) {
		if (!LASPoint_IsValid(p)) {
			p = LASReader_GetNextPoint(reader);
			continue;
		}
        
		ptr[n*3] = LASPoint_GetX(p);
		ptr[n*3+1] = LASPoint_GetY(p);
		ptr[n*3+2] = LASPoint_GetZ(p);
		n++;
		p = LASReader_GetNextPoint(reader);
	}
}

/* ---------------------------------------------------------------------------------- */
int  conditional_xyz ( mxArray *plhs[], LASReaderH reader, LASHeaderH header,
		double ang, int classif, int intens, int nRet, int srcID, int got_R,
		double west, double east, double south, double north, double z_min, double z_max) {
	/* Get the XYZ triplets that escape excluding clausules */

	unsigned int n = 0, n_alloc = 100000, nPoints;
	int	got_Z = FALSE, first_only, last_only;
	double	*ptr, *tmp, x, y, z, MinZ, MaxZ;
	LASPointH p = NULL;
 
	last_only  = (nRet > 9  ) ? TRUE : FALSE;
	first_only = (nRet == 1 ) ? TRUE : FALSE;
	if (z_min > -1000) {
		got_Z = TRUE;
		MinZ = LASHeader_GetMinZ(header);
		MaxZ = LASHeader_GetMaxZ(header);
	}

	ptr = (double *)mxMalloc((mwSize)(n_alloc * 3 * sizeof(double)));

	nPoints = LASHeader_GetPointRecordsCount(header); 
	p = LASReader_GetNextPoint(reader);
 
	while (p) {
		if (!LASPoint_IsValid(p)) {
			p = LASReader_GetNextPoint(reader);
			continue;
		}
        
		if (classif && LASPoint_GetClassification(p) != classif) goto fim;
		if (srcID && LASPoint_GetPointSourceId(p) != srcID) goto fim;
		if (intens && LASPoint_GetIntensity(p) < intens) goto fim;
		if (last_only && LASPoint_GetReturnNumber(p) != LASPoint_GetNumberOfReturns(p)) goto fim;
		if (first_only && LASPoint_GetReturnNumber(p) != 1) goto fim;
		if (ang && (LASPoint_GetScanAngleRank(p) > ang || LASPoint_GetScanAngleRank(p) < -ang)) goto fim;

		x = LASPoint_GetX(p);		y = LASPoint_GetY(p);
		if (got_R && (x < west || x > east || y < south || y > north)) goto fim; 
		z = LASPoint_GetZ(p);
		if (got_Z && (z < MinZ || z > MaxZ)) goto fim; 

		ptr[n*3] = x;
		ptr[n*3+1] = y;
		ptr[n*3+2] = z;
		n++;

		if (n >= n_alloc) {
			n_alloc += 50000;
			ptr = mxRealloc((void *)ptr, (mwSize)(n_alloc * 3 * sizeof(double)));
		}
fim:
		p = LASReader_GetNextPoint(reader);
	}

	if (n) {
		n--;
			
		ptr = mxRealloc((void *)ptr, (mwSize)(n * 3 * sizeof(double)));

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
void  get_classification_list ( mxArray *plhs[],  LASReaderH reader) {
	/* Scan the LAS file for classification and create a vector list of the different ones */
	double	this_classif, last_classif = 0;	/* Use doubles because they are few and makes life easier in ML */
	double	*class_list = NULL;	
	int	nClass = 0, onList, i;
	LASPointH p = LASReader_GetNextPoint(reader);

	class_list = (double *)mxMalloc(32 * sizeof(double));
	class_list[0] = 1;

	while (p) {
		if (!LASPoint_IsValid(p)) {
			p = LASReader_GetNextPoint(reader);
			continue;
		}

		this_classif = LASPoint_GetClassification(p);
		if ( this_classif != last_classif) { 
			onList = FALSE;
			for (i = 0; i < nClass; i++) {	/* Check if we have a new classification index */
				if (this_classif == class_list[i]) onList = TRUE;
			}
			if (!onList) {		/* Got a new classification. Store it */
				class_list[++nClass] = this_classif;
				last_classif = this_classif;
			}
		}

		p = LASReader_GetNextPoint(reader);
	}

	plhs[0] = mxCreateDoubleMatrix(1, (nClass + 1), mxREAL);
	memcpy(mxGetPr(plhs[0]), class_list, (nClass + 1) * sizeof(double));
	mxFree((void *) class_list);
}

/* ---------------------------------------------------------------------------------- */
void  get_ID_list ( mxArray *plhs[],  LASReaderH reader) {
	/* Scan the LAS file for Source ID and create a vector list of the different ones */
	double	this_ID, last_ID = -1;	/* Use doubles because they are few and makes life easier in ML */
	double	*ID_list = NULL;	
	int	nID = 0, onList, i;
	LASPointH p = LASReader_GetNextPoint(reader);

	ID_list = (double *)mxMalloc(4096 * sizeof(double));	/* 4096 should be enough no? */
	ID_list[0] = 0;

	while (p) {
		if (!LASPoint_IsValid(p)) {
			p = LASReader_GetNextPoint(reader);
			continue;
		}

		this_ID = LASPoint_GetPointSourceId(p);
		if ( this_ID != last_ID) { 
			onList = FALSE;
			for (i = 0; i < nID; i++) {	/* Check if we have a new classification index */
				if (this_ID == ID_list[i]) onList = TRUE;
			}
			if (!onList) {		/* Got a new classification. Store it */
				ID_list[nID++] = this_ID;
				last_ID = this_ID;
			}
		}

		p = LASReader_GetNextPoint(reader);
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
void print_header(LASHeaderH header, int bSkipVLR) {

    char *pszSignature = NULL;
    char *pszProjectId = NULL;
    char *pszSystemId = NULL;
    char *pszSoftwareId = NULL;
    /*char *pszProj4 = NULL;*/
    
    char *pszVLRUser = NULL;
    char *pszVLRDescription = NULL;
    uint16_t nVLRLength = 0;
    uint16_t nVLRRecordId = 0;
    
    LASVLRH pVLR = NULL;
    LASSRSH pSRS = NULL;
    uint32_t nVLR = 0;
    int i = 0;

#ifdef HAVE_GEOTIFF
    const GTIF* pGTIF = NULL;
#else
    const void* pGTIF = NULL;
#endif    

    pszSignature = LASHeader_GetFileSignature(header);
    pszProjectId = LASHeader_GetProjectId(header);
    pszSystemId = LASHeader_GetSystemId(header);
    pszSoftwareId = LASHeader_GetSoftwareId(header);
    
    pSRS = LASHeader_GetSRS(header);
    /*pszProj4 = LASSRS_GetProj4(pSRS);*/
    pGTIF = LASSRS_GetGTIF(pSRS);
    
    nVLR = LASHeader_GetRecordsCount(header);
 
    mexPrintf("\n---------------------------------------------------------\n");
    mexPrintf("  Header Summary\n");
    mexPrintf("---------------------------------------------------------\n");

    if (strcmp(pszSignature,"LASF") !=0) {
        LASError_Print("File signature is not 'LASF'... aborting");
        return;
    }
    mexPrintf("  Version:                    %d.%d\n", 
                    LASHeader_GetVersionMajor(header), LASHeader_GetVersionMinor(header));

    mexPrintf("  Source ID:                  %d\n", LASHeader_GetFileSourceId(header) ) ;
    mexPrintf("  Reserved:                   %d\n", LASHeader_GetReserved(header) );
    mexPrintf("  Project ID/GUID:           '%s'\n", pszProjectId);
    mexPrintf("  System Identifier:         '%s'\n", pszSystemId);
    mexPrintf("  Generating Software:       '%s'\n", pszSoftwareId);
    mexPrintf("  File Creation Day/Year:    %d/%d\n", 
                    LASHeader_GetCreationDOY(header), LASHeader_GetCreationYear(header));

    mexPrintf("  Header Size                %d\n", LASHeader_GetHeaderSize(header));
    mexPrintf("  Offset to Point Data       %d\n", LASHeader_GetDataOffset(header));
    mexPrintf("  Number Var. Length Records %d\n", LASHeader_GetRecordsCount(header));
    mexPrintf("  Point Data Format          %d\n", LASHeader_GetDataFormatId(header));
    mexPrintf("  Point Data Record Length   %d\n", LASHeader_GetDataRecordLength(header));
    mexPrintf("  Number of Point Records    %d\n", LASHeader_GetPointRecordsCount(header));
    mexPrintf("  Number of Points by Return %d %d %d %d %d\n", 
                    LASHeader_GetPointRecordsByReturnCount(header, 0), 
                    LASHeader_GetPointRecordsByReturnCount(header, 1), 
                    LASHeader_GetPointRecordsByReturnCount(header, 2), 
                    LASHeader_GetPointRecordsByReturnCount(header, 3), 
                    LASHeader_GetPointRecordsByReturnCount(header, 4));

    mexPrintf("  Scale Factor X Y Z         %.6g %.6g %.6g\n", 
                    LASHeader_GetScaleX(header), 
                    LASHeader_GetScaleY(header),
                    LASHeader_GetScaleZ(header));

    mexPrintf("  Offset X Y Z               %.6f %.6f %.6f\n", 
                    LASHeader_GetOffsetX(header), 
                    LASHeader_GetOffsetY(header), 
                    LASHeader_GetOffsetZ(header));

    mexPrintf("  Min X Y Z                  %.6f %.6f %.6f\n",
                    LASHeader_GetMinX(header), 
                    LASHeader_GetMinY(header), 
                    LASHeader_GetMinZ(header));

    mexPrintf("  Max X Y Z                  %.6f %.6f %.6f\n", 
                    LASHeader_GetMaxX(header), 
                    LASHeader_GetMaxY(header), 
                    LASHeader_GetMaxZ(header));
    
    /*mexPrintf(" Spatial Reference           %s\n", pszProj4);*/
#ifdef HAVE_LIBGEOTIFF
    if (pGTIF) GTIFPrint((GTIF*)pGTIF, 0, 0);
#endif
    if (nVLR && !bSkipVLR) {
        
    mexPrintf("\n---------------------------------------------------------\n");
    mexPrintf("  VLR Summary\n");
    mexPrintf("---------------------------------------------------------\n");

        for (i = 0; i < (int)nVLR; i++) {
            pVLR = LASHeader_GetVLR(header, i);

            if (LASError_GetLastErrorNum()) {
                LASError_Print("Unable to fetch VLR");
                return;
            }
            
            pszVLRUser = LASVLR_GetUserId(pVLR);
            pszVLRDescription = LASVLR_GetDescription(pVLR);
            nVLRLength = LASVLR_GetRecordLength(pVLR);
            nVLRRecordId = LASVLR_GetRecordId(pVLR);
            
            mexPrintf("   User: '%s' - Description: '%s'\n", pszVLRUser, pszVLRDescription);
            mexPrintf("   ID: %d Length: %d\n\n", nVLRRecordId, nVLRLength);
            
            LASVLR_Destroy(pVLR);
            pVLR = NULL;
            
            LASString_Free(pszVLRUser);
            LASString_Free(pszVLRDescription);
        }
        
    }
    LASString_Free(pszSignature);
    LASString_Free(pszProjectId);
    LASString_Free(pszSystemId);
    LASString_Free(pszSoftwareId);
    /*LASString_Free(pszProj4);*/
}
