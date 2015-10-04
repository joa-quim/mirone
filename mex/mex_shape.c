/*================================================================= *
 * MEX_SHAPE.C
 *     Gateway routine to interface with shapelib library.
 *
 *     This file evolved from shpdump.c as provided by the shapelib
 *     distribution.  Some flotsam and jetsam from that file are
 *     still present.
 *
 * The calling syntax is:
 *     [s, t] = mex_shapefile ( shapefile );
 *
 *     "s" is a structure array with at least three fields, "mx_data" and
 *     "my_data", which contain the vertices x and y data for each
 *     respective structure element, as well as the shape type.
 *
 *     "type" is a string defining the shape type.  This can be
 *     "Polygon"
 *
 * PARAMETERS:
 * Input:
 *   shapefile:
 *      character filename of "*.shp" file
 * Output:
 *   shape:
 *      structure array of shapefiles with their associated attributes.
 *   type:
 *      type of shapefile that was read.
 *
 * In case of an error, an exception is thrown.
 *
 *=================================================================*/
/* $Revision: 1.5 If POINT or POINTZ output struct with vectors instead of a vector os structs - JL 24-02-10 */
/* $Revision: 1.4 First BoundingBox contains the ensemble extent - JL 28-01-07 */
/* $Revision: 1.3 BoundingBox per element - JL 4-10-06 */
/*                Insert NaNs between rings on polylines */
/* $Revision: 1.2 Outputs the BoundingBox - JL 14-3-06 */
/* $Revision: 1.1 $ */

#include "mex.h"
#include <string.h>
#include "shapefil.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) { 

	/* Pointer to temporary matlab array */
	const mxArray *mx;
	mxArray *m_shape_type;

	int	nShapeType, nEntities, i, j, k, iPart, nFields, nStructElem, isPoint, firstInPoints = 1;
	int	num_dbf_fields, num_dbf_records; /* number of DBF attributes, records */
	int	buflen;		/* length of input shapefile name */
	int	status;		/* success or failure */
	char	*shapefile;	/* holder for input shapefile name */
	const	char *pszPartType = "";

	/* Not used. */
	double adfMinBound[4], adfMaxBound[4];

	/* pointer to the shapefile */
	SHPHandle	hSHP;
	SHPObject	*psShape;
	DBFHandle	dbh;	/* handle for DBF file */

	/* This structure will hold a description of each field in the DBF file. */
	typedef struct DBF_Field_Descriptor {
		char          pszFieldName[12];
		DBFFieldType  field_type;
	} DBF_Field_Descriptor;

	DBF_Field_Descriptor *dbf_field;

	/* stores individual values from the DBF.  */
	int	dbf_integer_val, dims[2], *p_parts_ptr, nNaNs, c, i_start, i_stop;
	char	*dbf_char_val, error_buffer[500];
	char	*fnames[100];  /* holds name of fields */
	mxArray *out_struct, *x_out, *y_out, *z_out, *bbox, *p_parts;
	double	*x_out_ptr, *y_out_ptr, *z_out_ptr, *bb_ptr, nan, dbf_double_val;
	size_t	sizebuf;

	/*  Initialize the dbf record.  */
	dbf_field = NULL;

	/* Check for proper number of arguments */
	if (nrhs != 1)
		mexErrMsgTxt("One input arguments are required."); 

	if (nlhs != 2)
		mexErrMsgTxt("Two output arguments required."); 

	/* Make sure the input is a proper string. */
	if ( mxIsChar(prhs[0]) != 1 )
		mexErrMsgTxt("Shapefile parameter must be a string\n" );

	if ( mxGetM(prhs[0]) != 1 )
		mexErrMsgTxt("Shapefile parameter must be a row vector, not a column string\n" );

	buflen = mxGetN(prhs[0]) + 1; 
	shapefile = mxCalloc( buflen, sizeof(char) );

	/* copy the string data from prhs[0] into a C string. */
	status = mxGetString( prhs[0], shapefile, buflen );
	if (status != 0)
		mexErrMsgTxt( "Not enough space for shapefile argument.\n" );

	/* -------------------------------------------------------------------- */
	/*      Open the passed shapefile.                                      */
	/* -------------------------------------------------------------------- */
	hSHP = SHPOpen( shapefile, "rb" );
	if (hSHP == NULL) {
		mexPrintf("Unable to open:%s\n", shapefile);
		mexErrMsgTxt("\n");
	}

	/* -------------------------------------------------------------------- */
	/*      Get the needed information about the shapefile.                 */
	/* -------------------------------------------------------------------- */
	SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

	/* Make sure that we can handle the type. */
	switch ( nShapeType ) {
		case SHPT_POINT:
			break;
		case SHPT_POINTZ:
			break;
		case SHPT_ARC:
			break;
		case SHPT_ARCZ:
			break;
		case SHPT_POLYGON:
			break;
		case SHPT_POLYGONZ:
			break;
		case SHPT_MULTIPOINT:		/* JL */
			break;
		default:
			sprintf ( error_buffer, "Unhandled shape code %d (%s)", nShapeType, SHPTypeName ( nShapeType ) );
	        	mexErrMsgTxt( error_buffer ); 
	}
    
	/* Create the output shape type parameter. */
	plhs[1] = mxCreateString ( SHPTypeName ( nShapeType ) );

	/* Open the DBF, and retrieve the number of fields and records */
	dbh = DBFOpen (shapefile, "rb");
	num_dbf_fields = DBFGetFieldCount ( dbh );
	num_dbf_records = DBFGetRecordCount ( dbh );

	/* Allocate space for a description of each record, and populate it.
	 * I allocate space for two extra "dummy" records that go in positions
	 * 0 and 1.  These I reserve for the xy data. */
	nFields = 3;
	if ( (nShapeType == SHPT_POLYGONZ) || (nShapeType == SHPT_ARCZ) || (nShapeType == SHPT_POINTZ) ) nFields++;
	dbf_field = (DBF_Field_Descriptor *) mxMalloc ( (num_dbf_fields+nFields) * sizeof ( DBF_Field_Descriptor ) );
	if ( dbf_field == NULL )
		mexErrMsgTxt("Memory allocation for DBF_Field_Descriptor failed."); 

	for ( j = 0; j < num_dbf_fields; ++j )
		dbf_field[j+nFields].field_type = DBFGetFieldInfo ( dbh, j, dbf_field[j+nFields].pszFieldName, NULL, NULL );

	fnames[0] = strdup ( "X" );
	fnames[1] = strdup ( "Y" );
	if ( (nShapeType == SHPT_POLYGONZ) || (nShapeType == SHPT_ARCZ) || (nShapeType == SHPT_POINTZ) ) {
		fnames[2] = strdup ( "Z" );
		fnames[3] = strdup ( "BoundingBox" );
	}
	else
		fnames[2] = strdup ( "BoundingBox" );

	for ( j = 0; j < num_dbf_fields; ++j )
		fnames[j+nFields] = strdup ( dbf_field[j+nFields].pszFieldName );

	/* To hold information on eventual polygons with rings */
	/*fnames[num_dbf_fields+3] = strdup ( "nParts" );*/
	/*fnames[num_dbf_fields+4] = strdup ( "PartsIndex" );*/

	/* Allocate space for the output structure. */
	isPoint = ( nShapeType == SHPT_POINT || nShapeType == SHPT_POINTZ ) ? 1: 0;
	nStructElem = ( nShapeType == SHPT_POINT || nShapeType == SHPT_POINTZ ) ? 1: nEntities;
	out_struct = mxCreateStructMatrix ( nStructElem, 1, nFields + num_dbf_fields, (const char **)fnames );

	/* create the BoundingBox */
	dims[0] = 4;		dims[1] = 2;
	bbox = mxCreateNumericArray ( 2, dims, mxDOUBLE_CLASS, mxREAL );
	bb_ptr = mxGetData ( bbox );
	for (i = 0; i < 4; i++) bb_ptr[i]   = adfMinBound[i];
	for (i = 0; i < 4; i++) bb_ptr[i+4] = adfMaxBound[i];
	mxSetField ( out_struct, 0, "BoundingBox", bbox );

	nan = mxGetNaN();
	/* -------------------------------------------------------------------- */
	/*	Skim over the list of shapes, printing all the vertices.	*/
	/* -------------------------------------------------------------------- */
	for( i = 0; i < nEntities; i++ ) {
		psShape = SHPReadObject( hSHP, i );

		/* Create the fields in this struct element.  */
		if ( !isPoint ) {
			nNaNs = psShape->nParts > 1 ? psShape->nParts : 0;
			dims[0] = psShape->nVertices + nNaNs;
			dims[1] = 1;
			x_out = mxCreateNumericArray ( 2, dims, mxDOUBLE_CLASS, mxREAL );
			x_out_ptr = mxGetData ( x_out );
			y_out = mxCreateNumericArray ( 2, dims, mxDOUBLE_CLASS, mxREAL );
			y_out_ptr = mxGetData ( y_out );
			if ( (nShapeType == SHPT_POLYGONZ) || (nShapeType == SHPT_POINTZ) || (nShapeType == SHPT_ARCZ)) {
				z_out = mxCreateNumericArray ( 2, dims, mxDOUBLE_CLASS, mxREAL );
				z_out_ptr = mxGetData ( z_out );
			}
		}
		else if (firstInPoints) {	/* Allocate all memory we'll need */
			x_out = mxCreateDoubleMatrix (nEntities, 1, mxREAL);
			x_out_ptr = mxGetPr ( x_out );
			y_out = mxCreateDoubleMatrix (nEntities, 1, mxREAL);
			y_out_ptr = mxGetPr ( y_out );
			if (nShapeType == SHPT_POINTZ) {
				z_out = mxCreateDoubleMatrix (nEntities, 1, mxREAL);
				z_out_ptr = mxGetPr ( z_out );
			}
			firstInPoints = 0;
		}

		if (!isPoint && psShape->nParts > 1) {
			for (k = c = 0; k < psShape->nParts; k++) {
				i_start = psShape->panPartStart[k];
				if (k < psShape->nParts - 1)
					i_stop = psShape->panPartStart[k+1];
				else
					i_stop = psShape->nVertices;

				if ( nShapeType == SHPT_POLYGONZ ) {
					for (j = i_start; j < i_stop; c++, j++) {
						x_out_ptr[c] = psShape->padfX[j];
						y_out_ptr[c] = psShape->padfY[j];
						z_out_ptr[c] = psShape->padfZ[j];
					}
					x_out_ptr[c] = nan;
					y_out_ptr[c] = nan;
					z_out_ptr[c] = nan;
				}
				else {
					for (j = i_start; j < i_stop; c++, j++) {
						x_out_ptr[c] = psShape->padfX[j];
						y_out_ptr[c] = psShape->padfY[j];
					}
					x_out_ptr[c] = nan;
					y_out_ptr[c] = nan;
				}
				c++;
			}
		}
		else if ( isPoint ) {
			x_out_ptr[i] = *psShape->padfX;
			y_out_ptr[i] = *psShape->padfY;
			if (nShapeType == SHPT_POINTZ) z_out_ptr[i] = *psShape->padfZ;
			if (i > 0) {
        			SHPDestroyObject( psShape );
				continue;
			}
		}
		else {		/* Just copy the vertices over. */
			sizebuf = mxGetElementSize ( x_out ) * psShape->nVertices;
			memcpy ( (void *) x_out_ptr, (void *) psShape->padfX, sizebuf );
			memcpy ( (void *) y_out_ptr, (void *) psShape->padfY, sizebuf );
			if ( (nShapeType == SHPT_POLYGONZ)  || (nShapeType == SHPT_ARCZ))
				memcpy ( (void *) z_out_ptr, (void *) psShape->padfZ, sizebuf );
		}

		mxSetField ( out_struct, i, "X", x_out );
		mxSetField ( out_struct, i, "Y", y_out );
		if ( (nShapeType == SHPT_POLYGONZ) || (nShapeType == SHPT_ARCZ) )
			mxSetField ( out_struct, i, "Z", z_out );

		bbox = mxCreateNumericMatrix ( 4, 2, mxDOUBLE_CLASS, mxREAL );
		bb_ptr = (double *)mxGetData ( bbox );
		bb_ptr[0] = psShape->dfXMin;		bb_ptr[1] = psShape->dfYMin;
		bb_ptr[2] = psShape->dfZMin;		bb_ptr[3] = psShape->dfMMin;
		bb_ptr[4] = psShape->dfXMax;		bb_ptr[5] = psShape->dfYMax;
		bb_ptr[6] = psShape->dfZMax;		bb_ptr[7] = psShape->dfMMax;
		if (i > 0)	/* First BB contains the ensemble extent */
			mxSetField ( out_struct, i, "BoundingBox", bbox );

		for ( j = 0; j < num_dbf_fields; ++j ) {
			switch ( dbf_field[j+nFields].field_type ) {
				case FTString:
					dbf_char_val = (char *) DBFReadStringAttribute ( dbh, i, j );
					mxSetField ( out_struct, i, dbf_field[j+nFields].pszFieldName, mxCreateString ( dbf_char_val ) );
					break;
				case FTDouble:
					dbf_double_val = DBFReadDoubleAttribute ( dbh, i, j );
					mxSetField ( out_struct, i, dbf_field[j+nFields].pszFieldName, mxCreateDoubleScalar ( dbf_double_val ) );
					break;
				case FTInteger:
				case FTLogical:
					dbf_integer_val = DBFReadIntegerAttribute ( dbh, i, j );
					dbf_double_val = dbf_integer_val;
					mxSetField ( out_struct, i, dbf_field[j+nFields].pszFieldName, mxCreateDoubleScalar ( dbf_double_val ) );
					break;
				default:
					sprintf(error_buffer, "Unhandled code %d, shape %d, record %d\n", dbf_field[j+nFields].field_type, i, j);
	        			mexErrMsgTxt("Unhandled code"); 
			}
		}
        	SHPDestroyObject( psShape );
	}

	if ( isPoint ) {	/* In this case we still need to "send the true data out" */
		mxSetField ( out_struct, 0, "X", x_out );
		mxSetField ( out_struct, 0, "Y", y_out );
		if (nShapeType == SHPT_POINTZ)
			mxSetField ( out_struct, 0, "Z", z_out );
	}

	/* Clean up, close up shop. */
	SHPClose( hSHP );
	DBFClose ( dbh );

	if ( dbf_field != NULL )
		mxFree ( (void *)dbf_field );

	plhs[0] = out_struct;
}
