/*--------------------------------------------------------------------
 *	$Id: ogrread.c 3516 2012-04-24 21:19:50Z j $
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

/* Program:	ogrread.c
 * Purpose:	matlab callable routine to read files supported by OGR
 * 		and dumping all vector data of that dataset.
 *
 * The calling syntax is:
 *	s = ogrread (vector_file);
 *
 *	"s" is a 2D or 3D structure array with fields:
 *
 *	Name:		A string holding the layer name.
 *	SRSWkt:		A string describing the reference system in the Well Known Format.
 *	SRSProj4:	A string describing the reference system as a Proj4 string
 *	BoundingBox:	The 2D dataset BoundingBox as a 2x2 matrix with Xmin/Xmax in first column and Y in second
 *	Type:		Geometry type. E.g. Point, Polygon or LineString
 *	X:		Column vector of doubles with the vector x-coordinates		
 *	Y:		Column vector of doubles with the vector y-coordinates		
 *	Z:		Same for z when vector is 3D, otherwise empty
 *	Islands:	2 columns matrix with start and ending indexes of the main Ring and its islands (if any).
 *			This only applies to Polygon geometries that have interior rings (islands).
 *	BB_geom:	Not currently assigned (would be the BoundingBox of each individual geometry)
 *	Att_number:	Number of attributes of a Feature
 *	Att_names:	Names of the attributes of a Feature
 *	Att_values:	Values of the attributes of a Feature as strings
 *	Att_types:	Because "Att_values" came as strings, this is a vector with the codes allowing
 *			the conversion into their original data types as returned by OGR_Fld_GetType(hField). 
 *			Thus if the code of element n is 2 (OFTReal) Att_values[n] can be converted to double with
 *			atof. See ogr_core.h for the list of codes and their meanings.
 *
 * Now, given the potential complexity of an OGR dataset the "s" structure can be either a 2D or 3D dimensional
 * struct array. On the simpler case of one single layer the structure is 2D. If more layers are present in the
 * dataset, each layer will be stored in each of the planes (pages) of the 3D array.
 * Also, besides the three simpler geometries (Point, LineString & Polygon) we can have also the more complex
 * combinations of MultiPoint|Line|Polygon or even the GeometryCollections. So each plane of the "s" array is
 * organized as in the following example:
 *
 *	F1 [G11 G12 ...   G1N]
 *	F2 |G21 G22 [] ... []|
 *	F3 |G31 [] ....... []|
 *	FM [.................]
 *
 * Each Gij element of the matrix is a "s" structure as describe above. The F1 .. FM are the features collection
 * in the Layer. Because one Feature can have one or many Geometries the array is defined as having as many columns
 * as the maximum number of Geometries in all Features. There are as many rows as Features in the Layer. However,
 * not all fields of the individual Gij are filled. Only those that contain real data. Also, to not waste space
 * and given that the Attributes of all Geometries of a Feature are equal, only the first column Gm1 elements
 * have assigned the "Att_names|values|types". That is, the others are defined but not filled. The same applies
 * to all other matrix elements represented as []. Again, they are defined (they have to be because the matrix
 * must be rectangular) but their members are not filled.
 *
 * One final note regarding the "Islands" struct element. Because we want to keep trace on paternity of interior
 * rings (islands) of a polygon, each one of those interior polygons is appended to the main polygon but separated
 * with one row of NaNs (2 or 3 NaNs depending if vector is 2d or 3D). The "Islands" element contains thus a Nx2
 * matrix with the indexes of the starting and ending positions of the N polygons that were once the Polygon and
 * its interior rings in the OGR model. For Polygons with no islands, "Islands" is an empty ([]) variable. 
 *
 * --------------------------------------------------------------------------------------------------------------
 * Revision 1.0  06/9/2011 Joaquim Luis
 */

#include "mex.h"
#include "ogr_srs_api.h"
#include "ogr_api.h"

char *mxStrdup(const char *s);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	int	i, j, iLayer, nEmptyLayers, nEmptyGeoms, nAttribs = 0, dims[3];
	int	region = 0, verbose = 0;
	int	*layers;		/* Array with layer numbers*/
	int	nLayers;		/* number of layers in dataset */
	char	**layerNames;		/* layers names */
	char	*fname;
	double	xmin, ymin, xmax, ymax;

	mxArray *out_struct, *mBBox;
	int	nDims, nFields, nFeature, nMaxFeatures, nMaxGeoms;
	char	*fnames[128];		/* holds the main struct field names */
	double	*bb_ptr;

	OGRDataSourceH hDS;
	OGRLayerH hLayer;
	OGRFeatureH hFeature;
	OGRFeatureDefnH hFeatureDefn;
	OGRGeometryH hGeom, hPolygon, poSpatialFilter;
	OGRSpatialReferenceH hSRS;
	OGREnvelope sEnvelop;
	OGRwkbGeometryType eType;

	xmin = ymin = xmax = ymax = 0.0;

	if (nrhs == 0)
		fname = mxStrdup("C:/SVN/ABW_adm0.shp");
	else
		fname = (char *)mxArrayToString(prhs[0]);

	OGRRegisterAll();

	hDS = OGROpen(fname, FALSE, NULL);	/* Open OGR Datasourse */
	if (hDS == NULL) {
		mexPrintf("Unable to open data source <%s>\n", fname);
		mexErrMsgTxt("\n");
	}

	nLayers = OGR_DS_GetLayerCount(hDS);	/* Get available layers */
	if (nLayers < 1)
		mexErrMsgTxt("No OGR layers available\n");

	layerNames = (char **)mxMalloc(nLayers * sizeof(char *));

	for (i = 0; i < nLayers; i++) {
		hLayer = OGR_DS_GetLayer(hDS, i);
		hFeatureDefn = OGR_L_GetLayerDefn(hLayer);
		layerNames[i] = mxStrdup((char *)OGR_FD_GetName(hFeatureDefn));
		if (verbose)
			mexPrintf("Layer name = %s\n", layerNames[i]);
	}

	if (region) {		/* If we have a sub-region request. */
		poSpatialFilter = OGR_G_CreateGeometry(wkbPolygon);
		hPolygon = OGR_G_CreateGeometry(wkbLinearRing);
		OGR_G_AddPoint(hPolygon, xmin, ymin, 0.0);
		OGR_G_AddPoint(hPolygon, xmin, ymax, 0.0);
		OGR_G_AddPoint(hPolygon, xmax, ymax, 0.0);
		OGR_G_AddPoint(hPolygon, xmax, ymin, 0.0);
		OGR_G_AddPoint(hPolygon, xmin, ymin, 0.0);
		OGR_G_AddGeometryDirectly(poSpatialFilter, hPolygon);
	}

	layers = (int *)mxMalloc(nLayers * sizeof(int));

	/* Get MAX number of features of all layers */
	for (i = j = nEmptyLayers = 0, nMaxFeatures = nMaxGeoms = 1; i < nLayers; i++) {
		hLayer = OGR_DS_GetLayer(hDS, i);

		if (region)
			OGR_L_SetSpatialFilter(hLayer, poSpatialFilter);

		nMaxFeatures = MAX(OGR_L_GetFeatureCount(hLayer, 1), nMaxFeatures);
		hFeature = OGR_L_GetNextFeature(hLayer);
		if (hFeature == NULL) {		/* Yes, this can happen. Probably on crapy files */
			nEmptyLayers++;
			continue;
		}
		layers[j++] = i;		/* Store indices of non-empty layers */
		hGeom = OGR_F_GetGeometryRef(hFeature);
		hFeatureDefn = OGR_L_GetLayerDefn(hLayer);
		eType = wkbFlatten(OGR_G_GetGeometryType(hGeom));
		if (eType != wkbPolygon)	/* For simple polygons, next would return only the number of interior rings */
			nMaxGeoms = MAX(OGR_G_GetGeometryCount(hGeom), nMaxGeoms);

		OGR_F_Destroy(hFeature);
	}
	if (nEmptyLayers) nLayers -= nEmptyLayers;

	nFields = 0;
	fnames[nFields++] = mxStrdup ("Name");
	fnames[nFields++] = mxStrdup ("SRSWkt");
	fnames[nFields++] = mxStrdup ("SRSProj4");
	fnames[nFields++] = mxStrdup ("BoundingBox");
	fnames[nFields++] = mxStrdup ("Type");
	fnames[nFields++] = mxStrdup ("X");
	fnames[nFields++] = mxStrdup ("Y");
	fnames[nFields++] = mxStrdup ("Z");
	fnames[nFields++] = mxStrdup ("Islands");
	fnames[nFields++] = mxStrdup ("BB_geom");		/* Not yet used */
	fnames[nFields++] = mxStrdup ("Att_number");
	fnames[nFields++] = mxStrdup ("Att_names");
	fnames[nFields++] = mxStrdup ("Att_values");
	fnames[nFields++] = mxStrdup ("Att_types");

	nDims = (nLayers > 1) ? 3 : 2;
	dims[0] = nMaxFeatures;		dims[1] = nMaxGeoms;
	dims[2] = nLayers;
	out_struct = mxCreateStructArray ( nDims, dims, nFields, (const char **)fnames );

	for (iLayer = nFeature = nEmptyGeoms = 0; iLayer < nLayers; iLayer++) {

		hLayer = OGR_DS_GetLayer(hDS, layers[iLayer]);
		OGR_L_ResetReading(hLayer);
		hFeatureDefn = OGR_L_GetLayerDefn(hLayer);

		hSRS = OGR_L_GetSpatialRef(hLayer);	/* Do not free it later */
		if (hSRS) {				/* Get Layer's SRS. */
			char	*pszWKT = NULL, *pszProj4 = NULL;
			mxArray	*mxPrjRef;
			if (OSRExportToProj4(hSRS, &pszProj4) == OGRERR_NONE) {
				mxPrjRef = mxCreateString (pszProj4);
				mxSetField (out_struct, 0*dims[0]*dims[1], "SRSProj4", mxPrjRef);
			}
			if (OSRExportToPrettyWkt(hSRS, &pszWKT, 1) == OGRERR_NONE) {
				mxPrjRef = mxCreateString (pszWKT);
				mxSetField (out_struct, 0*dims[0]*dims[1], "SRSWkt", mxPrjRef);
			}
		}

		/* Get this layer BoundingBox as two column vectors of X and Y respectively. */
		mBBox = mxCreateDoubleMatrix (2, 2, mxREAL);
		bb_ptr = mxGetPr (mBBox);
		if ((OGR_L_GetExtent(hLayer, &sEnvelop, 1)) == OGRERR_NONE) {
			bb_ptr[0] = sEnvelop.MinX;		bb_ptr[1] = sEnvelop.MaxX;
			bb_ptr[2] = sEnvelop.MinY;		bb_ptr[3] = sEnvelop.MaxY;
		}
		else {
			bb_ptr[0] = bb_ptr[2] = -mxGetInf();
			bb_ptr[1] = bb_ptr[3] =  mxGetInf();
		}
		mxSetField (out_struct, 0*dims[0]*dims[1], "BoundingBox", mBBox);
		mxSetField (out_struct, 0*dims[0]*dims[1], "name", mxCreateString(layerNames[iLayer]));

		nAttribs = OGR_FD_GetFieldCount(hFeatureDefn);

		if (verbose)
			mexPrintf("Importing %d features from layer <%s>\n", OGR_L_GetFeatureCount(hLayer, 1), layerNames[iLayer]);
		while ((hFeature = OGR_L_GetNextFeature(hLayer)) != NULL) {
			hGeom = OGR_F_GetGeometryRef(hFeature);
			eType = wkbFlatten(OGR_G_GetGeometryType(hGeom));
			if (hGeom != NULL)
				get_data(out_struct, hFeature, hFeatureDefn, hGeom, iLayer, nFeature, 
					nLayers, nAttribs, nMaxFeatures, 0);
			else
				nEmptyGeoms++;
			nFeature++;		/* Counter to number of features in this layer */

			OGR_F_Destroy(hFeature);
		}

		if (verbose && nEmptyGeoms > 0)
			mexPrintf("There were %d empty geometries in this layer\n", nEmptyGeoms);
	}

	OGR_DS_Destroy(hDS);

	if (nlhs > 0)
		plhs[0] = out_struct;
}

int get_data(mxArray *out_struct, OGRFeatureH hFeature, OGRFeatureDefnH hFeatureDefn, OGRGeometryH hGeom, 
		int iLayer, int nFeature, int nLayers, int nAttribs, int nMaxFeatures, int recursionLevel) {
	int	dims[2], is3D, i, j, jj, k, np, nPtsBase, nRings, indStruct, nGeoms, do_recursion;
	const char	**att_names, **att_values;  	/* To hold names and values of fields */
	OGRwkbGeometryType eType;
	OGRGeometryH hRing, hRingIsland;
	OGRFieldDefnH hField;
	mxArray *x_out, *y_out, *z_out, *bbox, *pTypes;
	double	*x_out_ptr, *y_out_ptr, *z_out_ptr, *bb_ptr, nan, *ptr_d;

	is3D = (OGR_G_GetCoordinateDimension(hGeom) > 2);
	eType = wkbFlatten(OGR_G_GetGeometryType(hGeom));
	do_recursion = (eType == wkbGeometryCollection || eType == wkbMultiPolygon ||	/* Find if we are going to do recursive calls */
			eType == wkbMultiLineString || eType == wkbMultiPoint);

	if (!do_recursion) {	/* Multis are break up into pieces by recursive calls, so no need to alloc anything */
		if (eType == wkbPoint || eType == wkbLineString) {
			np = OGR_G_GetPointCount(hGeom);
		}
		else {
			hRing = OGR_G_GetGeometryRef(hGeom, 0);
			np = OGR_G_GetPointCount(hRing);
			nRings = OGR_G_GetGeometryCount(hGeom);
			for (i = 1; i < nRings; i++) {		/* If we have islands (interior rings) */
				hRingIsland = OGR_G_GetGeometryRef(hGeom, i);
				np += OGR_G_GetPointCount(hRingIsland);
			}
			np += (nRings - 1);			/* To account for NaNs separating islands from outer ring */
		}
		dims[0] = np;		dims[1] = 1;
		x_out = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);
		y_out = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);
		x_out_ptr = mxGetData (x_out);
		y_out_ptr = mxGetData (y_out);
		if (is3D) {
			z_out = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);
			z_out_ptr = mxGetData (z_out);
		}
	}

	nGeoms = OGR_G_GetGeometryCount(hGeom);
	if ( (nGeoms == 0) && (eType == wkbPoint || eType == wkbLineString) )	/* These geometries will silently return 0 */
		nGeoms = 1;
	else if ( (nGeoms > 1) && (eType == wkbPolygon) )	/* Other geometries are Islands and those are dealt separately */
		nGeoms = 1;
	else if (nGeoms == 0) {
		mexPrintf("Screammm: No Geometries in this Feature\n");
		return(-1);
	}

	for (j = 0; j < nGeoms; j++) {		/* Loop over the number of geometries in this feature */

		jj = ((eType == wkbPolygon) ? 0 : j) + recursionLevel;	/* Islands don't count for matrix index (they are appended) */
		indStruct = (nFeature + jj * nMaxFeatures) * nLayers;

		if (eType == wkbPoint || eType == wkbLineString) {
			for (i = 0; i < np; i++) {
				x_out_ptr[i] = OGR_G_GetX(hGeom, i);
				y_out_ptr[i] = OGR_G_GetY(hGeom, i);
			}
			if (is3D) {
				for (i = 0; i < np; i++)
					z_out_ptr[i] = OGR_G_GetZ(hGeom, i);
			}
			if (eType == wkbPoint)
				mxSetField(out_struct, indStruct, "Type", mxCreateString("Point"));
			else
				mxSetField(out_struct, indStruct, "Type", mxCreateString("LineString"));
	
		}
		else if (eType == wkbPolygon) {
			nPtsBase = OGR_G_GetPointCount(hRing);	/* Need to ask it again because prev value may have eventual islands*/
			for (i = 0; i < nPtsBase; i++) {
				x_out_ptr[i] = OGR_G_GetX(hRing, i);
				y_out_ptr[i] = OGR_G_GetY(hRing, i);
			}
			if (is3D) {
				for (i = 0; i < nPtsBase; i++)
					z_out_ptr[i] = OGR_G_GetZ(hRing, i);
			}

			if (nRings > 1) {		/* Deal with the Islands */
				int	c, cz, nPtsRing, dimsI[2], *pi;
				double	nan;
				mxArray	*mxIslands;
				nan = mxGetNaN();
				c = nPtsBase;
				dimsI[0] = nRings;	dimsI[1] = 2;
				mxIslands = mxCreateNumericArray (2, dimsI, mxINT32_CLASS, mxREAL);
				pi = (int *)mxGetData(mxIslands);
				pi[0] = 0;
				for (k = 1; k < nRings; k++) {			/* Loop over islands (interior rings) */
					hRingIsland = OGR_G_GetGeometryRef(hGeom, k);
					nPtsRing = OGR_G_GetPointCount(hRingIsland);
					x_out_ptr[c] = y_out_ptr[c] = nan;
					if (is3D) z_out_ptr[c] = nan;
					c++;	cz = c;
					pi[k] = c;
					for (i = 0; i < nPtsRing; c++, i++) {	/* Loop over number of points in this island */
						x_out_ptr[c] = OGR_G_GetX(hRingIsland, i);
						y_out_ptr[c] = OGR_G_GetY(hRingIsland, i);
					}
					if (is3D) {
						for (i = 0; i < nPtsRing; cz++, i++)
							z_out_ptr[cz] = OGR_G_GetZ(hRingIsland, i);
					}
				}
				/* We still have to fill the second column of Islands, but only now we have the necessary info */
				for (k = 0; k < nRings - 1; k++)
					pi[nRings + k] = pi[k+1] - 2;
				pi[2*nRings - 1] = c - 1;		/* Last element was not assigned in the loop above */
				mxSetField(out_struct, indStruct, "Islands", mxIslands);
			}
			mxSetField(out_struct, indStruct, "Type", mxCreateString("Polygon"));
		}
		else if (do_recursion) {
			/* When we reach here it's because the current Geometry is of the Multi<something> type.
			 * The way we deal with it is to decompose it in its individual simple geometries, e.g.
			 * Polygon and call this function recursively until all basic geometries, controlled by
			 * the main for loop above [for (j = 0; j < nGeoms; j++)], are processed. */
			int	r;
	
			hRing = OGR_G_GetGeometryRef(hGeom, j);
			r = get_data(out_struct, hFeature, hFeatureDefn, hRing, iLayer, nFeature, nLayers, nAttribs, nMaxFeatures, j);
			if (r)
				mexPrintf("Unable to get data from element of a Multi<something>\n");
			continue;	/* We are done here */
		}
		else {
			mexPrintf("Unforeseen case -> unknown geometry type\n");
			return(-1);
		}

		mxSetField (out_struct, indStruct, "X", x_out);
		mxSetField (out_struct, indStruct, "Y", y_out);
		if (is3D) 
			mxSetField (out_struct, indStruct, "Z", z_out);

		if ((j + recursionLevel) == 0) {
			/* Only first column element is set with number of attributes (also called fields by ogr) */
			mxSetField(out_struct, indStruct, "Att_number", mxCreateDoubleScalar((double)nAttribs));

			att_names  = (const char **)mxCalloc((size_t)nAttribs, sizeof(char *));
			att_values = (const char **)mxCalloc((size_t)nAttribs, sizeof(char *));
			dims[0] = nAttribs;		dims[1] = 1;
			pTypes = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
			ptr_d = mxGetPr(pTypes);
			for (i = 0; i < nAttribs; i++) {
				hField        = OGR_FD_GetFieldDefn(hFeatureDefn, i);
				att_names[i]  = mxStrdup(OGR_Fld_GetNameRef(hField));
				att_values[i] = mxStrdup(OGR_F_GetFieldAsString(hFeature, i));
				ptr_d[i]      = OGR_Fld_GetType(hField);
			}
			mxSetField(out_struct, indStruct, "Att_names",  mxCreateCharMatrixFromStrings(nAttribs, att_names));
			mxSetField(out_struct, indStruct, "Att_values", mxCreateCharMatrixFromStrings(nAttribs, att_values));
			mxSetField(out_struct, indStruct, "Att_types",  pTypes);
		}
		else
			mxSetField(out_struct, indStruct, "Att_number", mxCreateDoubleScalar(0));
	}

	return(0);
}

char *mxStrdup(const char *s) {
    char *buf;

    if (s == NULL) {
	buf = mxMalloc(1);
	buf[0] = '\0';
    }
    else {
	buf = mxMalloc(strlen(s) + 1);
	strcpy(buf, s);
    }
    
    return buf;
}
