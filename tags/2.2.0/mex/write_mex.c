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

#include "mex.h"
#include "mwsize.h"

void stl_write_ascii(FILE *fp, float *v, mwSize *f, int nFacet, int nv);
void stl_write_binary(FILE *fp, float *v, mwSize *f, int nFacet, int nv);
void xyz_write_ascii(FILE *fp, float *v, int nv);
void triang_norm (float *Va, float *Vb, float *Vc, float *n);

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	float	*v;
	int	nv, nf, do_STL = 1, do_XYZ = 0, do_ascii = 1;
	mwSize	*f;
	char	*file, *type;
	FILE	*fp;
  
	/* Check for proper number of arguments */    
	if (nrhs == 0) {
		mexPrintf("Usage:  write_mex(V, F, fname [, type[, mode]])\n");
		mexPrintf("        Where F,V like in a patch object are:\n");
		mexPrintf("        	V -> Vertices; A Mx3 array of SINGLEs\n");
		mexPrintf("        	F -> Faces;    A Nx3 array of INT32\n");
		mexPrintf("        TYPE can be STL (Stereo Litographic) [Default] or XYX (plain x,y,z)\n");
		mexPrintf("        MODE is either ASC (save ascii) [Default] or BIN (binary)\n");
		mexPrintf("        NOTE: when selecting a MODE the TYPE must be given as well.\n");
		return;
	}
 
	if (nrhs < 3) 
		mexErrMsgIdAndTxt("write_mex:WrongNumberOfInputs", "Minimum is 3 -> write_mex(v, f, fname)");

	if (!mxIsSingle( prhs[0]))
		mexErrMsgIdAndTxt("write_mex:BadInputType", "The Vertices array must be of type SINGLE");

	if (!mxIsInt32( prhs[1]))
		mexErrMsgIdAndTxt("write_mex:BadInputType", "The Faces array must be of type INT32");

	v = (float *)mxGetData( prhs[0] );
	f = (mwSize *)mxGetData( prhs[1] );
	file = (char *)mxArrayToString( prhs[2] );

	if (nrhs >= 4) {		/* What type of format */
		type = (char *)mxArrayToString( prhs[3] );
		if (!strcmp(type, "STL") || !strcmp(type, "stl"))
			do_STL = 1;
		else if (!strcmp(type, "XYZ") || !strcmp(type, "xyz"))
			{do_XYZ = 1;	do_STL = 0;}
	}
	if (nrhs == 5) {		/* ASCII or Binary */
		type = (char *)mxArrayToString( prhs[4] );
		if (!strncmp(type, "ASC", 3) || !strncmp(type, "asc", 3))
			do_ascii = 1;
		else if (!strncmp(type, "BIN", 3) || !strncmp(type, "bin", 3))
			do_ascii = 0;
	}

	nv = mxGetM( prhs[0] );
	nf = mxGetM( prhs[1] );

	/* Open the file */
	if (do_ascii) {
		if ((fp = fopen(file, "w")) == NULL)
			mexErrMsgIdAndTxt("write_formated:Error", "Couldn't open file for writing");

		if (do_STL)
			stl_write_ascii(fp, v, f, nf, nv);
		else
			xyz_write_ascii(fp, v, nv);
	}
	else {
		if ((fp = fopen(file, "wb")) == NULL)
			mexErrMsgIdAndTxt("write_formated:Error", "Couldn't open file for writing");

		if (do_STL)
			stl_write_binary(fp, v, f, nf, nv);
	}

	fclose(fp);

}

void stl_write_ascii(FILE *fp, float *v, mwSize *f, int nFacet, int nv) {
	int	i, i1, i2, i3;
	float norm[3], Va[3], Vb[3], Vc[3];
 
	fprintf(fp, "solid  Created By Mirone\n");
  
	for (i = 0; i < nFacet; i++) {
		i1 = f[i];		i2 = f[i + nFacet];		i3 = f[i + 2*nFacet];
		Va[0] = v[i1-1];	Vb[0] = v[i2-1];		Vc[0] = v[i3-1];
		Va[1] = v[i1 + nv-1];	Vb[1] = v[i2 + nv-1];		Vc[1] = v[i3 + nv-1];
		Va[2] = v[i1 + 2*nv-1];	Vb[2] = v[i2 + 2*nv-1];		Vc[2] = v[i3 + 2*nv-1];
	
		triang_norm (Va, Vb, Vc, norm);

		fprintf(fp, "\tfacet normal %.9g %.9g %.9g\n", norm[0], norm[1], norm[2]);
		fprintf(fp, "\t\touter loop\n");
		fprintf(fp, "\t\t\tvertex %.12g %.12g %.12g\n", Va[0], Va[1], Va[2]);
		fprintf(fp, "\t\t\tvertex %.12g %.12g %.12g\n", Vb[0], Vb[1], Vb[2]);
		fprintf(fp, "\t\t\tvertex %.12g %.12g %.12g\n", Vc[0], Vc[1], Vc[2]);
		fprintf(fp, "\t\tendloop\n");
		fprintf(fp, "\tendfacet\n");
	}
  
	fprintf(fp, "endsolid  Created By Mirone\n");
}


void stl_write_binary(FILE *fp, float *v, mwSize *f, int nFacet, int nv) {
	char *label = "Created by Mirone";
	short int i0 = 0;
	int	i, i1, i2, i3;
	float norm[3], Va[3], Vb[3], Vc[3];

	fprintf(fp, "%s", label);
	for (i = strlen(label); i < 80; i++) putc(0, fp);

	fseek(fp, 80, SEEK_SET);

	fwrite ((void *)&nFacet, sizeof(int), 1, fp);
  
	for (i = 0; i < nFacet; i++) {
		i1 = f[i];		i2 = f[i + nFacet];		i3 = f[i + 2*nFacet];
		Va[0] = v[i1-1];	Vb[0] = v[i2-1];		Vc[0] = v[i3-1];
		Va[1] = v[i1 + nv-1];	Vb[1] = v[i2 + nv-1];		Vc[1] = v[i3 + nv-1];
		Va[2] = v[i1 + 2*nv-1];	Vb[2] = v[i2 + 2*nv-1];		Vc[2] = v[i3 + 2*nv-1];
	
		triang_norm (Va, Vb, Vc, norm);
		fwrite ((void *)norm, sizeof(float), 3, fp);
		fwrite ((void *)Va, sizeof(float), 3, fp);
		fwrite ((void *)Vb, sizeof(float), 3, fp);
		fwrite ((void *)Vc, sizeof(float), 3, fp);
		fwrite ((void *)&i0, sizeof(short int), 1, fp);
	}
}

void xyz_write_ascii(FILE *fp, float *v, int nv) {
	int	i;
 
	for (i = 0; i < nv; i++)
		fprintf(fp, "%.9g %.9g %.9g\n", v[i], v[i+nv], v[i+2*nv]);
}

void triang_norm (float *Va, float *Vb, float *Vc, float *n) {
	/* Computes the unit normal to trianglular facet */
	double dv1[3], dv2[3], dv3[3], mod;

	dv1[0] = Va[0] - Vb[0];
	dv1[1] = Va[1] - Vb[1];
	dv1[2] = Va[2] - Vb[2];

	dv2[0] = Vb[0] - Vc[0];
	dv2[1] = Vb[1] - Vc[1];
	dv2[2] = Vb[2] - Vc[2];

	dv3[0] = dv1[1]*dv2[2] - dv1[2]*dv2[1];
	dv3[1] = dv1[2]*dv2[0] - dv1[0]*dv2[2];
	dv3[2] = dv1[0]*dv2[1] - dv1[1]*dv2[0];

	mod = sqrt(dv3[0]*dv3[0] + dv3[1]*dv3[1] + dv3[2]*dv3[2]);
	n[0] = (float)(dv3[0] / mod);
	n[1] = (float)(dv3[1] / mod);
	n[2] = (float)(dv3[2] / mod);
}
