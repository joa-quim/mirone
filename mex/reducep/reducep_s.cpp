/* Copyright 1993-2007 The MathWorks, Inc. */

/*
  Source file for reducep MEX file
  This code is based on "QSlim Simplification Software"
  http://almond.srv.cs.cmu.edu/afs/cs/user/garland/www/quadrics/qslim.html
*/

static char rcsid[] = "$Revision: 1.8.4.7 $";

#include <math.h>
#include <stdlib.h>

#include "AdjModel.h"
#include "decimate.h"
#include "avars.h"
#include "mex.h"
#include "../mwsize.h"

void reduce_init(Model *M0, float *v, mwSize *f, mwSize nvin, mwSize nfin, int verbose) {
	mwSize i, initialVertCount, initialEdgeCount, initialFaceCount;
    
	will_constrain_boundaries = true;
	boundary_constraint_weight = 1000;
	placement_policy = PLACE_ENDPOINTS;

	// create vertices
	Vec3 vec;
	for (i = 0; i < nvin; i++) {    
		vec[X] = v[i];
		vec[Y] = v[i+nvin];
		vec[Z] = v[i+2*nvin];
		M0->in_Vertex(vec);
	}
    
	// create faces
	for (i = 0; i < nfin; i++)
		M0->in_Face((int)(f[i]), (int)(f[i+nfin]), (int)(f[i+2*nfin]));
    
	M0->bounds.complete();
    
	initialVertCount = M0->vertCount();
	initialEdgeCount = M0->edgeCount();
	initialFaceCount = M0->faceCount();
    
	// Get rid of degenerate faces
	int fkill=0;
	for(i=0; i<M0->faceCount(); i++)
		if( !M0->face(i)->plane().isValid() )
			M0->killFace(M0->face(i));
    
	M0->removeDegeneracy(M0->allFaces());

	// Get rid of unused vertices
	int vkill = 0;
	for(i=0; i<M0->vertCount(); i++) {
		if( M0->vertex(i)->edgeUses().length() == 0 ) {
			M0->vertex(i)->kill();
			vkill++;
		}
	}
    
	if (verbose) 
		printf("REDUCEPATCH: Cleaned Input: Vertices=%d, Faces=%d\n",
			M0->validVertCount-vkill,  M0->validFaceCount-fkill);	
}

void reduce_execute( Model *M0, int face_target, int verbose, float *vOut, int *fOut, 
			mwSize *nv, mwSize *nf, mwSize nvout, mwSize nfout) {
	mwSize i, j, numOrigFaces, verboseIncrement;
	mwSize uniqVerts = 0;

	decimate_init(*M0, pair_selection_tolerance);

	numOrigFaces = M0->validFaceCount;
	verboseIncrement = (mwSize) ceil((numOrigFaces-face_target)/50.0);
    
	if (verbose) printf("REDUCEPATCH: Working");

	while( M0->validFaceCount > face_target && M0->validFaceCount > 0
		/*&& decimate_min_error() < error_tolerance*/ ) {
		if (verbose && (numOrigFaces-M0->validFaceCount)%verboseIncrement==0) printf(".");
		decimate_contract(*M0);
	}
	if (verbose) printf("\n");
    
	for(i=0, j=0; i<M0->vertCount(); i++) {
		if( M0->vertex(i)->isValid() ) {
			M0->vertex(i)->tempID = ++uniqVerts;
			const Vertex& v = *M0->vertex(i);
			vOut[j]           = (float)v[X];
			vOut[j + nvout]   = (float)v[Y];
			vOut[j + 2*nvout] = (float)v[Z];
			j++;
		}
		else
			M0->vertex(i)->tempID = -1;
	}
    
	*nv = j;
    
	if (verbose) printf("REDUCEPATCH: Done: Vertices=%d, ",j);
     
	for(i=0, j=0; i<M0->faceCount(); i++) {
		if( M0->face(i)->isValid() ) {
			Face *f = M0->face(i);
			Vertex *v0 = (Vertex *)f->vertex(0);
			Vertex *v1 = (Vertex *)f->vertex(1);
			Vertex *v2 = (Vertex *)f->vertex(2);

			fOut[j]           = (int)v0->tempID;
			fOut[j + nfout]   = (int)v1->tempID;
			fOut[j + 2*nfout] = (int)v2->tempID;
			j++;
		}
	}

	*nf = j;

	if (verbose) printf("Faces=%d\n",j);
     
	decimate_cleanup();
}


void reduce( int target, float *v, mwSize *f, mwSize nvin, mwSize nfin, int verbose, mwSize nvout, 
		mwSize nfout, float *vOut, int *fOut, mwSize *nv, mwSize *nf) {
	Model M0;

	reduce_init(&M0, v, f, nvin, nfin, verbose);
	reduce_execute(&M0, target, verbose, vOut, fOut, nv, nf, nvout, nfout);
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
	float	*v, *vOut;
	mwSize	ndims, nvin, nfin, nvout, nfout, nv, nf, *f;
	int	target, verbose, *fOut;
  
	/* Check for proper number of arguments */    
	if (nrhs != 4) 
		mexErrMsgIdAndTxt("reducep:WrongNumberOfInputs",
                        "reducep requires 4 input arguments. [v f] = reducep(v, f, targetFaces, verbose)");
	else if (nlhs != 2) 
		mexErrMsgIdAndTxt("reducep:WrongNumberOfOutputs",
                        "reducep requires 2 output arguments. [v f] = reducep(v, f, targetFaces, verbose)");
    
    
	/* Check the dimensions of data It must be 2 */
    
	ndims = mxGetNumberOfDimensions(prhs[0]);
	if (ndims != 2) 
		mexErrMsgIdAndTxt("reducep:WrongNumberOfDimensions", "reducep requires 2 dimensional data");
	ndims = mxGetNumberOfDimensions(prhs[1]);
	if (ndims != 2) 
		mexErrMsgIdAndTxt("reducep:WrongNumberOfDimensions", "reducep requires 2 dimensional data");

	v = (float *)mxGetData( prhs[0] );
	f = (mwSize *)mxGetData( prhs[1] );
	target = (int)(*(mxGetPr( prhs[2] )));
	verbose = (int)(*(mxGetPr( prhs[3] )));

	nvin = mxGetM( prhs[0] );
	nfin = mxGetM( prhs[1] );
    
	nfout = target;
	nvout = min(nvin, 3*nfout);

	/* Create matrices for the return arguments */
	plhs[0] = mxCreateNumericMatrix (nvout,3,mxSINGLE_CLASS,mxREAL);
	plhs[1] = mxCreateNumericMatrix (nfout,3,mxINT32_CLASS,mxREAL);
    
	/* Assign pointers to the various parameters */
	vOut = (float *)mxGetData( plhs[0] );
	fOut = (int *)mxGetData( plhs[1] );

	reduce(target, v, f, nvin, nfin, verbose, nvout, nfout, vOut, fOut, &nv, &nf);

	memcpy( vOut+nv*1, vOut+nvout*1, nv*sizeof(float) );
	memcpy( vOut+nv*2, vOut+nvout*2, nv*sizeof(float) );
	mxSetM( plhs[0], nv);

	memcpy( fOut+nf*1, fOut+nfout*1, nf*sizeof(int) );
	memcpy( fOut+nf*2, fOut+nfout*2, nf*sizeof(int) );
	mxSetM( plhs[1], nf);
}
