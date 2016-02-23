/* 	$Id$ 
 *
 * Attempt to provide popen for matlab, to allow reading from a process.
 * 
 * 2004-09-28 dpwe@ee.columbia.edu  after PlayOn
 * 
 * Based on popenr.c from Daniel Ellis FEX contribution 13851
 * but changed in several ways. Namely output data type is the same as required.
 * No byte swapping.
 *
 * J. Luis 23-Feb-2016
 *
	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are
	met:

	    * Redistributions of source code must retain the above copyright
	      notice, this list of conditions and the following disclaimer.
	    * Redistributions in binary form must reproduce the above copyright
	      notice, this list of conditions and the following disclaimer in
	      the documentation and/or other materials provided with the distribution
	    * Neither the name of the Columbia University nor the names
	      of its contributors may be used to endorse or promote products derived
	      from this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
	POSSIBILITY OF SUCH DAMAGE.

 */
 
#include    <stdio.h>
#include    <math.h>
#include    <ctype.h>
#include    "mex.h"

#ifdef WIN32
#	define popen _popen
#	define pclose _pclose
#endif

#define FILETABSZ 32
static FILE *filetab[FILETABSZ];
static int filetabix = 0;

int findfreetab() {
	/* find an open slot in the file table */
	int i;
	for (i = 0; i < filetabix; i++) {
		if (filetab[i] == NULL)		/* NULL entries are currently unused */
			return i;
	}
	if (filetabix < FILETABSZ) {
		i = filetabix;
		/* initialize it */
		filetab[i] = NULL;
		filetabix++;
		return i;
	}
	return -1;			/* out of space */
}

/* ---------------------------------------------------------------------- */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int     i, ix, nRows, nCols, nBands, ngot, sz = 1, dims[3];
	size_t  nPts;
	double *pd;
	mxClassID classid = mxUINT8_CLASS;
	FILE *f = NULL;
	mxArray *rslt = NULL;

	if (nrhs < 1) {
		mexPrintf("popenr  Y = popenr(X[,N[,F]])  Open and read an external process\n");
		mexPrintf("        When X is a string, that string is executed as a new process\n");
		mexPrintf("        and Y is returned as a handle (integer) to that stream.\n");
		mexPrintf("        Subsequent calls to popenr(Y,N) with that handle return the next N\n");
		mexPrintf("        values read from the standard ouptut of converted to Matlab values\n");
		mexPrintf("        according to the format string F (default: char).\n");
		mexPrintf("        If N is a row vector with [nRows nCols nBands] and type 'char' then data is\n");
		mexPrintf("        internaly transposed to Matlab orientation and output is a [MxNxnBands] array.\n");
		mexPrintf("        A call with N set to -1 means to close the stream.\n");
		return;
	}
	/* look at the data */
	/* Check to be sure input argument is a string. */
	if ((mxIsChar(prhs[0]))) {
		/* first argument is string - opening a new command */
		char *cmd;
		int tabix = findfreetab();
		
		if (tabix < 0) mexErrMsgTxt("Out of file table slots.");

		cmd = mxArrayToString(prhs[0]);
		f = popen(cmd, "r");
		mxFree(cmd);
		if (f == NULL) mexErrMsgTxt("Error running external command.");

		/* else have a new command path - save the handle */
		filetab[tabix] = f;
		/* return the index */
		if (nlhs > 0) {
			mxArray *rslt = mxCreateDoubleMatrix(1,1, mxREAL);
			plhs[0] = rslt;
			pd = mxGetPr(rslt);
			*pd = (double)tabix;
		}
		return;
	}

	if (nrhs < 2) mexErrMsgTxt("apparently accessing handle, but no N argument");
	
	/* get the handle */
	if (mxGetN(prhs[0]) == 0) mexErrMsgTxt("handle argument is empty");

	ix = (int)*mxGetPr(prhs[0]);
	if (ix < filetabix) f = filetab[ix];
	if (f == NULL) mexErrMsgTxt("invalid handle");

	/* how many items required? */
	if (mxGetN(prhs[1]) == 0) mexErrMsgTxt("length argument is empty");

	nBands = 1;
	pd = mxGetPr(prhs[1]);
	if (mxGetN(prhs[1]) == 1) {
		nRows = (int)pd[0];		nCols = 1;
	}
	else if (mxGetN(prhs[1]) == 2){
		nRows = (int)pd[0];		nCols = (int)pd[1];
	}
	else if (mxGetN(prhs[1]) == 3){
		nRows = (int)pd[0];		nCols = (int)pd[1];		nBands = (int)pd[2];
	}
	else
		mexErrMsgTxt("More than 3 dimensional data is not supported");
	
	/* maybe close */
	if (nRows < 0) {
		pclose(f);
		filetab[ix] = NULL;
		return;
	}

	/* what is the format? */
	if (nrhs > 2) {
		char *fmtstr = NULL;

		if (!mxIsChar(prhs[2])) mexErrMsgTxt("format arg must be a string");

		fmtstr = mxArrayToString(prhs[2]);
		if (!strcmp(fmtstr, "uint8") || !strcmp(fmtstr, "char")) {
			classid = mxUINT8_CLASS;	sz = 1;
		}
		else if (!strcmp(fmtstr, "int8")) {
			classid = mxINT8_CLASS;		sz = 1;
		}
		else if (!strcmp(fmtstr, "int16")) {
			classid = mxINT16_CLASS;	sz = 2;
		}
		else if (!strcmp(fmtstr, "uint16")) {
			classid = mxUINT16_CLASS;	sz = 2;
		}
		else if (!strcmp(fmtstr, "int32")) {
			classid = mxINT32_CLASS;	sz = 4;
		}
		else if (!strcmp(fmtstr, "uint32")) {
			classid = mxUINT32_CLASS;	sz = 4;
		}
		else if (!strcmp(fmtstr, "float") || !strcmp(fmtstr, "single")) {
			classid = mxSINGLE_CLASS;	sz = 4;
		}
		else if (!strcmp(fmtstr, "double")) {
			classid = mxDOUBLE_CLASS;	sz = 8;
		}
		else
			mexErrMsgTxt("unrecognized format");
	}

	/* do the read */
	dims[0] = nRows;		dims[1] = nCols;		dims[2] = nBands;
	rslt = mxCreateNumericArray((nBands == 1 ? 2 : 3), dims, classid, mxREAL);
	nPts = (size_t)nRows * nCols * nBands;
	if (nCols == 1)
		ngot = fread(mxGetData(rslt), sz, nPts, f);
	else {
		int k, row, col, band;
		size_t nXY = (size_t)nRows * nCols;
		void *pv = mxGetData(rslt);
		if (classid == mxUINT8_CLASS) {
			unsigned char *tmp = mxMalloc(nCols * nBands * sz);
			unsigned char *p = (unsigned char *)pv;
			for (row = 0; row < nRows; row++) {
				fread(tmp, sz, nCols * nBands, f);		/* Read a row of nCols by nBands bytes of data */
				for (col = 0; col < nCols; col++) {
					for (band = 0; band < nBands; band++) {
						p[row + col*nRows + band*nXY] = tmp[band + col*nBands];
					}
				}
			}
			mxFree(tmp);
		}
	}


	if (ngot < nPts) {
		void *pv = mxGetData(rslt);
		mxRealloc(pv, ngot * sz);
	}

	plhs[0] = rslt;
}
