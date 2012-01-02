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
 * Create an array of zeros or any other value compatible with the selected class
 *
 * Class can be any of the  classes supported in Matlab, but only for REALS
 * That is, no COMPLEX arrays.
 *
 * This MEX file mimics and extends the behavior of both 'zeros' and 'ones'
 * functions of R14. Furthermore, it can be used in compiled code using the
 * true compiler of the R13 and previous releases.
 *
 * USAGE:
 *
 * alloc_mex(m, n,...,classname) or alloc_mex([m,n,...],classname) is an m-by-n-by-...
 * array of zeros of data type classname. classname is a string specifying the data
 * type of the output. classname can have the following values: 'double', 'single'
 * OR 'float', 'logical', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32',
 * 'int64', or 'uint64'
 *
 * alloc_mex(..., fill_value) fills the array with 'fill_values'. This means that
 * the MATLAB 'logical' function behavior can also be reprudeced
 *
 * EXAMPLES:
 *
 * x = alloc_mex(2,3,'int8');		-> Same as MATLAB R14 x = zeros(2,3,'int8'); 
 * x = alloc_mex(2,3,'int8',1);		-> Same as MATLAB R14 x = ones(2,3,'int8'); 
 * x = alloc_mex([2 3],'logical',1);	-> Same as MATLAB     x = logical([2 3],'true'); 
 *
 *
 * AUTHOR Joaquim Luis (jluis@ualg.pt)
 * 	  22-Feb-2006 
 *
 * Revision  1  25/10/2008 JL	Seting NaNs as 'fill_values' was returning zeros instead
 *
 */

#include "mex.h"
#include "string.h"

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double  *arg_dims, val = 0;
	int     i, argc, n_arg_char = 0, k = 0, n_dims, *dims, n_pts;
	char			*outI8, *classe, tmpI8;
	unsigned char		*outUI8, tmpUI8;
	short int		*outI16, tmpI16;
	unsigned short int	*outUI16, tmpUI16;
	int			*outI32, tmpI32;
	unsigned int		*outUI32, tmpUI32;
	long			*outI64, tmpI64;
	unsigned long		*outUI64, tmpUI64;
	float			*outF32, tmpF32;
	double			*outF64, tmpF64;
  
  	/*  check for proper number of arguments */
	if(nlhs != 1 || nrhs < 1) {
		mexPrintf ("usage: out = alloc_mex(m,n,p,...,class);\n");
		mexPrintf (" 	   out = alloc_mex([m n p ...],class);\n");
		mexPrintf (" 	   out = alloc_mex(...,fill_value);\n");
		mexPrintf (" 	 - where 'class' is any of:\n");
		mexPrintf (" 	    'logical' 'uint8', 'int8', 'uint16', 'int16', 'uint32'\n");
		mexPrintf (" 	    'int32', 'uint64', 'int64', 'float', 'single', 'double'\n");
		mexPrintf (" 	 - fill_value is a value used to fill the array (default = 0)\n");
		return;
	}

	argc = nrhs;
	for (i = 0; i < nrhs; i++) {		/* Check input to find how many arguments are of type char */
		if(mxIsChar(prhs[i])) {
			argc--;
			n_arg_char++;	/* Number of arguments that have a char type */
			k = i;
		}
	}
  

	if (n_arg_char > 1)
		mexErrMsgTxt("ALLOC_MEX ERROR: Class string can occur only once.\n");

	dims = mxCalloc(50,sizeof(int));

	n_dims = mxGetN(prhs[0]);	/* See if first arg is a scalar or a vector */
	arg_dims = mxGetPr(prhs[0]);
	if (n_dims > 1) {		/* Dimensions given in the [m n p ...] format */
		for (i = 0; i < n_dims; i++) {
			dims[i] = (int)arg_dims[i];
			if (dims[i] < 0) dims[i] = 0;
		}
	}
	else if (n_dims == 1 && k > 1) {
		n_dims = k;
		n_pts = 1;
		for (i = 0; i < k; i++) {
			dims[i] = (int)*mxGetPr(prhs[i]);
			if (dims[i] < 0) dims[i] = 0;
			n_pts *= dims[i];
		}
	}
	else {				/* We have a 2D square array */
		n_dims = 2;
		dims[0] = (int)arg_dims[0];
		dims[1] = (int)arg_dims[0];
		if (dims[0] < 0) dims[0] = 0;
		if (dims[1] < 0) dims[1] = 0;
		n_pts = dims[0] * dims[1];
	}

	if (n_arg_char == 1 && nrhs > k+1) {	/* Fill value was provided */
		val = *mxGetPr(prhs[nrhs-1]);
	}

	if (n_arg_char == 1) {		/* Find out in which data type was given the input array */
		classe = mxArrayToString(prhs[k]);
		if (!strncmp(classe,"logical",7)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxLOGICAL_CLASS, mxREAL);
			outUI8 = mxGetData(plhs[0]);
			if (val != 0) val = 1;		/* Here we only accept 0 or 1 */
			tmpUI8 = (unsigned char)val;
			for (i = 0; i < n_pts; i++) outUI8[i] = tmpUI8;
		}
		else if (!strncmp(classe,"int8",4)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxINT8_CLASS, mxREAL);
			outI8 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpI8 = (char)val;
				for (i = 0; i < n_pts; i++) outI8[i] = tmpI8;
			}
		}
		else if (!strncmp(classe,"uint8",5)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxUINT8_CLASS, mxREAL);
			outUI8 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpUI8 = (unsigned char)val;
				for (i = 0; i < n_pts; i++) outUI8[i] = tmpUI8;
			}
		}
		else if (!strncmp(classe,"int16",5)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxINT16_CLASS, mxREAL);
			outI16 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpI16 = (short int)val;
				for (i = 0; i < n_pts; i++) outI16[i] = tmpI16;
			}
		}
		else if (!strncmp(classe,"uint16",6)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxUINT16_CLASS, mxREAL);
			outUI16 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpUI16 = (unsigned short int)val;
				for (i = 0; i < n_pts; i++) outUI16[i] = tmpUI16;
			}
		}
		else if (!strncmp(classe,"uint32",6)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxUINT32_CLASS, mxREAL);
			outUI32 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpUI32 = (unsigned int)val;
				for (i = 0; i < n_pts; i++) outUI32[i] = tmpUI32;
			}
		}
		else if (!strncmp(classe,"int32",5)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxINT32_CLASS, mxREAL);
			outI32 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpI32 = (int)val;
				for (i = 0; i < n_pts; i++) outI32[i] = tmpI32;
			}
		}
		else if (!strncmp(classe,"uint64",6)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxUINT64_CLASS, mxREAL);
			outUI64 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpUI64 = (unsigned long int)val;
				for (i = 0; i < n_pts; i++) outUI64[i] = tmpUI64;
			}
		}
		else if (!strncmp(classe,"int64",5)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxINT64_CLASS, mxREAL);
			outI64 = mxGetData(plhs[0]);
			if (val != 0) {
				tmpI64 = (long int)val;
				for (i = 0; i < n_pts; i++) outI64[i] = tmpI64;
			}
		}
		else if (!strncmp(classe,"float",5) || !strncmp(classe,"single",6)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxSINGLE_CLASS, mxREAL);
			outF32 = (float *)mxGetData(plhs[0]);
			if (mxIsNaN(val) || val != 0) {
				tmpF32 = (float)val;
				for (i = 0; i < n_pts; i++) outF32[i] = tmpF32;
			}
		}
		else if (!strncmp(classe,"double",6)) {
			plhs[0] = mxCreateNumericArray (n_dims,dims,mxDOUBLE_CLASS, mxREAL);
			outF64 = mxGetPr(plhs[0]);
			if (mxIsNaN(val) || val != 0) {
				tmpF64 = (double)val;
				for (i = 0; i < n_pts; i++) outF64[i] = tmpF64;
			}
		}
		else
			mexErrMsgTxt("ALLOC_MEX ERROR: unknown class type\n");
	}
	else {		/* Default to doubles */
		plhs[0] = mxCreateNumericArray (n_dims,dims,mxDOUBLE_CLASS, mxREAL);
		outF64 = mxGetData(plhs[0]);
		if (val != 0) {
			tmpF64 = (double)val;
			for (i = 0; i < n_pts; i++) outF64[i] = tmpF64;
		}
	}

	mxFree(dims);
}
