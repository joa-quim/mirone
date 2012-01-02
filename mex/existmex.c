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

/*=================================================================
 * existmex.c
 *
 * existmex works like the Matlab 'exist()' but can be used with compiled
 * codes, wich is not the case with exist(variable,'var') since it is, 
 * guess what, yes BUGGED. It always returns 0
 *
 *================================================================*/

#include "mex.h"

void mexFunction( int nargout, mxArray *pargout[], int nargin, const mxArray *pargin[] ) {
	mxArray *rhs[2], *x;
	const char *var;
	int   status;

	if (!mxIsChar(pargin[0]))
    		mexErrMsgTxt("Function 'existmex' input arg(s) must be char strings.\n");

	if (nargin == 2) {
		if (!mxIsChar(pargin[1]))
    			mexErrMsgTxt("Function 'existmex' second arg must be a char string.\n");

		rhs[1] = mxCreateString( (char *)mxArrayToString(pargin[1]) );
	}
	else
		rhs[1] = mxCreateString("var");

  
	var = (char *)mxArrayToString(pargin[0]);

	rhs[0] = mxCreateString(var);
	status = mexCallMATLAB(1,&x,2, rhs, "exist");
	pargout[0] = x;
    
	mxFree(var);
}
