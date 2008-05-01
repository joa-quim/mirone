/*=================================================================
 * existmex.c
 *
 * existmex works like the Matlab 'exist()' but can be used with compiled
 * codes, wich is not the case with exist(variable,'var') since it is, 
 * guess what, yes BUGGED. It always returns 0
 *
 *      Coffeeright (c) 2002-2007 by J. Luis
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
