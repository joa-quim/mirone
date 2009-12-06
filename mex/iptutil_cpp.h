//
// Copyright 1993-2005 The MathWorks, Inc.
// $Revision: 1.1.8.3 $
//

//Utilities for C++ IPT mex functions

#include "mwsize.h"		// For R13 mainly (or only)

#ifndef IPTUTIL_CPP_H
#define IPTUTIL_CPP_H

//////////////////////////////////////////////////////////////////////////////
// MACROS
//////////////////////////////////////////////////////////////////////////////

#define DUMMY 0 //Part of the workaround for MS VC++ limitations

//Useful in formulating an error message for uninitialized class variables
#define ERR_STRING(a,b) "Class variable " a " must be" \
" properly initialized prior to calling " b ".\n"

#endif // IPTUTIL_CPP_H
