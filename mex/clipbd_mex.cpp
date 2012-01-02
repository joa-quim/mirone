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
 * This is very crude program to interface with the Windows clipboard.
 * Currently, it should only be used to put a numeric matrix string
 * as prepared by the MATLAB mat2clip .m function. Hopefully, it will
 * evolve to recognize other clipboard formats but meanwhile that's 
 * what it can do. Quite useful though, to paste directly into Excel.
 *
 *
 * Most of the code comes from this tutorial.
 * http://www.codeproject.com/KB/clipboard/archerclipboard1.aspx 
 *
 * To compile (but adjust to your own path and compiler) with VC6
 * mex clipbd_mex.cpp C:\programs\VisualStudio\VC98\Lib\USER32.LIB 
 * C:\programs\VisualStudio\VC98\Lib\GDI32.LIB -DWIN32 -O 
 *
 *
 * Date:	07-Sep-2008
 *--------------------------------------------------------------------*/

#include <windows.h>
#include "mex.h"

struct clip_matFormat {
	int nx, ny;
	float *z_s;
	double *z_d;
	double *head;
};


/* interface between MATLAB and the C function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* declare variables */
	int ns, n2, m, n, put = 0;
	char *str;
	char *pchData;
	HGLOBAL hClipboardData;

	if ( !OpenClipboard(NULL) )
		mexErrMsgTxt("CLIPBD ERROR: Error while opening clipboard.\n");

	put = nrhs;

	EmptyClipboard();

	if (put && mxIsChar(prhs[0])) {
		str = (char *)mxArrayToString(prhs[0]);
		ns = strlen(str) + 1;

		hClipboardData = GlobalAlloc(GMEM_DDESHARE, ns);

		pchData = (char*)GlobalLock(hClipboardData);
		  
		strcpy(pchData, str);
		  
		GlobalUnlock(hClipboardData);
		  
		SetClipboardData(CF_TEXT,hClipboardData);
		  
	}
	else if (put && mxIsNumeric(prhs[0])) {
		UINT format = RegisterClipboardFormat("000_FMT");
		int is_single = 0, is_double = 0;
		clip_matFormat data; 

		if (mxIsSingle(prhs[0]))
			is_single = 1;

		else if (mxIsDouble(prhs[0]))
			is_double = 1;

		data.ny = mxGetM(prhs[0]);	data.nx = mxGetN(prhs[0]);
		if (nrhs == 2) {
			if ( mxGetM(prhs[1]) == 9 || mxGetN(prhs[1]) == 9 )
				data.head = mxGetPr(prhs[1]);
		}

		if (is_single)
			data.z_s = (float *)mxGetData(prhs[0]);
		else
			data.z_d = mxGetPr(prhs[0]);

		hClipboardData = GlobalAlloc(GMEM_DDESHARE, sizeof(clip_matFormat));
		clip_matFormat *buffer = (clip_matFormat *)GlobalLock(hClipboardData);

		/* put the data into that memory */
		*buffer = data;

		/* Put it on the clipboard */
		GlobalUnlock(hClipboardData);
		SetClipboardData(format,hClipboardData);
		mexPrintf("FFSSS %d\n;", format);
		mexPrintf("FDS %d\n", IsClipboardFormatAvailable(format));
	}
	else if (!put) {
		UINT format = RegisterClipboardFormat("000_FMT");

		/* No idea why it is not working. The 'format' is no
		   longer recognized as a registered one .*/
		mexPrintf("FDS %d\n", IsClipboardFormatAvailable(format));
		if ( !IsClipboardFormatAvailable(format) )
			mexErrMsgTxt("fd-se.\n");

		clip_matFormat data; 
		HANDLE hData = GetClipboardData(format);
		clip_matFormat *buffer = (clip_matFormat *)GlobalLock(hData);

		/* make a local copy */
		data = *buffer;

		mexPrintf("mmm %d\n", (buffer->ny));

		GlobalUnlock( hData );
	}

	CloseClipboard();
}
