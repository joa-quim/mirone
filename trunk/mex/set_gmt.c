/*--------------------------------------------------------------------
 *	$Id$
 *
 *	Copyright (c) 2004-2013 by J. Luis
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

/*--------------------------------------------------------------------
 *    The Mirone/GMT-system:	@(#)set_gmt.c
 *
 *	test what GMT instalation we have. The output is a structure with fields:
 *	version - '0' if no GMT was found. Else '4'
 *	high    - 'y' if the high resolution shore file was found
 *	full    - 'y' if the full resolution shore file was found
 *	intermediate - 'y' if the intermediate resolution shore file was found
 *	low     - 'y' if the low resolution shore file was found
 *	crude   - 'y' if the crude resolution shore file was found
 *
 *	Joaquim Luis	14-Juin-2005
 *	 
 *		25/06/07 J Luis, If called with 2 argsin prepend first one to the PATH (second is ignored)
 *		03/04/07 J Luis, Now sets GMTHOME as GMT_SHAREDIR when need use pre 4.2.0 coastlines
 *		22/03/07 J Luis, Change of strategy to deal with changes introduced in 4.2.0
 *				 - Renamed to set_gmt
 *				 - It is now based on transmission of the GMT_USERDIR path as input
 *				   this allow us to use our minimalist share files and thus render
 *				   Mirone totally independent of a GMT installation.
 *				 - Next GMT_SHAREDIR is searched (ver 4.2), if not found try to set 
 *				   it from the also searched GMT_HOME. This way, the search for the
 *				   coastlines is expected to work both for 4.2 and older versions. 
 *				   When no GMT is present one must use the coastlines.conf file
 *
 *		08/11/06 J Luis, Defaulting HOME to "C:\" was bugged
 *		29/10/06 J Luis, When called with one string input, use that string as argument
 *				 to putenv (use this because R13 doesn't have setenv). 
 *		28/10/06 J Luis, Rewrite to not depend on GMT dlls. 
 *				 Also uses the coastline files independence of GMT
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *
 *--------------------------------------------------------------------*/

#include "mex.h"
#include <string.h>
#include <stdio.h>

#define CNULL		((char *)NULL)
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#define R_OK 04
#define F_OK 00

#ifdef WIN32	/* Start of Windows setup */
#include <io.h>
#define DIR_DELIM '\\'	/* Backslash as directory delimiter */
#endif		/* End of Windows setup */

#ifndef DIR_DELIM
#define DIR_DELIM '/'
#endif

char *shore_getpathname (char *stem, char *path, char *GMT_USERDIR, char *GMT_SHAREDIR);
char *getsharepath (const char *subdir, const char *stem, const char *suffix, char *path, char *GMT_USERDIR, char *GMT_SHAREDIR);
void chop (char *string);

/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	char *this, *fieldnames[6] = {"version", "full", "high", "intermediate", "low", "crude"};
	char file[64];
	char *GMT_SHAREDIR = CNULL;
	char *GMT_HOMEDIR = CNULL;
	char *GMT_USERDIR = CNULL;
	char path[BUFSIZ], *pato, *papato; 
	int	status;
	mxArray *mxStr, *info_struct;

	if (nrhs >= 1 && mxIsChar(prhs[0])) {
		char	*envString; 
		envString = (char *)mxArrayToString(prhs[0]);
		if (nrhs == 2 && nlhs == 0) {		/* The set PATH case */
			this = getenv ("PATH");
			pato = (char *) mxCalloc ((size_t)(strlen(this) + strlen(envString) + 2), (size_t)1);
			strcpy (pato, envString);
			strcat (pato, this);
			if (status = putenv(pato))
				mexPrintf("SET_GMT: Failure to set the PATH environmental variable\n %s\n", pato);
			mxFree(pato);

		}
		else if (nrhs == 2 && nlhs == 1) {	/* Return the contents of the 'envstring' env var (for debug) */
			if ((this = getenv (envString)) != CNULL)
				mxStr = mxCreateString(this);
			else
				mxStr = mxCreateString("");

			plhs[0] = mxStr;
			return;
		}
		else if (nrhs == 1 && (status = putenv(envString)))
			mexPrintf("SET_GMT: Failure to set the environmental variable\n %s\n", envString);

		if (nlhs == 0)
			return;
	}

	if (nlhs == 0 && nrhs == 0) {
		mexPrintf("Usage: info = set_gmt;\nReturns the info structure with information on GMT version and costlines\n");
		mexPrintf("info = set_gmt('envstring');\nDo the same as above and set 'envstring' in the environment list\n");
		mexPrintf("       set_gmt('envstring');\nJust set 'envstring' in the environment list and return.\n");
		mexPrintf("       set_gmt('envstring',whatever);\nPrepends 'envstring' to the PATH.\n");
		mexPrintf(" env = set_gmt('envstring',whatever);\nGets the contents of 'envstring' (for debugging mostly).\n");
		return;
	}

	/* Structure field names. If we have an official GMT installation, only full & high 
	   will be tested (4-6 will be empty). Otherwise, we'll check also for the lower costlines. */
	info_struct = mxCreateStructMatrix (1, 1, 6, (const char **)fieldnames );


	if ((this = getenv ("GMT_SHAREDIR")) != CNULL) {	/* We have a 4.2 or greater version */
		GMT_SHAREDIR = (char *) mxCalloc ((size_t)(strlen (this) + 1), (size_t)1);
		strcpy (GMT_SHAREDIR, this);
		mxStr = mxCreateString("4");
	}
	else if ((this = getenv ("GMT5_SHAREDIR")) != CNULL) {	/* We have a 5.0 or greater version */
		char *sdir;
		GMT_SHAREDIR = (char *) mxCalloc ((size_t)(strlen (this) + 1), (size_t)1);
		strcpy (GMT_SHAREDIR, this);
		mxStr = mxCreateString("5");
		/* Now cheat GMT4 reseting GMT_SHAREDIR with the GMT5_SHAREDIR's value */
		sdir = (char *) mxCalloc ((size_t)(strlen (this) + 1), (size_t)1);
		sprintf (sdir, "GMT_SHAREDIR=%s", GMT_SHAREDIR);
		if (status = putenv(sdir))
			mexPrintf("SET_GMT: Failure trick to local reset the GMT_SHAREDIR env variable\n %s\n", sdir);
	}
	else if ((this = getenv ("GMTHOME")) != CNULL) {	/* We have a pre 4.2 version */
		char *sdir;
		GMT_SHAREDIR = (char *) mxCalloc ((size_t)(strlen (this) + 7), (size_t)1);
		sdir = (char *) mxCalloc ((size_t)(strlen (this) + 20), (size_t)1);
		/* GMT_SHAREDIR is used in this MEX to report what coastlines we have */
		sprintf (GMT_SHAREDIR, "%s%c%s", this, DIR_DELIM, "share");
		/* sdir will be used by our gmt dlls (& shoredump) */
		sprintf (sdir, "GMT_SHAREDIR=%s%c%s", this, DIR_DELIM, "share");
		if (status = putenv(sdir))
			mexPrintf("SET_GMT: Failure to set the sharedir environmental variable\n %s\n", sdir);
		mxStr = mxCreateString("4");
	}
	else						/* No GMT in sight */
		mxStr = mxCreateString("0");

	mxSetField(info_struct, 0, fieldnames[0], mxStr);


	/* -------------------------------------------------------------------- */
	if ((this = getenv ("HOME")) != CNULL) {		/* HOME was set */
		GMT_HOMEDIR = (char *) mxCalloc ((size_t)(strlen (this) + 1), (size_t)1);
 		strcpy (GMT_HOMEDIR, this);
	}
	else {
#ifdef WIN32
		/* Set HOME to C:\ under Windows */
		GMT_HOMEDIR = (char *) mxCalloc (4, (size_t)1);
		sprintf (GMT_HOMEDIR, "C:%c", DIR_DELIM);
#else
		mexPrintf ("SET_GMT: Warning, could not determine home directory!\n");
#endif
	}

	if ((this = getenv ("GMT_USERDIR")) != CNULL) {		/* GMT_USERDIR was set */
		GMT_USERDIR = (char *) mxCalloc ((size_t)(strlen (this) + 1), (size_t)1);
		strcpy (GMT_USERDIR, this);
	}
	else {				/* Use default pato for GMT_USERDIR (~/.gmt) */
		GMT_USERDIR = (char *) mxCalloc ((size_t)(strlen (GMT_HOMEDIR) + 6), (size_t)1);
		sprintf (GMT_USERDIR, "%s%c%s", GMT_HOMEDIR, DIR_DELIM, ".gmt");
	}
	if (access(GMT_USERDIR,R_OK)) GMT_USERDIR = CNULL;


	/* -------------------------------------------------------------------- */
	/* See if we have the Full definition files */
	sprintf (file, "binned_GSHHS_%c", 'f');
       	if (shore_getpathname (file, path, GMT_USERDIR, GMT_SHAREDIR))
		mxStr = mxCreateString("y");
	else
		mxStr = mxCreateString("n");
	mxSetField(info_struct, 0, fieldnames[1], mxStr);

	/* See if we have the High definition files */
	sprintf (file, "binned_GSHHS_%c", 'h');
	if (shore_getpathname (file, path, GMT_USERDIR, GMT_SHAREDIR))
		mxStr = mxCreateString("y");
	else
		mxStr = mxCreateString("n");
	mxSetField(info_struct, 0, fieldnames[2], mxStr);

	sprintf (file, "binned_GSHHS_%c", 'i');	/* Intermediate */
	if (shore_getpathname (file, path, GMT_USERDIR, GMT_SHAREDIR))
		mxStr = mxCreateString("y");
	else
		mxStr = mxCreateString("n");
	mxSetField(info_struct, 0, fieldnames[3], mxStr);

	sprintf (file, "binned_GSHHS_%c", 'l');	/* Intermediate */
	if (shore_getpathname (file, path, GMT_USERDIR, GMT_SHAREDIR))
		mxStr = mxCreateString("y");
	else
		mxStr = mxCreateString("n");
	mxSetField(info_struct, 0, fieldnames[4], mxStr);

	sprintf (file, "binned_GSHHS_%c", 'c');	/* Intermediate */
	if (shore_getpathname (file, path, GMT_USERDIR, GMT_SHAREDIR))
		mxStr = mxCreateString("y");
	else
		mxStr = mxCreateString("n");
	mxSetField(info_struct, 0, fieldnames[5], mxStr);


	plhs[0] = info_struct;
}

void chop (char *string) {
	/* Chops off any CR or LF at end of string and ensures it is null-terminated */
	int i, n;
	if (!string) return;	/* NULL pointer */
	if ((n = strlen (string)) == 0) return;	/* Empty string */
	for (i = n - 1; i >= 0 && (string[i] == '\n' || string[i] == '\r'); i--);
	i++;
	if (i >= 0) string[i] = '\0';	/* Terminate string */
}


char *getsharepath (const char *subdir, const char *stem, const char *suffix, char *path, char *GMT_USERDIR, char *GMT_SHAREDIR) {
	/* stem is the name of the file, e.g., GMT_CPT.lis
	 * subdir is an optional subdirectory name in the $GMT_SHAREDIR directory.
	 * suffix is an optional suffix to append to name
	 * path is the full path to the file in question
	 * Returns full pathname if a workable path was found
	 * Looks for file stem in current directory, $GMT_USERDIR (default ~/.gmt) and $GMT_SHAREDIR[/subdir]
	 */

	/* First look in the current working directory */

	sprintf (path, "%s%s", stem, suffix);
	if (!access (path, R_OK)) return (path);	/* Yes, found it in current directory */

	/* Do not continue when full pathname is given */

#ifdef WIN32
	if (stem[0] == '\\' || stem[1] == ':') return (NULL);
#else
	if (stem[0] == '/') return (NULL);
#endif

	/* Not found, see if there is a file in the user's GMT_USERDIR (~/.gmt) directory */

	if (GMT_USERDIR) {
		sprintf (path, "%s%c%s%s", GMT_USERDIR, DIR_DELIM, stem, suffix);
		if (!access (path, R_OK)) return (path);
	}

	/* Try to get file from $GMT_SHAREDIR/subdir */

	if (subdir) {
		sprintf (path, "%s%c%s%c%s%s", GMT_SHAREDIR, DIR_DELIM, subdir, DIR_DELIM, stem, suffix);
		if (!access (path, R_OK)) return (path);
	}

	/* Finally try file in $GMT_SHAREDIR (for backward compatibility) */

	sprintf (path, "%s%c%s%s", GMT_SHAREDIR, DIR_DELIM, stem, suffix);
	if (!access (path, R_OK)) return (path);

	return (NULL);	/* No file found, give up */
}

char *shore_getpathname (char *stem, char *path, char *GMT_USERDIR, char *GMT_SHAREDIR) {
	/* Prepends the appropriate directory to the file name
	 * and returns path if file is readable, NULL otherwise */
	 
	FILE *fp = NULL;
	char dir[BUFSIZ];

	/* This is the order of checking:
	 * 1. Is there a file coastline.conf in current directory, GMT_USERDIR or GMT_SHAREDIR[/coast]?
	 *    If so, use its information
	 * 2. Look in current directory, GMT_USERDIR or GMT_SHAREDIR[/coast] for file "name".
	 */
	 
	/* 1. First check for coastline.conf */
	
	if (getsharepath ("coast", "coastline", ".conf", path, GMT_USERDIR, GMT_SHAREDIR)) {

		/* We get here if coastline.conf exists - search among its directories for the named file */

		fp = fopen (path, "r");
		while (fgets (dir, BUFSIZ, fp)) {	/* Loop over all input lines until found or done */
			if (dir[0] == '#' || dir[0] == '\n') continue;	/* Comment or blank */
			chop (dir);		/* Chop off LF or CR/LF */
			sprintf (path, "%s%c%s%s", dir, DIR_DELIM, stem, ".cdf");
			if (!access (path, R_OK)) {
				fclose (fp);
				return (path);
			}
		}
		fclose (fp);
	}
	
	/* 2. Then check for the named file itself, but now (2013) we may have both .nc or .cdf */

	if (getsharepath ("coast", stem, ".cdf", path, GMT_USERDIR, GMT_SHAREDIR)) return (path);
	if (getsharepath ("coast", stem, ".nc",  path, GMT_USERDIR, GMT_SHAREDIR)) return (path);

	return (NULL);
}
