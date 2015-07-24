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

/*--------------------------------------------------------------------
 *    The Mirone/GMT-system:	@(#)test_gmt.c
 *
 *	test what GMT instalation we have. The output is a structure with fields:
 *	version - '0' if no GMT was found. If we have a GMT than the last 3 fields are
 *		      empty because it is assumed that those files exist under the GMT dir
 *	high    - 'y' if the high resolution shore file was found
 *	full    - 'y' if the full resolution shore file was found
 *	intermediate - 'y' if the intermediate resolution shore file was found
 *	low     - 'y' if the low resolution shore file was found
 *	crude   - 'y' if the crude resolution shore file was found
 *
 *	Joaquim Luis	14-Juin-2005
 *	 
 *		04/06/06 J Luis, Updated to compile with version 4.1.3
 *		28/10/06 J Luis, Rewrite to not depend on GMT dlls. 
 *				 Also uses the coastline files independence of GMT
 *		29/10/06 J Luis, When called with one string input, use that string as argument
 *				 to putenv (use this because R13 doesn't have setenv). 
 *		08/11/06 J Luis, Defaulting HOME to "C:\" was bugged
 *
 *--------------------------------------------------------------------*/

#include "mex.h"
#include <string.h>

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
#define DIR_DELIM '\\'	/* Backslash as directory delimiter */
#include <io.h>
#endif		/* End of Windows setup */

#ifndef DIR_DELIM
#define DIR_DELIM '/'
#endif


int getpathname (char *name, char *path, char *GMTHOME, char *GMT_HOMEDIR, char *GMT_USERDIR);
int shore_conffile (char *name, char *dir, char *path);
void chop (char *string);

/* Matlab Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	char *this, *fieldnames[6] = {"version", "full", "high", "intermediate", "low", "crude"};
	char file[64];
	char *GMTHOME = CNULL;
	char *GMT_HOMEDIR = CNULL;
	char *GMT_USERDIR = CNULL;
	char path[BUFSIZ]; 
	mxArray *mxStr, *info_struct;

	if (nrhs == 1 && mxIsChar(prhs[0])) {
		char	*envString; 
		int	status;
		envString = (char *)mxArrayToString(prhs[0]);
		if (status = putenv(envString))
			mexPrintf("TEST_GMT: Failure to set the environmental variable\n %s\n", envString);

		if (nlhs == 0)
			return;
	}

	if (nlhs == 0 && nrhs == 0) {
		mexPrintf("Usage: info = test_gmt;\nReturns the info structure with information on GMT version and costlines\n");
		mexPrintf("info = test_gmt('envstring');\nDo the same as above and set 'envstring' in the environment list\n");
		mexPrintf("       test_gmt('envstring');\nJust set 'envstring' in the environment list and return.\n");
		return;
	}

	/* Structure field names. If we have an official GMT installation, only full & high 
	   will be tested (4-6 will be empty). Otherwise, we'll check also for the lower costlines. */
	info_struct = mxCreateStructMatrix (1, 1, 6, (const char **)fieldnames );

	if ((this = getenv ("GMTHOME")) == CNULL) {	/* No official GMT in sight */
		mxStr = mxCreateString("0");
	}
	else {
		GMTHOME = (char *) mxCalloc ((size_t)(strlen (this) + 1), (size_t)1);
		strcpy (GMTHOME, this);
		mxStr = mxCreateString("4");
	}
	mxSetField(info_struct, 0, fieldnames[0], mxStr);

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
		mexPrintf ("Warning: Could not determine home directory!\n");
#endif
	}

	if ((this = getenv ("GMT_USERDIR")) != CNULL) {		/* GMT_USERDIR was set */
		GMT_USERDIR = (char *) mxCalloc ((size_t)(strlen (this) + 1), (size_t)1);
		strcpy (GMT_USERDIR, this);
	}
	else if (GMT_HOMEDIR) {			/* Use default pato for GMT_USERDIR (~/.gmt) */
		GMT_USERDIR = (char *) mxCalloc ((size_t)(strlen (GMT_HOMEDIR) + 6), (size_t)1);
		sprintf (GMT_USERDIR, "%s%c%s", GMT_HOMEDIR, DIR_DELIM, ".gmt");
	}
	if (access(GMT_USERDIR,R_OK)) GMT_USERDIR = CNULL;

	/* See if we have the Full definition files */
	sprintf (file, "binned_GSHHS_%c.cdf", 'f');
       	if (getpathname (file, path, GMTHOME, GMT_HOMEDIR, GMT_USERDIR))
		mxStr = mxCreateString("y");
	else
		mxStr = mxCreateString("n");
	mxSetField(info_struct, 0, fieldnames[1], mxStr);

	/* See if we have the High definition files */
	sprintf (file, "binned_GSHHS_%c.cdf", 'h');
       	if (getpathname (file, path, GMTHOME, GMT_HOMEDIR, GMT_USERDIR))
		mxStr = mxCreateString("y");
	else
		mxStr = mxCreateString("n");
	mxSetField(info_struct, 0, fieldnames[2], mxStr);

	if (!GMTHOME) {
		/* An official GMT installation does not exist. So check if we have the lower costlines */
		sprintf (file, "binned_GSHHS_%c.cdf", 'i');	/* Intermediate */
       		if (getpathname (file, path, GMTHOME, GMT_HOMEDIR, GMT_USERDIR))
			mxStr = mxCreateString("y");
		else
			mxStr = mxCreateString("n");
		mxSetField(info_struct, 0, fieldnames[3], mxStr);

		sprintf (file, "binned_GSHHS_%c.cdf", 'l');	/* Intermediate */
       		if (getpathname (file, path, GMTHOME, GMT_HOMEDIR, GMT_USERDIR))
			mxStr = mxCreateString("y");
		else
			mxStr = mxCreateString("n");
		mxSetField(info_struct, 0, fieldnames[4], mxStr);

		sprintf (file, "binned_GSHHS_%c.cdf", 'c');	/* Intermediate */
       		if (getpathname (file, path, GMTHOME, GMT_HOMEDIR, GMT_USERDIR))
			mxStr = mxCreateString("y");
		else
			mxStr = mxCreateString("n");
		mxSetField(info_struct, 0, fieldnames[5], mxStr);

	}

	plhs[0] = info_struct;
}

int getpathname (char *name, char *path, char *GMTHOME, char *GMT_HOMEDIR, char *GMT_USERDIR) {
	/* Prepends the appropriate directory to the file name and returns TRUE if file is readable. */
	int found = FALSE;
	char dir[BUFSIZ];

	/* This is the order of checking:
	 * 1. Is there a GMT_USERDIR with a coastline.conf in it? If so use its information
	 * 2. Look in GMTHOME/share/coast
	 * 3. Look in GMTHOME/share (backward check)
	 * 4. Look for GMTHOME/share/coastline.conf and use its information
	 * 5. Give up
	 */
	 
	/* 1. First check the $GMT_USERDIR directory */
	
	if (GMT_USERDIR) {		/* See if we have a GMT_USERDIR with a coastline.conf */
		sprintf (dir, "%s%ccoastline.conf", GMT_USERDIR, DIR_DELIM);
		found = shore_conffile (name, dir, path);
		if (found) return (TRUE);
	}
	
	/* 2. Then check the $GMTHOME/share/coast directory */

	sprintf (path, "%s%cshare%ccoast%c%s", GMTHOME, DIR_DELIM, DIR_DELIM, DIR_DELIM, name);
	if (!access (path, R_OK)) return (TRUE);	/* File exists and is readable, return with name */

	/* File was not readable.  Now check if it exists */

	if (!access (path, F_OK))  { /* Kicks in if file is there, meaning it has the wrong permissions */
		mexPrintf ("%s: Error: does not have permission to open %s!\n", "test_gmt", path);
		return (0);
	}
	
	/* 3. Nothing in share/coast; do a backwards-compatible check in the $GMTHOME/share directory */

	sprintf (path, "%s%cshare%c%s", GMTHOME, DIR_DELIM, DIR_DELIM, name);
	if (!access (path, R_OK)) return (TRUE);	/* File exists and is readable, return with name */

	/* File was not readable.  Now check if it exists */

	if (!access (path, F_OK))  { /* Kicks in if file is there, meaning it has the wrong permissions */
		mexPrintf ("%s: Error: does not have permission to open %s!\n", "test_gmt", path);
		return (0);
	}

	/* 4. File is not there.  Thus, we check if a coastline.conf file exists
	 * It is not an error if we cannot find the named file, only if it is found
	 * but cannot be read due to permission problems */

	sprintf (dir, "%s%cshare%ccoastline.conf", GMTHOME, DIR_DELIM, DIR_DELIM);
	found = shore_conffile (name, dir, path);
	
	return (found);
}

int shore_conffile (char *name, char *dir, char *path) {
	int found = FALSE;
	FILE *fp;
	
	/* Given the dir of a coastline.conf file, look to see if it can be found/read,
	 * and if so return the directory given.  */
	 
	if (!access (dir, F_OK))  { /* File exists... */
		if (access (dir, R_OK)) {	/* ...but cannot be read */
			mexPrintf ("%s: Error: does not have permission to open %s!\n", "test_gmt", dir);
			return (0);
		}
	}
	else 	/* There is no coastline.conf file to use; we're out of luck */
		return (FALSE);


	/* We get here if coastline.conf exists - search among its directories for the named file */
	if ((fp = fopen (dir, "r")) == NULL) {	/* This shouldn't be necessary, but cannot hurt */
		mexPrintf ("%s: Error: Cannot open configuration file %s\n", "test_gmt", dir);
		return (0);
	}

	found = FALSE;
	while (!found && fgets (dir, BUFSIZ, fp)) {	/* Loop over all input lines until found or done */
		if (dir[0] == '#' || dir[0] == '\n') continue;	/* Comment or blank */

		chop (dir);		/* Chop off LF or CR/LF */
		sprintf (path, "%s%c%s", dir, DIR_DELIM, name);
		if (!access (path, F_OK)) {	/* TRUE if file exists */
			if (!access (path, R_OK)) 	/* TRUE if file is readable */
				found = TRUE;
			else {
				mexPrintf ("%s: Error: does not have permission to open %s!\n", "test_gmt", path);
				return (0);
			}
		}
	}

	fclose (fp);
	
	return (found);
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

