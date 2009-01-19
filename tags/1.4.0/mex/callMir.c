/*--------------------------------------------------------------------
 *    The Mirone system:	@(#)callMir.c
 *
 *	Helper function to call mirone.exe but setting the path before the call
 *	This function replaces the mirone.bat batch in order to allow file
 *	associations. In order that the association works a environmental var
 *	MIRONE_HOME pointing to the mirone installation dir must exist. That is 
 *	so because when we double click a file the current dir is that of the
 *	file and not of the program, which leaves on the state of not knowing
 *	where we are.
 *
 *	Joaquim Luis	10-July-2007
 *	 
 *
 *--------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define CNULL		((char *)NULL)
#ifdef _WIN32	/* Start of Windows setup */
#define R_OK 04
#endif		/* End of Windows setup */

int main(int argc, char **argv) {

	char *path, *pato, *papato, *cd; 
	int	status, size_cd;

	size_cd = 128;

	path = getenv ("PATH");

	cd = (char *)malloc((size_t)size_cd);
	cd = getcwd (cd, (size_t)size_cd);
	papato = (char *) calloc ((size_t)(strlen(cd) + 12), (size_t)1);
	strcpy (papato, cd);
	strcat (papato, "\\mirone.exe");

	if (!access (papato, R_OK)) {		/* Found, so we are at home */
		pato = (char *) calloc ((size_t)(strlen(path) + 2*strlen(cd) + 18), (size_t)1);
	}
	else {		/* Called via file association. See if we have a MIRONE_HOME */
		cd = getenv ("MIRONE_HOME");
		pato = (char *) calloc ((size_t)(strlen(path) + 2*strlen(cd) + 18), (size_t)1);
	}
	strcpy (pato, "PATH=");
	strcat (pato, cd);
	strcat (pato, ";");
	strcat (pato, cd);
	strcat (pato, "\\bin\\win32;");
	strcat (pato, path);


	if (status = putenv(pato))
		fprintf(stderr, "callMir: Failure to set the environmental variable\n %s\n", pato);
	free((void *)pato);

	if (argc == 1)		/* Call a blank mirone window */
		system("mirone.exe");
	else {			/* Call with a file as argument */
		papato = (char *) calloc ((size_t)(size_cd), (size_t)1);
		strcpy (papato, "mirone.exe ");
		strcat (papato, argv[1]);
		system (papato);
		free((void *)papato);
	}

	return(1);
}
