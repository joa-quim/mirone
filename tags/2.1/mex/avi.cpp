//#include "avi.h"
//#include "avistandalone.h"

/* 
 * avi.cpp
 *
 * This MEX-file is an interface to the Microsoft AVIFILE routines.  
 * There are three modes of operation for this MEX-file; 'open',
 * 'addframe' and 'close'.  The syntax for calling this MEX-file is as follows:
 *
 * ID = AVI('open',FILENAME) initializes the interface to the AVIFILE routines
 *		and returns ID, a unique integer ID corresponding to the open file
 *		FILENAME.  The MATLAB M-file AIVFILE should be used to call this routine.
 *
 *  AVI('addframe',FRAME,BITMAPINFO,FRAMENUM,FPS,QUALITY,ID, STREAMNAME) adds FRAME 
 *		number FRAMENUM to the stream in the AVI file represented by ID.  The 
 *		BITMAPINFO is the bitmapheader structure of the AVIFILE object.  
 *		The FPS (frames per second) and QUALITY parameters are required by the
 *		AVIFILE routines. STREAMNAME is a string describing the video stream.
 *      The MATLAB M-file ADDFRAME should be used to call this routine.
 *
 *  AVI('close',ID) finishes writing the AVI file represented by ID.  This will 
 *		close the stream and file handles. This routine should be called by
 *		the MATLAB M-file AVICLOSE.
 *
 *	aviutils.c contains all the routines called by AVI.  
 */

/* 
 * Copyright 1984-2001 The MathWorks, Inc. 
 * $Revision: 1.3 $  $Date: 2003/03/05 06:47:58 $
*/

#include <windows.h>
#include <vfw.h>
#include "mex.h"

#ifdef __cplusplus
   extern "C"
   {
#endif

void aviopen(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void addframe(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void aviclose(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

extern CRITICAL_SECTION AVI_NODE_ACCESS_SECTION;

typedef struct node {
	int number;              /* Unique identifier */
	PAVIFILE pfile;          /* File handle */
	PAVISTREAM psCompressed; /* Stream handle */
	struct node *next;
	struct node *previous;
    CRITICAL_SECTION section; /* Used to prevent multiple thread entries. */
} NodeType;

void cleanList(void);
void InitNode(void);
int addNodetoList(PAVIFILE,PAVISTREAM);
NodeType *FindNodeInList(int);
void deleteNodeFromList(int);
int openFile(const char *filename);
void closeFile(int identifier);
void exitAVI(void);
char ** getVidCodecs(void);
int getNumVidCodecs(void);
void writeFrame(int identifier, char *StreamName, int FrameNumber, unsigned char *frameData, char* compression, long width, long height, int bitcount, int ImageSize, int ClrUsed, unsigned char *colormap, int Quality, int FramesPerSecond, int KeyFrameRate);

#ifdef __cplusplus
    }   /* extern "C" */
#endif

extern "C"
{

static char rcsid[] = "$Id: avi.cpp,v 1.3 2003/03/05 06:47:58 cportal Exp $";

static int mexAtExitIsSet = 0;
static bool GLOBAL_SECTION_INITIALIZED = false;

static void exitMex() {
    try
    {
        exitAVI();
    }catch (char *warnmsg)
    {
        mexWarnMsgTxt(warnmsg);
    }

    // Destroy the global critical section.
    if (GLOBAL_SECTION_INITIALIZED){
        DeleteCriticalSection(&AVI_NODE_ACCESS_SECTION);
        GLOBAL_SECTION_INITIALIZED = false;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


    char *mode;
    void (*opMode)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

    if( (nrhs < 2) && (nrhs > 8) )
        mexErrMsgTxt("Invalid number of input arguments.");

    // Initialize the global critical section.
    if (!GLOBAL_SECTION_INITIALIZED){
        InitializeCriticalSection(&AVI_NODE_ACCESS_SECTION);
        GLOBAL_SECTION_INITIALIZED = true;
    }

    if(!mexAtExitIsSet) {
            InitNode();
            mexAtExit(exitMex);
            mexAtExitIsSet = 1;
        }

    if ( mxIsChar(prhs[0]) ) 
        mode = mxArrayToString(prhs[0]);
    else
        mexErrMsgTxt("First input to AVI must be a string.");

    if (mode == NULL)
        mexErrMsgTxt("Memory allocation failed.");
	
    if(!strcmp(mode,"open"))
    	opMode = aviopen;
    else if(!strcmp(mode,"addframe"))
    	opMode = addframe;
	else if(!strcmp(mode,"close"))
    	opMode = aviclose;
	else
		mexErrMsgTxt("Unrecognized mode for AVI.");

	mxFree(mode);
		
	(*opMode)(nlhs, plhs, nrhs, prhs);
}

/*
 * AVIOPEN opens the AVI file and returns a unique file ID in plhs
*/

void aviopen(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    char *filename;
    int id;
    
    /* Validate inputs */
    if (nrhs != 2)
        mexErrMsgTxt("Invalid number of inputs to AVIOPEN.");
    
    if ( !mxIsChar(prhs[1]) )
        mexErrMsgTxt("The first input to AVIOPEN must be the filename.");

    filename = mxArrayToString(prhs[1]);
    if( filename == NULL )
        mexErrMsgTxt("Memory allocation failed.");

    try
    {
        id = openFile(filename);
    }
    catch(char *errmsg)
    {
        mexErrMsgTxt(errmsg);
    }
    
    /* Return a unique number */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    *((double *) mxGetPr(plhs[0])) = id;
    mxFree(filename);
}

/* 
 * ADDFRAME  Adds a frame to the video stream in the avifile opened by aviopen. 
 *	 If this is the first stream, the stream is created, otherwise it is 
 *	 appended to the end.  
 */
 
void addframe(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	uint8_T KeyFrameRate, *frame, *colormap;
	mxArray *mxWidth, *mxHeight;
        long width, height;
	char *compression, *StreamName;
	int FrameNumber, identifier, ImageSize, FramesPerSecond, bitcount, ClrUsed, Quality;


	/* Validate Inputs */
	if (!mxIsUint8(prhs[1]))
		mexErrMsgTxt("Frame must be a uint8 matrix.");
	if(!mxIsStruct(prhs[2]))
		mexErrMsgTxt("Third input to addframe must be a structure representing BITMAPINFO.");
	if(!mxIsNumeric(prhs[3]))
		mexErrMsgTxt("Fourth input to addframe must be the frame number.");
	if(!mxIsNumeric(prhs[4]))
		mexErrMsgTxt("Fifth input to addframe must be the frames per second.");
	if(!mxIsNumeric(prhs[5]))
		mexErrMsgTxt("Sixth input to addframe must be the Quality.");
	if(!mxIsNumeric(prhs[6]))
		mexErrMsgTxt("Seventh input to addframe must be the Id number.");
	if(!mxIsChar(prhs[7]))
		mexErrMsgTxt("Eighth input to addframe must be a stream name.");
	if(!mxIsNumeric(prhs[8]))
		mexErrMsgTxt("Last input to addframe must be the key frame rate.");

	/* Find the stream handle for this particular file. */
	identifier = (int) *mxGetPr(prhs[6]);

	compression = mxArrayToString(mxGetField(prhs[2],0,"biCompression"));
	if (compression == NULL)
		mexErrMsgTxt("Failed to allocate memory.");

	frame = (uint8_T *)mxGetData(prhs[1]);
	FrameNumber = (int) *mxGetPr(prhs[3]);

	KeyFrameRate = (int) *mxGetPr(prhs[8]);
        
	mxWidth = mxGetField(prhs[2],0,"biWidth");
	mxHeight = mxGetField(prhs[2],0,"biHeight");
        width = (long) *mxGetPr(mxWidth);
        height = (long) *mxGetPr(mxHeight);
        bitcount = (int) *mxGetPr(mxGetField(prhs[2],0,"biBitCount"));
        ImageSize = (int) *mxGetPr(mxGetField(prhs[2],0,"biSizeImage"));
        ClrUsed = (int) *mxGetPr(mxGetField(prhs[2],0,"biClrUsed"));
        colormap = (uint8_T *)mxGetData(mxGetField(prhs[2],0,"Colormap"));
        StreamName = mxArrayToString(prhs[7]);
        Quality = (int) *mxGetPr(prhs[5]);
	if (StreamName == NULL)
            mexErrMsgTxt("Memory allocation failed.");

        FramesPerSecond = (int) *mxGetPr(prhs[4]);

        try {
            writeFrame(identifier, StreamName, FrameNumber, frame, compression, width, height, bitcount, ImageSize, ClrUsed, colormap, Quality, FramesPerSecond, KeyFrameRate);
	}
        catch(char *errmsg) {
            mexErrMsgTxt(errmsg);
	}

	mxFree(compression);
	mxFree(StreamName);
}

/* 
 * AVICLOSE close the stream and the file then remove this information from
 *	the list.
 */

void aviclose(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int identifier;

	/*Need to find the stream handle for this particular file. */
	identifier = (int) *mxGetPr(prhs[1]);
        try {
            closeFile(identifier);
	}
        catch(char *errmsg) {
            mexErrMsgTxt(errmsg);
	}
}


static NodeType OpenFiles;

/* InitNode 
 *	Initialize the list.
 */
void InitNode(void) {
	OpenFiles.number = -1;
	OpenFiles.pfile = NULL;
	OpenFiles.psCompressed = NULL;
	OpenFiles.next = NULL;
	OpenFiles.previous = NULL;
}

/* addNodetoList
 *	Insert a new node at the begninning of the OpenFiles list.  
 *
 *	Inputs:  pfile        - AVI file identifier
			 psCompressed - Video stream identifier
 *   
 *	Outputs: A unique identification number
 */
int addNodetoList(PAVIFILE pfile, PAVISTREAM psCompressed) {
	NodeType *NewNode = NULL;
	NodeType *ListHead;
	ListHead = &OpenFiles;

	NewNode = (NodeType *) malloc(sizeof(NodeType));
	if (NewNode != NULL) {
            /* mexMakeMemoryPersistent((void *) NewNode); */
		NewNode->number = (ListHead->next == NULL)?1:ListHead->next->number+1;
		NewNode->pfile = pfile;
		NewNode->psCompressed = psCompressed;
		NewNode->next = NULL;
		NewNode->previous = NULL;
		InitializeCriticalSection(&NewNode->section);
	}
    
	else
        	throw("Out of memory in AVI MEX-file");
    
	/* Add node to begining of list */
	NewNode ->next = ListHead->next;
	NewNode->previous = ListHead;
	ListHead->next = NewNode;
	if (NewNode->next != NULL)
	        NewNode->next->previous = NewNode;

	return(NewNode->number);
}

/* FindNodeInList
 *	Find a node in the OpenFiles list given a unique identification number.
 *
 *	Inputs:  number - a unique number identifying the node
 *	Outputs: pointer to the node found
 */
NodeType *FindNodeInList(int number) {
	NodeType *current;
	int found = 0;

	current =  &OpenFiles;
	current = current->next;
	while (current != NULL) {
		if(current->number == number) {
			found = 1;
			break;
		}
		current = current->next;
	}
	if (found == 0)
		current = NULL;

	return(current);
}

/* deleteNodeFromList
 *	Remove a node from the OpenFiles list given a unique identification number.
 *	
 *	Inputs: number - a unique number identifying the node.
 *  Outputs: none
 */

void deleteNodeFromList(int number) {
	NodeType *MatchNode;
	MatchNode = FindNodeInList(number);

	if(MatchNode != NULL) {
		MatchNode->previous->next = MatchNode->next;
		if (MatchNode->next != NULL)
			MatchNode->next->previous = MatchNode->previous;

        	DeleteCriticalSection(&MatchNode->section);

		free((void *)MatchNode);

	}else{
		throw("Unable to find node in list");
    }
}	

/* cleanList
	For each node in the OpenFiles list, close the stream and the file. Also close
 *  the AVIFILE interface. 
 */

void cleanList(void) {
	NodeType *ListHead;
	NodeType *current;

	ListHead = &OpenFiles;

	while (ListHead->next != NULL) {
		current = ListHead->next;

		if (current->psCompressed) {
			AVIStreamClose(current->psCompressed);
			current->psCompressed = NULL;
		}
		if (current->pfile) {
			AVIFileClose(current->pfile);
			current->pfile = NULL;
		}

		ListHead->next = current->next;
        	DeleteCriticalSection(&current->section);
		free((void *)current);
	}
	     
     AVIFileExit();
 }
		

void exitAVI(void) {
	NodeType *ListHead;
	ListHead = &OpenFiles;

	if (ListHead ->next != NULL) {
		cleanList();
		throw("Closing all open AVI files. It's no longer possible to write to any previously open AVI files.");
	}
}

/// This global section allows exported fcns to access nodes safely.
CRITICAL_SECTION AVI_NODE_ACCESS_SECTION;

/* openFile opens the AVI file and returns a unique file ID */

int openFile(const char *filename) {
  PAVIFILE pfile;
  HRESULT hr;
  int id;
  
  AVIFileInit();
  hr = AVIFileOpen(&pfile,		    /* returned file pointer */
                   filename,		            /* file name */
                   OF_WRITE | OF_CREATE,	    /* mode to open file with */
                   0L);

  if (hr != AVIERR_OK)
      throw("Failed to open file.");
  
  /* Keep track of stream handles */
  id = addNodetoList(pfile,NULL);

  /* Return a unique number */
  return id;
}

/*
 * writeFrame writes a frame to the AVI
 */

void writeFrame(int identifier, char *StreamName, int FrameNumber, unsigned char *frameData, char* compression, long width, long height, int bitcount, int ImageSize, int ClrUsed, unsigned char *colormap, int Quality, int FramesPerSecond, int KeyFrameRate)
{
    HRESULT hr;
    PBITMAPINFO bi;
    BITMAPINFOHEADER bih;
    AVICOMPRESSOPTIONS opts;
    AVISTREAMINFO strhdr;
    
    FOURCC CompCode;

    NodeType *handle;
    PAVIFILE pfile;
    PAVISTREAM ps, psCompressed;

    unsigned int i, j;
    
    // Enter global section.
    EnterCriticalSection(&AVI_NODE_ACCESS_SECTION);
    handle = FindNodeInList(identifier);
    if(handle == NULL){
        // Leave global section.
        LeaveCriticalSection(&AVI_NODE_ACCESS_SECTION);
        throw("Failed to find the open file handle.");
    }

    // Overlapping critical sections prevent 
    // someone from closing the handle on us 
    // right after we accessed the handle.
    EnterCriticalSection(&handle->section);
    LeaveCriticalSection(&AVI_NODE_ACCESS_SECTION);

    psCompressed = handle->psCompressed;
    pfile = handle->pfile;
    if (FrameNumber == 0){        
        /* Populate bih structure */
        bih.biSize = sizeof(BITMAPINFOHEADER);
        bih.biWidth = width;
        bih.biHeight = abs(height);
        bih.biPlanes = 1;
        bih.biBitCount = bitcount; 
        bih.biCompression = 0; /* 0 For uncompressed */
        bih.biSizeImage = ImageSize;
        
        bih.biXPelsPerMeter = 0;
        bih.biYPelsPerMeter = 0;
        bih.biClrUsed = ClrUsed; 
        bih.biClrImportant = 0;
        
        bi = (PBITMAPINFO) malloc(sizeof(BITMAPINFOHEADER) + 
            sizeof(RGBQUAD) * bih.biClrUsed);
        
        bi->bmiHeader = bih;
        
        if( colormap != NULL ) {
            j=0;

            for(i=0;i<bih.biClrUsed; i++) {
                bi->bmiColors[i].rgbBlue = colormap[j++];
                bi->bmiColors[i].rgbGreen = colormap[j++];
                bi->bmiColors[i].rgbRed = colormap[j++];
                bi->bmiColors[i].rgbReserved = colormap[j++];
            }
        }
        
        memset(&strhdr, 0, sizeof(strhdr));
        strhdr.fccType                = streamtypeVIDEO;/* stream type */
        strhdr.fccHandler             = 0;
        strhdr.dwScale                = 100;
        strhdr.dwRate                 = FramesPerSecond;		
        strhdr.dwSuggestedBufferSize  = ImageSize;
        
        strcpy(strhdr.szName, StreamName);
        
        SetRect(&strhdr.rcFrame, 0, 0,		    /* rectangle for stream */
            (int) width, abs((int) height));
        
        /* 	For first frame, create the stream.	*/
        hr = AVIFileCreateStream(pfile,		    /* file pointer */
            &ps,		    /* returned stream pointer */
            &strhdr);	    /* stream header */
        if (hr != AVIERR_OK){
            LeaveCriticalSection(&handle->section);
            throw("Failed to create video stream.");
        }
        
        memset(&opts, 0, sizeof(opts));
        
        CompCode = mmioFOURCC(compression[0],compression[1],compression[2],compression[3]);        
        if (CompCode ==0)				/* Uncompressed AVI */
        {
            opts.fccType = streamtypeVIDEO;
            opts.fccHandler = CompCode;
            opts.dwKeyFrameEvery = 0;
            opts.dwQuality = 0;
            opts.dwFlags = 0;
            opts.lpFormat = bi;
            opts.cbFormat = bi->bmiHeader.biSize + 
                bi->bmiHeader.biClrUsed * sizeof(RGBQUAD);
            
        }
        else							/* Compressed AVI */
        {
            opts.fccType = 0; /* streamtypeVIDEO doesn't seem to work for all compressors. */
            opts.fccHandler = CompCode;
            opts.dwKeyFrameEvery = KeyFrameRate;
            opts.dwQuality = Quality;
            opts.dwFlags = AVICOMPRESSF_KEYFRAMES;
            opts.lpFormat = bi;
            opts.cbFormat = bi->bmiHeader.biSize + bi->bmiHeader.biClrUsed * sizeof(RGBQUAD);
        }
        
        hr = AVIMakeCompressedStream(&psCompressed, ps, &opts, NULL);
        if (hr == AVIERR_NOCOMPRESSOR) {		
            if (ps) {
                AVIStreamClose(ps);
                ps = NULL;
            }
            LeaveCriticalSection(&handle->section);
            throw("Can not locate compressor.");
        }
        else if (hr == AVIERR_MEMORY) {	
            if (ps) {
                AVIStreamClose(ps);
                ps = NULL;
            }
            LeaveCriticalSection(&handle->section);
            throw("Not enough memory available for compressor.");
        }
        else if (hr == AVIERR_UNSUPPORTED) {	
            if (ps) {
                AVIStreamClose(ps);
                ps = NULL;
            }
            LeaveCriticalSection(&handle->section);
            throw("Compressor can not compress this data type.");
        }
        else if (hr != AVIERR_OK) {	
            if (ps) {
                AVIStreamClose(ps);
                ps = NULL;
            }
            LeaveCriticalSection(&handle->section);
            throw("Failed to make compressed stream. Unable to determine reason.");
        }
        
        /* Assign the stream to the list */
        handle->psCompressed = psCompressed;
        if (ps) {
            AVIStreamClose(ps);
            ps = NULL;
        }
        
        hr = AVIStreamSetFormat(psCompressed, 0,
            bi,	    /* stream format */
            bi->bmiHeader.biSize + bi->bmiHeader.biClrUsed * sizeof(RGBQUAD));
        if (hr != AVIERR_OK){
            LeaveCriticalSection(&handle->section);
            throw("Failed to set stream format.");
        }
        
        free(bi);

    } /* End of code required for first frame */
      
    // Write the frame out.
    hr = AVIStreamWrite(psCompressed,	/* stream pointer */
                        FrameNumber,				/* time of this frame */
                        1,				/* number to write */
                        frameData,/* pointer to data */
                        ImageSize,	/* size of this frame */
                        0,			 /* flags.... */
                        NULL, NULL);
    if (hr != AVIERR_OK){
        LeaveCriticalSection(&handle->section);
        throw("Failed to write stream data.");
    }

    // Done writing.
    LeaveCriticalSection(&handle->section);
}

/*
 * closeFile closes the AVI file
 */

void closeFile(int identifier) {
    PAVIFILE pfile;
    PAVISTREAM psCompressed;
    NodeType *handle;
    
    // Enter global section.
    EnterCriticalSection(&AVI_NODE_ACCESS_SECTION);

    /*Need to find the stream handle for this particular file. */
    handle = FindNodeInList(identifier);    
    if(handle == NULL){
        // Leave global section.
        LeaveCriticalSection(&AVI_NODE_ACCESS_SECTION);
        throw("The file is not open.");
    }
    
    // Overlapping critical sections prevent 
    // someone from closing the handle on us 
    // right after we accessed the handle.
    EnterCriticalSection(&handle->section);

    psCompressed = handle->psCompressed;
    pfile = handle->pfile;
    
    if (psCompressed) {
        AVIStreamClose(psCompressed);
        handle->psCompressed = NULL;
    }
    
    if (pfile) {
        AVIFileClose(pfile);
        handle->pfile = NULL;
    }
    
    LeaveCriticalSection(&handle->section);
    deleteNodeFromList(identifier);
    LeaveCriticalSection(&AVI_NODE_ACCESS_SECTION);
}

/*
 * Returns the number of video codecs on the user's machine
 */

int getNumVidCodecs() {
    HKEY hBase;
    HKEY hKey;
    LONG lret;
    DWORD numValues;
    const char * name = NULL;
    const char * path = NULL;
    int counter, inner_counter = 0;
    hBase = HKEY_LOCAL_MACHINE;

    /*
     * This doesn't work on Win98.  Not much we can do since they store
     * this information in the system.ini file rather than the registry.
     */
    path = "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion\\Drivers32";
    lret = RegOpenKeyEx(hBase, path, 0, KEY_QUERY_VALUE, &hKey);

    if (lret == ERROR_SUCCESS) {
        LONG lret = ERROR_SUCCESS;
        
        /*
         * Also use this to get the number of name and value of the keys...
         */
        RegQueryInfoKey(hKey, NULL, NULL, NULL, NULL, NULL, NULL, &numValues, NULL, NULL, NULL, NULL);

        for (counter = 0; counter < (int) numValues; counter++) {
            CHAR name[80];
            DWORD numChars = 80;

            lret = RegEnumValue(hKey, counter, name, (LPDWORD) &numChars, NULL, NULL, NULL, NULL);
            if (lret == ERROR_SUCCESS) {

                /*
                 * Get all the ones following "vidc." or "VIDC."
                 */
                if (name != NULL) {
                    if ((strncmp((LPTSTR) name, "vidc.", 5) == 0) || (strncmp((LPTSTR) name, "VIDC.", 5) == 0))
                        inner_counter++;

                }else
                {
                    /*
                    printf("Error with call to RegEnumValue\n");
                    */
                }
            }else {
                LPVOID lpMsgBuf;
                if (!FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                    FORMAT_MESSAGE_IGNORE_INSERTS, NULL, GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                    (LPTSTR) &lpMsgBuf, 0, NULL))
                {
                }
                /*
                    printf("RegEnumValue failed: %s\n", (LPCTSTR) lpMsgBuf);
                */
            }
        }
    } 
    RegCloseKey(hKey);
    return inner_counter;
}

/*
 * Gets a list of the FourCCs in the user's Windows registry
 * which begin with vidc. or VIDC.
 */
char ** getVidCodecs(void) {
    HKEY hBase;
    HKEY hKey;
    LONG lret;
    DWORD numValues;
    const char * name = NULL;
    const char * path = NULL;
    char **codecList;
    int counter, inner_counter = 0;
    hBase = HKEY_LOCAL_MACHINE;

    /*
     * This doesn't work on Win98.  Not much we can do since they store
     * this information in the system.ini file rather than the registry.
     */
    path = "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion\\Drivers32";
    lret = RegOpenKeyEx(hBase, path, 0, KEY_QUERY_VALUE, &hKey);

    codecList = (char **) malloc(getNumVidCodecs() * sizeof(char *));
    
    if (lret == ERROR_SUCCESS) {
        LONG lret = ERROR_SUCCESS;

        /*
         * Also use this to get the number of name and value of the keys...
         */
        RegQueryInfoKey(hKey, NULL, NULL, NULL, NULL, NULL, NULL, &numValues, NULL, NULL, NULL, NULL);

        for (counter = 0; counter < (int) numValues; counter++) {
            CHAR name[80];
            DWORD numChars = 80;

            lret = RegEnumValue(hKey, counter, name, (LPDWORD) &numChars, NULL, NULL, NULL, NULL);
            if (lret == ERROR_SUCCESS) {
                
                /*
                 * Get all the ones following "vidc." or "VIDC."
                 */
                if (name != NULL) {
                    if ((strncmp((LPTSTR) name, "vidc.", 5) == 0) || (strncmp((LPTSTR) name, "VIDC.", 5) == 0)) {
                        /* 
                         * Then we have a video codec, so store it on the list
                         */
                        codecList[inner_counter] = (char *) malloc(4 * sizeof(char) + 1);
                        strncpy(codecList[inner_counter], ((LPTSTR) name) + 5, 4);
                        codecList[inner_counter][4] = 0;
                        inner_counter++;
                    }
                    
                }else
                {
                    /*
                    printf("Error with call to RegEnumValue\n");
                    */
                }
            }else {
                LPVOID lpMsgBuf;
                if (!FormatMessage( FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                    FORMAT_MESSAGE_IGNORE_INSERTS, NULL, GetLastError(),
                    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPTSTR) &lpMsgBuf, 0, NULL))
                {
                }
                /*
                    printf("RegEnumValue failed: %s\n", (LPCTSTR) lpMsgBuf);
                */
            }
        }
    } 
    RegCloseKey(hKey);
    return codecList;
}

}
