// WindowAPI.c
// Set window state using Windows API
// WindowAPI(FigureHandle, Command, [Value])
// INPUT:
//   FigureHandle:  Matlab Handle of a figure with enabled 'Visible' property.
//                  For all commands except 'FullScreen' and 'Clip' the Windows
//                  handle (e.g. obtained by the 'GetHWnd' command) is accepted
//                  also.
//   Command:       String, not case sensitive, first 3 characters matter:
//     'TopMost':   Keep the window on top of all other windows. This does not
//                  enable the modal state of the window.
//     'NoTopMost': Disable the topmost state.
//     'Front':     Lift the window on top of all others, flash the icon in the
//                  taskbar. The window is not made active. This may interrupt
//                  the current work of the user.
//     'Minimize':  Minimize the window.
//     'Restore':   Restore the original size before minimization or
//                  maximization.
//     'Maximize':  Maximize: Full screen with visible Windows taskbar and
//                  menubar of the figure, but no border around the figure.
//     'XMax', 'YMax': Maximize the window horizontally or vertically.
//     'Position':  Set inner position of the figure realtive to the current
//                  monitor, or to the monitor with index defined as 4th input.
//                  The position value can be:
//                    DOUBLE vector [Left, Top, Right, Bottom]: pixel units
//                            relative to monitor.
//                    'work': Full monitor position without taskbar and sidebar.
//                    'full': Full monitor position. Using this you see only the
//                            figure's contents without title, border and
//                            menubar and without the taskbar.
//     'OuterPosition': As 'Position', but sets the outer position.
//     'ToScreen':  Move window completely to the nearest monitor. See MOVEGUI.
//     'Flash':     Short flashing the window border or the taskbar icon.
//     'Alpha':     WindowAPI(FigH, 'Alpha', A): Set the transparency level of
//                  the complete window from A=0.0 (invisible) to A=1.0 (opaque).
//                  WindowAPI(FigH, 'Alpha', A, [R,G,B]): The pixels with the
//                  color [R,G,B] are not drawn such that the pixels behind the
//                  window are visible. [R,G,B] must be integers in the range
//                  from 0 to 255. See NOTES.
//     'Opaque':    Release memory needed for Alpha blending to save resources.
//     'Clip':      Clip the window border or a specified rectangle. 3rd input:
//                    TRUE:  Clip the figure border, "splash screen".
//                    FALSE: Show the full window.
//                    [X, Y, Width, Height]: Rectangle in coordinates relative
//                           to figure as [X, Y, Width, Height] measured from
//                           bottom left in pixels.
//     'LockCursor': Keep cursor inside a rectangle. 3rd input:
//                    1, TRUE: Limit to figure,
//                    [X, Y, Width, Height]: Rectangle in pixel units relative
//                             to figure, DOUBLE vector.
//                    0, FALSE or omitted: free the cursor.
//                  NOTE: This does not lock the topmost window. You can
//                  activate the command window using Alt-Tab and free the
//                  cursor by:  WindowAPI('UnlockCursor').
//     'SetFocus':  Set the keyboard focus to the specified figure. Actually
//                  "figure(FigHandle)" should do this according to Matlab's
//                  documentation, but it doesn't from version 6.5 to 2009a or
//                  higher.
//
// GET INFORMATION:
// Reply = WindowAPI(FigureHandle, Command)
//     'GetStatus': Current window statue: 'maximized', 'minimized', 'restored'.
//     'GetHWnd':   Get the OS handle of the figure as UINT64 value. Most
//                  commands of WindowAPI are faster using this handle, but the
//                  MATLAB handle is needed if the inner figure position is used
//                  e.g. in 'Position', 'Clip' or 'LockCursor'.
//                  NOTE: The HWnd handle changes if the visibility of a figure
//                        is disabled!
//     'Monitor':   Information about monitor with the largest overlap to the
//                  figure. Struct with fields:
//                    FullPosition: [X, Y, W, H] monitor size.
//                    WorkPosition: [X, Y, W, H] size without taskbar / sidebar.
//                    FigureOnScreen: LOGICAL flag, TRUE if any part of the
//                                  figure overlaps with this monitor. Without
//                                  overlap the nearest monitor is replied.
//                    isPrimaryMonitor: This monitor is the primary monitor.
//     'Position', 'OuterPosition': Using these commands with 2 inputs only
//                  reply the size relative to the current monitor as struct
//                  with the fields 'Position' and 'MonitorIndex'.
//
// CONTROL FEATURES:
//   FormerStatus = WindowAPI('topmost', Status): If Status is 'off', the
//      TopMost property is not set in future calls. This keeps Matlab in the
//      background e.g. during a long test, even if a function asks WindowAPI
//      to set a figure as top-most window.
//   WindowAPI('UnlockCursor'): Frees a cursor locked to a rectangle.
//
// NOTES:
//   This function calls Windows-API functions => No Linux, no MacOS - sorry!
//   Suggestions for Unix implementations are very appreciated!
//
//   Enabling the Alpha blending can cause a short black flashing of the
//   figure. It might  be nicer to create the figure outside the visible screen
//   area, enable the Alpha blending and move the window to the desired position
//   afterwards.
//
//   To find the window handle of the OS, the figure title is modified for some
//   milliseconds. This works without Java and for all known Matlab versions.
//   It fails, if two Matlab sessions call this function at the same time for
//   figures with the same title - a very unlikely constellation.
//
//   Alpha blending does not work reliably with the OpenGL renderer on my
//   computer. Setting the FIGURE's WVisual property to '07' ("RGB 16
//   bits(05 06 05 00) zdepth 16, Hardware Accelerated, Opengl, Double Buffered,
//   Window") helps most of the times. I suggest using the Painters or ZBuffer
//   renderer for Alpha blending.
//
//   The StencilRGB value of the command Alpha does not exactly equal the RGB
//   value set by Matlab: The Matlab colors [88/255,0,0] to [95/255,0,0] are
//   matching the StencilRGB=[90,0,0] on my laptop. I suggest to use the RGB
//   colors with either 0 or 255 as components.
//
// EXAMPLES: See WindowAPI.m and demo_WindowAPI.m.
//
// COMPILE:
// Sorry: Windows only.
// Not compatible with LCC compiler shipped with Matlab. MSVC compiler works.
// Activate the compiler on demand: mex -setup
//   mex -O WindowAPI.c
// You can download the precompiled MEX file: http://www.n-simon.de/mex
// The case-insensitive string comparison depends on the compiler:
//   strnicmp, strncmpi, _strnicmp, ...
// Define the macro STRINGCMPI in the source code on demand.
// Run the unit-test uTest_WindowAPI to check validity.
//
// Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
//         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
// Assumed Compatibility: higher Matlab versions
// Author: Jan Simon, Heidelberg, (C) 2008-2011 matlab.THISYEAR(a)nMINUSsimon.de
//
// See in the FEX:
// ShowWindow, Matthew Simoneau:
//   http://www.mathworks.com/matlabcentral/fileexchange/3407
// Window Manipulation, Phil Goddard:
//   http://www.mathworks.com/matlabcentral/fileexchange/3434
// api_showwindow, Mihai Moldovan:
//   http://www.mathworks.com/matlabcentral/fileexchange/2041
// maxfig, Mihai Moldovan:
//   http://www.mathworks.com/matlabcentral/fileexchange/6913
// setFigTransparency, Yair Altman:
//   http://www.mathworks.com/matlabcentral/fileexchange/30583
// FigureManagement (multi-monitor setup), Mirko Hrovat:
//   http://www.mathworks.com/matlabcentral/fileexchange/12607

/*
% $JRev: R-w V:067 Sum:TwYOiG8iAZqi Date:04-Nov-2011 09:35:52 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_WindowAPI $
% $File: Tools\Mex\Source\WindowAPI.c $
% History:
% 023: 10-Sep-2008 11:13, Only 3 characters of [Command] matter.
% 034: 13-Sep-2010 23:24, [GetStatus] command, nicer temporary name.
% 041: 11-May-2011 22:15, Alpha, GetHWnd. SetWindowLong->SetWindowLongPtr.
%      SetWindowLongPtr supports Windows 64bit.
%      Renamed: OSWindowFcn -> WindowAPI.
% 051: 24-May-2011 20:34, FullScreen: improved placement, not topmost.
%      Marc Lalancette found a more consistent method to consider the window
%      border for the fullscreen mode. Thanks!
% 057: 22-Jun-2011 10:59, New commands: Clip, Monitor, Position, OuterPosition.
%      'Position' and 'OuterPosition' replace the former 'Screen' and
%      'FullScreen' commands. 'ToScreen', 'SetFocus'.
% 067: 04-Nov-2011 09:25, LockCursor, tested under WinXP/32, Win7/64.
*/

// Headers: --------------------------------------------------------------------
#if !defined(__WINDOWS__) && !defined(_WIN32) && !defined(_WIN64)
#error Implemented for Windows only!
#endif

#include <windows.h>

// MSVC 2008 needs a defined COMPILE_MULTIMON_STUBS to include "multimon.h",
// but including "WinUser.h" works without also.
// #define COMPILE_MULTIMON_STUBS
//#include "multimon.h"
#include "WinUser.h"

#include <string.h>
#include <shellapi.h>
#include <math.h>
#include "mex.h"

// Compiler dependent settings: ------------------------------------------------
#if defined(__LCC__)
#  ifndef WS_EX_LAYERED
// We could define the missing constants for LCC2.4 (shipped with Matlab), but
// the functions SetLayeredWindowAttributes, Set/GetWindowLongPtr are included
// in the libs. But LCC 3.8 downloaded from the net works!
#    error  Not compatible with LCC 2.4 - sorry
#    define WS_EX_LAYERED  0x00080000  // unused dummy!
#    define LWA_ALPHA      0x00000002  // unused dummy!
#  endif
#  define STRINGCMPI strnicmp
#  define uint64_T unsigned long long
#elif defined(__BORLANDC__)
#  define STRINGCMPI strnicmp
#else
// Different names for case-insensitive string comparison:
#  define STRINGCMPI _strnicmp
#endif

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_ID    "JSimon:WindowAPI:"
#define ERR_HEAD  "*** WindowAPI[mex]: "
#define WARN_HEAD "### WindowAPI[mex]: "

// Maximal string length of input command:
#define Command_LEN 16   // Longest command: "OuterPosition"
#define Param_LEN 4      // "on" or "off"
#define Area_LEN 5       // "work" or "full"

// Enabled properties:
static int allow_TopMost = 1;  // Allow setting window to topmost

// Prototypes: -----------------------------------------------------------------
void GetPosCorrect(const double Fig_Handle, int PosCorrect[4]);
HWND myGetHWnd(const double Fig_Handle);
mxArray *myMonitorInfo(HWND hWnd);
void myWindowRgn(HWND hWnd, const double Fig_Handle, const mxArray *Value);
void myAlpha(HWND hWnd, const mxArray *InAlpha, const mxArray *InStencilRGB);
void mySetPosition(HWND hWnd, const double Fig_Handle, const BOOL inner,
                   const mxArray *Value, const mxArray *MonitorIndex);
void myToScreen(HWND hWnd);
mxArray *myGetPosition(HWND hWnd, const double Fig_Handle, const BOOL inner);
void myMaximizeXY(HWND hWnd, BOOL horz);
void myClipCursor(HWND hWnd, const double Fig_Handle, const mxArray *Param);

// Multi-monitor management: ---------------------------------------------------
typedef struct tagGETMONITOR {
   HMONITOR target;
   int      index;
   BOOL     found;
} GETMONITOR, *LPGETMONITOR;

int  myGetIndexFromMonitor(HMONITOR hMonitor);
BOOL CALLBACK myMonitorMatchProc(HMONITOR hMonitor, HDC hdc,
                                LPRECT lprc, LPARAM dwData);

HMONITOR myGetMonitorFromIndex(HWND hWnd, int MonitorIndex);
BOOL CALLBACK myIndexMatchProc(HMONITOR hMonitor, HDC hdcMonitor,
                               LPRECT lprcMonitor, LPARAM dwData);


// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Inputs:
  double Fig_Handle;
  char   Command[Command_LEN], Param[Param_LEN];
  BOOL   orphanedInput;
  
  // Output:
  mxArray *Reply = NULL;
  
  // Computations:
  HWND hWnd;
  
  // Check inputs:
  if (nrhs < 1 || nrhs > 4) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                       ERR_HEAD "1 to 4 inputs accepted.");
     
  } else if (mxIsEmpty(prhs[0])) {
     if (nrhs != 0) {
        mexErrMsgIdAndTxt(ERR_ID   "EmptyInput",
                          ERR_HEAD "No output for empty handle.");
     }
     return;  // Nothing to do for empty input
     
  } else if (mxIsChar(prhs[0])) {  // Meta commands: ---------------------------
     // Set control flags to disable setting a window to topmost or in front of
     // other programs:
     mxGetString(prhs[0], Command, Command_LEN);
     if (STRINGCMPI(Command, "topmost", 7) == 0) {
        // Reply former status:
        if (allow_TopMost) {
           plhs[0] = mxCreateString("on");
        } else {
           plhs[0] = mxCreateString("off");
        }
        
        // Set the new value:
        if (nrhs >= 2) {
           mxGetString(prhs[1], Param, Param_LEN);
           if (STRINGCMPI(Param, "on", 2) == 0) {
              allow_TopMost = 1;
           } else {
              allow_TopMost = 0;
           }
        }
        
     } else if (STRINGCMPI(Command, "unlockcursor", 6) == 0) {
        ClipCursor(NULL);
        
     } else {
        mexErrMsgIdAndTxt(ERR_ID   "BadCommand",
                          ERR_HEAD "Unknown command: [%s].", Command);
     }
     return;
          
  } else if (nrhs < 2) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                       ERR_HEAD "Setting properties needs 2 to 4 inputs.");
     
  } else if (!mxIsChar(prhs[1])) {  // Standard: WindowAPI(H, Command, ...)
     mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput2",
                       ERR_HEAD "2nd input must be a string.");
     
  } else if (mxGetNumberOfElements(prhs[0]) != 1) {
     mexErrMsgIdAndTxt(ERR_ID   "BadInput1",
                       ERR_HEAD "1st input must be a scalar figure handle.");
  }
  
  // Need at least 3 characters to identify the command:
  if (mxGetNumberOfElements(prhs[1]) < 3) {
     mexErrMsgIdAndTxt(ERR_ID   "BadInput2Size",
                       ERR_HEAD "[Command] needs at least 3 characters.");
  }
  
  // Get the figure handle from the 1st input: ---------------------------------
  if (mxIsDouble(prhs[0])) {  // Empty matrix excluded before
     Fig_Handle = mxGetScalar(prhs[0]);
     hWnd       = myGetHWnd(Fig_Handle);
     
  } else if (mxIsUint64(prhs[0])) {
     Fig_Handle = -1.0;  // Invalid Matlab handle
     hWnd       = *(HWND *) mxGetData(prhs[0]);
     if (!IsWindow(hWnd)) {
        mexErrMsgIdAndTxt(ERR_ID "BadHWND",
                         ERR_HEAD "HWND does not point to an existing window.");
     }
     
  } else {
     mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput1",
                       ERR_HEAD "1st input must be a figure handle.");
  }
  
  // Get command (longest valid command string: "OuterPosition": 14 chars):
  mxGetString(prhs[1], Command, Command_LEN);
  orphanedInput = (BOOL) (nrhs >= 3);
  
  // The actions: ==============================================================
  if (STRINGCMPI(Command, "TopMost", 3) == 0) {
     if (allow_TopMost) {
       SetWindowPos(hWnd, HWND_TOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
     } else {
       SetWindowPos(hWnd, HWND_TOP,     0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
     }

  } else if (STRINGCMPI(Command, "NoTopMost", 5) == 0) {
     SetWindowPos(hWnd, HWND_NOTOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);

  } else if (STRINGCMPI(Command, "Maximize", 3) == 0) {
     // Full screen with visible Windows taskbar and menubar of the figure,
     // but no border around the figure, figure gets on top:
     ShowWindow(hWnd, SW_MAXIMIZE);
     
  } else if (STRINGCMPI(Command, "Minimize", 3) == 0) {
     ShowWindow(hWnd, SW_MINIMIZE);

  } else if (STRINGCMPI(Command, "Restore", 7) == 0) {
     ShowWindow(hWnd, SW_RESTORE);

  } else if (STRINGCMPI(Command, "XMax", 4) == 0) {
     myMaximizeXY(hWnd, TRUE);
     
  } else if (STRINGCMPI(Command, "YMax", 4) == 0) {
     myMaximizeXY(hWnd, FALSE);

  } else if (STRINGCMPI(Command, "Flash", 5) == 0) {
     FlashWindow(hWnd, TRUE);
     
  } else if (STRINGCMPI(Command, "SetFocus", 8) == 0) {
     SetFocus(hWnd);
     
  } else if (STRINGCMPI(Command, "Position", 3) == 0) {
     // Set inner figure position without Matlab's auto-resize.
     orphanedInput = FALSE;
     switch (nrhs) {
        case 4:  mySetPosition(hWnd, Fig_Handle, TRUE, prhs[2], prhs[3]);
                 break;
        case 3:  mySetPosition(hWnd, Fig_Handle, TRUE, prhs[2], NULL);
                 break;
        case 2:  Reply = myGetPosition(hWnd, Fig_Handle, TRUE);
                 break;
        default: mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                                   ERR_HEAD "[Position] needs 2 to 4 inputs.");
     }
     
  } else if (STRINGCMPI(Command, "OuterPosition", 5) == 0) {    // -------------
     // Set outer figure position without Matlab's auto-resize:
     orphanedInput = FALSE;
     switch (nrhs) {
        case 4:  mySetPosition(hWnd, Fig_Handle, FALSE, prhs[2], prhs[3]);
                 break;
        case 3:  mySetPosition(hWnd, Fig_Handle, FALSE, prhs[2], NULL);
                 break;
        case 2:  Reply = myGetPosition(hWnd, Fig_Handle, FALSE);
                 break;
        default: mexErrMsgIdAndTxt(ERR_ID "BadNInput",
                               ERR_HEAD "[OuterPosition] needs 2 to 4 inputs.");
     }
     
  } else if (STRINGCMPI(Command, "ToScreen", 8) == 0) {    // ------------------
     myToScreen(hWnd);
     
  } else if (STRINGCMPI(Command, "Front", 5) == 0) {  // -----------------------
     // Raise window even in front of other programs, but do not activate it if
     // another program has the focus - then only taskbar icon flashes.
     if (allow_TopMost) {
        // Set window on top of all others, but disable the TOPMOST flag again:
        SetWindowPos(hWnd, HWND_TOPMOST,   0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
        SetForegroundWindow(hWnd);
        SetWindowPos(hWnd, HWND_NOTOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
     } else {  // TOPMOST disabled: At least catch attention by flashing:
        SetForegroundWindow(hWnd);
     }
     
     // BringWindowToTop(hWnd) triggers a resize from Matlab
     
  } else if (STRINGCMPI(Command, "GetStatus", 9) == 0) {  // -------------------
     // Reply current status as string:
     if (IsZoomed(hWnd)) {
        Reply = mxCreateString("maximized");
     } else if (IsIconic(hWnd)) {
        Reply = mxCreateString("minimized");
     } else {
        Reply = mxCreateString("restored");
     }

  } else if (STRINGCMPI(Command, "GetHWnd", 7) == 0) {  // ---------------------
     // Reply the handle as UINT64:
     Reply = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
     *(uint64_T *) mxGetData(Reply) = (uint64_T) hWnd;

  } else if (STRINGCMPI(Command, "Alpha", 5) == 0) {    // ---------------------
     // Set aplha blending and stencil color:
     orphanedInput = (nrhs > 4);
     if (nrhs == 3) {
        myAlpha(hWnd, prhs[2], NULL);
     } else {
        myAlpha(hWnd, prhs[2], prhs[3]);
     }
     
  } else if (STRINGCMPI(Command, "Opaque", 6) == 0) {  // ----------------------
     // Set opaque again: Remove WS_EX_LAYERED from this window styles to save
     // memory:
     SetWindowLongPtr(hWnd, GWL_EXSTYLE,
                      GetWindowLongPtr(hWnd, GWL_EXSTYLE) & ~WS_EX_LAYERED);
     
     // Update components:
     RedrawWindow(hWnd, NULL, NULL,
                  RDW_ERASE | RDW_INVALIDATE | RDW_FRAME | RDW_ALLCHILDREN);
     mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
  
  } else if (STRINGCMPI(Command, "Monitor", 7) == 0)  {  // --------------------
     // Get monitor info for the monitor of this window:
     Reply = myMonitorInfo(hWnd);

  } else if (STRINGCMPI(Command, "Clip", 4) == 0)  {  // -----------------------
     // Clip rectangular subregion of this window:
     if (nrhs == 2) {
        myWindowRgn(hWnd, Fig_Handle, NULL);
     } else if (nrhs >= 3) {
        orphanedInput = (BOOL) (nrhs > 3);
        myWindowRgn(hWnd, Fig_Handle, prhs[2]);
     }
     
  } else if (STRINGCMPI(Command, "LockCursor", 4) == 0) {  // ------------------
     if (nrhs == 2) {
        myClipCursor(hWnd, Fig_Handle, NULL);
     } else {
        myClipCursor(hWnd, Fig_Handle, prhs[2]);
        orphanedInput = (BOOL) (nrhs > 3);
     }
     
  } else {  // -----------------------------------------------------------------
     mexErrMsgIdAndTxt(ERR_ID   "BadCommand",
                       ERR_HEAD "Unknown command; [%s]", Command);
  }
    
  if (orphanedInput) {  // Mention unused inputs:
     mexWarnMsgIdAndTxt(ERR_ID    "OrphanedInput",
                        WARN_HEAD "Some inputs are not used.");
  }
  
  // Create output: ------------------------------------------------------------
  if (Reply != NULL) {
     if (nlhs > 1) {
       mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                         ERR_HEAD "Only 1 output allowed.");
     }
     plhs[0] = Reply;
     
  } else if (nlhs > 0) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                       ERR_HEAD "No output for this command.");
  }
  
  return;
  
  // TODO:
  // Not exactly matching the context, because it does not belong to a figure:
  // } else if (STRINGCMPI(Command, "ScrollBarExtent", 6) == 0)  {
  //   // Get size of horizontal and vertical scrollbars - I've never seen that
  //   // these size differ.
  //   Reply = mxCreateDoubleMatrix(1, 2, mxREAL);
  //   Dp    = mxGetPr(Reply);
  //   Dp[0] = (double) GetSystemMetrics(SM_CYHSCROLL);
  //   Dp[1] = (double) GetSystemMetrics(SM_CXVSCROLL);
}

// *****************************************************************************
HWND myGetHWnd(const double Fig_Handle) {
  // Get Windows window handle from Matlab's figure handle.
  // Method: A magic string is appended temporaily to the name of the figure.
  // Then the WindowAPI function FindWindow gets the HWnd handle and the
  // original name is restored.
   
  mxArray       *Fig_Name, *Fig_NumberTitle, *Tmp_Name;
  mwSize        Name_Len;
  const mxArray *Fig_NameP, *Fig_NumberTitleP;
  const char    *MagicKey = "       .";   // Arbitrary magic key, 8 chars + /0
  char          *Tmp_NameC;
  
  // Properties of the figure:
  const char *Name_Str        = "Name",
             *NumberTitle_Str = "NumberTitle";
  
  // Concerning Windows functions:
  HWND hWnd;
  
  // Store the original figure Name and NumberTitle:
  Fig_NameP        = mexGet(Fig_Handle, Name_Str);
  Fig_NumberTitleP = mexGet(Fig_Handle, NumberTitle_Str);
  if (Fig_NameP == NULL || Fig_NumberTitleP == NULL) {
     // Most likely the number is not a (figure) handle!
     mexErrMsgIdAndTxt(ERR_ID   "BadInput1",
                       ERR_HEAD "Invalid figure handle.");
  }
  
  // Get the current figure name and status of NumberTitle:
  Fig_Name        = mxDuplicateArray(Fig_NameP);
  Fig_NumberTitle = mxDuplicateArray(Fig_NumberTitleP);
  
  // Create temporary window title:
  Name_Len  = mxGetNumberOfElements(Fig_Name);
  Tmp_NameC = (char *) mxMalloc(Name_Len + 9 * sizeof(mxChar));
  if (Tmp_NameC == NULL) {
     mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                       ERR_HEAD "No memory for C-string.");
  }
  mxGetString(Fig_Name, Tmp_NameC, Name_Len + 9);
  memcpy(Tmp_NameC + Name_Len, MagicKey, 9);
  Tmp_Name = mxCreateString(Tmp_NameC);
  
  // Set name to temporary string, because searching a window would fail if
  // figure names are not unique:
  if (mexSet(Fig_Handle, Name_Str,        Tmp_Name) ||
      mexSet(Fig_Handle, NumberTitle_Str, mxCreateString("off"))) {
     // NOTE: No idea why mexSet could fail here - but a check is cheap.
     mexErrMsgIdAndTxt(ERR_ID   "BadInput2Size",
                       ERR_HEAD "Cannot set figure properties.");
  }
  
  // Matlab 7 needs a DRAWNOW to activate changes:
  mexCallMATLAB(0, NULL, 0, NULL, "drawnow");

  // ==========================================
  // Until here the code is Linux compatible...
  // ==========================================
  
  // Get Windows-handle of the figure:
  hWnd = FindWindow(NULL, Tmp_NameC);
  
  // Restore original figure name:
  mexSet(Fig_Handle, Name_Str,        Fig_Name);
  mexSet(Fig_Handle, NumberTitle_Str, Fig_NumberTitle);
  mexCallMATLAB(0, NULL, 0, NULL, "drawnow");

  // Free memory as soon as possible:
  mxDestroyArray(Tmp_Name);
  mxDestroyArray(Fig_Name);
  mxDestroyArray(Fig_NumberTitle);
  mxFree(Tmp_NameC);
  
  if (hWnd == 0) {
     mexErrMsgIdAndTxt(ERR_ID   "CannotGetHWND",
                       ERR_HEAD "Cannot get HWND for figure handle.");
  }
     
  return hWnd;
}

// *****************************************************************************
void GetPosCorrect(const double Fig_Handle, int PosCorrect[4]) {
  // Get difference between inner and outer position. This depends on the
  // visibility of menus and toolbars and the dimensions of the window border
  // and the title bar. Because these values can be changed by the user, the
  // difference must be measured dynamically.
  // 'OuterPosition' is not documented but works since Matlab 6.5.
  // Multi-monitor proof: Difference do not depend on monitor.

  double *InnerPos, *OuterPos;
  
  if (Fig_Handle == -1.0) {
     mexErrMsgIdAndTxt(ERR_ID   "NoMatlabHandle",
                       ERR_HEAD "Matlab figure handle needed as input.");
  }
  
  // mexGet replies just a handle to the mxArrays, so no mxFree is needed:
  InnerPos = mxGetPr(mexGet(Fig_Handle, "Position"));
  OuterPos = mxGetPr(mexGet(Fig_Handle, "OuterPosition"));  // % $Undocumented$
  
  PosCorrect[0] = (int) (OuterPos[0] - InnerPos[0]);
  PosCorrect[1] = (int) (OuterPos[1] - InnerPos[1]);
  PosCorrect[2] = (int) (OuterPos[2] - InnerPos[2]);
  PosCorrect[3] = (int) (OuterPos[3] - InnerPos[3]);
  
  return;
}

// *****************************************************************************
void myAlpha(HWND hWnd, const mxArray *Alpha_In, const mxArray *StencilRGB_In) {
  // Set Alpha blending of the complete window and StencilRGB value of not drawn
  // pixels.
   
  double   Alpha    = 0.8;        // default
  DWORD    Flags    = LWA_ALPHA;  // default: no stencil color
  COLORREF RGBValue = 0;
  double   *dColor;
  uint8_T  *u8Color;
  mxArray  *RGBArray[1];
  
  // Get 3rd input [Alpha]:
  if (mxIsNumeric(Alpha_In) && mxGetNumberOfElements(Alpha_In) == 1) {
     Alpha = mxGetScalar(Alpha_In);
     if (Alpha < 0.0 || Alpha > 1.0) {
        mexErrMsgIdAndTxt(ERR_ID   "BadAlpha",
                          ERR_HEAD "[Alpha] out or range: 0.0 <= Alpha <= 1.0");
     }
  } else {
     mexErrMsgIdAndTxt(ERR_ID   "BadAlhpa",
                       ERR_HEAD "[Alpha] must be a numeric scalar.");
  }
     
  // Get 4th input [StencilRGB]:
  if (StencilRGB_In != NULL) {
     Flags |= LWA_COLORKEY;
     if (!mxIsNumeric(StencilRGB_In) ||
         mxGetNumberOfElements(StencilRGB_In) != 3) {
        mexErrMsgIdAndTxt(ERR_ID "BadAlhpa",
                     ERR_HEAD "[StencilRGB] must be a numeric [1 x 3] vector.");
     }
     
     // Get values in the UINT8 range from all numeric types:
     if (mxIsUint8(StencilRGB_In)) {
        u8Color  = (uint8_T *) mxGetData(StencilRGB_In);
        RGBValue = RGB(u8Color[0], u8Color[1], u8Color[2]);
     } else if (mxIsDouble(StencilRGB_In)) {
        dColor   = mxGetPr(StencilRGB_In);
        RGBValue = RGB((BYTE) dColor[0], (BYTE) dColor[1], (BYTE) dColor[2]);
     } else {
        mexCallMATLAB(1, RGBArray, 1, (mxArray **) &StencilRGB_In, "double");
        dColor   = mxGetPr(RGBArray[0]);
        RGBValue = RGB((BYTE) dColor[0], (BYTE) dColor[1], (BYTE) dColor[2]);
     }
  }
  
  // Set WS_EX_LAYERED on this window:
  SetWindowLongPtr(hWnd, GWL_EXSTYLE,
                   GetWindowLongPtr(hWnd, GWL_EXSTYLE) | WS_EX_LAYERED);
  
  // Ask the window and its children to repaint:
  SetLayeredWindowAttributes(hWnd, RGBValue, (BYTE) (255 * Alpha), Flags);
  
  // Update components:
  mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
  RedrawWindow(hWnd, NULL, NULL,
               RDW_ERASE | RDW_INVALIDATE | RDW_FRAME | RDW_ALLCHILDREN);
  
  return;
}

// *****************************************************************************
void mySetPosition(HWND hWnd, const double Fig_Handle, const BOOL inner,
                   const mxArray *Pos_In, const mxArray *Monitor_In) {
  // Set figure position relative to current or specified monitor.
  // Equivalent to "set(Fig_Handle, 'Position', Value)", but does not trigger
  // Matlab's auto resizing if the area exceeds the screen area. Two strings
  // are accepted as Value also: 'full' is the full monitor size, 'work'
  // excludes the taskbar.
  // Multi-monitor proof: Window position relative to current monitor.
  
  RECT        Rect;
  int         PosCorrect[4], Left, Top, Width, Height;
  double      *p, Index_d;
  char        Area[Area_LEN];
  HMONITOR    hMonitor;
  MONITORINFO mInfo;
  
  // Get monitor of figure or specified one:
  if (Monitor_In == NULL) {
     // Use the nearest if the figure does not overlap with any monitor:
     hMonitor = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONEAREST);
  } else {
     if (!mxIsNumeric(Monitor_In) || mxGetNumberOfElements(Monitor_In) != 1) {
        mexErrMsgIdAndTxt(ERR_ID   "BadMonitorIndex",
                          ERR_HEAD "[MonitorIndex] must be a scalar integer.");
     }
     Index_d = mxGetScalar(Monitor_In);
     if (Index_d < 1.0 || Index_d != floor(Index_d)) {
        mexErrMsgIdAndTxt(ERR_ID   "BadMonitorIndex",
                          ERR_HEAD "[MonitorIndex] must be a scalar integer.");
     }
     
     hMonitor = myGetMonitorFromIndex(hWnd, (int) Index_d);
  }
  
  // Get monitor info:
  mInfo.cbSize = sizeof(mInfo);                  // This is a MONITORINFO struct
  if (GetMonitorInfo(hMonitor, &mInfo) == 0)  {  // Failing is unlikely!
     mexErrMsgIdAndTxt(ERR_ID   "GetMonitorInfo",
                       ERR_HEAD "GetMonitorInfo failed [%d].", GetLastError());
  }
  
  if (mxIsChar(Pos_In)) {
     // Multi-monitor proof: Window position relative to current monitor:
     mxGetString(Pos_In, Area, Area_LEN);
     if (STRINGCMPI(Area, "work", 4) == 0)
        Rect = mInfo.rcWork;
     else if (STRINGCMPI(Area, "full", 4) == 0)
        Rect = mInfo.rcMonitor;
     else
        mexErrMsgIdAndTxt(ERR_ID   "BadNInput", ERR_HEAD "Unknown [Position] string.");
     
     // RECT to [X, Y, Width, Height] relative to virtual screen:
     Left   = Rect.left;
     Top    = Rect.top;
     Width  = Rect.right  - Left;
     Height = Rect.bottom - Top;
     
  } else if (mxIsDouble(Pos_In) && mxGetNumberOfElements(Pos_In) == 4) {
     // Convert Matlab position [X, Y, Width, Height] measured from bottom
     // left to monitor position measured from top left:
     Rect   = mInfo.rcMonitor;   // Relative to full monitor rect
     p      = mxGetPr(Pos_In);
     Left   = Rect.left   + (int) p[0] - 1;
     Top    = Rect.bottom - ((int) p[1] + (int) p[3] - 1);
     Width  = (int) p[2];
     Height = (int) p[3];
     
  } else {
     mexErrMsgIdAndTxt(ERR_ID   "BadInput3",
                       ERR_HEAD "Bad value for [Position] command.");
  }
    
  // Offset between inner and outer position:
  if (inner) {
     GetPosCorrect(Fig_Handle, PosCorrect);
     Left   += PosCorrect[0];
     Top    -= PosCorrect[3] + PosCorrect[1];
     Width  += PosCorrect[2];       // Value is negative
     Height += PosCorrect[3];       // Value is negative
  }
  
  // Set the adjusted window:
  SetWindowPos(hWnd, HWND_TOP, Left, Top, Width, Height,
               SWP_SHOWWINDOW | SWP_NOSENDCHANGING |
               SWP_NOACTIVATE | SWP_NOZORDER);
  
  // Let Matlab update the display:
  mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
  
  return;
}

// *****************************************************************************
mxArray *myGetPosition(HWND hWnd, const double Fig_Handle, const BOOL inner) {
  // Get figure position relative to full current monitor.
  // Multi-monitor proof: Window position relative to current monitor.
  
  RECT        Mon_Rect, Win_Rect;
  int         PosCorrect[4];
  double      *p;
  HMONITOR    hMonitor;
  MONITORINFO mInfo;
  mxArray     *Pos, *MonitorIndex, *Reply;
  const char  *Field[2] = {"Position", "MonitorIndex"};
  
  // Use the nearest if the figure does not overlap with any monitor:
  hMonitor = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONEAREST);
  
  // Get monitor info:
  mInfo.cbSize = sizeof(mInfo);                  // This is a MONITORINFO struct
  if (GetMonitorInfo(hMonitor, &mInfo) == 0)  {  // Failing is unlikely!
     mexErrMsgIdAndTxt(ERR_ID   "GetMonitorInfo",
                       ERR_HEAD "GetMonitorInfo failed [%d].", GetLastError());
  }
  
  // Create output:
  Pos = mxCreateDoubleMatrix(1, 4, mxREAL);
  p   = mxGetPr(Pos);
  
  // RECT to [X, Y, Width, Height] relative to current window:
  Mon_Rect = mInfo.rcMonitor;
  GetWindowRect(hWnd, &Win_Rect);
  
  p[0] = Win_Rect.left   - Mon_Rect.left + 1;    // X from the left
  p[1] = Mon_Rect.bottom - Win_Rect.bottom + 1;  // Y measured from bottom
  p[2] = Win_Rect.right  - Win_Rect.left;        // Width
  p[3] = Win_Rect.bottom - Win_Rect.top;         // Height
  
  // Offset between inner and outer position:
  if (inner) {
     GetPosCorrect(Fig_Handle, PosCorrect);
     p[0] -= PosCorrect[0];
     p[1] -= PosCorrect[1];
     p[2] -= PosCorrect[2];       // Value is negative
     p[3] -= PosCorrect[3];       // Value is negative
  }
  
  // Get monitor index - shortcut for primary monitor:
  if (mInfo.dwFlags & MONITORINFOF_PRIMARY)
     MonitorIndex = mxCreateDoubleScalar(1);
  else
     MonitorIndex = mxCreateDoubleScalar(myGetIndexFromMonitor(hMonitor));

  
  // Create output struct:
  Reply = mxCreateStructMatrix(1, 1, 2, Field);
  mxSetFieldByNumber(Reply, 0, 0, Pos);
  mxSetFieldByNumber(Reply, 0, 1, MonitorIndex);
    
  return Reply;
}

// *****************************************************************************
void myMaximizeXY(HWND hWnd, BOOL horz) {
  // Maximize a window horizontally or vertically only.
  // Multi-monitor proof: Changes are made relative to current monitor.
  
  HMONITOR    hMonitor;
  MONITORINFO mInfo;
  RECT        monitorRect, windowRect;
  int         Left, Top, Width, Height;
  
  // Use the nearest if the figure does not overlap with any monitor:
  hMonitor = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONEAREST);
  
  mInfo.cbSize = sizeof(mInfo);
  if (GetMonitorInfo(hMonitor, &mInfo) == 0)  {  // Failing is unlikely...
     mexErrMsgIdAndTxt(ERR_ID   "GetMonitorInfo",
                       ERR_HEAD "GetMonitorInfo failed [%d].", GetLastError());
  }
  
  // Get work rect of monitor and rect of window:
  monitorRect = mInfo.rcWork;
  GetWindowRect(hWnd, &windowRect);

  if (horz) {  // XMax:
     Left   = monitorRect.left;
     Top    = windowRect.top;
     Width  = monitorRect.right - monitorRect.left;
     Height = windowRect.bottom - Top;
  } else {     // YMax:
     Left   = windowRect.left;
     Top    = monitorRect.top;
     Width  = windowRect.right - Left;
     Height = monitorRect.bottom - monitorRect.top;
  }
  
  // Set the adjusted window:
  SetWindowPos(hWnd, HWND_TOP,
               Left, Top, Width, Height,
               SWP_SHOWWINDOW | SWP_NOSENDCHANGING |
               SWP_NOACTIVATE | SWP_NOZORDER);
  
  return;
}

// *****************************************************************************
void myToScreen(HWND hWnd) {
  // Move a window to the nearest monitor.
  // Multi-monitor proof: Changes are made relative to current monitor.
  
  HMONITOR    hMonitor;
  MONITORINFO mInfo;
  RECT        Mon_Rect, Win_Rect;
  int         Left, Right, Top, Bottom;
  BOOL        maximize = FALSE;
  
  // Use the nearest if the figure does not overlap with any monitor:
  hMonitor = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONEAREST);
  
  mInfo.cbSize = sizeof(mInfo);
  if (GetMonitorInfo(hMonitor, &mInfo) == 0)  {  // Failing is unlikely...
     mexErrMsgIdAndTxt(ERR_ID   "GetMonitorInfo",
                       ERR_HEAD "GetMonitorInfo failed [%d].", GetLastError());
  }
  
  // Not working for maximized window:
  if (IsZoomed(hWnd)) {
     ShowWindow(hWnd, SW_RESTORE);
     maximize = TRUE;
  }
  
  // Get work rect of monitor and rect of window:
  Mon_Rect = mInfo.rcWork;
  GetWindowRect(hWnd, &Win_Rect);
  Left   = Win_Rect.left;
  Top    = Win_Rect.top;
  Right  = Win_Rect.right;
  Bottom = Win_Rect.bottom;
  
  // Move the window completely to the screen. If the window is larger than the
  // work area, the top left corner is visible:
  if (Right > Mon_Rect.right)
     Left -= Right - Mon_Rect.right;

  if (Left < Mon_Rect.left)
     Left = Mon_Rect.left;
  
  if (Bottom > Mon_Rect.bottom)
     Top -= Bottom - Mon_Rect.bottom;

  if (Top < Mon_Rect.top)
     Top = Mon_Rect.top;
    
  // Set the adjusted window:
  SetWindowPos(hWnd, HWND_TOP, Left, Top, 0, 0,
               SWP_SHOWWINDOW | SWP_NOSENDCHANGING |
               SWP_NOACTIVATE | SWP_NOZORDER | SWP_NOSIZE);
  
  if (maximize)
     ShowWindow(hWnd, SW_MAXIMIZE);
  
  return;
}

// *****************************************************************************
void myWindowRgn(HWND hWnd, const double Fig_Handle, const mxArray *Value) {
  // Show only the part of the window inside a RECT, which is defined relative
  // to the current window position. If Value is omitted (then it is NULL here),
  // the border and titlebar is clip.
  // Multi-monitor proof: Coordinates are relative to window.
   
  RECT   Rect;
  HRGN   WinRgn;
  int    PosCorrect[4], Left, Top, Right, Bottom;
  BOOL   Enable;
  mwSize lenValue;
  double *p;
  
  // Get difference between inner and outer position and window's RECT:
  GetPosCorrect(Fig_Handle, PosCorrect);
  GetWindowRect(hWnd, &Rect);
  
  // Get Value:
  if (Value == NULL) {              // WindowAPI(FigH, 'clip'):
     Enable = TRUE;
     Left   = -PosCorrect[0];
     Top    = PosCorrect[3] + PosCorrect[1];
     Right  = Rect.right    - Rect.left + PosCorrect[0];
     Bottom = Rect.bottom   - Rect.top  + PosCorrect[1];
     
  } else {
     lenValue = mxGetNumberOfElements(Value);
     if ((mxIsLogical(Value) || mxIsNumeric(Value)) && lenValue == 1) {
        // Value is TRUE, FALSE, 0 or any number:
        Enable = (mxGetScalar(Value) != 0.0);
        if (Enable) {
           Left   = -PosCorrect[0];
           Top    = PosCorrect[3] + PosCorrect[1];
           Right  = Rect.right    - Rect.left + PosCorrect[0];
           Bottom = Rect.bottom   - Rect.top  + PosCorrect[1];
        }
        
     } else if (mxIsDouble(Value) && lenValue == 4) {
        Enable = TRUE;
        p      = mxGetPr(Value);
        Left   = ((int) p[0]) - PosCorrect[0] - 1;
        Right  = Left + (int) p[2];
        Bottom = Rect.bottom - Rect.top - ((int) p[1]) + PosCorrect[1] + 1;
        Top    = Bottom - (int) p[3];
        
        // Reject negative width or height:
        if (Left > Right || Bottom < Top)
           mexErrMsgIdAndTxt(ERR_ID  "InvalidPosition",
                             ERR_HEAD "(Clip): Invalid position.");
        
     } else {
        mexErrMsgIdAndTxt(ERR_ID  "BadValueLength",
                  ERR_HEAD "(Clip): [Value] must be a flag or [1 x 4] double.");
     }
  }
  
  // Set or remove the mask:
  if (Enable) {
     WinRgn = CreateRectRgn(Left, Top, Right, Bottom);
     SetWindowRgn(hWnd, WinRgn, TRUE);
     DeleteObject(WinRgn);
  }
  else      // Remove the clip region:
     SetWindowRgn(hWnd, NULL, TRUE);
  
  return;
}
        
// *****************************************************************************
mxArray *myMonitorInfo(HWND hWnd) {
  // Get info for the monitor with largest overlap or nearest to the figure.
  // For the secondary monitor FullPosition equals WorkPosition.
  // Multi-monitor proof.

  HMONITOR    hMonitor;
  MONITORINFO mInfo;
  RECT        Rect;
  double      MonitorIndex = 1.0;
  mxLogical   isOnScreen = TRUE;
  const char  *Field[4] = {"Position", "WorkPosition",
                           "FigureOnScreen", "MonitorIndex"};
  mxArray     *FullPos, *WorkPos, *Reply;
  double      *p;
  int         PrimaryHeight = GetSystemMetrics(SM_CYSCREEN);
 
  // Get monitor with the largest overlapping area:
  hMonitor = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONULL);
  
  if (hMonitor == NULL)  {  // No overlap with any visible monitor:
     hMonitor   = MonitorFromWindow(hWnd, MONITOR_DEFAULTTONEAREST);
     isOnScreen = FALSE;
  }
  
  // Get the info for this monitor:
  mInfo.cbSize = sizeof(mInfo);
  if (GetMonitorInfo(hMonitor, &mInfo) == 0)  {  // Failing is unlikely...
     mexErrMsgIdAndTxt(ERR_ID   "MonitorInfo",
                       ERR_HEAD "GetMonitorInfo failed [%d].", GetLastError());
  }
  
  // Get monitor index - shortcut for primary monitor:
  if ((mInfo.dwFlags & MONITORINFOF_PRIMARY) == 0) {
     MonitorIndex = (double) myGetIndexFromMonitor(hMonitor);
  }

  // Create output values:
  FullPos = mxCreateDoubleMatrix(1, 4, mxREAL);
  WorkPos = mxCreateDoubleMatrix(1, 4, mxREAL);
  
  // Copy RECT values [Left,Top,Right,Bottom] measured from top left:
  // WorkRect = mxCreateDoubleMatrix(1, 4, mxREAL);
  //p    = mxGetPr(WorkRect);   // Monitor rect without taskbar and sidebar
  //Rect = mInfo.rcWork;
  //p[0] = Rect.left;
  //p[1] = Rect.top;
  //p[2] = Rect.right;
  //p[3] = Rect.bottom;
  
  // FullRect = mxCreateDoubleMatrix(1, 4, mxREAL);
  //p    = mxGetPr(FullRect);   // Full monitor rect
  //Rect = mInfo.rcMonitor;
  //p[0] = Rect.left;
  //p[1] = Rect.top;
  //p[2] = Rect.right;
  //p[3] = Rect.bottom;
  
  // MATLAB measures the secondary monitors from the bottom of the primary one.
  // This has a distinct level of logic, but is not intuitive.
  
  // Monitor position as [X, Y, Width, Height] measured from bottom left:
  Rect = mInfo.rcMonitor;
  p    = mxGetPr(FullPos);
  p[0] = Rect.left + 1;
  p[1] = PrimaryHeight - Rect.bottom + 1;  // Y measured from bottom or primary
  p[2] = Rect.right    - Rect.left;        // Width
  p[3] = Rect.bottom   - Rect.top;         // Height
  
  // Monitor position as [X, Y, Width, Height] measured from bottom left:
  Rect = mInfo.rcWork;
  p    = mxGetPr(WorkPos);
  p[0] = Rect.left + 1;
  p[1] = PrimaryHeight - Rect.bottom + 1;  // Y measured from bottom
  p[2] = Rect.right    - Rect.left;        // Width
  p[3] = Rect.bottom   - Rect.top;         // Height

  // Create output struct:
  Reply = mxCreateStructMatrix(1, 1, 4, Field);
  mxSetFieldByNumber(Reply, 0, 0, FullPos);
  mxSetFieldByNumber(Reply, 0, 1, WorkPos);
  mxSetFieldByNumber(Reply, 0, 2, mxCreateLogicalScalar(isOnScreen));
  mxSetFieldByNumber(Reply, 0, 3, mxCreateDoubleScalar(MonitorIndex));
  
  return Reply;
}

// *****************************************************************************
int myGetIndexFromMonitor(HMONITOR hMonitor) {
  // Get the index of the specified monitor.
  // I expect this function to find any monitor inside the virtual screen.
  // But I do not have experiences with multi-monitor setups and I hope that
  // this is not confused by non-display devices.
  
  GETMONITOR getMonitor;
  getMonitor.target = hMonitor;
  getMonitor.index  = 0;
  getMonitor.found  = FALSE;
  
  // Is it helpful to limit the search to the virtual screen coordinates???
  // int WINAPI GetSystemMetrics(int nIndex);
  // SM_XVIRTUALSCREEN  SM_CXVIRTUALSCREEN    X Width
  // SM_YVIRTUALSCREEN  SM_CYVIRTUALSCREEN    Y Height
  // Number of display devices: GetSystemMetrics(SM_CMONITORS)
  
  EnumDisplayMonitors(NULL, NULL, myMonitorMatchProc, (LPARAM) &getMonitor);
  
  // Return the monitor index of 1 as primary monitor as fallback:
  if (getMonitor.found) {
     return getMonitor.index;
  }

  // This should never happen!
  mexWarnMsgIdAndTxt(ERR_ID    "MonitorNotFound",
                     WARN_HEAD "Cannot find monitor.");
  return 0;
}

// *****************************************************************************
BOOL CALLBACK myMonitorMatchProc(HMONITOR hMonitor, HDC hdcMonitor,
                                 LPRECT lprcMonitor, LPARAM dwData) {
  // Increase the counter and stop on matching hMonitor.
   
  LPGETMONITOR getMonitorP = (LPGETMONITOR) dwData;
  (getMonitorP->index)++;
  
  if (hMonitor == getMonitorP->target) {
     getMonitorP->found = TRUE;
     return FALSE;      // Stop enumeration
  }
  
  return TRUE;
}

// *****************************************************************************
HMONITOR myGetMonitorFromIndex(HWND hWnd, int MonitorIndex) {
  // Get the monitor with specified index. Use the primary monitor as fallback.
  POINT      Origin = {0, 0};
  GETMONITOR getMonitor;
  
  getMonitor.target = (HMONITOR) 0;   // Dummy value
  getMonitor.index  = MonitorIndex;   // Counter
  getMonitor.found  = FALSE;          // Flag
  
  // Search the matching monitor:
  EnumDisplayMonitors(NULL, NULL, myIndexMatchProc, (LPARAM) &getMonitor);
  
  if (getMonitor.found) {
     return getMonitor.target;
  }
  
  // Reply the primary monitor as fallback:
  return MonitorFromPoint(Origin, MONITOR_DEFAULTTOPRIMARY);
}

// *****************************************************************************
BOOL CALLBACK myIndexMatchProc(HMONITOR hMonitor, HDC hdcMonitor,
                               LPRECT lprcMonitor, LPARAM dwData) {
  // Decrease the counter and stop at zero.
   
  LPGETMONITOR getMonitorP = (LPGETMONITOR) dwData;
  
  if (--(getMonitorP->index) == 0) {
     getMonitorP->found  = TRUE;
     getMonitorP->target = hMonitor;
     return FALSE;      // Stop enumeration
  }
  
  return TRUE;
}

// *****************************************************************************
void myClipCursor(HWND hWnd, const double Fig_Handle, const mxArray *Param) {
  // Restrict cursor motion to window rectangle.
  // Called "LockCursor" instead of the Windows API 'ClipCursor' to reduce the
  // confusion with 'ClipRegion'.
   
  RECT   Rect, *lpRect = NULL;
  int    PosCorrect[4];
  double *p;
  
  if (Param != NULL) {
     switch (mxGetNumberOfElements(Param)) {
        case 0:
           break;
           
        case 1:  // WindowAPI(FigH, 'LockCursor', 1):
           if (!mxIsNumeric(Param) && !mxIsLogical(Param)) {
              mexWarnMsgIdAndTxt(ERR_ID "BadParam",
                      WARN_HEAD
                      "(LockCursor): Parameter must be numerical or LOGICAL.");
           }
           
           // Clip cursor to full window if input is neither FALSE nor 0:
           if (mxGetScalar(Param) != 0.0) {
              GetWindowRect(hWnd, &Rect);
              lpRect = &Rect;
           }
           break;
           
        case 4:  // WindowAPI(FigH, 'LockCursor', [X, Y, Width, Height]):
           // Accept a DOUBLE vector only:
           if (!mxIsDouble(Param)) {
              mexWarnMsgIdAndTxt(ERR_ID "BadParam",
                      WARN_HEAD
                      "(LockCursor): Rect must be a [1 x 4] DOUBLE vector.");
           }
           p = mxGetPr(Param);
           
           // Get window rect:
           GetPosCorrect(Fig_Handle, PosCorrect);
           GetWindowRect(hWnd, &Rect);
           
           // Calculate input area relative to window rect:
           Rect.left   += (int) p[0]  - PosCorrect[0] - 1;
           Rect.bottom -= (int) p[1]  - PosCorrect[1] - 1;
           Rect.right   = Rect.left   + (int) p[2];
           Rect.top     = Rect.bottom - (int) p[3];
           lpRect       = &Rect;
           break;
           
        default:
           mexWarnMsgIdAndTxt(ERR_ID    "BadRect",
                   WARN_HEAD
                   "(LockCursor): Parameter must be scalar or [1 x 4] vector");
     }
  }
  
  // Apply the clipping of the cursor:
  ClipCursor(lpRect);
  
  return;
}
