//////////////////////////////////////////////////////////////////////////////
// Copyright 1993-2007 The MathWorks, Inc.
// $Revision: 1.1.6.4 $
//////////////////////////////////////////////////////////////////////////////

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include <string.h>
#include "mex.h"
#include "mwsize.h"
#include "iptutil_cpp.h"

//Values used for marking the starting and boundary pixels.
#define START_DOUBLE     -2
#define BOUNDARY_DOUBLE  -3

#define START_UINT8       2
#define BOUNDARY_UINT8    3

#define INITIAL_SCRATCH_LENGTH 200

#define INVALID -10

typedef enum {
    counterclockwise = 0,
    clockwise = 1
} tdir_T;

class Boundaries {
  private:   

    //Member variables
    //////////////////
    mwSize fOrigNumRows;
    mwSize fNumCols;
    mwSize fNumRows;
    int fConn;

    //Variables used by the tracing algorithm
    mwSignedIndex  fOffsets[8];
    mwSignedIndex  fVOffsets[8];
    int     *fNextDirectionLut;
    int     *fNextSearchDirectionLut;
    int      fNextSearchDir;
    mwSize   fScratchLength;
    mwSize  *fScratch;
    int      fStartMarker;
    int      fBoundaryMarker;

    //Methods
    /////////
    void initTraceLUTs(mxClassID classID, tdir_T    dir       = clockwise);

    //The defaults for this method are used by findBoundaries()
    void setNextSearchDirection(uint8_T *bwImage    = NULL, 
                                mwSize   idx        = 0,
                                int      firstStep  = 0,
                                tdir_T   dir        = clockwise);

    bool isBoundaryPixel(uint8_T *bwImage, mwSize idx);
    inline mwSize calculateIdxIntoOrig(mwSize idx);
    void copyCoordsToBuf(mwSize numPixels,double *buf);
    double calculateNumRegions(double *matrix, mwSize numElements);
    mxArray *createOutputArray(double numRegions);
    void expandScratchSpace(void);

    //Templetized Methods

    //////////////////////////////////////////////////////////////////////////
    //Adds a border of zeros to the input image so that we will not have 
    //to worry about going out of bounds while tracing objects 
    //////////////////////////////////////////////////////////////////////////
    template<typename _T>
        _T *padBuffer(_T *buff) {
            _T *paddedBuf = (_T *)mxCalloc(fNumRows*fNumCols, sizeof(_T));
            for(mwSize i=0; i < fNumCols-2; i++)
                memcpy(&paddedBuf[(i+1)*fNumRows+1], &buff[i*(fNumRows-2)], (fNumRows-2)*sizeof(_T));

            return(paddedBuf);
        }

    //////////////////////////////////////////////////////////////////////////
    // 
    //////////////////////////////////////////////////////////////////////////
    template<typename _T>
    bool isRegionBoundaryPixel(double *origLabelMatrix, _T *pLabelMatrix, 
                         mwSize idx, int currentLabel) {
        mxAssert( fVOffsets[0] != INVALID, "fVOffsets[] must be calculated");
        
        mwSize idxIntoOrig;
        bool isBoundaryPix = false;
        
        //First, make sure that it's not a background pixel, otherwise it's
        //an automatic failure.
        if(pLabelMatrix[idx]) {
            //make sure that pixel has the same label as the currentLabel,
            //otherwise it's an automatic failure.
            idxIntoOrig = calculateIdxIntoOrig(idx);
            if (origLabelMatrix[idxIntoOrig] == currentLabel)
                isBoundaryPix = true;
        }
        return(isBoundaryPix);
    }

    //////////////////////////////////////////////////////////////////////////
    //This method traces a single contour of a region surrounded by zeros.  It
    //takes a label matrix, linear index to the initial border pixel belonging
    //to the object that's going to be traced, and it needs to be
    //pre-configured by invoking initializeTracer routine. It returns mxArray
    //containing X-Y coordinates of border pixels.
    //////////////////////////////////////////////////////////////////////////
    template<typename _T>
        mxArray *traceBoundary(_T *bwImage, mwSize idx, mwSize maxNumPixels=MWSIZE_MAX) {
            mxAssert( (fConn != -1),
                      ERR_STRING("Boundaries::fConn","traceBoundary()"));

            mxAssert( (fStartMarker != -1),
                      ERR_STRING("Boundaries::fStartMarker", "traceBoundary()"));

            mxAssert( (fBoundaryMarker != -1),
                      ERR_STRING("Boundaries::fBoundaryMarker", "traceBoundary()"));
            

            //Initialize loop variables
            fScratch[0]           = idx;
            bwImage[idx]          = fStartMarker;
            bool done             = false;
            mwSize  numPixels     = 1;
            mwSize  currentPixel  = idx;
            int  nextSearchDir    = fNextSearchDir;
            int  initDepartureDir = -1;
            
            while(!done) {
                // Find the next boundary pixel.
                int  direction      = nextSearchDir;
                bool foundNextPixel = false;

                for(int k=0; k < fConn; k++) {
                    //Try to locate the next pixel in the chain
                    mwSize neighbor = currentPixel + fOffsets[direction];
                    
                    if(bwImage[neighbor]) {
                        //Found the next boundary pixel.
                        if(bwImage[currentPixel] == fStartMarker && initDepartureDir == -1) {
                            //We are making the initial departure from the
                            //starting pixel.
                            initDepartureDir = direction;
                        }
                        else if(bwImage[currentPixel] == fStartMarker && initDepartureDir == direction) {
                            //We are about to retrace our path.
                            //That means we're done.
                            done = true;
                            foundNextPixel = true;
                            break;
                        }
                        
                        //Take the next step along the boundary.
                        nextSearchDir = fNextSearchDirectionLut[direction];
                        foundNextPixel = true;
                        
                        if(fScratchLength <= numPixels)
                            expandScratchSpace();
                        
                        //First use numPixels as an index into scratch array,
                        //then increment it
                        fScratch[numPixels] = neighbor;
                        numPixels++;

                        //note to self: numPixels at this point will
                        //be at least 2;
                        if(numPixels == maxNumPixels) {
                            done = true;
                            break;
                        }
                        
                        if(bwImage[neighbor] != fStartMarker)
                            bwImage[neighbor] = fBoundaryMarker;
                        
                        currentPixel = neighbor;
                        break;
                    }
                    
                    direction = fNextDirectionLut[direction];
                }
                
                if (!foundNextPixel) {
                    //If there is no next neighbor, the region must 
                    //just have a single pixel.
                    numPixels = 2;
                    fScratch[1] = fScratch[0];
                    done = true;
                }
            }

            mxArray *boundary = mxCreateDoubleMatrix(numPixels,2,mxREAL);
            //Stuff the array with proper data
            copyCoordsToBuf(numPixels,(double *)mxGetData(boundary));
            
            return boundary;
        }    

    //////////////////////////////////////////////////////////////////////////
    //This method traces a single contour of a region based on its label (more
    //work than traceBoundary). It takes the original label matrix, the padded
    //label matrix, a linear index to the initial border pixel belonging to
    //the object that's going to be traced, and the current pixel label. This
    //method needs to be pre-configured by invoking initializeTracer
    //routine. It returns mxArray containing X-Y coordinates of border pixels.
    //
    // traceRegionBoundary() and traceBoundary() share a lot of code but it is
    // tricky to remove the common code into a separate method. Combining both
    // methods into one is prettier but it reduces the speed of
    // bwtraceboundary.m and bwboundaries.m in my sandbox. Revisit this issue
    // once this code is in Aimage.  (isimon 17 Nov 2003)
    //////////////////////////////////////////////////////////////////////////
    template<typename _T>
        mxArray *traceRegionBoundary(double *origLabelMatrix,_T *bwImage, mwSize idx, int currentLabel) {
            mxAssert( (fConn != -1),
                      ERR_STRING("Boundaries::fConn","traceBoundary()"));

            mxAssert( (fStartMarker != -1),
                      ERR_STRING("Boundaries::fStartMarker", "traceBoundary()"));

            mxAssert( (fBoundaryMarker != -1),
                      ERR_STRING("Boundaries::fBoundaryMarker", "traceBoundary()"));


            //Initialize loop variables
            fScratch[0]           = idx;
            bwImage[idx]          = fStartMarker;
            bool done             = false;
            mwSize  numPixels        = 1;
            mwSize  currentPixel     = idx;
            int  nextSearchDir    = fNextSearchDir;
            int  initDepartureDir = -1;

            while(!done) {
                // Find the next boundary pixel.
                int  direction      = nextSearchDir;
                bool foundNextPixel = false;

                for(int k=0; k < fConn; k++) {
                    //Try to locate the next pixel in the chain
                    mwSize neighbor = currentPixel + fOffsets[direction];
                    
                    if(isRegionBoundaryPixel(origLabelMatrix,bwImage, neighbor,currentLabel)) {
                        //Found the next boundary pixel.
                        if(bwImage[currentPixel] == fStartMarker && initDepartureDir == -1) {
                            //We are making the initial departure from the
                            //starting pixel.
                            initDepartureDir = direction;
                        }
                        else if(bwImage[currentPixel] == fStartMarker && initDepartureDir == direction) {
                            //We are about to retrace our path.
                            //That means we're done.
                            done = true;
                            foundNextPixel = true;
                            break;
                        }
                        
                        //Take the next step along the boundary.
                        nextSearchDir = fNextSearchDirectionLut[direction];
                        foundNextPixel = true;
                        
                        if(fScratchLength <= numPixels)
                            expandScratchSpace();
                        
                        //First use numPixels as an index into scratch array,
                        //then increment it
                        fScratch[numPixels] = neighbor;
                        numPixels++;
                        
                        if(bwImage[neighbor] != fStartMarker)
                            bwImage[neighbor] = fBoundaryMarker;
                        
                        currentPixel = neighbor;
                        break;
                    }
                    direction = fNextDirectionLut[direction];
                }
                
                if (!foundNextPixel)
                {
                    //If there is no next neighbor, the region must 
                    //just have a single pixel.
                    numPixels = 2;
                    fScratch[1] = fScratch[0];
                    done = true;
                }
            }

            mxArray *boundary = mxCreateDoubleMatrix(numPixels,2,mxREAL);

            //Stuff the array with proper data
            copyCoordsToBuf(numPixels,(double *)mxGetData(boundary));
            return boundary;
        }    

 public:
    
    Boundaries();
    virtual ~Boundaries();
    
    mxArray *findRegionBoundaries(double *matrix);
    mxArray *findBoundaries(double *matrix);
    mxArray *traceBoundary(uint8_T *bwImage, mwSize startRow, 
                           mwSize startCol, int firstStep, 
                           tdir_T dir, mwSize maxNumPoints);
    
    //Accessor methods
    //////////////////
    void setNumCols(mwSize numCols);
    void setNumRows(mwSize numRows);
    void setConnectivity(int conn);
    void setOrigLabelMatrix(const mxArray *matrix);
    void setOrigNumRows(mwSize num);
};

#endif
