/* Copyright 1993-2007 The MathWorks, Inc. */

/*
 * First helper MEX-file for BWLABEL
 * 1. Computes run-length encoding for input binary image
 * 2. Assigns initial labels to runs while recording equivalences
 *
 * Syntax:
 *        [SR,ER,C,LABELS,I,J] = BWLABEL1(BW,MODE)
 *
 * BW is a 2-D uint8 array.  MODE is a scalar, either 8 or 4, indicating
 * the connectedness of the foreground of BW.
 *
 * SR, ER, C, and LABELS are column vectors whose length is the number
 * of runs.  SR contains the starting row for each run.  ER contains the
 * ending row for each run.  C contains the column for each run.  LABELS
 * contains the initial labels determined for each run.  I and J contain
 * labels equivalence information.  For example, if I(4) = 10 and J(4) = 20,
 * that means that labels 10 and 20 are equivalent.
 * 
 */

#include "mex.h"
#include "mwsize.h"

/*
 * Make sure input and output arguments are correct.
 */
void ValidateInputs(int nrhs, const mxArray *prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("Images:bwlabel1:wrongNumInputs",
                          "BWLABEL1 requires two input arguments.");
    }    
    if ((mxGetScalar(prhs[1]) != 4.0) && (mxGetScalar(prhs[1]) != 8.0)) {
        mexErrMsgIdAndTxt("Images:bwlabel1:badConn",
                          "%s",
                          "Function bwlabel1 expected its second input argument, CONN, to be either 4 or 8.");
    }
}

/*
 * Scan the input array, counting the number of runs present.
 */
mwSize NumberOfRuns(mxLogical *in, mwSize M, mwSize N) {
    mwSize result = 0;
    mwSize col;
    mxLogical *pr;
    mwSize k;
    mwSize offset;

    if ( (M != 0) && (N != 0) ) {
        for (col = 0; col < N; col++) {
            offset = col*M;
            pr = in + offset;
            
            if (pr[0] != 0)
                result++;
            
            for (k = 1; k < M; k++) {
                if ((pr[k] != 0) && (pr[k-1] == 0))
                    result++;
            }
        }
    }
    
    return(result);
}

/*
 * Scan the array, recording start row, end row, and column information
 * about each run.  The calling function must allocate sufficient space 
 * for sr, er, and c.
 */
void FillRunVectors(mxLogical *in, mwSize M, mwSize N, double *sr, double *er, double *c) {
    mwSize col;
    mwSize offset;
    mxLogical *pr;
    mwSize runCounter = 0;
    mwSize k;
        
    for (col = 0; col < N; col++) {
        offset = col*M;
        pr = in + offset;
        
        k = 0;
        while (k < M) {
            /* Look for the next run. */
            while ((k < M) && (pr[k] == 0))
                k++;
            
            if ((k < M) && (pr[k] != 0)) {
                c[runCounter] = (double) (col + 1);
                sr[runCounter] = (double) (k + 1);
                while ((k < M) && (pr[k] != 0))
                    k++;
                er[runCounter] = (double) k;
                runCounter++;
            }
        }
    }
}

typedef struct equivnode {
    mwSize rowIndex;
    mwSize colIndex;
    struct equivnode *next;
}
EquivNode;

typedef struct equivlist {
    EquivNode *head;
    mwSize length;
}
EquivList;

/*
 * Initialize equivalence linked list
 */
void InitEquivalenceList(EquivList *list) {
    list->head = NULL;
    list->length = 0;
}

/*
 * Allocate and add new node to head of equivalence list
 */
void PrependEquivalenceToList(EquivList *list, double rowIndex, double colIndex) {
    EquivNode *newNode = NULL;
    
    newNode = mxCalloc(1, sizeof(*newNode));
    newNode->next = list->head;
    newNode->rowIndex = (mwSize) rowIndex;
    newNode->colIndex = (mwSize) colIndex;
    list->head = newNode;
    list->length++;
}

/*
 * Remove and free the first node in the equivalence list
 */
void DeleteHeadNode(EquivList *list) {
    EquivNode *killNode;
    
    if (list->head != NULL) {
        killNode = list->head;
        list->head = killNode->next;
        mxFree((void *) killNode);
        list->length--;
    }
}

/*
 * Free the equivalence list
 */
void DestroyList(EquivList *list) {
    while (list->head != NULL)
        DeleteHeadNode(list);
}

/*
 * Scan the runs, assigning labels and remembering equivalences.
 * This function allocates space for the row and column equivalence index
 * vectors rowEquivalences and colEquivalences.
 */
void FirstPass(mwSize numRuns, int mode, double *sr, double *er, double *c, 
               double *labels, double **rowEquivalences, double **colEquivalences,
               mwSize *numEquivalences)
{
    mwSize k;
    mwSignedIndex p;
    mwSize offset;
    double currentColumn = 0;
    mwSize nextLabel = 1;
    mwSignedIndex firstRunOnPreviousColumn = -1;
    mwSignedIndex lastRunOnPreviousColumn = -1;
    mwSignedIndex firstRunOnThisColumn = -1;
    EquivList equivList;
    EquivNode *node;

    InitEquivalenceList(&equivList);

    if (mode == 8)
        /* This value is used in the overlap test below. */
        offset = 1;
    else
        offset = 0;
    
    for (k = 0; k < numRuns; k++) {
        /* Process k-th run */
        
        if (c[k] == (currentColumn + 1)) {
            /* We are starting a new column adjacent to previous column */
            firstRunOnPreviousColumn = firstRunOnThisColumn;
            firstRunOnThisColumn = k;
            lastRunOnPreviousColumn = k-1;
            currentColumn = c[k];
        }
        else if (c[k] > (currentColumn + 1)) {
            /* We are starting a new column not adjacent to previous column */
            firstRunOnPreviousColumn = -1;
            lastRunOnPreviousColumn = -1;
            firstRunOnThisColumn = k;
            currentColumn = c[k];
        }
        else {
            /* Not changing columns; nothing to do here */
        }
        
        if (firstRunOnPreviousColumn >= 0) {
            /*
             * Look for overlaps on previous column
             */
            
            p = firstRunOnPreviousColumn;
            while ((p <= lastRunOnPreviousColumn) && (sr[p] <= (er[k] + offset))) {
                if ((er[k] >= (sr[p]-offset)) && (sr[k] <= (er[p]+offset))) {
                    /* 
                     * We've got an overlap; it's 4-connected or 8-connected
                     * depending on the value of offset.
                     */
                    if (labels[k] == 0) {
                        /* 
                         * This run hasn't yet been labeled;
                         * copy over the overlapping run's label
                         */
                        labels[k] = labels[p];
                    }
                    else {
                        if (labels[k] != labels[p]) {
                            /* This run and the overlapping run
                             * have been labeled with different
                             * labels.  Remember the equivalence.
                             */
                            PrependEquivalenceToList(&equivList, labels[k], labels[p]);
                        }
                        else {
                            /* This run and the overlapping run
                             * have been labeled with the same label;
                             * nothing to do here.
                             */
                        }
                        
                    }
                }
                p++;
            }
        }
        
        if (labels[k] == 0) {
            /*
             * This run hasn't yet been labeled because we
             * didn't find any overlapping runs.  Label it
             * with a new label.
             */
            labels[k] = (double) nextLabel;
            nextLabel++;
        }
    }
    
    *numEquivalences = equivList.length;
    if (*numEquivalences > 0) {
        *rowEquivalences = (double *) mxCalloc(equivList.length, sizeof(double));
        *colEquivalences = (double *) mxCalloc(equivList.length, sizeof(double));

        /* 
         * Traverse the equivalence list, recording indices in the
         * output arrays.
         */
        k = 0;
        node = equivList.head;
        while (node != NULL) {
            (*rowEquivalences)[k] = (double) node->rowIndex;
            (*colEquivalences)[k] = (double) node->colIndex;
            k++;
            node = node->next;
        }
    }
    else {
        *rowEquivalences = NULL;
        *colEquivalences = NULL;
    }
    
    DestroyList(&equivList);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *sr;
    double *er;
    double *c;
    double *labels;
    double *rowEquivalences;
    double *colEquivalences;
    double *pr;
    mwSize numEquivalences;
    mwSize numRuns;
    int mode;
    mwSize k;
    const mxArray *BW;
    
    ValidateInputs(nrhs, prhs);
    BW = prhs[0];
    mode = (int) mxGetScalar(prhs[1]);

    numRuns = NumberOfRuns((mxLogical *) mxGetData(BW), mxGetM(BW), mxGetN(BW));
    sr = mxCalloc(numRuns, sizeof(*sr));
    er = mxCalloc(numRuns, sizeof(*er));
    c = mxCalloc(numRuns, sizeof(*c));
    labels = mxCalloc(numRuns, sizeof(*labels));
        
    FillRunVectors((mxLogical *) mxGetData(BW), mxGetM(BW), mxGetN(BW),
                   sr, er, c);

    FirstPass(numRuns, (int) mode, sr, er, c, labels, &rowEquivalences, 
              &colEquivalences, &numEquivalences);

    /* First output argument */
    plhs[0] = mxCreateDoubleMatrix(numRuns, 1, mxREAL);
    pr = (double *) mxGetData(plhs[0]);
    for (k = 0; k < numRuns; k++)
        pr[k] = sr[k];
    
    /* Second output argument */
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(numRuns, 1, mxREAL);
        pr = (double *) mxGetData(plhs[1]);
        for (k = 0; k < numRuns; k++)
            pr[k] = er[k];
    }

    /* Third output argument */
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleMatrix(numRuns, 1, mxREAL);
        pr = (double *) mxGetData(plhs[2]);
        for (k = 0; k < numRuns; k++)
            pr[k] = c[k];
    }

    /* Fourth output argument */
    if (nlhs > 3) {
        plhs[3] = mxCreateDoubleMatrix(numRuns, 1, mxREAL);
        pr = (double *) mxGetData(plhs[3]);
        for (k = 0; k < numRuns; k++)
            pr[k] = labels[k];
    }

    /* Fifth output argument */
    if (nlhs > 4) {
        plhs[4] = mxCreateDoubleMatrix(numEquivalences, 1, mxREAL);
        pr = (double *) mxGetData(plhs[4]);
        for (k = 0; k < numEquivalences; k++)
            pr[k] = rowEquivalences[k];
    }

    /* Sixth output argument */
    if (nlhs > 5) {
        plhs[5] = mxCreateDoubleMatrix(numEquivalences, 1, mxREAL);
        pr = (double *) mxGetData(plhs[5]);
        for (k = 0; k < numEquivalences; k++)
            pr[k] = colEquivalences[k];
    }

    mxFree(sr);
    mxFree(er);
    mxFree(c);
    mxFree(labels);
    mxFree(rowEquivalences);
    mxFree(colEquivalences);
}
