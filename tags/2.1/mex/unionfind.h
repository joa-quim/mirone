/* 
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.4.4.2 $
 */

#ifndef UNIONFIND_H
#define UNIONFIND_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mex.h"
#include "mwsize.h"

typedef struct ufType {
    mwSize *id;
    mwSize *sz;
    mwSize allocated_length;
    mwSize num_nodes;
    mwSize num_sets;
    int finalized;
}
ufType;

extern ufType *uf_init(mwSize hint);
void uf_new_node(ufType *uf);
void uf_new_pair(ufType *uf, mwSize p, mwSize q);
mwSize uf_renumber(ufType *uf, mwSize first);
mwSize uf_query_set(ufType *uf, mwSize p);
void uf_destroy(ufType *uf);

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* UNIONFIND_H */
