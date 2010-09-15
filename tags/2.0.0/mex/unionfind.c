/* 
 * Copyright 1993-2007 The MathWorks, Inc.
 * $Revision: 1.4.4.2 $
 */

#include "unionfind.h"
#include "mex.h"
#include "mwsize.h"

#define DEFAULT_ALLOCATED_LENGTH 16
#define REALLOC_FACTOR 2

ufType *uf_init(mwSize hint) {
    ufType *uf = NULL;

    if (hint == 0)
        hint = DEFAULT_ALLOCATED_LENGTH;

    uf = (ufType *) mxCalloc(1, sizeof(*uf));
    uf->id = (mwSize *) mxCalloc(hint, sizeof(*(uf->id)));
    uf->sz = (mwSize *) mxCalloc(hint, sizeof(*(uf->sz)));
    uf->allocated_length = hint;
    uf->num_nodes = 0;
    uf->num_sets = 0;
    uf->finalized = false;

    return uf;
}

void uf_new_node(ufType *uf) {
    if (uf->num_nodes >= uf->allocated_length) {
        uf->allocated_length *= REALLOC_FACTOR;
        uf->id = (mwSize *) mxRealloc(uf->id, uf->allocated_length * sizeof(mwSize));
        uf->sz = (mwSize *) mxRealloc(uf->sz, uf->allocated_length * sizeof(mwSize));
    }

    uf->num_nodes++;
    uf->id[uf->num_nodes - 1] = uf->num_nodes - 1;
    uf->sz[uf->num_nodes - 1] = 1;
}

mwSize uf_find(ufType *uf, mwSize p) {
    mwSize i;
    mwSize t;
    mwSize *id = uf->id;

    mxAssert(p < uf->num_nodes,"");
    
    for (i = p; i != id[i]; i = id[i]) {
        t = i;
        i = id[id[t]];
        id[t] = i;
    }

    return i;
}

static void uf_union(ufType *uf, mwSize p, mwSize q) {
    if (uf->sz[p] < uf->sz[q]) {
        uf->id[p] = q;
        uf->sz[q] += uf->sz[p];
    }
    else
    {
        uf->id[q] = p;
        uf->sz[p] += uf->sz[q];
    }
}

void uf_new_pair(ufType *uf, mwSize p, mwSize q) {
    mwSize i;
    mwSize j;

    mxAssert(p >= 0,"");
    mxAssert(q >= 0,"");
    mxAssert(uf,"");
    mxAssert(! uf->finalized,"");

    i = uf_find(uf, p);
    j = uf_find(uf, q);
    if (i != j)
	uf_union(uf, i, j);
}

mwSize uf_renumber(ufType *uf, mwSize first) {
    mwSize k;
    mwSize counter = first;

    uf->finalized = true;

    for (k = 0; k < uf->num_nodes; k++) {
        if (uf->id[k] == k)
            uf->sz[k] = counter++;
    }

    uf->num_sets = counter - first;

    return uf->num_sets;
}

mwSize uf_query_set(ufType *uf, mwSize p) {
    mwSize k;

    mxAssert(uf->finalized,"");

    k = uf_find(uf, p);
    return uf->sz[k];
}

void uf_destroy(ufType *uf) {
    mxAssert(uf, "");

    mxFree(uf->id);
    mxFree(uf->sz);
}
