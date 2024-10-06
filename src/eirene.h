#ifndef EIRENE_H
#define EIRENE_H

#include <stdint.h>
#include "indsz.hpp"

#include <stdio.h>
#include <stdlib.h>

extern "C" 
{

typedef struct MorsePair {
    indsz_t lowInd;
    indsz_t highInd;
} MorsePair;

typedef struct CPairArr {
    size_t len;
    MorsePair* data;
} CPairArr;

typedef struct U32Arr {
    size_t len;
    indsz_t* data;
} U32Arr;

typedef struct CycleRepArr {
    size_t len;
    U32Arr* data;
} CycleRepsArr;

typedef struct CSparseBool {
    U32Arr rowVals;
    U32Arr colPtrs;
    size_t numRows;
} CSparseBool;

typedef struct CMorseReducedResult
{
    CSparseBool leftFactor;   
    U32Arr tid;
    CPairArr pairArr;
    CycleRepArr cycleReps;
} CMorseReducedResult;

typedef struct CMorseReducedResultArr 
{
    size_t len;
    CMorseReducedResult* data;
} CMorseReducedResultArr;

inline void freeCycleRepArr(CycleRepArr arr)
{
    for (size_t i = 0; i < arr.len; ++i) {
        free(arr.data[i].data);
    }
}

inline void freeSparseBool(CSparseBool sps)
{
    free(sps.rowVals.data);
    free(sps.colPtrs.data);
}

inline void freeMorseReducedResultArr(CMorseReducedResultArr res) 
{
    for (size_t i = 0; i < res.len; ++i) {
        free(res.data[i].tid.data);
        free(res.data[i].pairArr.data);
        freeCycleRepArr(res.data[i].cycleReps);
        freeSparseBool(res.data[i].leftFactor);
    }
}

CMorseReducedResultArr cellularPersistence(CSparseBool* chainComplex, double* filtVals, size_t filtLen, indsz_t* dimVec, size_t numDim);
}
#endif
