from libcpp.vector cimport vector

cdef extern from "eirene.cpp":
        cdef cppclass SparseBool:
                vector[long] rowVals
                vector[long] colPtrs
                long numRows
        
        ctypedef struct MorsePair:
                long lowInd
                long highInd

        ctypedef struct PyreneResult:
                vector[vector[unsigned long long]] tid
                vector[vector[vector[unsigned long long]]] cycleReps
                vector[vector[double]] birthDeath

        PyreneResult pyVripsPersistence(double* dMat, long numPoints, double maxRad)
        PyreneResult pyCellularPersistence(unsigned long long* rowVals, unsigned long long numVals, unsigned long long* colPtrs, 
                unsigned long long numCols, unsigned long long numRows, double* filtVals, size_t filtLen, unsigned long long* dimVec, size_t numDim)
        PyreneResult pyRelativeCellularPersistence(unsigned long long* rowVals, unsigned long long numVals, unsigned long long* colPtrs, 
                unsigned long long numCols, unsigned long long numRows, double* filtVals, size_t filtLen, unsigned long long* dimVec, size_t numDim, unsigned long long* openComp, unsigned long long compLen)