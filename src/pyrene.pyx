# distutils: language = c++

cimport numpy as np
cimport src.pyrene as pyrene
import cython
import math
from scipy.sparse import csc_array

@cython.boundscheck(False)
@cython.wraparound(False)
def vripsPersistence(np.ndarray[double, ndim=2, mode="c"] dMat not None, long numPoints):
	cdef double* dPtr = &dMat[0, 0]
	return pyrene.pyVripsPersistence(dPtr, numPoints, math.inf)

@cython.boundscheck(False)
@cython.wraparound(False)
def relativePersistence(np.ndarray[unsigned long long, ndim=1, mode="c"] colInds not None, np.ndarray[unsigned long long, ndim=1, mode="c"] rowVals not None, numCells, np.ndarray[double, ndim=1, mode="c"] filtVals not None, np.ndarray[unsigned long long, ndim=1, mode="c"] dimVec not None, np.ndarray[unsigned long long, ndim=1, mode="c"] openComp not None):
	cdef unsigned long long* indPtr = &colInds[0]
	cdef unsigned long long* rowValsPtr = &rowVals[0]
	cdef unsigned long long* dimVecPtr = &dimVec[0]
	cdef unsigned long long* openCompPtr = &openComp[0]
	cdef double* filtValsPtr = &filtVals[0]
	return pyrene.pyRelativeCellularPersistence(rowValsPtr, len(rowVals), indPtr, len(colInds) - 1, numCells, filtValsPtr, len(filtVals), dimVecPtr, len(dimVec), openCompPtr, len(openComp))

@cython.boundscheck(False)
@cython.wraparound(False)
def cellularPersistence(np.ndarray[unsigned long long, ndim=1, mode="c"] colInds not None, np.ndarray[unsigned long long, ndim=1, mode="c"] rowVals not None, numCells, np.ndarray[double, ndim=1, mode="c"] filtVals not None, np.ndarray[unsigned long long, ndim=1, mode="c"] dimVec not None):
	print('Taking addresses')
	cdef unsigned long long* indPtr = &colInds[0]
	cdef unsigned long long* rowValsPtr = &rowVals[0]
	cdef unsigned long long* dimVecPtr = &dimVec[0]
	cdef double* filtValsPtr = &filtVals[0]
	print('Calling extension')
	return pyrene.pyCellularPersistence(rowValsPtr, len(rowVals), indPtr, len(colInds) - 1, numCells, filtValsPtr, len(filtVals), dimVecPtr, len(dimVec))