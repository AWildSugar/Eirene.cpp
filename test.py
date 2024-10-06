import numpy as np
import pyrene
from scipy.sparse import csc_array

# a vr-comlex from a point cloud
a = np.array([[1,0,1,0],
         [1,1,0,0],
         [1,0,1,0],
         [0,0,1,1]])

b = a.reshape(a.shape[0], 1, a.shape[1])

distmat = np.sqrt(np.einsum('ijk, ijk->ij', a-b, a-b))

vripsRes = pyrene.vripsPersistence(distmat, a.shape[0])
print(vripsRes)


# relative persistence on an nsphere 
numCells = 6

denseArr = np.zeros((numCells, numCells))
denseArr[0, 2] = 1
denseArr[0, 3] = 1
denseArr[1, 2] = 1
denseArr[1, 3] = 1
denseArr[2, 4] = 1
denseArr[2, 5] = 1
denseArr[3, 4] = 1
denseArr[3, 5] = 1

sparseComp = csc_array(denseArr, dtype=np.ulonglong)

filtVals = np.array([0.01 * (i + 1) for i in range(numCells)], dtype=np.double)
dimVec = np.array([2, 2, 2], dtype=np.ulonglong)
openComp = np.array([1, 2], dtype=np.ulonglong)
indPtr = np.array(sparseComp.indptr, dtype=np.ulonglong)
rowInds = np.array(sparseComp.indices, dtype=np.ulonglong)

cellRes = pyrene.relativePersistence(indPtr, rowInds, numCells, filtVals, dimVec, openComp)
print(cellRes)
