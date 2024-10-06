#pragma once

#include <variant>
#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cstdlib>
#include <optional>
#include <fmt/core.h>

#include "sparsebool.hpp"
#include "util.hpp"
#include "eirene.h"

namespace eirene
{

using PairVec = std::vector<MorsePair>;

typedef struct 
{
    spla::SparseBool leftFactor;
    std::vector<indsz_t> tid;
    PairVec pairVec;
    std::vector<std::vector<indsz_t>> cycleReps;
    std::vector<double> birthDeath; // even size
} MorseReducedResult;

typedef struct {
    std::vector<std::vector<indsz_t>> tid;
    std::vector<std::vector<std::vector<indsz_t>>> cycleReps;
    std::vector<std::vector<double>> birthDeath; // even size
} PyreneResult;

struct Complex;

struct ComplexInclusion 
{
    /*
     * Contains the information of how a subcomplex is included in it's 'parent'.
     * Cell i in the subcomplex maps to incl[i] in the parent.
     */

    std::shared_ptr<eirene::Complex> parent;
    std::vector<indsz_t> incl;
};

struct Complex
{
    /*
     *   A persistence-filtered (cellular) chain complex
     *   splitCells:  A vector sparse matrix, representing the boundaries of the CW complex. Entry [i,j] is true iff cell i is
     *                  codimension-1 coface of cell j. A row represents all the cells of which i is a (co-dim 1)
     *                  coface, and a column represents all the the cells which are (co-dim 1) faces of cell j.
     *   splitFiltInds: Vector of filtration values for the cells, in each dimension. Each value is an index into uniqFiltVals.
     *                  Smaller indices represent larger birth times.
     *   dimPatt:   Vector whose ith entry is the number of cells with dimension i ('Euler Vector')
    */

    std::vector<spla::SparseBool> splitCells;
    std::vector<indsz_t> dimPatt;

    std::vector<std::vector<indsz_t>> splitFiltInds;
    std::vector<double> uniqFiltVals;

    std::optional<ComplexInclusion> parent = std::nullopt;
};

int toComplex(const spla::SparseBool& adj, const std::vector<double>& filts,
                                           const std::vector<indsz_t>& dimVec,  
                                           eirene::Complex& comp, const bool check);

std::vector<MorseReducedResult> cellularPersistence(const eirene::Complex& comp);
PyreneResult pyCellularPersistence(indsz_t* rowVals, indsz_t numVals, indsz_t* colPtrs, indsz_t numCols, indsz_t numRows, 
                                            double* filtVals, size_t filtLen, indsz_t* dimVec, size_t numDim);

std::vector<MorseReducedResult> vripsPersistence(double* dMat, indsz_t numPoints, double maxRad = INFINITY);
PyreneResult pyVripsPersistence(double* dMat, indsz_t numPoints, double maxRad = INFINITY);

// Returns the star in complex of the cells represented in global indexes in globalCellInds, also represented as vec of global inds.
// Assumes globalCellInds is sorted.
std::vector<uint64_t> openStar(const std::vector<uint64_t>& globalCellInds, const eirene::Complex& complex); 

// Returns the complement of the open subset of complex, given as a argument open, relative to complex, 
// which is a complex in it's own right (i.e. Alexandrov closed)
std::vector<uint64_t> complement(const std::vector<uint64_t>& open, const eirene::Complex& complex);

// Returns the relative cellular persistence of the complex subcomplex represented by subsetInds in the complex.
// From paper "Local homology of abstract simplical complexes":
// For an open subset U ⊆ X of an abstract simplicial complex, the local homology at U is Hk(X, X\U).
std::vector<eirene::MorseReducedResult> relativeCellularPersistence(const eirene::Complex& comp, const std::vector<uint64_t>& openInds);
PyreneResult pyRelativeCellularPersistence(indsz_t* rowVals, indsz_t numVals, indsz_t* colPtrs, indsz_t numCols, indsz_t numRows, 
                                            double* filtVals, size_t filtLen, indsz_t* dimVec, size_t numDim, indsz_t* openComp, indsz_t compLen);

}   // namespace eirene

/*
 *  C types -> C++ types
 */

inline spla::SparseBool sparseBoolConvertInv(const CSparseBool* sBool) 
{
    auto result = spla::SparseBool();
    result.rowVals = std::vector<indsz_t>(sBool->rowVals.data, sBool->rowVals.data + sBool->rowVals.len);
    result.colPtrs = std::vector<indsz_t>(sBool->colPtrs.data, sBool->colPtrs.data + sBool->colPtrs.len);
    result.numRows = sBool->numRows;

    return result;
}

/*
 *  C++ Types -> C Types

inline U32Arr u32VecToArr(const std::vector<indsz_t>& vec)
{
    const size_t len  = vec.size();
    indsz_t* arrData = (indsz_t*)malloc(len * sizeof(indsz_t));

    std::copy(vec.data(), vec.data(), arrData);

    return U32Arr{
        .len  = len,
        .data = arrData,
    };
}

inline CPairArr pairVecConvert(const eirene::PairVec& vec)
{
    const size_t len   = vec.size();
    MorsePair* arrData = (MorsePair*)malloc(len * sizeof(MorsePair));

    std::transform(vec.begin(), vec.end(), arrData,
                   [](const MorsePair& p) { return MorsePair{.lowInd = p.lowInd, .highInd = p.highInd}; });

    return CPairArr{
        .len  = len,
        .data = arrData,
    };
}

inline CycleRepArr cycleRepsConvert(const std::vector<std::vector<indsz_t>>& vec)
{
    const size_t len = vec.size();
    U32Arr* arrData  = (U32Arr*)malloc(len * sizeof(U32Arr));

    std::transform(vec.begin(), vec.end(), arrData, [](const std::vector<indsz_t>& p) { return u32VecToArr(p); });

    return CycleRepArr{
        .len  = len,
        .data = arrData,
    };
}

inline CSparseBool sparseBoolConvert(const spla::SparseBool& mat)
{
    return CSparseBool{
        .rowVals = u32VecToArr(mat.rowVals),
        .colPtrs = u32VecToArr(mat.colPtrs),
        .numRows = mat.numRows,
    };
}

*/