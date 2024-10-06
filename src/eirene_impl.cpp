#include <algorithm>
#include <cassert>
#include <cmath>
#include <array>
#include <numeric>
#include <iostream>
#include <fstream>
#include <limits>

#include "eirene_impl.hpp"

using spla::SparseBool;
using std::tuple;
using std::vector;

namespace
{
template <typename T>
inline void print_vector(const std::vector<T> vec)
{
    if (vec.size() == 0)
    {
        std::cerr << "[]";
        return;
    }
    std::cerr << "[";

    for (indsz_t i = 0; i < vec.size() - 1; ++i) std::cerr << vec[i] << ", ";

    std::cerr << vec.back() << "]";
}

inline vector<indsz_t> getColFiltTransitionPtr(const vector<indsz_t>& colFilt)
{
    /*
     *  ColFilt is a sorted list of integers, this returns a vector whose ith element is the
     *  beginning of the ith string of constant values in colFilt.
     */

    vector<indsz_t> toRet(0);

    if (colFilt.size() == 0) return toRet;

    toRet.reserve(colFilt.size() / 2);   // random guess
    toRet.push_back(0);
    indsz_t currVal = colFilt[0];

    for (indsz_t i = 1; i < colFilt.size(); ++i)
    {
        if (colFilt[i] != currVal)
        {
            toRet.push_back(i);
            currVal = colFilt[i];
        }
    }

    toRet.push_back(colFilt.size());
    return toRet;
}

vector<indsz_t> integersinsameorderbycolumn2(const vector<indsz_t>& coBdryColPtrs, const vector<indsz_t>& columnFiltPtr)
{
    const auto numCols = columnFiltPtr.size() - 1;

    auto colWiseSum = std::vector<indsz_t>(coBdryColPtrs.size() - 1);
    for (indsz_t i = 0; i < coBdryColPtrs.size() - 1; ++i) {
        colWiseSum[i] = coBdryColPtrs[i + 1] - coBdryColPtrs[i];
    }

    const indsz_t maxColWiseSum = *std::max_element(colWiseSum.begin(), colWiseSum.end());
    auto x = std::vector<indsz_t>(maxColWiseSum + 1);
    std::fill(x.begin(), x.end(), 0);
    auto z = std::vector<indsz_t>(colWiseSum.size());

    for (indsz_t j = 0; j < numCols; ++j) {
        if (columnFiltPtr[j] == columnFiltPtr[j + 1]) continue;

        for (indsz_t i = columnFiltPtr[j]; i < columnFiltPtr[j + 1]; ++i) {
            x[colWiseSum[i]] += 1;
        }

        auto maxv = colWiseSum[columnFiltPtr[j]];
        auto minv = maxv;

        for (indsz_t i = columnFiltPtr[j] + 1; i < columnFiltPtr[j + 1]; ++i) {
            if (colWiseSum[i] > maxv) maxv = colWiseSum[i];
            else if (colWiseSum[i] < minv) minv = colWiseSum[i];
        }

        auto prevsum = columnFiltPtr[j];
        for (indsz_t i = minv; i <= maxv; ++i) {
            auto sum = prevsum + x[i];
            x[i] = prevsum;
            prevsum = sum;
        }

        for (indsz_t i = columnFiltPtr[j]; i < columnFiltPtr[j + 1]; ++i) {
            auto u = colWiseSum[i];
            z[i] = x[u];
            x[u] += 1;
        }

        for (indsz_t i = minv; i <= maxv; ++i) {
            x[i] = 0;
        }
    }
    return z;
}

vector<indsz_t> filtSubordColSortByMass(const vector<indsz_t>& coBdryColPtrs, const vector<indsz_t>& columnFiltPtr)
{
    /*
        Given that ColumnFiltPtr is a sorted list of indices
        Returns a permutation vector such that:
           1) For every i < numLowCells, the range columnFiltPtr[i]:columnFiltPtr[i + 1] maps to itself under the
           permutation
           2) For every i < numLowCells, the indices in columnFiltPtr[i]:columnFiltPtr[i + 1] correspond to columns
           in coBdryPtrs in an order such that the number of elements in the indexed columns is nondecreasing.

        See integersinsameorderbycolumn2 in Eirene.jl
    */

    // doing std::sort() on each column's indices with the sorting function being coBdryColPtrs[i + 1]
    // -coBdryColPtrs[i]. Also just sorting the whole thing ignoring criticalLowCount
    // Gotta make it fast sometime (radixsort)

    const indsz_t numLowCells = coBdryColPtrs.size() - 1;
    vector<indsz_t> toRet(numLowCells);
    std::iota(toRet.begin(), toRet.end(), 0);

    const auto sortFunc = [&coBdryColPtrs](indsz_t colIndA, indsz_t colIndB)
    {
        indsz_t aNumElem = coBdryColPtrs[colIndA + 1] - coBdryColPtrs[colIndA];
        indsz_t bNumElem = coBdryColPtrs[colIndB + 1] - coBdryColPtrs[colIndB];

        return aNumElem < bNumElem;
    };

    for (indsz_t colFiltValInd = 0; colFiltValInd < columnFiltPtr.size(); ++colFiltValInd)
    {
        indsz_t endDist = std::min(columnFiltPtr[colFiltValInd + 1], numLowCells);
        auto start = toRet.begin() + columnFiltPtr[colFiltValInd], end = toRet.begin() + endDist;
        std::sort(start, end, sortFunc);

        if (endDist == numLowCells) break;
    }

    return toRet;
}

void getPairs(spla::SparseBool& coBdry, const vector<indsz_t>& rowFilt, const vector<indsz_t>& colFilt,
              eirene::PairVec& pairVec)
{
    /*
     * Finds pairs for the bi-filtration
     * These define an acyclic f-relation on the complex, as defined in MFACPH.
     */

    const indsz_t criticalLowCount = coBdry.numCols(), criticalHighCount = coBdry.numRows;
    vector<indsz_t>&coBdryRowVals = coBdry.rowVals, &coBdryColPtrs = coBdry.colPtrs;
    vector<indsz_t> earliestBornCoFace(criticalLowCount, 0);
    vector<indsz_t> numFaces(criticalHighCount, 0);
    pairVec.clear();

    /*
     *  For each low dim. cell (coBdry col), if it has cofaces, find the
     *  one with the largest filtration index (corresp. to earliest birth time). Also as we go through the
     *  cofaces (rows), aggregate the sum of all values in each row (number of cells with that row's cell as a coface)
     *  Remember: columns are ordered s.t. cells born latest come first.
     *  Thus, we try and match the low cells born latest with the high cells born earliest.
     */

    for (indsz_t lowCellInd = 0; lowCellInd < criticalLowCount; ++lowCellInd)
    {
        indsz_t highestCoFaceInd = coBdryColPtrs[lowCellInd];
        // if (highestCoFaceInd >= coBdryRowVals.size()) continue;
        if (highestCoFaceInd >= coBdryColPtrs[lowCellInd + 1]) continue;
        const auto firstRow = coBdryRowVals[highestCoFaceInd];
        numFaces[firstRow]++;

        indsz_t highestFiltVal = rowFilt[firstRow];

        for (indsz_t coFaceInd = coBdryColPtrs[lowCellInd] + 1; coFaceInd < coBdryColPtrs[lowCellInd + 1]; ++coFaceInd)
        {
            indsz_t highCellInd  = coBdryRowVals[coFaceInd];
            indsz_t otherFiltVal = rowFilt[highCellInd];

            ++numFaces[highCellInd];

            if (otherFiltVal > highestFiltVal)
            {
                highestCoFaceInd = coFaceInd;
                highestFiltVal   = otherFiltVal;
            }
        }

        earliestBornCoFace[lowCellInd] = highestCoFaceInd;
    }

    vector<indsz_t> columnFiltPtr        = getColFiltTransitionPtr(colFilt);
    // vector<indsz_t> colWiseSumLinearized = filtSubordColSortByMass(coBdryColPtrs, columnFiltPtr);
    vector<indsz_t> colWiseSumLinearized = integersinsameorderbycolumn2(coBdryColPtrs, columnFiltPtr);
    vector<indsz_t> massSortedColsInds(criticalLowCount);
    std::iota(massSortedColsInds.begin(), massSortedColsInds.end(), 0);
    // perm::reindexInplace(colWiseSumLinearized, massSortedColsInds);
    for (indsz_t i = 0; i < colWiseSumLinearized.size(); ++i) {
        massSortedColsInds[colWiseSumLinearized[i]] = i;
    }
    // massSortedColsInds = perm::reindex(colWiseSumLinearized, massSortedColsInds);

    // As we make pairs (low, hi), once a col is associated to a pair, all rows it has 1's in (cofaces)
    // not eligible to be paired this time around
    // In other words, when pairing a new row, it cannot have a 1 in an already paired column,
    // In other words, if we reindex the paired submatrix rows/cols by pairs' hiInd/lowInds after this,
    // it's upper triangular with 1s on the diagonal
    vector<bool> nCoveredSupp(criticalHighCount, true);
    indsz_t nCoveredCount = 0;

    // Perm is such that first cols after perm have small coface degree,
    // but all remain in reverse filt. order
    for (indsz_t colPermIndex = 0; colPermIndex < criticalLowCount; ++colPermIndex)
    {
        const indsz_t colIndex = massSortedColsInds[colPermIndex];
        // Thus, this is nondecreasing within each subcomplex/filtration value
        const indsz_t numCoFaces = coBdryColPtrs[colIndex + 1] - coBdryColPtrs[colIndex];

        if (numCoFaces == 0) continue;

        const indsz_t largeFiltFaceInd   = earliestBornCoFace[colIndex];      // index in sparse bdry
        indsz_t largeCoFace              = coBdryRowVals[largeFiltFaceInd];   // row index of coface
        indsz_t largeCoFaceDegree        = numFaces[largeCoFace];             // num low cells with coface as coface
        const indsz_t largeCoFaceFiltVal = rowFilt[largeCoFace];



        // If largeCoFace is incident to a face that has already been paired (its 'covered')
        if (!nCoveredSupp[largeCoFace])
        {
            // ... try and find one that hasn't.
            for (indsz_t otherCoFaceInd = largeFiltFaceInd + 1; otherCoFaceInd < coBdryColPtrs[colIndex + 1];
                 ++otherCoFaceInd)
            {
                indsz_t otherCoFace = coBdryRowVals[otherCoFaceInd], otherFiltVal = rowFilt[otherCoFace],
                        otherDegree = numFaces[otherCoFace];

                // Trying to pick a coface with the earliest possible birth time that has the smallest number of
                // incident faces and that has not yet had an incident face paired
                if (otherFiltVal == largeCoFaceFiltVal && otherDegree <= largeCoFaceDegree && !nCoveredSupp[largeCoFace] && nCoveredSupp[otherCoFace])
                {
                    largeCoFaceDegree = otherDegree;
                    largeCoFace       = otherCoFace;
                }
            }
        }

        // If we found one
        if (nCoveredSupp[largeCoFace])  {
            pairVec.push_back(MorsePair{colIndex, largeCoFace});
        }

        // Set all of the lowcell collIndex's cofaces as covered
        for (indsz_t colElem = coBdryColPtrs[colIndex]; colElem < coBdryColPtrs[colIndex + 1]; ++colElem)
        {
            nCoveredCount += nCoveredSupp[coBdryRowVals[colElem]];
            nCoveredSupp[coBdryRowVals[colElem]] = false;
        }

        // If we've already touched on every coface, not gonna add any more pairs, so we can just break
        if (nCoveredCount == criticalHighCount) break;
    }
}

}   // anonymous namespace

// sparse linear algebra
namespace spla
{
struct F2Result
{
    /*
     * Vectors of this struct are used in sparse F2 matrix operations.
     * Allows merging three length "numRows" arrays into one.
     * Makes using it a little more confusing but simplifies cache locality and
     * reduces # of vectors created.
     */
    indsz_t colInd;
    indsz_t rowInd;
    bool indicator;
};

inline void F2MultAlternative(const indsz_t rowInd, const indsz_t colInd, indsz_t& prodRowCounter,
                              vector<F2Result>& colResult)
{
    /*
     *  Forms the computation of the inner loop of many sparse F2 algorithms.
     *  "Alternative" because we either add the row as being potentially nonzero (set indicator to true and
     *  add one to prodRowCounter), or toggle the indicated (bc. F2).
     */
    if (colResult[rowInd].colInd != (colInd + 1))
    {
        colResult[rowInd].colInd         = colInd + 1;
        colResult[rowInd].indicator      = true;
        colResult[prodRowCounter].rowInd = rowInd;
        prodRowCounter++;
    }
    else
        colResult[rowInd].indicator = !colResult[rowInd].indicator;
}

inline void F2SilentPremult(const indsz_t startInd, const indsz_t endInd, const indsz_t colInd,
                            const vector<indsz_t>& rowVals, vector<F2Result>& colResult)
{
    /*
     *  A "Silent" sparse matrix has implicit 1s on the diagonal.
     *  When doing operations on these matrices, this gets called in the beginning of
     *  each loop computing a column of the output matrix.
     */

    for (indsz_t valInd = startInd; valInd < endInd; valInd++)
    {
        indsz_t rowInd = rowVals[valInd];

        colResult[valInd - startInd].rowInd = rowInd;
        colResult[rowInd].colInd            = colInd + 1;
        colResult[rowInd].indicator         = true;
    }
}

inline void F2DoubleSilentPremult(const indsz_t startInd, const indsz_t endInd, const indsz_t colInd,
                                  const vector<indsz_t>& rowVals, vector<F2Result>& colResult)
{
    /*
     *  A "Silent" sparse matrix has implicit 1s on the diagonal.
     *  When doing operations on these matrices, this gets called in the beginning of
     *  each loop computing a column of the output matrix.
     */

    for (indsz_t valInd = startInd; valInd < endInd; valInd++)
    {
        indsz_t rowInd = rowVals[valInd];

        colResult[valInd - startInd].rowInd = rowInd;
        colResult[rowInd].colInd            = colInd + 1;
        colResult[rowInd].indicator         = true;
    }

    colResult[endInd - startInd].rowInd = colInd;
    colResult[colInd].colInd            = colInd + 1;
    colResult[colInd].indicator         = true;
}

auto F2UpTriInvSilentOut(const spla::SparseBool& A)
{
    /*
     * Given a matrix A which is upper triangular and invertible, this returns a "silent Identity"
     * inverse of a matrix. Remember that an upper triangular matrix is invertible iff
     * every diagonal entry is nonzero and that the inverse of an upper triangular matrix is
     * always upper triangular.
     *
     * Silent here is in the sense that every invertible UT matrix can be written as B' = B + I.
     * Silent refers to only storing the B part.
     *
     * Does not assume row values inside individual column segments of aColPtr are sorted!
     */

    const indsz_t aNumCols = A.numCols(), aNumRows = A.numRows;

    spla::SparseBool aInv{aNumCols, aNumRows};

    aInv.rowVals.reserve(A.rowVals.size());
    vector<F2Result> colResult(aNumCols);

    for (indsz_t aColInd = 0; aColInd < aNumCols; ++aColInd)
    {
        // If the column has only one nonzero element, it must be the (aColInd, aColInd) element since A is
        // assumed to be invertible. Then, the inverse's column must also have only the diagonal entry as nonzero,
        // this can be proven by looking at the expansion of (A^(-1) * A)_ij combined with the assumptions.
        if (A.colPtrs[aColInd] + 1 == A.colPtrs[aColInd + 1]) aInv.colPtrs[aColInd + 1] = aInv.colPtrs[aColInd];

        // Elif the column has two vals, so that the high dimension cell has two cofaces...
        else if (A.colPtrs[aColInd] + 2 == A.colPtrs[aColInd + 1])
        {
            indsz_t k = 0;

            // If the row index of the element is < column index (i.e. it's in the upper tri part),
            // set k to that row index. Otherwise, set k to the row index of the next one.
            // Since A is invertible, one will be on the diagonal. So picking the off diagonal one.
            if (A.rowVals[A.colPtrs[aColInd]] < aColInd)
                k = A.rowVals[A.colPtrs[aColInd]];
            else
                k = A.rowVals[A.colPtrs[aColInd] + 1];

            // Append row indices of elements in column k
            aInv.rowVals.insert(aInv.rowVals.end(), aInv.rowVals.begin() + aInv.colPtrs[k],
                                aInv.rowVals.begin() + aInv.colPtrs[k + 1]);

            aInv.rowVals.push_back(k);
            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
        }

        // General case: basically, do the two element case for more elements.
        else
        {
            F2SilentPremult(A.colPtrs[aColInd], A.colPtrs[aColInd + 1], aColInd, A.rowVals, colResult);

            colResult[aColInd].indicator = false;

            indsz_t newRowCounter = A.colPtrs[aColInd + 1] - A.colPtrs[aColInd];

            for (indsz_t colElemInd = A.colPtrs[aColInd]; colElemInd < A.colPtrs[aColInd + 1]; ++colElemInd)
            {
                indsz_t elemRowVal = A.rowVals[colElemInd];

                // For each off-diagonal value. Unsorted rows again.
                if (elemRowVal < aColInd)
                {
                    for (indsz_t aInvColElInd = aInv.colPtrs[elemRowVal]; aInvColElInd < aInv.colPtrs[elemRowVal + 1];
                         ++aInvColElInd)
                    {
                        indsz_t k = aInv.rowVals[aInvColElInd];

                        F2MultAlternative(k, aColInd, newRowCounter, colResult);
                    }
                }
            }

            aInv.rowVals.reserve(aInv.rowVals.size() + newRowCounter);

            for (indsz_t newRowInd = 0; newRowInd < newRowCounter; ++newRowInd)
            {
                indsz_t rowInd = colResult[newRowInd].rowInd;

                if (colResult[rowInd].indicator) aInv.rowVals.push_back(rowInd);
            }

            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
            // std::sort(aInv.rowVals.begin() + aInv.colPtrs[aColInd], aInv.rowVals.begin() + aInv.colPtrs[aColInd + 1]);
        }
    }

    return aInv;
}

auto F2UpTriInvSilentInOut(const spla::SparseBool& A)
{
    /*
     * Given an upper triangular invertible matrix represented as a silent matrix,
     * meaning if we write A' = A + I, this takes A as input and returns B, where
     * B' * A' + I and B' = B + I.
     *
     * Suggest looking at F2UpTriInvSilentOut first.
     *
     * This one does assume that A's colptr column segments are sorted.
     */

    const indsz_t aNumCols = A.numCols();

    spla::SparseBool aInv{aNumCols, aNumCols};

    aInv.rowVals.reserve(A.rowVals.size());
    vector<F2Result> colResult(aNumCols);

    for (indsz_t aColInd = 0; aColInd < aNumCols; ++aColInd)
    {
        // If the column has no nonzero elements, the inverse has none (silent in/out).
        if (A.colPtrs[aColInd] == A.colPtrs[aColInd + 1])
            aInv.colPtrs[aColInd + 1] = aInv.colPtrs[aColInd];

        else if (A.colPtrs[aColInd] + 1 == A.colPtrs[aColInd + 1])
        {
            indsz_t k = A.rowVals[A.colPtrs[aColInd]];

            aInv.rowVals.insert(aInv.rowVals.end(), aInv.rowVals.begin() + aInv.colPtrs[k],
                                aInv.rowVals.begin() + aInv.colPtrs[k + 1]);

            aInv.rowVals.push_back(k);
            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
        }

        else
        {
            F2SilentPremult(A.colPtrs[aColInd], A.colPtrs[aColInd + 1], aColInd, A.rowVals, colResult);

            indsz_t newRowCounter = A.colPtrs[aColInd + 1] - A.colPtrs[aColInd];

            for (indsz_t colElemInd = A.colPtrs[aColInd]; colElemInd < A.colPtrs[aColInd + 1]; ++colElemInd)
            {
                indsz_t elemRowVal = A.rowVals[colElemInd];

                for (indsz_t aInvColElInd = aInv.colPtrs[elemRowVal]; aInvColElInd < aInv.colPtrs[elemRowVal + 1];
                     ++aInvColElInd)
                {
                    indsz_t k = aInv.rowVals[aInvColElInd];

                    F2MultAlternative(k, aColInd, newRowCounter, colResult);
                }
            }

            aInv.rowVals.reserve(aInv.rowVals.size() + newRowCounter);

            for (indsz_t newRowInd = 0; newRowInd < newRowCounter; ++newRowInd)
            {
                indsz_t rowInd = colResult[newRowInd].rowInd;

                if (colResult[rowInd].indicator) aInv.rowVals.push_back(rowInd);
            }

            aInv.colPtrs[aColInd + 1] = aInv.rowVals.size();
        }
    }

    return aInv;
}

auto matMultF2(const spla::SparseBool& A, const spla::SparseBool& B)
{
    /*
     * Straight up F2 matmult. No silent ins or outs.
     */

    const size_t aNumRows = A.numRows, bNumCols = B.numCols();

    spla::SparseBool prodMat{aNumRows, bNumCols};

    size_t sizeGuess = std::min(aNumRows * bNumCols, static_cast<size_t>(A.numNonZero() + B.numNonZero()));
    prodMat.rowVals.reserve(sizeGuess);

    vector<F2Result> colResult(aNumRows);

    for (indsz_t prodColInd = 0; prodColInd < bNumCols; ++prodColInd)
    {
        indsz_t newRowCounter = 0;

        for (indsz_t bValInd = B.colPtrs[prodColInd]; bValInd < B.colPtrs[prodColInd + 1]; ++bValInd)
        {
            indsz_t bRowInd = B.rowVals[bValInd];

            for (indsz_t aValInd = A.colPtrs[bRowInd]; aValInd < A.colPtrs[bRowInd + 1]; ++aValInd)
            {
                indsz_t aRowInd = A.rowVals[aValInd];

                F2MultAlternative(aRowInd, prodColInd, newRowCounter, colResult);
            }
        }

        prodMat.colPtrs[prodColInd + 1] = prodMat.colPtrs[prodColInd];

        if (newRowCounter == 0) continue;

        for (indsz_t j = 0; j < newRowCounter; ++j)
        {
            indsz_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) prodMat.rowVals.push_back(rowVal);
        }

        prodMat.colPtrs[prodColInd + 1] = prodMat.rowVals.size();
    }

    return prodMat;
}

auto siliLeftMatMultF2(const spla::SparseBool& A, const spla::SparseBool& B)
{
    /*
     *  A * B = C over F2 for sparse matrices, with A "sili".
     *
     *  Silent means: implicit 1s on the diagonal.
     *
     *  Output is not sili (1's on diagonal can't be guaranteed in general!).
     *
     *  c_ij = \Sigma_k a_ik * b_kj
     */

    const indsz_t aNumRows = A.numRows, bNumCols = B.numCols();

    spla::SparseBool prodMat{aNumRows, bNumCols};
    prodMat.rowVals.reserve(A.colPtrs.back() + B.colPtrs.back());

    vector<F2Result> colResult(aNumRows);

    for (indsz_t bColInd = 0; bColInd < bNumCols; ++bColInd)
    {
        indsz_t prodNonzeroCounter = B.colPtrs[bColInd + 1] - B.colPtrs[bColInd];

        F2SilentPremult(B.colPtrs[bColInd], B.colPtrs[bColInd + 1], bColInd, B.rowVals, colResult);

        for (indsz_t bValInd = B.colPtrs[bColInd]; bValInd < B.colPtrs[bColInd + 1]; ++bValInd)
        {
            indsz_t bRowInd = B.rowVals[bValInd];

            for (indsz_t aValInd = A.colPtrs[bRowInd]; aValInd < A.colPtrs[bRowInd + 1]; aValInd++)
            {
                // A has nonzero entry at (aRowInd, bRowInd)
                indsz_t aRowInd = A.rowVals[aValInd];
                F2MultAlternative(aRowInd, bColInd, prodNonzeroCounter, colResult);
            }
        }

        prodMat.colPtrs[bColInd + 1] = prodMat.colPtrs[bColInd];

        if (prodNonzeroCounter == 0) continue;

        for (indsz_t j = 0; j < prodNonzeroCounter; ++j)
        {
            indsz_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) prodMat.rowVals.push_back(rowVal);
        }

        prodMat.colPtrs[bColInd + 1] = prodMat.rowVals.size();
        const indsz_t colStart = prodMat.colPtrs[bColInd], colEnd = prodMat.colPtrs[bColInd + 1];
        std::sort(prodMat.rowVals.begin() + colStart, prodMat.rowVals.begin() + colEnd);
    }

    return prodMat;
}

auto fullSiliMatMultF2(const spla::SparseBool& A, const spla::SparseBool& B)
{
    const indsz_t aNumRows = A.numRows, bNumCols = B.numCols();

    spla::SparseBool prodMat{aNumRows, bNumCols};
    prodMat.rowVals.reserve(A.colPtrs.back() + B.colPtrs.back());

    vector<F2Result> colResult(aNumRows);

    for (indsz_t bColInd = 0; bColInd < bNumCols; ++bColInd)
    {
        indsz_t prodNonzeroCounter = 1 + B.colPtrs[bColInd + 1] - B.colPtrs[bColInd];

        F2DoubleSilentPremult(B.colPtrs[bColInd], B.colPtrs[bColInd + 1], bColInd, B.rowVals, colResult);

        for (indsz_t bValInd = B.colPtrs[bColInd]; bValInd < B.colPtrs[bColInd + 1]; ++bValInd)
        {
            indsz_t bRowInd = B.rowVals[bValInd];

            for (indsz_t aValInd = A.colPtrs[bRowInd]; aValInd < A.colPtrs[bRowInd + 1]; aValInd++)
            {
                // A has nonzero entry at (aRowInd, bRowInd)
                indsz_t aRowInd = A.rowVals[aValInd];
                F2MultAlternative(aRowInd, bColInd, prodNonzeroCounter, colResult);
            }
        }

        prodMat.colPtrs[bColInd + 1] = prodMat.colPtrs[bColInd];

        if (prodNonzeroCounter == 0) continue;

        for (indsz_t j = 0; j < prodNonzeroCounter; ++j)
        {
            indsz_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) prodMat.rowVals.push_back(rowVal);
        }

        prodMat.colPtrs[bColInd + 1] = prodMat.rowVals.size();
    }

    return prodMat;
}

auto blockProdSum(const spla::SparseBool& D, const spla::SparseBool& C, const spla::SparseBool& B,
                  spla::SparseBool& result)
{
    /*
     *  Calculates S = D + C*B, stores the result in result, overwriting it.
     *  (Only) B is Silent!!!
     */

    bool validSizes = (D.numRows == C.numRows) && (C.numCols() == B.numRows) && (B.numCols() == D.numCols());
    if (!validSizes) { fmt::print(stderr, "Invalid sizes in blockProdSum"); }

    const indsz_t dNumCols = D.numCols(), dNumRows = D.numRows;

    result.colPtrs.resize(dNumCols + 1);
    result.rowVals.clear();
    std::fill(result.colPtrs.begin(), result.colPtrs.end(), 0);
    result.numRows = D.numRows;

    vector<F2Result> colResult(dNumRows);

    for (indsz_t resultColInd = 0; resultColInd < dNumCols; ++resultColInd)
    {
        // Compensate for B being sili and D's contribution to the columns.
        F2SilentPremult(D.colPtrs[resultColInd], D.colPtrs[resultColInd + 1], resultColInd, D.rowVals, colResult);

        indsz_t newRowsCounter = D.colPtrs[resultColInd + 1] - D.colPtrs[resultColInd];

        // B might have more columns than D; if we're not there yet, calculate the contribution of C*B to prod.
        if (B.colPtrs[resultColInd] < B.colPtrs[dNumCols])
        {
            for (indsz_t bElemInd = B.colPtrs[resultColInd]; bElemInd < B.colPtrs[resultColInd + 1]; ++bElemInd)
            {
                indsz_t bRowInd = B.rowVals[bElemInd];

                for (indsz_t cElemInd = C.colPtrs[bRowInd]; cElemInd < C.colPtrs[bRowInd + 1]; ++cElemInd)
                {
                    indsz_t cRowVal = C.rowVals[cElemInd];

                    F2MultAlternative(cRowVal, resultColInd, newRowsCounter, colResult);
                }
            }
        }

        result.colPtrs[resultColInd + 1] = result.colPtrs[resultColInd];

        if (newRowsCounter == 0) continue;

        for (indsz_t j = 0; j < newRowsCounter; ++j)
        {
            indsz_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) result.rowVals.push_back(rowVal);
        }

        result.colPtrs[resultColInd + 1] = result.rowVals.size();
    }
}

auto blockProdSumSilentIColsLeft(const spla::SparseBool& D, const spla::SparseBool& C, const spla::SparseBool& B,
                                 spla::SparseBool& result, const vector<indsz_t>& col2SilentI)
{
    /*
     *  Calculates S = D + C*B, stores the result in result, overwriting it.
     *  Only C is silent.
     *  col2SilentI determines where the silent I entries are in each column of C.
     */

    const indsz_t dNumRows = D.numRows, dNumCols = D.numCols();

    result.colPtrs.resize(D.numCols() + 1);
    std::fill(result.colPtrs.begin(), result.colPtrs.end(), 0);
    result.rowVals.clear();
    result.numRows = D.numRows;

    vector<F2Result> colResult(dNumRows + 1);

    for (indsz_t resultColInd = 0; resultColInd < dNumCols; ++resultColInd)
    {
        F2SilentPremult(D.colPtrs[resultColInd], D.colPtrs[resultColInd + 1], resultColInd, D.rowVals, colResult);

        indsz_t newRowsCounter = D.colPtrs[resultColInd + 1] - D.colPtrs[resultColInd];

        // This loop handles the silent contribution with the col2SilentI permutation
        for (indsz_t bElemInd = B.colPtrs[resultColInd]; bElemInd < B.colPtrs[resultColInd + 1]; ++bElemInd)
        {
            indsz_t bRowInd = B.rowVals[bElemInd];
            // B's row ind determines which columns of C is being added
            indsz_t iRowVal = col2SilentI[bRowInd];

            F2MultAlternative(iRowVal, resultColInd, newRowsCounter, colResult);
        }

        for (indsz_t bElemInd = B.colPtrs[resultColInd]; bElemInd < B.colPtrs[resultColInd + 1]; ++bElemInd)
        {
            indsz_t bRowInd = B.rowVals[bElemInd];

            for (indsz_t cElemInd = C.colPtrs[bRowInd]; cElemInd < C.colPtrs[bRowInd + 1]; ++cElemInd)
            {
                indsz_t cRowVal = C.rowVals[cElemInd];

                F2MultAlternative(cRowVal, resultColInd, newRowsCounter, colResult);
            }
        }

        result.colPtrs[resultColInd + 1] = result.colPtrs[resultColInd];

        if (newRowsCounter == 0) continue;

        for (indsz_t j = 0; j < newRowsCounter; ++j)
        {
            indsz_t rowVal = colResult[j].rowInd;

            if (colResult[rowVal].indicator) result.rowVals.push_back(rowVal);
        }

        result.colPtrs[resultColInd + 1] = result.rowVals.size();
    }
}

auto stackedSubMat(const spla::SparseBool& blockMat, const vector<indsz_t>& rows1, const vector<indsz_t>& rows2,
                   const vector<indsz_t>& cols)
{
    /*
        Returns the sparse submatrices whose elements are in columns cols, and rows rows 1/2.

        The row values in the submatrices are sorted and the actual row indices are relative to the
        order in rows1.

        Rows1 .> Rows2?
    */

    const indsz_t numCols = cols.size(), numRows = blockMat.numRows;
    const vector<indsz_t>&colPtrs = blockMat.colPtrs, &rowVals = blockMat.rowVals;

    auto suppCol1 = vector<bool>(numRows, false);
    auto suppCol2 = vector<bool>(numRows, false);
    for (auto r : rows1) suppCol1[r] = true;
    for (auto r : rows2) suppCol2[r] = true;
    indsz_t nz1 = 0;
    indsz_t nz2 = 0;

    for (indsz_t selCol : cols) {
        for (indsz_t i = colPtrs[selCol]; i < colPtrs[selCol + 1]; ++i) {
            const indsz_t rowVal = rowVals[i];
            if (suppCol1[rowVal]) nz1++;
            else if (suppCol2[rowVal]) nz2++;
        }
    }


    // Create the submatrix rows/cols:
    spla::SparseBool topBlock{rows1.size(), numCols}, bottomBlock{rows2.size(), numCols};
    topBlock.rowVals.reserve(nz1);
    bottomBlock.rowVals.reserve(nz2);

    for (indsz_t col = 0; col < cols.size(); ++col)   // colInd : cols)
    {
        indsz_t colInd = cols[col];

        for (indsz_t elemInd = colPtrs[colInd]; elemInd < colPtrs[colInd + 1]; ++elemInd)
        {
            indsz_t rowVal    = rowVals[elemInd];

            if (suppCol1[rowVal])
                topBlock.rowVals.push_back(rowVal);
            else if (suppCol2[rowVal])
                bottomBlock.rowVals.push_back(rowVal);
        }

        topBlock.colPtrs[col + 1]    = topBlock.rowVals.size();
        bottomBlock.colPtrs[col + 1] = bottomBlock.rowVals.size();

        /*
        std::sort(topBlock.rowVals.begin() + topBlock.colPtrs[col],
                  topBlock.rowVals.begin() + topBlock.colPtrs[col + 1]);
        std::sort(bottomBlock.rowVals.begin() + bottomBlock.colPtrs[col],
                  bottomBlock.rowVals.begin() + bottomBlock.colPtrs[col + 1]);
        */
    }

    return std::make_tuple(topBlock, bottomBlock);
}

inline void appendSparseToBack(spla::SparseBool& A, const indsz_t lastACol, spla::SparseBool& B, const indsz_t lastBCol)
{
    /*
     *  Appends the first 0:lastBCol columns of B to the back of the first 0:lastACol columns of A.
     */
    const indsz_t lastAValInd = A.colPtrs[lastACol], lastBValInd = B.colPtrs[lastBCol];

    A.rowVals.resize(lastAValInd);
    A.colPtrs.resize(lastACol);
    B.rowVals.resize(lastBValInd);
    B.colPtrs.resize(lastBCol + 1);

    std::transform(B.colPtrs.begin(), B.colPtrs.end(), B.colPtrs.begin(),
                   [lastAValInd](indsz_t val) { return val + lastAValInd; });

    A.colPtrs.insert(A.colPtrs.end(), B.colPtrs.begin(), B.colPtrs.end());
    A.rowVals.insert(A.rowVals.end(), B.rowVals.begin(), B.rowVals.end());
}

vector<indsz_t> findDownStreamInUpperTriangMat(const spla::SparseBool& coBdry, const vector<indsz_t>& unpairedCols,
                                               const eirene::PairVec& pairVec)
{
    /*
     *  Returns ("keeps") pairs whose columns have a 1 in the same row as currently unpaired columns.
     *  The benefit of the algorithm seems to be that an unkept pairs saves significant computation one could skip this
     *  and keep all pairs, always, and the algorithm would still be correct.
     */

    const indsz_t numPairs = pairVec.size(), numRows = coBdry.numRows;
    const auto &coBdryColPtrs = coBdry.colPtrs, &coBdryRowVals = coBdry.rowVals;

    if (numPairs == 0) return vector<indsz_t>{};

    // Find out what rows are nonzero in the unpairedCols. Faster algo?
    vector<bool> criticalCellRowSupp(numRows);

    for (indsz_t colInd : unpairedCols)
        for (indsz_t elemInd = coBdryColPtrs[colInd]; elemInd < coBdryColPtrs[colInd + 1]; ++elemInd)
            criticalCellRowSupp[coBdryRowVals[elemInd]] = true;

    // unpairedColRowSupp contains indexes of rows with 1's in critical cols.
    vector<indsz_t> unpairedColRowSupp(0);
    unpairedColRowSupp.reserve(numRows);

    for (indsz_t i = 0; i < numRows; ++i)
        if (criticalCellRowSupp[i]) unpairedColRowSupp.push_back(i);

    // rowIndToPairInd[rowIndex] = pairIndex
    vector<indsz_t> rowIndToPairInd = vector<indsz_t>(numRows);
    for (indsz_t i = 0; i < numPairs; ++i) rowIndToPairInd[pairVec[i].highInd] = i;

    // given a row index, returns true iff a pair is supported by that row
    vector<bool> pairedRowInds(numRows, false);
    for (indsz_t i = 0; i < numPairs; ++i) pairedRowInds[pairVec[i].highInd] = true;

    // for each row with a 1 in an unpaired cell's col, if it's a paired row,
    // set dSS[pairInd] = true
    vector<bool> downStreamSupp(pairVec.size(), false);
    for (indsz_t rowInd : unpairedColRowSupp)
        if (pairedRowInds[rowInd]) downStreamSupp[rowIndToPairInd[rowInd]] = true;

    // The effect of all the following seems to be:
    // 'keep' any pairs whose cols contain a 1 in the same row as any paired col which has a 1 in the
    // same row as an unpaired col. In other words:
    // "keep any pairs whose low-dim cells have as coface a paired high-dim cell which has an unpaired face"
    for (int32_t pairInd = numPairs - 1; pairInd >= 0; --pairInd)
    {
        // For each rowval of the paired column...
        const indsz_t pairCol = pairVec[pairInd].lowInd;

        for (indsz_t elemInd = coBdryColPtrs[pairCol]; elemInd < coBdryColPtrs[pairCol + 1]; ++elemInd)
        {
            indsz_t rowVal = coBdryRowVals[elemInd];

            // if a pair is supported on the row AND that pair's row is in dSS <==>
            // there's a 1 in the row in the critical cell part of the matrix
            if (pairedRowInds[rowVal] && downStreamSupp[rowIndToPairInd[rowVal]])
            {
                // then go through all of the elements in the column again
                for (indsz_t elemInd = coBdryColPtrs[pairCol]; elemInd < coBdryColPtrs[pairCol + 1]; ++elemInd)
                {
                    indsz_t rowVal = coBdryRowVals[elemInd];

                    // And mark any paired cofaces as having a critical face cell
                    if (pairedRowInds[rowVal]) downStreamSupp[rowIndToPairInd[rowVal]] = true;
                }

                break;
            }
        }
    }

    // Finally, return the indexes of pairs set true in dSS
    const indsz_t downStreamNum = std::accumulate(downStreamSupp.begin(), downStreamSupp.end(), 0);

    vector<indsz_t> downStreamPairInds(downStreamNum);
    indsz_t counter = 0;

    // i iterates over pairs
    for (indsz_t i = 0; counter < downStreamNum; ++i)
        if (downStreamSupp[i]) downStreamPairInds[counter++] = i;

    return downStreamPairInds;
}

void insertCols(const spla::SparseBool& A, spla::SparseBool& B, const vector<indsz_t>& columnInds,
                const indsz_t destInd)
{
    /*
        Copies columns with indices in columnInds from sparse matrix A to B,
        starting at column destInd, and overwriting the existing columns in B.
    */
    indsz_t numNewElem = 0;
    for (indsz_t ind : columnInds) numNewElem += A.colPtrs[ind + 1] - A.colPtrs[ind];

    B.rowVals.resize(B.colPtrs[destInd] + numNewElem);

    if (B.colPtrs.size() < destInd + columnInds.size() + 1) B.colPtrs.resize(destInd + columnInds.size() + 1);

    B.colPtrs[0] = 0;

    for (indsz_t newCol = 0; newCol < columnInds.size(); ++newCol)
    {
        indsz_t newColInd = destInd + newCol, oldColInd = columnInds[newCol];
        B.colPtrs[newColInd + 1] = B.colPtrs[newColInd] + A.colPtrs[oldColInd + 1] - A.colPtrs[oldColInd];
        std::copy(A.rowVals.begin() + A.colPtrs[oldColInd], A.rowVals.begin() + A.colPtrs[oldColInd + 1],
                  B.rowVals.begin() + B.colPtrs[newColInd]);
    }
}

inline static auto copySubMatrix(const spla::SparseBool& mat, const vector<indsz_t>& colIndices)
{
    /*
        Returns a pair of vectors, representing the sparse submatrix [:] x [colIndices[:]).
    */

    indsz_t numValAcc = 0;
    for (indsz_t colInd = 0; colInd < colIndices.size(); ++colInd)
        numValAcc += mat.colPtrs[colIndices[colInd] + 1] - mat.colPtrs[colIndices[colInd]];

    const indsz_t numCol = colIndices.size(), numVal = numValAcc;

    spla::SparseBool subMat{mat.numRows, numCol};
    subMat.rowVals.resize(numVal);

    for (indsz_t colInd = 0; colInd < numCol; ++colInd)
    {
        indsz_t origInd        = colIndices[colInd];
        indsz_t elemBeforeNext = mat.colPtrs[origInd + 1];
        indsz_t elemBeforeCurr = mat.colPtrs[origInd];

        subMat.colPtrs[colInd + 1] = subMat.colPtrs[colInd] + elemBeforeNext - elemBeforeCurr;

        // pointer arithmetic
        std::copy(mat.rowVals.data() + elemBeforeCurr, mat.rowVals.data() + elemBeforeNext,
                  subMat.rowVals.data() + subMat.colPtrs[colInd]);
    }

    return subMat;
}

void schurIt(spla::SparseBool& coBdry, vector<indsz_t>& criticalLow, vector<indsz_t>& criticalHigh,
             eirene::PairVec& pairVec, eirene::PairVec& excisedPairs, const vector<indsz_t>& unpairedRows,
             const vector<indsz_t>& unpairedCols, spla::SparseBool& criticalPart, spla::SparseBool& pairedPart)
{
    /*
     * Does an iteration of Schur complement-based reduction.
     * Specifically, tries to do a block pivot clearing operation, removing as many pairs as it can
     * at once from the co-boundary/complex and updates all the boundary and index related structres.
     * Section 4 of "Matroid Filtrations and Computational Persistent Homology" [MFCPH] is vital for this part.
     */

    const indsz_t numPairs = pairVec.size();

    vector<indsz_t> indsToCopy(numPairs);
    std::transform(pairVec.rbegin(), pairVec.rend(), indsToCopy.begin(),
                   [](const MorsePair& pair) { return pair.lowInd; });

    // Copies the boundary info for new pairs from criticalPart to pairedPart.
    insertCols(criticalPart, pairedPart, indsToCopy, excisedPairs.size());

    // note the reindexing to low/hi cells inds from col/row inds by passing thru criticalLow/Hi
    std::transform(pairVec.rbegin(), pairVec.rend(), std::back_inserter(excisedPairs),
                   [&criticalHigh, &criticalLow](const MorsePair& pair) {
                       return MorsePair{criticalLow[pair.lowInd], criticalHigh[pair.highInd]};
                   });

    auto pairsKept = findDownStreamInUpperTriangMat(coBdry, unpairedCols, pairVec);

    const indsz_t numKept = pairsKept.size();
    auto pairedRows = vector<indsz_t>(numKept), pairedCols = vector<indsz_t>(numKept);

    for (indsz_t i = 0; i < numKept; ++i)
    {
        indsz_t keptIndex = pairsKept[i];
        pairedRows[i]     = pairVec[keptIndex].highInd;
        pairedCols[i]     = pairVec[keptIndex].lowInd;
    }

    /*
    pairedCols unpairedCols
    _____ _______

    [ X .. - Y ..]  |
    [ ..   - ..  ]  | pairedRows
    [----------- ]
    [ Z .. - D ..]  |
    [ ..   - ..  ]  | unpairedRows

    Note pairedRows/Cols contains only those that were kept by downstream,
    the rest disappear without further computation.

    X is upper triangular and invertible by construction.
    (In particular, see definition of f-relation in section 4).
    */

    auto [X, Z] = stackedSubMat(coBdry, pairedRows, unpairedRows, pairedCols);
    auto [Y, D] = stackedSubMat(coBdry, pairedRows, unpairedRows, unpairedCols);

    // criticalPart holds that part of reduced cobdry that represents the images of the critical cells...
    auto L = copySubMatrix(criticalPart, pairedCols);     // totalRows x pairedCols
    auto R = copySubMatrix(criticalPart, unpairedCols);   // totalRows x unpairedCols

    vector<indsz_t> translator(coBdry.numRows);
    for (indsz_t i = 0; i < pairedRows.size(); ++i) {
        translator[pairedRows[i]] = i;
    }

    for (indsz_t i = 0; i < X.rowVals.size(); ++i) {
        X.rowVals[i] = translator[X.rowVals[i]];
    }
    for (indsz_t i = 0; i < Y.rowVals.size(); ++i) {
        Y.rowVals[i] = translator[Y.rowVals[i]];
    }
    

    for (indsz_t i = 0; i < unpairedRows.size(); ++i) {
        translator[unpairedRows[i]] = i;
    }
    for (indsz_t i = 0; i < Z.rowVals.size(); ++i) {
        Z.rowVals[i] = translator[Z.rowVals[i]];
    }
    for (indsz_t i = 0; i < D.rowVals.size(); ++i) {
        D.rowVals[i] = translator[D.rowVals[i]];
    }

    // Multiply aInv x Y, dimensions are (pairedCols x pairedRows) x (pairedRows x unpairedCols),
    // so results in E of size pairedCols x unpairedCols.
    // X upper triangular because of getPairs, see comments there.
    auto Xi = F2UpTriInvSilentOut(X);

    auto E  = siliLeftMatMultF2(Xi, Y);

    // criticalLow hasn't been updated to reflect that pairings were made and pairs excised yet! That is done
    // at the end of this function, and so translator takes a kept pair index and returns the cell index of the
    // low cell in the pair
    for (indsz_t j = 0; j < numKept; ++j) translator[j] = criticalLow[pairedCols[j]];

    // Result same size as D : (unpairedRows x unpairedCols)
    // !This is an inplace operation on coBdry!
    // Computes D + Z * X^{-1} * Y, i.e. the Schur Complement.
    // After this coBdry is the remaining unreduced part of the coBdry.
    blockProdSum(D, Z, E, coBdry);

    // Updates criticalPart to compensate for removed pairs.
    // R is not treated as silent here, enough though L is (and they come from cols. of the same matrix).
    // Similar to above, does R + L*E = R + L * X^{-1} * Y = * in [MFCPH]
    // Result has size (totalRows x unpairedCols)
    blockProdSumSilentIColsLeft(R, L, E, criticalPart, translator);

    // Updates the row/col to high/low cell mapping.
    for (indsz_t i = 0; i < unpairedRows.size(); ++i) criticalHigh[i] = criticalHigh[unpairedRows[i]];
    criticalHigh.resize(unpairedRows.size());

    for (indsz_t i = 0; i < unpairedCols.size(); ++i) criticalLow[i] = criticalLow[unpairedCols[i]];
    criticalLow.resize(unpairedCols.size());
}

std::tuple<spla::SparseBool, eirene::PairVec, std::vector<indsz_t>> morseLU(spla::SparseBool& coBdry, eirene::PairVec& pairVec, const vector<indsz_t>& highFiltTemp,
             const vector<indsz_t>& lowFiltTemp)
{
    /*
        Here bdryRow/Col is a sparse representation of the coboundary.
        Thus, the rows correspond to higher dimensional cells, and cols to the lower.

        criticalLow contains the indices of lower dimensional cells that were unpaired in the last iteration.
        criticalHigh contains the indices of all high dimensional cells.
        Both of these indexes are defined relative to the unpaired cells passed into this algorithm, as
        opposed to the global indexes in the complex, and the cells they refer to in lowlabels (not passed in)
        are in sorted order w.r.t the filtration.

        Then, unpairedRows/Cols defined below indicate the row/cols that correspond to these critical high/low
        cells. Thus it ends up partially representing the permutation that reduces the matrix.

        As the matrix is reduced, it('s nontrivial part/representation) shrinks, and so do criticalLow/Hi and
       unpairedRows/Cols.

        The difference between criticalLow/Hi and unpairedLo/Hi is that one tracks the row/col inds and the other the
        cell inds.
    */

    vector<indsz_t>& coBdryColPtrs = coBdry.colPtrs;

    indsz_t criticalLowCount = coBdry.numCols(), criticalHighCount = coBdry.numRows;

    vector<indsz_t> criticalHigh(criticalHighCount), criticalLow(criticalLowCount);
    std::iota(criticalLow.begin(), criticalLow.end(), 0);
    std::iota(criticalHigh.begin(), criticalHigh.end(), 0);

    // Set up the indexing so that the paired rows and cols come first, followed by the critical cells in
    // unpairedrows/cols
    indsz_t maxNumPairs = std::min(criticalLowCount, criticalHighCount), numPairs = pairVec.size();

    pairVec.reserve(maxNumPairs);

    // Note excisedpairs is indexed by cell while pairVec is indexed by col/row index.
    auto excisedPairs = eirene::PairVec();
    excisedPairs.reserve(maxNumPairs);

    auto unpairedRows = vector<indsz_t>(criticalHighCount - numPairs),
         unpairedCols = vector<indsz_t>(criticalLowCount - numPairs);

    std::iota(unpairedRows.begin(), unpairedRows.end(), numPairs);
    std::iota(unpairedCols.begin(), unpairedCols.end(), numPairs);

    // criticalPart is related to (is, I think?) the nonidentity part of R, i.e. X^-1*Y, in [MFCPH].
    // pairedPart represents the paired part of the boundary, which is nonsingular.
    // At the end of this routine, criticalPart's cols are inserted at the end of pairedPart, which is returned.
    // Both updated it schurIt.
    spla::SparseBool criticalPart{criticalHighCount, criticalLowCount}, pairedPart{criticalHighCount, criticalLowCount};

    if (pairVec.size())
    {
        schurIt(coBdry, criticalLow, criticalHigh, pairVec, excisedPairs, unpairedRows, unpairedCols, criticalPart,
                pairedPart);
        criticalLowCount  = unpairedCols.size();
        criticalHighCount = unpairedRows.size();
    }

    auto rowFilt = vector<indsz_t>(unpairedRows.size()), colFilt = vector<indsz_t>(unpairedCols.size());

    // While critical cells with nontrivial coboundary exist, or maybe it's while critical cells whose whose cobdry
    // intersects unpaired higher dimensional cells?
    vector<indsz_t> halfPairVec(numPairs);

    while (coBdryColPtrs[criticalLowCount] > 0)
    {
        rowFilt.resize(criticalHigh.size());
        colFilt.resize(criticalLow.size());
        std::transform(criticalHigh.begin(), criticalHigh.end(), rowFilt.begin(),
                       [&highFiltTemp](const indsz_t cellInd) { return highFiltTemp[cellInd]; });
        std::transform(criticalLow.begin(), criticalLow.end(), colFilt.begin(),
                       [&lowFiltTemp](const indsz_t cellInd) { return lowFiltTemp[cellInd]; });

        // Find new pairs to excise from complex (mutates pairVec)
        getPairs(coBdry, rowFilt, colFilt, pairVec);

        // Set the unpaired Rows and Columns to those that are not in the list of paired row/columns.
        halfPairVec.resize(pairVec.size());

        std::transform(pairVec.begin(), pairVec.end(), halfPairVec.begin(),
                       [](const MorsePair pair) { return pair.highInd; });
        unpairedRows = uniqueIntsNotInUnsortedUpTo(halfPairVec, criticalHigh.size());

        std::transform(pairVec.begin(), pairVec.end(), halfPairVec.begin(),
                       [](const MorsePair pair) { return pair.lowInd; });
        unpairedCols = uniqueIntsNotInUnsortedUpTo(halfPairVec, criticalLow.size());

        criticalLowCount  = unpairedCols.size();
        criticalHighCount = unpairedRows.size();

        schurIt(coBdry, criticalLow, criticalHigh, pairVec, excisedPairs, unpairedRows, unpairedCols, criticalPart,
                pairedPart);
    }

    // After the loop, criticalPart has size (totalRows x unpairedCols)
    // pairedPart has size (totalRows x pairedCols)

    // Append criticalPart to the back of pairedPart
    appendSparseToBack(pairedPart, excisedPairs.size(), criticalPart, criticalLowCount);

    // The first (rank of M) elements of tlabels index the nonzero columns of the
    // reduced matrix, corresponding to paired cells.
    std::vector<indsz_t> tlabels(0);
    tlabels.reserve(excisedPairs.size() + criticalLow.size());
    std::transform(excisedPairs.begin(), excisedPairs.end(), std::back_inserter(tlabels),
                   [](const MorsePair& pair) { return pair.lowInd; });
    tlabels.insert(tlabels.end(), criticalLow.begin(), criticalLow.end());

    return std::make_tuple(pairedPart, excisedPairs, tlabels);
}

spla::SparseBool transposeSparseSubmatrix(const spla::SparseBool& sparseMat, const vector<indsz_t>& rows,
                                          const vector<indsz_t>& cols)
{
    /*
        Returns the sparse representation of the transpose of the submatrix given delineated
        by the arguments.

        The outputs' rows are in the order given in cols. Does not assume rows, cols are sorted.

        Notation: input matrix has dimension numRows x (colPtrs.size() - 1)
                  output (transposed) submatrix has dimensions cols.size() x rows.size() = subCols x subRows

        TODO: Bad algo?
    */
    const vector<indsz_t>&rowVals = sparseMat.rowVals, &colPtrs = sparseMat.colPtrs;
    const indsz_t subRows = rows.size();

    spla::SparseBool transMat(cols.size(), rows.size());

    // Compute the column counts of transposed submatrix:
    for (indsz_t elemInd = 0; elemInd < rowVals.size(); ++elemInd)
    {
        auto subRowPos = std::find(rows.begin(), rows.end(), rowVals[elemInd]);
        if (subRowPos == rows.end()) continue;

        indsz_t newRowInd = std::distance(rows.begin(), subRowPos);
        for (indsz_t j = newRowInd + 1; j < subRows + 1; ++j) transMat.colPtrs[j]++;
    }

    const indsz_t numNotZero = transMat.colPtrs.back();
    transMat.rowVals.reserve(numNotZero);

    // Sort the entries into their final resting spots (bad algo)
    // Should this not reindex the rows as well?
    for (indsz_t rowInd : rows)   // for each of the rows (transposed cols)
        // go through the cols (transposed rows) in the order given
        for (indsz_t colInd = 0; colInd < cols.size(); ++colInd)
        {
            const auto col = cols[colInd];
            for (indsz_t colSt = colPtrs[col]; colSt < colPtrs[col + 1]; ++colSt)
                if (rowVals[colSt] == rowInd) transMat.rowVals.push_back(colInd);
        }

    return transMat;
}

auto transposeSparseMatrix(const spla::SparseBool& sparseMat)
{
    /*
     *   Returns the sparse representation of the transpose of the whole sparse matrix.
     *   TODO: Optimize above algo and this one.
     */
    const vector<indsz_t>&rowVals = sparseMat.rowVals, &colPtrs = sparseMat.colPtrs;

    spla::SparseBool transMat(colPtrs.size() - 1, sparseMat.numRows);

    // Compute the column counts of transposed
    for (indsz_t elemInd = 0; elemInd < rowVals.size(); ++elemInd)
    {
        auto subRowPos = rowVals[elemInd];

        for (indsz_t j = subRowPos + 1; j < sparseMat.numRows + 1; ++j) transMat.colPtrs[j]++;
    }

    const indsz_t numNotZero = transMat.colPtrs.back();
    transMat.rowVals.reserve(numNotZero);

    // Sort the entries into their final resting spots (bad algo)
    // Should this not reindex the rows as well?
    for (indsz_t rowInd = 0; rowInd < sparseMat.numRows; ++rowInd)
        for (indsz_t colInd = 0; colInd < sparseMat.colPtrs.size() - 1; ++colInd)
            for (indsz_t colSt = colPtrs[colInd]; colSt < colPtrs[colInd + 1]; ++colSt)
                if (rowVals[colSt] == rowInd) transMat.rowVals.push_back(colInd);

    return transMat;
}

auto maxNonsingularBlockRelative(const spla::SparseBool& mat, const eirene::PairVec& pairVec, const vector<indsz_t>& tid,
                         const indsz_t numLowerPairs, const uint32_t dim, const eirene::Subcomplex& closedComp)
{
    const auto clInds = closedComp.cellInds;
    const auto clFilt = closedComp.dimFilt;
    const auto openDimCard = clFilt[dim] - clFilt[dim - 1];

    const indsz_t numLowCells = mat.numRows - openDimCard, numPairs = pairVec.size(), numCriticalLow = numLowCells - numLowerPairs;

    if (numPairs == 0) return spla::SparseBool(tid.size(), numCriticalLow);

    vector<indsz_t> emptyVec{};
    vector<indsz_t> highPairInds(numPairs);
    vector<indsz_t> translator(mat.numRows);
    for (indsz_t i = 0; i < tid.size(); ++i) {
        translator[tid[i]] = i;
    }

    std::transform(pairVec.begin(), pairVec.end(), highPairInds.begin(), [](const auto& pair) { return pair.highInd; });

    // A has size tid.size() x numCriticalLow
    auto [A, Z] = stackedSubMat(mat, tid, emptyVec, highPairInds);

    // Columns numPairs:numCriticalLow are 0
    A.colPtrs.resize(numCriticalLow + 1);
    std::fill(A.colPtrs.begin() + numPairs, A.colPtrs.end(), A.numNonZero());

    for (indsz_t i = 0; i < A.rowVals.size(); ++i) {
        A.rowVals[i] = translator[A.rowVals[i]];
    }

    return A;
}

auto maxNonsingularBlock(const spla::SparseBool& mat, const eirene::PairVec& pairVec, const vector<indsz_t>& tid,
                         const indsz_t numLowerPairs)
{
    /*
     *  Returns the submatrix of mat, (in this case, the original boundary operator), with rows in tid and
     *  columns in highPairInds, in the order they are in in those vectors.
     *  Since tid is itself formed from the low pair indices, this basically just takes the submatrix
     *  defined by the rows/cols from low/high pair indices, respectively.
     */
    const indsz_t numLowCells = mat.numRows, numPairs = pairVec.size(), numCriticalLow = numLowCells - numLowerPairs;

    if (numPairs == 0) return spla::SparseBool(tid.size(), numCriticalLow);

    vector<indsz_t> emptyVec{};
    vector<indsz_t> highPairInds(numPairs);
    vector<indsz_t> translator(numLowCells);
    for (indsz_t i = 0; i < tid.size(); ++i) {
        translator[tid[i]] = i;
    }

    std::transform(pairVec.begin(), pairVec.end(), highPairInds.begin(), [](const auto& pair) { return pair.highInd; });

    // A has size tid.size() x numCriticalLow
    auto [A, Z] = stackedSubMat(mat, tid, emptyVec, highPairInds);

    // Columns numPairs:numCriticalLow are 0
    A.colPtrs.resize(numCriticalLow + 1);
    std::fill(A.colPtrs.begin() + numPairs, A.colPtrs.end(), A.numNonZero());

    for (indsz_t i = 0; i < A.rowVals.size(); ++i) {
        A.rowVals[i] = translator[A.rowVals[i]];
    }

    return A;
}

}   // namespace spla

namespace eirene
{

auto getCycleRep(const spla::SparseBool& lowBoundary, const std::array<SparseBool, 3>& factors,
                 const std::array<SparseBool, 3>& lowFactors, const PairVec& lowPairs, const vector<indsz_t>& tid,
                 const vector<indsz_t>& lowTid, const indsz_t twoBelowCells, const indsz_t cycleNumber)
{
    // Calculates a cycle rep of dimension (dim - 1)
    const auto& [L, R, Li] = factors;

    // Step one: get the cycleNumber-th column from Li, pass it though tid reindexing, back into complex indexing.
    // Note, row values in columns of Li correspond to cells in (dim - 1).
    const indsz_t numChainSummands = Li.colPtrs[cycleNumber + 1] - Li.colPtrs[cycleNumber];

    vector<indsz_t> summands(numChainSummands + 1);

    for (indsz_t i = 0; i < numChainSummands; ++i)
    {
        auto LiInd  = Li.colPtrs[cycleNumber] + i;
        summands[i] = tid[Li.rowVals[LiInd]];
    }

    // Because matrix is sili, the diagonal value is supressed.
    summands.back() = tid[cycleNumber];

    // Remember in the algo, at each dimension, only the critical cells from the dimension below were fed in.
    // i.e. the rows/cols corresponding to paired cells from the lower dim were removed, specifically when
    // transposing. Similarly, now the row values/low cell inds in the column we got from Li correspond to critical
    // cells, and so we need to calculate what the cycle corresponds to in the full, non-reduced complex. Unless,
    // of course, there were no 'twobelowcells', i.e., dim == 1.
    if (twoBelowCells == 0) { return summands; }

    // Pass the cell through the (dim - 1) -> (dim - 2)
    auto support                    = vector<bool>(twoBelowCells, false);
    const auto& [lowL, lowR, lowLi] = lowFactors;

    for (auto j : summands)
    {
        for (indsz_t k = lowBoundary.colPtrs[j]; k < lowBoundary.colPtrs[j + 1]; ++k)
        {
            indsz_t i  = lowBoundary.rowVals[k];
            support[i] = !support[i];
        }
    }

    // Assemble the image into a (column) vector.
    // Represents the boundary of the cells in summands
    indsz_t suppCard = std::accumulate(support.begin(), support.end(), 0);
    auto brv         = vector<indsz_t>(0);
    brv.reserve(suppCard);

    // Pass from global indexes back into tid-indexing (rows/cols in reduced matrix).
    // This is the reverse direction of that done above.
    // We're talking in (dim - 2) cells now
    for (indsz_t lowTidIndex = 0; lowTidIndex < lowTid.size(); ++lowTidIndex)
    {
        indsz_t lowIndex = lowTid[lowTidIndex];
        if (support[lowIndex]) brv.push_back(lowTidIndex);
    }

    auto bcp = vector<indsz_t>{0, static_cast<indsz_t>(brv.size())};
    auto b   = spla::SparseBool{brv, bcp, lowL.numCols()};

    // Now b is the aforementioned boundary as a vector, with cells indexed by tid.
    // Hit it with RL - we already used the boundary above, so result is
    // "the fundamental cycle matrix w.r.t to doubly minimal basis"
    // Honestly, I have only a vague idea what that means, the map RL is closely related to the
    // discrete gradient flow from Morse Theory for Cell Complexes
    b = spla::siliLeftMatMultF2(lowL, b);
    b = spla::siliLeftMatMultF2(lowR, b);

    // lowPairs is pairing (dim - 2) x (dim - 1)
    vector<indsz_t> lowToHigh(twoBelowCells, 0);
    // lowToHigh[pairVec[i].lowInd] = pairVec[i].highInd
    for (indsz_t i = 0; i < lowPairs.size(); ++i) lowToHigh[lowPairs[i].lowInd] = lowPairs[i].highInd;

    // "Nonzero entries might lie in nonbasis rows"
    // We can be sure that b's nonzero rowvals correspond to paired cells because it's the output of RL above
    for (indsz_t i = 0; i < b.rowVals.size(); ++i) b.rowVals[i] = lowToHigh[lowTid[b.rowVals[i]]];

    //
    b.rowVals.insert(b.rowVals.end(), summands.begin(), summands.end());

    return b.rowVals;
}

auto toTidIndexing(const std::vector<indsz_t>& indexes, const std::vector<indsz_t>& tid, std::vector<indsz_t>& helper,
                   const indsz_t numCells)
{
    helper.resize(numCells);
    for (indsz_t i = 0; i < tid.size(); ++i) helper[tid[i]] = i;
    return perm::reindex(indexes, helper);
};

void extractHomology(const Complex& complex, vector<eirene::MorseReducedResult>& reducedVec, const std::optional<eirene::Subcomplex>& closedCompO)
{
    /*
     *  This stuff happens in getcycle and unpack! in eirene.jl
     *
     *  After doing the Morse LU transform on all the boundary matrices, calculate homology.
     *  Section 5 of MFCPH.
     *
     *  Greg's notes in Eirene.jl:
     *  -	L, Li, and R are all indexed by tid, just like (trv,tcp)
     *  - 	L (not Li) is the transpose of (trv,tcp)
     *   - 	up to permutation and deletion of some zero rows, RLM is the fundamental
     *       cycle matrix of the d-dimensional boundary operator with respect
     *       to the doubly-minimal basis we have constructed, where M the submatrix of
     *       the boundary with rows indexed by tid.
     *
     * # of pairs dim - 1 x dim == rank of boundary operator from dim -> (dim - 1), contained in reducedVec[dim]
     * The first (rank of boundary) elements of tid index the nonzero columns of the reduced matrix
     *
     * Note, to pass from complex-relative indexes into tid-indexes, we do
     *    compInd |-> (i : tid[i] == compInd)
     * Thus, to go the other way, we do
     *    tidInd |-> tid[tidInd]
     */


    const auto& dimPatt = complex.dimPatt;
    const auto& splitComplex = complex.splitCells;
    const auto& sortedSplitFilt = complex.splitFiltInds;

    vector<std::array<spla::SparseBool, 3>> dimFactors(numberOfDim(dimPatt));
    const indsz_t topDim = numberOfDim(dimPatt);
    vector<indsz_t> lowTranslator(0);

    for (indsz_t dim = 1; dim < numberOfDim(dimPatt); ++dim)
    {
        const indsz_t numLowCells = splitComplex[dim - 1].numCols();

        // leftFactor was highCellNum x numCriticalLow
        SparseBool Lt{reducedVec[dim].leftFactor};

        const auto& tid = reducedVec[dim].tid;

        Lt.rowVals = toTidIndexing(Lt.rowVals, tid, lowTranslator, numLowCells);
        auto Lit   = spla::F2UpTriInvSilentInOut(Lt);
        auto Li = spla::transposeSparseMatrix(Lit);

        // for the last dim, only want to calculate Li
        if (dim == topDim)
        {
            if (!(closedCompO == std::nullopt)) 
            dimFactors[dim][2] = Li;
            continue;
        }

        SparseBool nonSingularReduced;
        if (closedCompO == std::nullopt) {
            // These are the non-singular columns of the reduced (co?)boundary matrix

            std::ofstream file;
            file.open("cpp_debug.txt", std::ofstream::out | std::ofstream::trunc);
            for (auto p: reducedVec[dim].pairVec) file << p.highInd + 1 << "," << p.lowInd + 1 << "\n";
            file.close();

            nonSingularReduced = spla::maxNonsingularBlock(splitComplex[dim], reducedVec[dim].pairVec, tid,
                                                                reducedVec[dim - 1].pairVec.size());
        } else {
            nonSingularReduced = spla::maxNonsingularBlockRelative(splitComplex[dim], reducedVec[dim].pairVec, tid,
                                                                reducedVec[dim - 1].pairVec.size(), dim, closedCompO.value());
        }

        auto numnlpl = splitComplex[dim - 1].numCols() - reducedVec[dim - 1].pairVec.size();
        /*
         * Note according to Lemma 7 of MFCPH, when restructed to the paired (= nonsingular)
         * portion of the complex, LAR = I, and so LA = R^(-1).
         */
        auto L = spla::transposeSparseMatrix(Lt);
        auto B = spla::siliLeftMatMultF2(L, nonSingularReduced);
        auto R = spla::F2UpTriInvSilentOut(B);

        dimFactors[dim] = std::array<SparseBool, 3>{L, R, Li};
    }

    const auto getCycleReps = [&](const indsz_t dim, const indsz_t twoBelowCells)
    {
        /*
         * Computes the cocycle representatives for one dimension.
         * Calling this with dim = n returns representatives of dimension n - 1.
         * Factors and tid are those from dimension n, lowFactors/lowTid are dim (n - 1).
         * lowBoundary is from (n - 1) -> (n - 2)
         */

        const auto& pairs    = reducedVec[dim].pairVec;
        const auto& lowPairs = reducedVec[dim - 1].pairVec;
        const auto& tid      = reducedVec[dim].tid;
        const auto& lowTid   = reducedVec[dim - 1].tid;

        vector<vector<indsz_t>> cycles;
        vector<double> birthDeath;
        cycles.reserve(tid.size());
        birthDeath.reserve(tid.size() * 2);

        // mortal cycles, from pairs
        for (indsz_t pairInd = 0; pairInd < pairs.size(); ++pairInd)
        {
            auto pair = pairs[pairInd];
            if (sortedSplitFilt[dim - 1][pair.lowInd] == sortedSplitFilt[dim][pair.highInd]) continue;

            cycles.push_back(getCycleRep(splitComplex[dim - 1], dimFactors[dim], dimFactors[dim - 1], lowPairs, tid,
                                         lowTid, twoBelowCells, pairInd));

            birthDeath.push_back(complex.uniqFiltVals[sortedSplitFilt[dim - 1][pair.lowInd]]);
            birthDeath.push_back(complex.uniqFiltVals[sortedSplitFilt[dim][pair.highInd]]);
        }

        // immortals / critical cells
        for (indsz_t cycleInd = pairs.size(); cycleInd < tid.size(); ++cycleInd)
        {
            cycles.push_back(getCycleRep(splitComplex[dim - 1], dimFactors[dim], dimFactors[dim - 1], lowPairs, tid,
                                         lowTid, twoBelowCells, cycleInd));

            birthDeath.push_back(complex.uniqFiltVals[sortedSplitFilt[dim - 1][tid[cycleInd]]]);
            birthDeath.push_back(INFINITY);
        }

        return std::tuple(cycles, birthDeath);
    };

    for (indsz_t dim = 1; dim < numberOfDim(dimPatt); ++dim)
    {
        const indsz_t twoBelowCells = dim > 1 ? splitComplex[dim - 2].numCols() : 0;
        auto [repVec, birthDeathVec] = getCycleReps(dim, twoBelowCells);

        reducedVec[dim].cycleReps = repVec;
        reducedVec[dim].birthDeath = birthDeathVec;
    }
}

vector<MorseReducedResult> performPersistence(const Complex& complex)
{
    /*
     *   Runs the persistence algorithm:
     *   First reduce all of the boundary operators with Discrete Morse Theory/Eirene algo.
     *   Then can extract barcodes and cycle reps from the reduced complex.
     */

    const vector<SparseBool>& splitBdry         = complex.splitCells;
    const vector<indsz_t>& dimPatt                 = complex.dimPatt;
    const vector<vector<indsz_t>>& sortedSplitFilt = complex.splitFiltInds;

    const indsz_t numDim = numberOfDim(dimPatt);

    vector<MorseReducedResult> resultVec(numDim + 1);
    vector<vector<indsz_t>> prePairs(numDim + 1, vector<indsz_t>{});

    // Calculating cohomology: start from low dimensions
    // Each iteration factorizes boundary from dim to (dim - 1)
    for (indsz_t dim = 1; dim < numDim; ++dim)
    {
        // Get the high cell indices from the pairs calculated in the lower dimension
        vector<indsz_t> lowbasisnames(resultVec[dim - 1].pairVec.size());
        std::transform(resultVec[dim - 1].pairVec.begin(), resultVec[dim - 1].pairVec.end(), lowbasisnames.begin(),
                       [](const MorsePair& pair) { return pair.highInd; });

        const SparseBool& dimBoundary = splitBdry[dim];
        indsz_t numHighCells = dimBoundary.numCols(), numLowCells = splitBdry[dim - 1].numCols();

        // lowlabels will hold the indices of cells of dimension (dim - 1) which do not appear in
        // lowbasisnames, i.e. were unpaired in any previous iterations.
        // lowlabels/hilabels will be used after reducing to map from the indices used in the LU transform back
        // to the correct global/cellular indexes.
        vector<indsz_t> lowlabels = uniqueIntsNotInUnsortedUpTo(lowbasisnames, numLowCells);

        // also, low labels will be in permuted order
        // sortedSplitFilt is indices into a vector of unique filtration values in REVERSE
        // sorted order, so a bigger index means the cell is born earlier
        vector<double> subFilt = vector<double>(lowlabels.size());
        std::transform(lowlabels.begin(), lowlabels.end(), subFilt.begin(),
                       [&](indsz_t ind) { return sortedSplitFilt[dim - 1][ind]; });

        vector<indsz_t> lowCellFiltPerm = getSortPerm(subFilt);

        perm::reindexInplace(lowCellFiltPerm, lowlabels);

        // Nothing to reduce. For example, when dim == numDim - 1 and the boundary operator is zero.
        if (dimBoundary.numNonZero() == 0)
        {
            resultVec[dim].tid = lowlabels;

            // the reduced coboundary is zero, but make sure it's the correct size...
            resultVec[dim].leftFactor.colPtrs.resize(lowlabels.size() + 1);
            resultVec[dim].leftFactor.numRows = dimBoundary.numCols();

            continue;
        }

        vector<indsz_t> highlabels(numHighCells);
        std::iota(highlabels.begin(), highlabels.end(), 0);

        // Transpose from boundary to coboundary. Note this also reindexes the
        // input's rows/output's columns according to lowlabels
        auto dimCoBdry = spla::transposeSparseSubmatrix(dimBoundary, lowlabels, highlabels);

        // Bunch of copies. Good spot for some memory mgmnt..
        auto highFiltTemp = vector<indsz_t>(highlabels.size());
        for (indsz_t i = 0; i < highlabels.size(); ++i) highFiltTemp[i] = sortedSplitFilt[dim][highlabels[i]];

        auto lowFiltTemp = vector<indsz_t>(lowlabels.size());
        for (indsz_t i = 0; i < lowlabels.size(); ++i) lowFiltTemp[i] = sortedSplitFilt[dim - 1][lowlabels[i]];

        PairVec pairs;

        // The first pairs.size elements of tlabels are the indices of low cells in pairs (correspond to
        // nonsingular/cleared part), rest are critical cells
        auto [leftFactor, excisedPairs, tlabels] = spla::morseLU(dimCoBdry, pairs, highFiltTemp, lowFiltTemp);

        // Do the aforementioned reindexing into complex-relative indices: pass everything through low/highlabels
        std::transform(excisedPairs.begin(), excisedPairs.end(), excisedPairs.begin(),
                       [&lowlabels, &highlabels](const MorsePair& pair) {
                           return MorsePair{lowlabels[pair.lowInd], highlabels[pair.highInd]};
                       });

        leftFactor.rowVals = perm::reindex(leftFactor.rowVals, lowlabels);
        resultVec[dim]     = eirene::MorseReducedResult{leftFactor,                      // trv and tcp
                                                    perm::reindex(tlabels, lowlabels),   // tid
                                                    excisedPairs,                        // plo/hi
                                                    vector<vector<indsz_t>>(0)};
    }

    extractHomology(complex, resultVec, std::nullopt);

    return resultVec;
}

}   // namespace eirene

namespace vrips
{

vector<double> offDiagMean(const vector<indsz_t>& indMat, indsz_t numPoints)
{
    // Mean of values in each row off the diagonal.
    vector<double> result = vector<double>(numPoints, 0.0f);

    for (indsz_t j = 0; j < numPoints; ++j)
    {
        double acc = 0.0;

        for (indsz_t k = 0; k < numPoints; ++k)
        {
            acc += k == j ? 0.0 : indMat[j * numPoints + k];
            result[j] = acc / (numPoints - 1);
        }
    }

    return result;
}

std::tuple<SparseBool, vector<indsz_t>> generate1Simplices(const vector<indsz_t>& indMat, indsz_t numPoints,
                                                           const vector<indsz_t>& indPerm)
{
    // Find number of pairs of points with finite (relevant) distances
    indsz_t finiteFiltCount = 0;

    for (indsz_t i = 0; i < numPoints; ++i)
        for (indsz_t j = i + 1; j < numPoints; ++j)
            if (indMat[indPerm[j] * numPoints + indPerm[i]] > 0) finiteFiltCount += 1;

    auto bdry     = SparseBool(numPoints, numPoints);
    auto& rowVals = bdry.rowVals;
    auto& colPtrs = bdry.colPtrs;
    auto filtVals = vector<indsz_t>();
    rowVals.reserve(finiteFiltCount);
    filtVals.reserve(finiteFiltCount);

    // In the 2D simplex (edge) case, for each point, look at each other point and create an edge at the filtration
    // index correspoding to the distane between them.
    for (indsz_t colInd = 0; colInd < numPoints; ++colInd)
    {
        colPtrs[colInd + 1] = colPtrs[colInd];
        for (indsz_t rowInd = colInd + 1; rowInd < numPoints; ++rowInd)
        {
            const indsz_t filtInd = indMat[indPerm[rowInd] * numPoints + indPerm[colInd]];
            // 0 corresponds to distance greater than the cutoff, filtration value Infinity
            if (filtInd > 0)
            {
                colPtrs[colInd + 1] += 1;
                rowVals.push_back(rowInd);
                filtVals.push_back(filtInd);
            }
        }
    }

    return std::tuple(bdry, filtVals);
}

vector<indsz_t> integersInSameOrder(const vector<indsz_t>& v)
{
    /*
     *  Returns the permutation z on {1,...,length(v)} such z[i]<z[j] iff either
     *  (a) v[i] < v[j], or
     *  (b) v[i] = v[j] and i < j
     *  So z[i] >= z[j] ==> (v[i] > v[j]) or (v[i] = v[j] and i >= j)
     *
     *  For example: iiso([1, 2, 3, 1, 3]) = [1, 3, 4, 2, 5] and
     *  ([1, 2, 3, 1, 3])[[1, 3, 4, 2, 5]] = [1, 3, 1, 2, 3]
     *  Or, applying the inverse: [1, 1, 2, 3, 3]
     */

    if (v.size() == 0) return vector<indsz_t>();

    const auto vLen = v.size();
    const auto vMax = *std::max_element(v.begin(), v.end());
    const auto vMin = (*std::min_element(v.begin(), v.end())) - 1;

    auto helper = vector<indsz_t>(vMax - vMin, 0);
    auto result = vector<indsz_t>(vLen, 0);

    for (indsz_t i = 0; i < vLen; ++i)
    {
        const auto ind = v[i] - vMin - 1;
        helper[ind] += 1;
    }

    uint64_t prevSum = 1;
    for (uint64_t i = 0; i < helper.size(); ++i)
    {
        const uint64_t sum = prevSum + helper[i];
        helper[i]          = prevSum;
        prevSum            = sum;
    }

    for (uint64_t i = 0; i < vLen; ++i)
    {
        const uint64_t u = v[i];
        result[i]        = helper[u - vMin - 1] - 1;
        helper[u - vMin - 1] += 1;
    }

    return result;
}

void faceUpdateHelper(indsz_t& faceCount, vector<indsz_t>& r, vector<indsz_t>& z, indsz_t k, indsz_t farFilt,
                      indsz_t stepSize, vector<indsz_t>& s, indsz_t i)
{
    if (faceCount >= r.size())
    {
        const auto toApp = vector<indsz_t>(stepSize, 0);
        r.insert(r.end(), toApp.begin(), toApp.end());
        z.insert(z.end(), toApp.begin(), toApp.end());
        s.insert(s.end(), toApp.begin(), toApp.end());
    }

    r[faceCount] = k;
    z[faceCount] = farFilt;
    s[faceCount] = i;

    faceCount += 1;
}

std::tuple<SparseBool, vector<indsz_t>, vector<indsz_t>> generate2Simplices(const SparseBool& oneBdry,
                                                                            const vector<indsz_t>& filtVals,
                                                                            const indsz_t numPoints)
{
    const auto& rowVals    = oneBdry.rowVals;
    const auto& colPtrs    = oneBdry.colPtrs;
    const indsz_t numVerts = colPtrs.size() - 1;
    const indsz_t numEdges = rowVals.size();
    const indsz_t stepSize = numPoints * numPoints;
    indsz_t faceCount      = 0;

    auto nzValColInds = vector<indsz_t>(numEdges, 0);

    for (uint64_t i = 0; i < numPoints; ++i)
        for (uint64_t j = colPtrs[i]; j < colPtrs[i + 1]; ++j) nzValColInds[j] = i;

    auto iso            = integersInSameOrder(rowVals);
    auto coBdryRowVals  = vector<indsz_t>(numEdges, 0);
    auto coBdryFiltVals = vector<indsz_t>(numEdges, 0);
    for (uint64_t i = 0; i < iso.size(); ++i)
    {
        coBdryRowVals[iso[i]]  = nzValColInds[i];
        coBdryFiltVals[iso[i]] = filtVals[i];
    }

    // Really just the transpose of the above, the "far faces" boundary
    auto coBdryColPtrs = vector<indsz_t>(numPoints + 1);
    for (const auto face : rowVals) coBdryColPtrs[face + 1] += 1;
    coBdryColPtrs[0] = 0;
    for (uint64_t i = 1; i < numPoints + 1; ++i) coBdryColPtrs[i] = coBdryColPtrs[i - 1] + coBdryColPtrs[i];

    auto aDist = vector<indsz_t>(numPoints, 0);
    auto iDist = vector<indsz_t>(numPoints, 0);

    auto twoBdryRowVals  = vector<indsz_t>(numEdges, 0);
    auto twoBdryFiltVals = vector<indsz_t>(numEdges, 0);
    auto s               = vector<indsz_t>(numEdges, 0);

    auto clawVec       = vector<indsz_t>(numPoints, 0);
    auto nCheckedEdges = vector<bool>(numEdges, true);

    // Going to do nested iteration over points: point a, point i, point j, etc.
    for (indsz_t aInd = 0; aInd < numPoints; ++aInd)
    {
        // we set aDist[cellInd] = persistence index distance between point a and point cellInd
        // persistence values of the edges in column aInd
        std::fill(aDist.begin(), aDist.end(), 0);
        for (indsz_t edgeInd = colPtrs[aInd]; edgeInd < colPtrs[aInd + 1]; ++edgeInd)
            aDist[rowVals[edgeInd]] = filtVals[edgeInd];

        // for each edge incident with a
        for (indsz_t aEdgeInd = colPtrs[aInd]; aEdgeInd < colPtrs[aInd + 1]; ++aEdgeInd)
        {
            // edge is between aInd and iInd, persistance index dai
            const auto iInd = rowVals[aEdgeInd];
            const auto dai  = filtVals[aEdgeInd];   // == aDist[iInd]?

            // iDist[row] = persistence index between point at iInd and point at row
            std::fill(iDist.begin(), iDist.end(), 0);
            for (indsz_t edgeInd = colPtrs[iInd]; edgeInd < colPtrs[iInd + 1]; ++edgeInd)
                iDist[rowVals[edgeInd]] = filtVals[edgeInd];

            for (indsz_t edgeInd = coBdryColPtrs[iInd]; edgeInd < coBdryColPtrs[iInd + 1]; ++edgeInd)
                iDist[coBdryRowVals[edgeInd]] = coBdryFiltVals[edgeInd];

            // for each edge incident with i
            for (indsz_t ijEdgeInd = colPtrs[iInd]; ijEdgeInd < colPtrs[iInd + 1]; ijEdgeInd++)
            {
                if (!nCheckedEdges[ijEdgeInd]) continue;
                auto jInd = rowVals[ijEdgeInd];
                auto dij  = filtVals[ijEdgeInd];

                // check if we should make a 2 simplex (a, i, j)
                // we can say we should maybe if the edge (i, j) is born after (a, i) and (a, j)
                // we can say we should not if the edge (i, j) is born after (a, i) and (a, j)
                if (!(dij <= dai && dij <= aDist[jInd])) continue;

                // only have to do this once per edge over the whole outer iteration
                nCheckedEdges[ijEdgeInd] = false;

                // over edges to those l with lInd < jInd
                std::fill(clawVec.begin(), clawVec.end(), 0);
                for (indsz_t jCoEdgeInd = coBdryColPtrs[jInd]; jCoEdgeInd < coBdryColPtrs[jInd + 1]; ++jCoEdgeInd)
                {
                    auto lInd = coBdryRowVals[jCoEdgeInd];
                    // sorted according to filtration: only consider l born after i
                    if (lInd >= iInd)
                        break;
                    else if (iDist[lInd] != 0)
                        clawVec[lInd] = std::min(iDist[lInd], coBdryFiltVals[jCoEdgeInd]);
                }

                // the case where jInd < lInd
                for (indsz_t lEdgeInd = colPtrs[jInd]; lEdgeInd < colPtrs[jInd + 1]; ++lEdgeInd)
                {
                    auto lInd = rowVals[lEdgeInd];
                    auto djl  = filtVals[lEdgeInd];
                    auto dal  = aDist[lInd];
                    auto dil  = iDist[lInd];

                    if (dal < dij && dal < djl && dal < dil)
                    {
                        auto dijl     = std::min(dij, std::min(dil, djl));
                        bool keepface = true;

                        // look at edges k to b with b born
                        for (indsz_t kCoEdgeInd = coBdryColPtrs[lInd]; kCoEdgeInd < coBdryColPtrs[lInd + 1];
                             ++kCoEdgeInd)
                        {
                            auto bInd = coBdryRowVals[kCoEdgeInd];
                            if (bInd >= iInd) break;
                            // keep the face if the max(dil, djl, djl) > dijk
                            else if (std::min(clawVec[bInd], coBdryFiltVals[kCoEdgeInd]) >= dijl)
                            {
                                keepface = false;
                                break;
                            }
                        }
                        if (keepface)
                            faceUpdateHelper(faceCount, twoBdryRowVals, twoBdryFiltVals, lEdgeInd, dijl, stepSize, s,
                                             iInd);
                    }
                }
            }
        }
    }

    indsz_t holdi;
    for (indsz_t edgeInd = 0; edgeInd < nCheckedEdges.size(); ++edgeInd)
    {
        if (!nCheckedEdges[edgeInd]) continue;
        indsz_t closeFaceInd = nzValColInds[edgeInd];
        indsz_t farFaceInd   = rowVals[edgeInd];
        indsz_t closeFarDist = filtVals[edgeInd];

        if (edgeInd == 0 || closeFaceInd != holdi)
        {
            std::fill(iDist.begin(), iDist.end(), 0);

            for (indsz_t edgeInd = colPtrs[closeFaceInd]; edgeInd < colPtrs[closeFaceInd + 1]; ++edgeInd)
                iDist[rowVals[edgeInd]] = filtVals[edgeInd];

            for (indsz_t edgeInd = coBdryColPtrs[closeFaceInd]; edgeInd < coBdryColPtrs[closeFaceInd + 1]; ++edgeInd)
                iDist[coBdryRowVals[edgeInd]] = coBdryFiltVals[edgeInd];

            holdi = closeFaceInd;
        }
        std::fill(clawVec.begin(), clawVec.begin() + closeFaceInd, 0);

        for (indsz_t lp = coBdryColPtrs[farFaceInd]; lp < coBdryColPtrs[farFaceInd + 1]; ++lp)
        {
            auto lInd = coBdryRowVals[lp];
            if (lInd >= closeFaceInd)
                break;
            else if (iDist[lInd] != 0)
                clawVec[lInd] = std::min(iDist[lInd], coBdryFiltVals[lp]);
        }

        // facsimile of above
        for (indsz_t kInd = colPtrs[farFaceInd]; kInd < colPtrs[farFaceInd + 1]; ++kInd)
        {
            auto k   = rowVals[kInd];
            auto dik = iDist[k];
            if (dik == 0) continue;
            auto djk      = filtVals[kInd];
            auto dijk     = std::min(closeFarDist, std::min(dik, djk));
            bool keepface = true;

            for (indsz_t bp = coBdryColPtrs[k]; bp < coBdryColPtrs[k + 1]; ++bp)
            {
                auto b = coBdryRowVals[bp];
                if (b >= closeFaceInd)
                    break;
                else if (std::min(clawVec[b], coBdryFiltVals[bp]) >= dijk)
                {
                    keepface = false;
                    break;
                }
            }
            if (keepface)
                faceUpdateHelper(faceCount, twoBdryRowVals, twoBdryFiltVals, kInd, dijk, stepSize, s, closeFaceInd);
        }
    }

    const auto num3faces = faceCount;
    twoBdryRowVals.erase(twoBdryRowVals.begin() + num3faces, twoBdryRowVals.end());
    twoBdryFiltVals.erase(twoBdryFiltVals.begin() + num3faces, twoBdryFiltVals.end());
    s.erase(s.begin() + num3faces, s.end());

    iso         = integersInSameOrder(s);
    auto isoInv = std::vector<indsz_t>(iso.size());
    for (uint64_t i = 0; i < iso.size(); ++i) isoInv[iso[i]] = i;

    std::copy(isoInv.begin(), isoInv.end(), iso.begin());
    perm::reindexInplace(iso, twoBdryRowVals);
    perm::reindexInplace(isoInv, twoBdryFiltVals);

    auto fv3 = vector<indsz_t>(numVerts + 1, 0);
    fv3[0]   = 0;

    for (auto faceInd : s) fv3[faceInd + 1] += 1;
    for (indsz_t i = 1; i < numVerts + 1; ++i) fv3[i] = fv3[i - 1] + fv3[i];

    // for each edge (faceInd) look at the corresponding 0-cell/face, if it's unpaired and the filtration values are the
    // same, make a pair: set prepairs[pairIndex] = faceInd/highIndex and so twoBdryRowVals[prePairs[parIndex]] =
    // twoBdryRowVals[faceInd] = lowCellInd
    indsz_t pairIndex   = 0;
    auto prepairs       = vector<indsz_t>(numEdges, 0);
    auto nonPairedEdges = vector<bool>(numEdges, true);

    uint64_t maxVal = 0;
    for (auto f : twoBdryRowVals) maxVal = std::max(maxVal, f);

    for (indsz_t faceInd = 0; faceInd < num3faces; ++faceInd)
    {
        auto edge = twoBdryRowVals[faceInd];
        if (nonPairedEdges[edge] && twoBdryFiltVals[faceInd] == filtVals[edge])
        {
            nonPairedEdges[edge] = false;
            prepairs[pairIndex]  = faceInd;
            pairIndex += 1;
        }
    }
    prepairs.erase(prepairs.begin() + pairIndex, prepairs.begin() + numEdges);

    auto toRet = SparseBool(twoBdryRowVals, fv3, num3faces);
    return std::tuple(toRet, twoBdryFiltVals, prepairs);
}

std::tuple<eirene::Complex, vector<indsz_t>, vector<vector<indsz_t>>> distMatToVRComplex(const vector<indsz_t>& indMat,
                                                                                         indsz_t numPoints,
                                                                                         indsz_t maxSimplexDim)
{
    vector<vector<indsz_t>> splitFiltInds(maxSimplexDim + 1);
    vector<vector<indsz_t>> prePairs(maxSimplexDim + 1);
    vector<SparseBool> farFaces(maxSimplexDim + 1);

    splitFiltInds[maxSimplexDim] = vector<indsz_t>();
    farFaces[maxSimplexDim]      = SparseBool();
    prePairs[maxSimplexDim]      = vector<indsz_t>();

    vector<double> w         = offDiagMean(indMat, numPoints);
    vector<indsz_t> sortPerm = getSortPerm(w, true);

    // the 0 bdry is diagonal - the matrix is sili
    auto zeroBdry    = SparseBool(numPoints, numPoints);
    zeroBdry.rowVals = vector<indsz_t>(numPoints, 0);
    std::iota(zeroBdry.rowVals.begin(), zeroBdry.rowVals.end(), 0);
    std::iota(zeroBdry.colPtrs.begin(), zeroBdry.colPtrs.end(), 0);
    farFaces[0] = zeroBdry;

    // set splitFiltInds[0] to the diagonal of the permuted index matrix
    splitFiltInds[0] = vector<indsz_t>(numPoints, 0);
    for (indsz_t i = 0; i < numPoints; ++i) splitFiltInds[0][i] = indMat[sortPerm[i] * numPoints + sortPerm[i]];
    prePairs[0] = vector<indsz_t>();

    // the 1-boundary is encoded by a sparse matrix with nonzero value at (i,j) when an edge exists for a finite
    // filtration value between (i, j) for i > j, the diagonal is 0
    const auto [b2, f2] = generate1Simplices(indMat, numPoints, sortPerm);

    farFaces[1]      = b2;
    splitFiltInds[1] = f2;
    prePairs[1]      = std::vector<indsz_t>();

    if (maxSimplexDim == 3)
    {
        const auto [b3, z3, pp3] = generate2Simplices(b2, f2, numPoints);
        farFaces[2]              = b3;
        splitFiltInds[2]         = z3;
        prePairs[2]              = pp3;
    }

    return std::tuple(
        eirene::Complex{
            .splitCells    = farFaces,
            .splitFiltInds = splitFiltInds,
        },
        sortPerm, prePairs);
}

struct DenseIndMat
{
    indsz_t rows;
    indsz_t cols;
    std::vector<indsz_t> vals;

    DenseIndMat(indsz_t rows, indsz_t cols) : rows(rows), cols(cols) { vals.resize(rows * cols); }

    indsz_t* at(indsz_t row, indsz_t col) { return &vals[cols * row + col]; }

    indsz_t get(indsz_t row, indsz_t col) const { return vals[cols * row + col]; }
};

DenseIndMat buildAllFromClose(const SparseBool& lowBdry, const SparseBool& highBdry, DenseIndMat& lclosefaces,
                              const vector<indsz_t>& selectedcolumnindices)
{
    const indsz_t m            = highBdry.numCols();
    const indsz_t numHighCells = highBdry.numRows;

    const indsz_t numSelected = selectedcolumnindices.size();
    const indsz_t rowDepth    = lclosefaces.rows;
    const indsz_t sd          = rowDepth + 1;

    auto hclosefaces = DenseIndMat(sd + 1, numSelected);

    if (numSelected == 0) return hclosefaces;
    auto rosettacol = vector<indsz_t>(*std::max_element(lowBdry.rowVals.begin(), lowBdry.rowVals.end()));
    auto columnSupp = vector<bool>(numHighCells, false);

    for (const auto v : selectedcolumnindices) columnSupp[v] = true;
    indsz_t columnmarker = 0;
    const auto& lColPtrs = lowBdry.colPtrs;
    const auto& lRowVals = lowBdry.rowVals;

    const auto& hColPtrs = highBdry.colPtrs;
    const auto& hRowVals = highBdry.rowVals;
    for (indsz_t i = 0; i < m; ++i)
    {
        for (indsz_t j = lColPtrs[i]; j < lColPtrs[i + 1]; ++j) { rosettacol[lRowVals[j] - 1] = j; }

        for (indsz_t j = hColPtrs[i]; j < hColPtrs[i + 1]; ++j)
        {
            if (columnSupp[j])
            {
                auto farface = hRowVals[j];
                for (indsz_t k = 0; k < rowDepth; ++k)
                {
                    *hclosefaces.at(k, columnmarker) = rosettacol[*lclosefaces.at(k, farface) - 1];
                }
                *hclosefaces.at(rowDepth, columnmarker) = rosettacol[lRowVals[farface] - 1];
                *hclosefaces.at(sd, columnmarker)       = farface;
                columnmarker += 1;
            }
        }
    }
    return hclosefaces;
}

DenseIndMat buildCloseFromClose(const SparseBool& lowBdry, const SparseBool& highBdry, const DenseIndMat& lclosefaces,
                                indsz_t facecard)
{
    const indsz_t m = highBdry.numCols();
    const indsz_t n = highBdry.rowVals.size();

    auto hclosefaces = DenseIndMat(facecard, n);
    if (n == 0) return hclosefaces;
    const indsz_t rowdepth = facecard - 1;
    auto rosettacol        = vector<indsz_t>(
        *std::max_element(lowBdry.rowVals.begin(), lowBdry.rowVals.end()));   // Array{Int64}(undef,maximum(lrowval))

    const auto& lColPtrs = lowBdry.colPtrs;
    const auto& lRowVals = lowBdry.rowVals;
    for (uint64_t i = 0; i < m; ++i)
    {
        for (uint64_t j = lColPtrs[i]; j < lColPtrs[i + 1]; ++j) { rosettacol[lRowVals[j]] = j; }

        for (uint64_t j = highBdry.colPtrs[i]; j < highBdry.colPtrs[i + 1]; ++j)
        {
            auto farface = highBdry.rowVals[j];
            // rowdepth is always 1
            for (uint64_t k = 0; k < rowdepth; ++k) { *hclosefaces.at(k, j) = rosettacol[lclosefaces.get(k, farface)]; }
            *hclosefaces.at(facecard, j) = rosettacol[lRowVals[farface]];
        }
    }

    return hclosefaces;
}

DenseIndMat buildAllFromFar(const eirene::Complex& farFaces, const vector<indsz_t>& colsInOrder, uint64_t dim)
{
    const indsz_t numPoints = farFaces.splitCells[0].numCols();

    if (dim == 0) return DenseIndMat(0, numPoints);

    const auto& oneBdry          = farFaces.splitCells[1];
    const indsz_t oneBdryNonzero = oneBdry.colPtrs.back();
    auto lclosefaces             = DenseIndMat(1, oneBdryNonzero);   // one col for each edge

    for (uint64_t i = 0; i < numPoints; ++i)
        for (uint64_t j = oneBdry.colPtrs[i]; j < oneBdry.colPtrs[i + 1]; ++j) *lclosefaces.at(0, j) = i;

    if (dim == 1)
    {
        auto toRet = DenseIndMat(2, colsInOrder.size());
        for (indsz_t i = 0; i < colsInOrder.size(); ++i)
        {
            *toRet.at(0, i) = *lclosefaces.at(0, colsInOrder[i]);
            *toRet.at(1, i) = farFaces.splitCells[dim].rowVals[colsInOrder[i]];
        }
        return toRet;
    }

    for (indsz_t i = 2; i < dim - 1; ++i)
    {
        lclosefaces = buildCloseFromClose(farFaces.splitCells[i - 1], farFaces.splitCells[i], lclosefaces, i - 1);
    }
    lclosefaces = buildAllFromClose(farFaces.splitCells[dim - 1], farFaces.splitCells[dim], lclosefaces, colsInOrder);
    return lclosefaces;
}

constexpr auto INDSZ_LIM = std::numeric_limits<indsz_t>::max();
SparseBool presparsefull2unsortedsparsetranspose(DenseIndMat& M, const vector<indsz_t>& row02row1translator,
                                                 const vector<indsz_t>& col02col1translator)
{
    /*
    Accepts an mxn integer array, M.  To M we implicitly associate an array N,
    as follows.
    If M has no nonzero entries, then N is the zeron-one array such that
    supp(N[:,i]) = M[:,i].  Here we regard M[:,i] as a set.  I *believe* one ignores
    duplicate entries, but have not checked. A similar interpretation holds when M
    has zero entries - one simply discards the zeros.

    Let S denote the support of N, r denote row02row1translator, and c denote
    col02col1translator.  This function returns the data specifying a sparse
    zero-one matrix whose support is
    { [c[j],r[i]] : [i,j] \in N and neither c[j] nor r[i] are zero }.
    */
    const indsz_t mRows = M.rows;
    const indsz_t mCols = M.cols;

    if (mCols == 0)
    {
        auto rowval1    = vector<indsz_t>(0);
        auto colptr1    = vector<indsz_t>(1);
        indsz_t maxElem = 0;
        for (auto v : row02row1translator)
            if (v < INDSZ_LIM) maxElem = std::max(v, maxElem);
        const auto maxVal = row02row1translator.size() == 0 ? 1 : 1 + maxElem;
        colptr1.resize(maxVal);
        std::fill(colptr1.begin(), colptr1.end(), 0);

        return SparseBool(rowval1, colptr1, mCols);
    }

    for (indsz_t i = 0; i < mRows; ++i)
    {
        for (indsz_t j = 0; j < mCols; ++j) { *M.at(i, j) = row02row1translator[*M.at(i, j)]; }
    }

    indsz_t m0 = 0;
    for (const auto i : row02row1translator)
    {
        if (i == INDSZ_LIM) continue;
        m0 = std::max(i, m0);
    }

    auto rowcounter = vector<indsz_t>(m0 + 1, 0);
    // allow for zero values
    for (auto k : M.vals)
    {
        if (k < INDSZ_LIM) rowcounter[k] += 1;
    }

    auto colPtr1 = vector<indsz_t>(m0 + 2, 0);
    for (indsz_t i = 0; i < m0 + 1; ++i) colPtr1[i + 1] = colPtr1[i] + rowcounter[i];
    auto rowVals1 = vector<indsz_t>(colPtr1.back());
    auto placer   = vector<indsz_t>(colPtr1);

    for (indsz_t j = 0; j < mCols; ++j)
    {
        auto jval = col02col1translator[j];
        for (indsz_t i = 0; i < mRows; ++i)
        {
            auto row = *M.at(i, j);
            if (row < INDSZ_LIM)
            {
                rowVals1[placer[row]] = jval;
                placer[row] += 1;
            }
        }
    }

    return SparseBool(rowVals1, colPtr1, mCols);
}

std::tuple<SparseBool, std::vector<indsz_t>, std::vector<indsz_t>> filteredMatrixFromFarFaces(
    const eirene::Complex& farFaces, vector<vector<indsz_t>>& prePairs, vector<indsz_t>& lowBasis, uint64_t dim)
{
    /*
    Returns a tuple (M, lowLabels, highLabels)
    lowLabels and highLabels are reindexed cell indices - lowlab reindexed s.t. the filtration indices are
    in ascending order, i.e. cells born last appear first.

    M is a sparse matrix representation of the dim - 1 -> dim coboundary reindexed by lowlab/highlab.
    */

    const uint64_t numHigh = farFaces.splitCells[dim].rowVals.size();
    const uint64_t numLow  = farFaces.splitCells[dim - 1].rowVals.size();

    const vector<indsz_t>& ppHigh = prePairs[dim];
    const uint64_t numPrePairs    = ppHigh.size();

    vector<indsz_t> ppLow(numPrePairs);
    for (uint64_t i = 0; i < numPrePairs; ++i) ppLow[i] = farFaces.splitCells[dim].rowVals[ppHigh[i]];

    // highBasis is the set of pre-paired cells in the high dimension
    // the pairs are (high, low) = (prePairs[pairIndex], bdryRowVals[prePairs[pairIndex]])
    vector<indsz_t> highBasis(prePairs[dim + 1].size());
    const auto& highHighRowVs = farFaces.splitCells[dim + 1].rowVals;
    for (uint64_t i = 0; i < prePairs[dim + 1].size(); ++i) highBasis[i] = highHighRowVs[prePairs[dim + 1][i]];
    const auto highLabelsInPointOrder = uniqueIntsNotInUnsortedUpTo(highBasis, numHigh);

    // here, high pairs come from the pre-pairing in the VR construction - ppLow/ppHigh -
    // and low pairs come from the morse reduction pairing - through lowBasis
    const auto numHighUnpaired = numHigh - highBasis.size();
    const auto numLowUnpaired  = numLow - lowBasis.size();

    // note these could intersect, but we make it unique after
    lowBasis.insert(lowBasis.end(), ppLow.begin(), ppLow.end());
    highBasis.insert(highBasis.end(), ppHigh.begin(), ppHigh.end());
    vector<indsz_t> nonPairLowLab  = uniqueIntsNotInUnsortedUpTo(lowBasis, numLow);
    vector<indsz_t> nonPairHighLab = uniqueIntsNotInUnsortedUpTo(highBasis, numHigh);

    auto highTranslator = vector<indsz_t>(numHighUnpaired, 0);
    auto lowTranslator  = vector<indsz_t>(numLow, INDSZ_LIM);
    for (uint64_t i = 0; i < numPrePairs; ++i) lowTranslator[ppLow[i]] = i;

    // Compute the order the non-paired cells should appear in the filtration (after the paired cells)
    auto npFilt  = vector<indsz_t>();
    auto npOrder = vector<indsz_t>();
    if (dim > 0)
    {
        const auto& lowFilts = farFaces.splitFiltInds[dim - 1];

        npFilt.resize(nonPairLowLab.size());
        for (uint64_t i = 0; i < nonPairLowLab.size(); ++i) npFilt[i] = lowFilts[nonPairLowLab[i]];

        npOrder = integersInSameOrder(npFilt);
        for (auto& v : npOrder) v += numPrePairs;
    }

    // set lowtranslator for the non-paired cells
    for (uint64_t i = 0; i < nonPairLowLab.size(); ++i) lowTranslator[nonPairLowLab[i]] = npOrder[i];

    // lowLabels, highLabels, then are a reindexing of the cells that were unpaired in the last dimension (when it was
    // high)
    auto lowLabels = vector<indsz_t>(numLowUnpaired);

    // npOrder is a permutation of [numPrePairs:numPrePairs + nonPairedLow.size()]
    for (uint64_t i = 0; i < numPrePairs; ++i) lowLabels[i] = ppLow[i];
    for (uint64_t i = 0; i < npOrder.size(); ++i) lowLabels[npOrder[i]] = nonPairLowLab[i];

    // start with a copy of the pre-paired high cells, add the non-MORSE-paired labeling
    auto highLabels = std::vector<indsz_t>(ppHigh);
    highLabels.insert(highLabels.end(), nonPairHighLab.begin(), nonPairHighLab.end());

    auto ppSupp = vector<bool>(numHigh);
    for (uint64_t i = 0; i < ppHigh.size(); ++i) { ppSupp[ppHigh[i]] = true; }

    // set highTranslator to
    uint64_t ppMarker  = 0;
    uint64_t nppMarker = numPrePairs;
    for (uint64_t i = 0; i < numHighUnpaired; ++i)
    {
        auto hig = highLabelsInPointOrder[i];
        if (ppSupp[hig])
        {
            highTranslator[i] = ppMarker;
            ppMarker += 1;
        }
        else
        {
            highTranslator[i] = nppMarker;
            nppMarker += 1;
        }
    }

    auto allfaces = buildAllFromFar(farFaces, highLabelsInPointOrder, dim);
    auto M        = presparsefull2unsortedsparsetranspose(allfaces, lowTranslator, highTranslator);

    return std::tuple(M, lowLabels, highLabels);
}

vector<eirene::MorseReducedResult> performVRPersistence(eirene::Complex farFaces, vector<vector<indsz_t>>& prePairs)
{
    const indsz_t numDim = 3;   // numberOfDim(dimPatt);

    vector<eirene::MorseReducedResult> resultVec(numDim + 1);
    // vector<vector<indsz_t>> prePairs(numDim + 1, vector<indsz_t>{});

    // Calculating cohomology: start from low dimensions
    // Each iteration factorizes boundary from dim to (dim - 1)
    for (indsz_t dim = 1; dim < numDim; ++dim)
    {
        if (dim > farFaces.splitCells.size()) continue;

        // Get the high cell indices from the pairs calculated in the lower dimension
        vector<indsz_t> lowbasisnames(resultVec[dim - 1].pairVec.size());
        std::transform(resultVec[dim - 1].pairVec.begin(), resultVec[dim - 1].pairVec.end(), lowbasisnames.begin(),
                       [](const MorsePair& pair) { return pair.highInd; });

        auto [M, lowLab, highLab] = filteredMatrixFromFarFaces(farFaces, prePairs, lowbasisnames, dim);
        auto lowLabTemp           = vector<indsz_t>(lowLab.size());    // convert(Array{Int64,1},1:length(lowlab))
        auto highLabTemp          = vector<indsz_t>(highLab.size());   // (Array{Int64,1},1:length(higlab))
        std::iota(lowLabTemp.begin(), lowLabTemp.end(), 0);
        std::iota(highLabTemp.begin(), highLabTemp.end(), 0);

        auto lowFiltTemp  = vector<indsz_t>(lowLab.size());    // convert(Array{Int64,1},1:length(lowlab))
        auto highFiltTemp = vector<indsz_t>(highLab.size());   // (Array{Int64,1},1:length(higlab))
        for (indsz_t i = 0; i < lowLab.size(); ++i) lowFiltTemp[i] = farFaces.splitFiltInds[dim - 1][lowLab[i]];
        for (indsz_t i = 0; i < highLab.size(); ++i) highFiltTemp[i] = farFaces.splitFiltInds[dim][highLab[i]];

        // pplow = convert(Array,length(prepairs[sd]):-1:1)
        // pphig = convert(Array,length(prepairs[sd]):-1:1)
        auto ppLow  = vector<indsz_t>(prePairs[dim].size());
        auto ppHigh = vector<indsz_t>(prePairs[dim].size());

        std::iota(ppLow.rbegin(), ppLow.rend(), 0);
        std::iota(ppHigh.rbegin(), ppHigh.rend(), 0);

        const auto numPP = prePairs[dim].size();
        auto pairs       = eirene::PairVec(numPP);

        for (indsz_t i = 0; i < numPP; ++i) { pairs[i] = MorsePair{.lowInd = numPP - i - 1, .highInd = numPP - i - 1}; }

        auto [leftFactor, excisedPairs, tlabels] = spla::morseLU(M, pairs, highFiltTemp, lowFiltTemp);

        std::transform(excisedPairs.begin(), excisedPairs.end(), excisedPairs.begin(),
                       [&lowLab, &highLab](const MorsePair& pair) {
                           return MorsePair{lowLab[pair.lowInd], highLab[pair.highInd]};
                       });

        leftFactor.rowVals = perm::reindex(leftFactor.rowVals, lowLab);
        resultVec[dim]     = eirene::MorseReducedResult{leftFactor,                       // trv and tcp
                                                    perm::reindex(tlabels, lowLab),   // tid
                                                    excisedPairs,                     // plo/hi
                                                    vector<vector<indsz_t>>(0)};
    }

    // extractHomology(resultVec, splitComplex, dimPatt, sortedSplitFilt, comp.uniqFiltVals);

    return resultVec;
}

}   // namespace vrips
