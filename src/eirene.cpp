#include "eirene.hpp"
#include "eirene_impl.hpp"

using std::tuple;
using std::vector;
using namespace eirene;

/*
 *   Functions that process the inital complex data.
 */

namespace
{
inline std::vector<indsz_t> computeDimensionPattern(const std::vector<indsz_t>& dimVec)
{
    /*
        Transforms the dimension vector into the dimension pattern:
        Dimension vector has number of cells with dimension i in entry i
        Dimension pattern has (# of cells with dimension < i) in entry i
        i.e. euler vector's prefix sums
    */
    auto dimPatt = std::vector<indsz_t>(dimVec.size() + 1, 0);
    dimPatt[0]   = 0;
    for (indsz_t i = 0; i < dimVec.size(); ++i) dimPatt[i + 1] = dimPatt[i] + dimVec[i];
    return dimPatt;
}

spla::SparseBool copySubMatrix(const spla::SparseBool A, indsz_t firstCol, indsz_t lastCol)
{
    /*
        Returns the submatrix consisting of columns [firstCol, lastCol)

        Obviously, lastCol >= firstCol.
    */
    const size_t numVal  = A.colPtrs[lastCol] - A.colPtrs[firstCol];
    const indsz_t numCol = lastCol - firstCol;

    auto rowCopies = vector<indsz_t>(numVal, 0);
    auto colCopies = vector<indsz_t>(lastCol - firstCol + 1, 0);

    for (indsz_t colInd = 0; colInd < numCol; ++colInd)
    {
        indsz_t origInd        = firstCol + colInd;
        indsz_t elemBeforeNext = A.colPtrs[origInd + 1];
        indsz_t elemBeforeCurr = A.colPtrs[origInd];

        colCopies[colInd + 1] = colCopies[colInd] + elemBeforeNext - elemBeforeCurr;

        // pointer arithmetic
        std::copy(A.rowVals.data() + elemBeforeCurr, A.rowVals.data() + elemBeforeNext,
                  rowCopies.data() + colCopies[colInd]);
    }

    return spla::SparseBool{rowCopies, colCopies, A.numRows};
}

auto checkFiltrationMonotone(const std::vector<spla::SparseBool>& splitComp,
                             const std::vector<std::vector<indsz_t>>& sortSplitFilt,
                             const std::vector<indsz_t>& dimPatt)

{
    bool ret = true;

    for (indsz_t dim = 0; dim < numberOfDim(dimPatt); ++dim)
    {
        const auto& dimAdj = splitComp[dim];

        for (indsz_t cell = 0; cell < dimAdj.numCols(); ++cell)
        {
            const indsz_t cellFiltVal = sortSplitFilt[dim][cell];

            for (indsz_t faceInd = dimAdj.colPtrs[cell]; faceInd < dimAdj.colPtrs[cell + 1]; ++faceInd)
            {
                indsz_t faceCell = dimAdj.rowVals[faceInd];

                // >= because they're reverse sorted, see the splitting function
                const bool flag = sortSplitFilt[dim - 1][faceCell] >= cellFiltVal;
                if (!flag)
                {
                    fmt::print(stderr,
                               "Error checking monotone, cell {} of dim {} (global index {}) has filtration value {} "
                               "but is face of cell {} of dim {} "
                               "(global index {}) with filtration value {}.\n",
                               faceCell, dim - 1, faceCell + dimPatt[dim - 1], sortSplitFilt[dim - 1][faceCell], cell,
                               dim, cell + dimPatt[dim], cellFiltVal);
                }
                ret &= flag;
            }

            if (!ret) return ret;
        }
    }

    return ret;
}

int splitByDim(const spla::SparseBool& comp, const vector<indsz_t>& dimPatt, vector<spla::SparseBool>& splitComplex,
               const bool check = true)
{
    /*
        Splits the information up by dimension.
        Checks that bdry operator is well graded if check = true;
    */

    // Largest dimension to compute
    const indsz_t topDim = numberOfDim(dimPatt);
    splitComplex.resize(topDim);

    for (indsz_t currDim = 0; currDim < topDim; ++currDim)
    {
        // The columns corresponding to cells with dimension currDim are those with indexes in [dimPatt[currDim] - 1,
        // dimPatt[currDim + 1] - 1)
        // topDim is greatest lower bound on 0 chain groups... (bdry_i is 0 for i >= topDim - 1)
        if (currDim < topDim - 1) splitComplex[currDim] = copySubMatrix(comp, dimPatt[currDim], dimPatt[currDim + 1]);

        // transform the row indexing for each dimension to be valid relative, i.e. start from 0 for each
        const indsz_t smallerCells        = (currDim == 0) ? 0 : dimPatt[currDim - 1],
                      smallerOrEqualCells = (currDim == topDim) ? dimPatt[currDim - 1] : dimPatt[currDim];

        for (indsz_t& rowVal : splitComplex[currDim].rowVals)
        {
            if (check && (rowVal < smallerCells || rowVal > smallerOrEqualCells))
            {
                fmt::print(stderr, "Chain boundary does not appear to be graded of dimension 1");
                return -1;
            }

            rowVal -= smallerCells;
        }

        splitComplex[currDim].numRows = smallerOrEqualCells - smallerCells;
    }

    return 0;
}

tuple<vector<vector<indsz_t>>, vector<double>> sortAndSplitFilts(const vector<double>& filtVals,
                                                                 const vector<indsz_t>& dimPatt)
{
    /*
        Returns a reverse sorted list of unique filtration values, and a vector of vector of indices
        into the (unique) filtration values, where we split the cells by dimension as we did with the complex
       information.
    */

    vector<vector<indsz_t>> splitIndices(numberOfDim(dimPatt));
    vector<double> uniqFilts(filtVals);

    // reverse sort, then make unique
    std::sort(uniqFilts.rbegin(), uniqFilts.rend());
    uniqFilts.erase(std::unique(uniqFilts.begin(), uniqFilts.end()), uniqFilts.end());

    // For each dim, fill a vector with the indices of the filtration value for each cell in the
    // reverse sorted list.
    indsz_t totalInd = 0;

    for (indsz_t dim = 0; dim < numberOfDim(dimPatt); ++dim)
    {
        bool overDim     = dim + 1 >= dimPatt.size();
        indsz_t numCells = overDim ? 0 : (dimPatt[dim + 1] - dimPatt[dim]);

        splitIndices[dim].resize(numCells);

        for (indsz_t cellInd = 0; cellInd < numCells; ++cellInd)
        {
            // Not using sorted property
            splitIndices[dim][cellInd] =
                std::distance(uniqFilts.begin(), std::find(uniqFilts.begin(), uniqFilts.end(), filtVals[totalInd]));
            totalInd++;
        }
    }

    return std::make_tuple(splitIndices, uniqFilts);
}

tuple<vector<indsz_t>, vector<double>> distMatToIndices(double* distMat, size_t numPoints)
{
    /*
     * Similar to sortAndSplitFilt, given a distance matrix return a list of unique distances in it and a matrix
     * indexing into that list. We don't reverse sort in this (VR) case
     */

    vector<indsz_t> indexMat(numPoints * numPoints, 0);
    vector<double> uniqFilts(distMat, distMat + numPoints * numPoints);

    // sort, then make unique
    std::sort(uniqFilts.rbegin(), uniqFilts.rend());
    uniqFilts.erase(std::unique(uniqFilts.begin(), uniqFilts.end()), uniqFilts.end());

    // For each dim, fill a vector with the indices of the filtration value for each cell in the
    // reverse sorted list.
    for (size_t pointInd = 0; pointInd < numPoints * numPoints; ++pointInd)
    {
        // Not using sorted property
        const auto ind =
            std::distance(uniqFilts.begin(), std::find(uniqFilts.begin(), uniqFilts.end(), distMat[pointInd]));
        indexMat[pointInd] = ind;
    }

    return std::make_tuple(indexMat, uniqFilts);
}
}   // namespace

int eirene::toComplex(const spla::SparseBool& adj, const std::vector<double>& filts, const std::vector<indsz_t>& dimVec,
                      eirene::Complex& comp, const bool check)
{
    /*
     * This function just splits the info by dimension and checks that it actually defines a complex.
     */

    if (dimVec.size() == 0)
    {
        fmt::print(stderr, "Eirene was passed a euler vector with size 0.");
        return -1;
    }

    auto dimPatt = computeDimensionPattern(dimVec);

    if (filts.size() != dimPatt.back())
    {
        fmt::print(stderr,
                   "The filtration vector should have a double value for every cell, but there "
                   "are {} vals in fv and {} cells according to dimPatt.\n",
                   filts.size(), dimPatt.back());
        return -1;
    }

    // SortedSplitFilt: vector of vector of indices, each represents the filt Val of a
    // cell - so sortedSplitFilt[i][j] represents the jth cell with dimension i.
    // The index is into uniqFilt, which is a reverse sorted list of unique filtration values.
    auto [sortSplitFilt, uniqFilt] = sortAndSplitFilts(filts, dimPatt);

    // split the complex by cell dimension: in original eirene, called 'segmenting'
    auto splitComp = std::vector<spla::SparseBool>();
    auto result    = splitByDim(adj, dimPatt, splitComp, check);
    if (result < 0) { return -1; }

    if (check && !checkFiltrationMonotone(splitComp, sortSplitFilt, dimPatt))
    {
        fmt::print(stderr, "Filtration is not monotone on the cell complex poset.\n");
        return -1;
    }
    comp = eirene::Complex{splitComp, dimPatt, sortSplitFilt, uniqFilt};
    return 0;
}

inline std::vector<uint64_t> subsetDimFilt(const std::vector<uint64_t> subsetInds, const std::vector<uint64_t> dimPatt)
{
    // Return the first element of subsetInds in each dimension
    // SubsetInds is assumed to be sorted
    std::vector<uint64_t> result;

    uint64_t subInd = 0;
    for (uint64_t dim = 0; dim < dimPatt.size(); ++dim)
    {
        while (subInd < subsetInds.size() && subsetInds[subInd] < dimPatt[dim]) subInd++;
        result.push_back(subInd);
    }

    return result;
}

std::vector<uint64_t> eirene::openStar(const std::vector<uint64_t>& subsetInds, const eirene::Complex& complex)
{
    /*
     * It's just the upward closure.
     * Assumed subsetInds is sorted.
     */

    auto cellsByDim                  = std::vector<std::vector<uint64_t>>(complex.dimPatt.size());
    std::vector<uint64_t> subDimFilt = subsetDimFilt(subsetInds, complex.dimPatt);

    const auto& dimPatt = complex.dimPatt;

    for (uint64_t dim = 0; dim < dimPatt.size() - 1; ++dim)
    {

        for (uint64_t subInd = subDimFilt[dim]; subInd < subDimFilt[dim + 1]; ++subInd)
        {
            cellsByDim[dim].push_back(subsetInds[subInd]);
        }

        if (dim == dimPatt.size() - 1) continue;
        const spla::SparseBool& bdry = complex.splitCells[dim + 1];
        for (const uint64_t cell : cellsByDim[dim])
        {
            const uint64_t relCellInd = cell - dimPatt[dim];

            for (uint64_t highInd = 0; highInd < bdry.numCols(); ++highInd)
                // the columns of the bdry represents the faces of the cells
                for (uint64_t faceInd = bdry.colPtrs[highInd]; faceInd < bdry.colPtrs[highInd + 1]; ++faceInd)
                    if (bdry.rowVals[faceInd] == relCellInd) cellsByDim[dim + 1].push_back(dimPatt[dim + 1] + highInd);
        }
    }

    std::vector<uint64_t> result;
    for (auto& dimCells : cellsByDim)
    {
        std::sort(dimCells.begin(), dimCells.end());
        dimCells.erase(std::unique(dimCells.begin(), dimCells.end()), dimCells.end());
        result.insert(result.end(), dimCells.begin(), dimCells.end());
    }

    return result;
}

std::vector<uint64_t> eirene::complement(const std::vector<uint64_t>& subset, const eirene::Complex& complex)
{
    const auto dimPatt = complex.dimPatt;

    auto result     = std::vector<uint64_t>();
    uint64_t subPtr = 0;
    for (uint64_t dim = 0; dim < complex.dimPatt.size() - 1; ++dim)
    {
        for (uint64_t cell = dimPatt[dim]; cell < dimPatt[dim + 1]; ++cell)
        {
            while (subPtr < subset.size() && subset[subPtr] <= cell)
            {
                // if the cell is in subset, dont add it
                if (subset[subPtr] == cell) { goto cell_top; }
                subPtr++;
            }

            result.push_back(cell);
        cell_top:
            continue;
        }
    }

    return result;
}

using std::vector;
using namespace eirene;

std::vector<eirene::MorseReducedResult> eirene::relativeCellularPersistence(const eirene::Complex& comp,
                                                                             const std::vector<uint64_t>& openInds)
{
    /*
     * Same as cellularPersistence, but relative to the open subset defined by openInds.
     * As if we contracted X/O.
     */

    const vector<spla::SparseBool>& splitComplex         = comp.splitCells;
    const vector<indsz_t>& dimPatt                 = comp.dimPatt;
    const vector<vector<indsz_t>>& sortedSplitFilt = comp.splitFiltInds;

    const indsz_t numDim = numberOfDim(dimPatt);

    vector<MorseReducedResult> resultVec(numDim + 1);
    vector<vector<indsz_t>> prePairs(numDim + 1, vector<indsz_t>{});

    const auto closedComplement = complement(openInds, comp);
    const auto cComplDimFilt    = subsetDimFilt(closedComplement, dimPatt);

    // Calculating cohomology: start from low dimensions
    // Each iteration factorizes boundary from dim to (dim - 1)
    for (indsz_t dim = 1; dim < numDim; ++dim)
    {
        // get the low (high below) cell indices in the complement
        auto lowCompLabels = vector<indsz_t>(closedComplement.begin() + cComplDimFilt[dim - 1],
                                             closedComplement.begin() + cComplDimFilt[dim]);
        for (auto& label : lowCompLabels) label -= dimPatt[dim - 1];

        // Get the high cell indices from the pairs calculated in the lower dimension
        vector<indsz_t> lowbasisnames(resultVec[dim - 1].pairVec.size());
        std::transform(resultVec[dim - 1].pairVec.begin(), resultVec[dim - 1].pairVec.end(), lowbasisnames.begin(),
                       [](const MorsePair& pair) { return pair.highInd; });

        // Add the lowCompLabels to the lowbasisnames
        lowbasisnames.insert(lowbasisnames.end(), lowCompLabels.begin(), lowCompLabels.end());

        const spla::SparseBool& dimBoundary = splitComplex[dim];
        indsz_t numHighCells = dimBoundary.numCols(), numLowCells = splitComplex[dim - 1].numCols();

        // lowlabels will hold the indices of cells of dimension (dim - 1) which do not appear in
        // lowbasisnames, i.e. were unpaired in any previous iterations (and not in the complement).
        vector<indsz_t> lowlabels = uniqueIntsNotInUnsortedUpTo(lowbasisnames, numLowCells);

        // also, low labels will be in permuted order
        // bigger index means the cell is born earlier
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

        auto highCompLabels = vector<indsz_t>(closedComplement.begin() + cComplDimFilt[dim],
                                              closedComplement.begin() + cComplDimFilt[dim + 1]);

        for (auto& label : highCompLabels) label -= dimPatt[dim];

        vector<indsz_t> highlabels = uniqueIntsNotInUnsortedUpTo(highCompLabels, numHighCells);

        // as above but for high cells
        subFilt = vector<double>(highlabels.size());
        std::transform(highlabels.begin(), highlabels.end(), subFilt.begin(),
                       [&](indsz_t ind) { return sortedSplitFilt[dim][ind]; });

        vector<indsz_t> highCellFiltPerm = getSortPerm(subFilt);
        perm::reindexInplace(highCellFiltPerm, highlabels);

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
        resultVec[dim]     = eirene::MorseReducedResult{leftFactor,                          // trv and tcp
                                                    perm::reindex(tlabels, lowlabels),   // tid
                                                    excisedPairs,                        // plo/hi
                                                    vector<vector<indsz_t>>(0)};
    }

    const auto complement = eirene::Subcomplex { .cellInds = closedComplement, .dimFilt = cComplDimFilt }; 
    eirene::extractHomology(comp, resultVec, complement);

    return resultVec;
}

PyreneResult eirene::pyRelativeCellularPersistence(indsz_t* rowVals, indsz_t numVals, indsz_t* colPtrs, indsz_t numCols, indsz_t numRows, 
                                            double* filtVals, size_t filtLen, indsz_t* dimVec, size_t numDim, indsz_t* openCompPtr, indsz_t compLen) {

    const CSparseBool cSparse = CSparseBool {
        .colPtrs = U32Arr {
            .data = colPtrs,
            .len = numCols + 1,
        },
        .rowVals = U32Arr {
            .data = rowVals,
            .len = numVals,
        },
        .numRows = numRows,
    };

    const spla::SparseBool cComp = sparseBoolConvertInv(&cSparse);
    const auto cFltVals          = std::vector<double>(filtVals, filtVals + filtLen);
    const auto cDimVec           = std::vector<indsz_t>(dimVec, dimVec + numDim);
    const auto subset          = std::vector<indsz_t>(openCompPtr, openCompPtr + compLen);

    eirene::Complex complex;
    eirene::toComplex(cComp, cFltVals, cDimVec, complex, true);

    const auto openComp = eirene::openStar(subset, complex);

    const vector<MorseReducedResult> res = eirene::relativeCellularPersistence(complex, openComp);
    auto pyResult = PyreneResult();
  
    for (const auto &dimRes : res) {
        pyResult.tid.push_back(dimRes.tid);
        pyResult.cycleReps.push_back(dimRes.cycleReps);
        pyResult.birthDeath.push_back(dimRes.birthDeath);
    }

    return pyResult;
}

vector<MorseReducedResult> eirene::vripsPersistence(double* dMat, indsz_t numPoints, double maxRad)
{
    auto colMaxs = std::vector<double>(numPoints, 0.0);
    for (uint64_t colInd = 0; colInd < numPoints; ++colInd)
    {
        colMaxs[colInd] = *std::max_element(dMat + numPoints * colInd, dMat + numPoints * (colInd + 1));
    }

    maxRad        = std::min(maxRad, *std::min_element(colMaxs.begin(), colMaxs.end()));
    auto distCopy = std::vector<double>();
    distCopy.reserve(numPoints * numPoints);
    std::copy(dMat, dMat + numPoints * numPoints, std::back_inserter(distCopy));

    for (uint64_t i = 0; i < numPoints * numPoints; ++i)
    {
        if (distCopy[i] > maxRad) distCopy[i] = INFINITY;
    }

    // uniqFilt is ocg2rad in eirene.jl
    auto [indMat, uniqFilt] = distMatToIndices(distCopy.data(), numPoints);

    auto [vrComplex, nvl2ovl, prePairs] = vrips::distMatToVRComplex(indMat, numPoints, 3);
    vrComplex.uniqFiltVals              = uniqFilt;

    return vrips::performVRPersistence(vrComplex, prePairs);
}

PyreneResult eirene::pyVripsPersistence(double* dMat, indsz_t numPoints, double maxRad) {
    const auto res = vripsPersistence(dMat, numPoints, maxRad);
    auto pyResult = PyreneResult();
  
    for (const auto &dimRes : res) {
        pyResult.tid.push_back(dimRes.tid);
        pyResult.cycleReps.push_back(dimRes.cycleReps);
    }

    return pyResult;
}

vector<MorseReducedResult> eirene::cellularPersistence(const eirene::Complex& comp)
{
    auto resultVec = eirene::performPersistence(comp);

    return resultVec;
}

PyreneResult eirene::pyCellularPersistence(indsz_t* rowVals, indsz_t numVals, indsz_t* colPtrs, indsz_t numCols, indsz_t numRows, 
                                            double* filtVals, size_t filtLen, indsz_t* dimVec, size_t numDim) {
    fmt::print(stderr, "In C++\n");
    const CSparseBool cSparse = CSparseBool {
        .colPtrs = U32Arr {
            .data = colPtrs,
            .len = numCols + 1,
        },
        .rowVals = U32Arr {
            .data = rowVals,
            .len = numVals,
        },
        .numRows = numRows,
    };

    const spla::SparseBool cComp = sparseBoolConvertInv(&cSparse);
    const auto cFltVals          = std::vector<double>(filtVals, filtVals + filtLen);
    const auto cDimVec           = std::vector<indsz_t>(dimVec, dimVec + numDim);

    fmt::print(stderr, "Making complex\n");
    eirene::Complex complex;
    eirene::toComplex(cComp, cFltVals, cDimVec, complex, true);

    fmt::print(stderr, "Performing persistence.\n");
    const vector<MorseReducedResult> res = eirene::performPersistence(complex);
    auto pyResult = PyreneResult();
  
    fmt::print(stderr, "Building pyresult\n");
    for (const auto &dimRes : res) {
        pyResult.tid.push_back(dimRes.tid);
        pyResult.cycleReps.push_back(dimRes.cycleReps);
        pyResult.birthDeath.push_back(dimRes.birthDeath);
    }

    fmt::print(stderr, "Returning\n");
    return pyResult;
}