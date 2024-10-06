
#include "util.hpp"

std::vector<indsz_t> perm::reindex(const std::vector<indsz_t>& permVec, const std::vector<indsz_t>& setVec)
{
    /*
     *  Applies the permutation/relabeling in the first argument to the "set" of indices in the second.
     *  Formally, consider the indices in setVec as a set S, and define a function
     *  Ψ : S -> S by Ψ(s) = permVec[s]. Then this returns the image of S under Ψ.
     */
    // indsz_t outSz = std::min(permVec.size(), setVec.size());
    indsz_t outSz = permVec.size();

    std::vector<indsz_t> result(outSz);

    for (indsz_t i = 0; i < outSz; ++i) result[i] = setVec[permVec[i]];

    return result;
}

void perm::reindexInplace(std::vector<indsz_t>& perm, std::vector<indsz_t>& vals)
{
    /*
     *   Applies the permutation specified by perm to vals, in-place.
     *   Overwrites both. No allocations. Faster?
     */
    for (indsz_t i = 0; i < perm.size(); ++i)
    {
        indsz_t holeInd = i;

        while (i != perm[holeInd])
        {
            indsz_t next = perm[holeInd];

            std::swap(vals[holeInd], vals[next]);
            perm[holeInd] = holeInd;
            holeInd       = next;
        }

        perm[holeInd] = holeInd;
    }
}


indsz_t numberOfDim(const std::vector<indsz_t>& dimPatt)
{
    /*
     * If the complex has maximum cell dimension of n, returns n + 1.
     */
    if (dimPatt.size() == 0) { fmt::print(stderr, "Error! dimPatt should never have size 0!"); }

    return dimPatt.size();
}

std::vector<indsz_t> uniqueIntsNotInUnsortedUpTo(const std::vector<indsz_t>& uniqIntVec, indsz_t endPt)
{
    /*
     *   Returns the integers from the interval [0, endPt) NOT present in
     *   uniqIntVec[0 : subVecEnd>0 ? subVecEnd : uniqIntVec.size()].
     *
     *   Note that as suggested, uniqIntVec may be unsorted std::vector of unique nonnnegative integers.
     *   This only works if endPt >= uniqIntVec.size(), and (x ∈ uniqIntVec) ⇒ (x < endPt)
     */

    const indsz_t vecSz = uniqIntVec.size();
    std::vector<indsz_t> complement{};

    // If vecSz is 0, act like every integer in [0, endPt) is not present,
    // so the complement is the whole interval
    if (vecSz == 0)
    {
        complement.resize(endPt);
        std::iota(complement.begin(), complement.end(), 0);
        return complement;
    }

    // Else, if vecSz == endPt, act like evry integer is present, so
    // the complement is nothing
    else if (vecSz == endPt) { return complement; }

    else if (endPt < vecSz)
    {
        // will error below
        fmt::print(stderr, "Error! endPt < size of passed vector in uniqueIntsNotInUnsortedUpTo!\n");
    }

    std::vector<bool> complementSupport(endPt, true);
    for (indsz_t i = 0; i < vecSz; ++i) complementSupport[uniqIntVec[i]] = false;

    complement.resize(endPt - vecSz);
    indsz_t marker = 0;

    for (indsz_t i = 0; i < endPt; ++i)
    {
        if (complementSupport[i])
        {
            complement[marker] = i;
            marker++;
        }
    }

    if (marker != (endPt - vecSz))
    {
        fmt::print(stderr, "Warning! uniqueIntsNotInUnsortedUpTo is not working as expected...");
        std::getchar();
    }

    return complement;
}