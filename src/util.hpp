#pragma once 

#include <vector>
#include <numeric>
#include "eirene.hpp"

namespace perm {

std::vector<indsz_t> reindex(const std::vector<indsz_t>& permVec, const std::vector<indsz_t>& setVec);
void reindexInplace(std::vector<indsz_t>& perm, std::vector<indsz_t>& vals);

}   // namespace perm

indsz_t numberOfDim(const std::vector<indsz_t>& dimPatt);
std::vector<indsz_t> uniqueIntsNotInUnsortedUpTo(const std::vector<indsz_t>& uniqIntVec, indsz_t endPt);

template <typename flt_t>
std::vector<indsz_t> getSortPerm(const std::vector<flt_t>& fltVec, const bool rev = false)
{
    /*
     * Returns the permutation that sortes fltVec, if rev then for reverse order.
     * In the use case for eirene, fltVec is guaranteed to not have duplicates.
     */
    auto temp = std::vector<std::tuple<flt_t, indsz_t>>(fltVec.size());

    for (indsz_t i = 0; i < fltVec.size(); ++i) temp[i] = std::make_tuple(fltVec[i], i);

    if (rev)
        std::sort(temp.rbegin(), temp.rend());   // reverse sort
    else
        std::sort(temp.begin(), temp.end());

    auto toRet = std::vector<indsz_t>(fltVec.size());
    std::transform(temp.begin(), temp.end(), toRet.begin(), [](auto tup) { return std::get<indsz_t>(tup); });

    return toRet;
}
