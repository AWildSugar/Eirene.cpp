#pragma once
#include "eirene.hpp"

namespace spla
{
spla::SparseBool transposeSparseSubmatrix(const spla::SparseBool& sparseMat, const std::vector<indsz_t>& rows,
                                          const std::vector<indsz_t>& cols);

std::tuple<spla::SparseBool, eirene::PairVec, std::vector<indsz_t>> morseLU(spla::SparseBool& coBdry,
                                                                            eirene::PairVec& pairVec,
                                                                            const std::vector<indsz_t>& highFiltTemp,
                                                                            const std::vector<indsz_t>& lowFiltTemp);

}   // namespace spla
namespace eirene
{

std::vector<eirene::MorseReducedResult> performPersistence(const eirene::Complex& comp);

struct Subcomplex
{
    std::vector<indsz_t> cellInds;
    std::vector<indsz_t> dimFilt;
};
void extractHomology(const Complex& complex, std::vector<eirene::MorseReducedResult>& reducedVec,
                     const std::optional<eirene::Subcomplex>& closedCompO);

}   // namespace eirene

namespace vrips
{
std::tuple<eirene::Complex, std::vector<indsz_t>, std::vector<std::vector<indsz_t>>> distMatToVRComplex(
    const std::vector<indsz_t>& indMat, indsz_t numPoints, indsz_t maxSimplexDim);

std::vector<eirene::MorseReducedResult> performVRPersistence(eirene::Complex farFaces,
                                                             std::vector<std::vector<indsz_t>>& prePairs);
}   // namespace vrips