#pragma once
#include <vector>
#include <cstdint>
#include <string>

#include "indsz.hpp"

namespace spla
{

struct SparseBool
{
    std::vector<indsz_t> rowVals;
    std::vector<indsz_t> colPtrs;
    size_t numRows;

    SparseBool() : SparseBool(0, 0) {}

    SparseBool(size_t numRows, size_t numCols) : numRows{numRows} { colPtrs.resize(numCols + 1); }

    SparseBool(std::vector<indsz_t> rowVals, std::vector<indsz_t> colPtrs, size_t numRows)
        : rowVals{rowVals}, colPtrs{colPtrs}, numRows{numRows}
    {
    }

    SparseBool(const SparseBool& other) : rowVals{other.rowVals}, colPtrs{other.colPtrs}, numRows{other.numRows} {}

    SparseBool& operator=(SparseBool other)
    {
        // Use vector::swap instead, moves mem, more efficient.
        std::swap(rowVals, other.rowVals);
        std::swap(colPtrs, other.colPtrs);
        std::swap(numRows, other.numRows);
        return *this;
    }

    bool at(indsz_t i, indsz_t j) const
    {
        for (indsz_t sInd = colPtrs[j]; sInd < colPtrs[j + 1]; ++sInd)
        {
            if (rowVals[sInd] > i) return false;
            if (rowVals[sInd] == i) return true;
        }

        return false;
    }

    indsz_t numCols() const { return colPtrs.size() - 1; }

    indsz_t numNonZero() const { return rowVals.size(); }

    void emplace_at(const bool toPlace, const indsz_t rowInd, const indsz_t colInd)
    {
        // Stitched in, not tested
        const indsz_t startInd = colPtrs[colInd], stopInd = colPtrs[colInd + 1];

        // If the column in nonempty, see if the location is already nonzero
        for (indsz_t i = startInd; i < stopInd && rowVals[i] <= rowInd; ++i)
        {
            // if it is, either erase it or overwrite it, depending
            if (rowVals[i] == rowInd)
            {
                if (!toPlace)
                {
                    rowVals.erase(rowVals.begin() + i);

                    for (i = colInd + 1; i < colPtrs.size(); ++i) this->colPtrs[i] -= 1;

                    return;
                }

                return;
            }
        }

        // if we're here, the location is zero
        if (!toPlace) return;

        // keep rowInds in increasing order for each column
        indsz_t insertInd = startInd;
        while (insertInd < stopInd && rowVals[insertInd] < rowInd) insertInd++;

        rowVals.insert(rowVals.begin() + insertInd, rowInd);

        for (indsz_t i = colInd + 1; i < colPtrs.size(); ++i) colPtrs[i]++;
    }
};

}   // namespace spla

inline bool sparseBoolAt(const spla::SparseBool& sb, indsz_t i, indsz_t j)
{
    for (indsz_t sInd = sb.colPtrs[j]; sInd < sb.colPtrs[j + 1]; ++sInd)
    {
        if (sb.rowVals[sInd] > i) return false;
        if (sb.rowVals[sInd] == i) return true;
    }

    return false;
}
inline std::string printSparseBool(const spla::SparseBool& obj)
{
    const indsz_t numRows = obj.numRows, numCols = obj.numCols();
    indsz_t toWrite = 0;

    std::string toRet{"\n"};

    for (indsz_t irow = 0; irow < numRows; ++irow)
    {
        for (indsz_t icol = 0; icol < numCols; ++icol)
        {
            indsz_t colStart = obj.colPtrs[icol], colEnd = obj.colPtrs[icol + 1];
            toWrite = 0;

            while (colStart < colEnd)
            {
                if (obj.rowVals[colStart] == irow)
                {
                    toWrite = 1;
                    break;
                }

                ++colStart;
            }

            toRet.push_back(toWrite ? '1' : '0');
        }
        toRet.push_back('\n');
    }

    toRet.push_back('\n');
    toRet.push_back('\n');
    return toRet;
}
