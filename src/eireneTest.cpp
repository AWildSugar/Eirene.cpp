#include "eirene.hpp"
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>

using std::vector;

auto testLine()
{
    spla::SparseBool testComplex(3, 3);
    testComplex.emplace_at(true, 0, 2);
    testComplex.emplace_at(true, 1, 2);

    printSparseBool(testComplex);

    std::vector<double> testFilt{0.01, 0.02, 0.03};
    std::vector<indsz_t> eulerVec{2, 1};

    eirene::Complex comp;
    eirene::toComplex(testComplex, testFilt, eulerVec, comp, true);

    return eirene::cellularPersistence(comp);
}

auto testTriangle()
{
    spla::SparseBool testComplex(6, 6);

    testComplex.emplace_at(true, 0, 3);
    testComplex.emplace_at(true, 1, 3);

    testComplex.emplace_at(true, 1, 4);
    testComplex.emplace_at(true, 2, 4);

    testComplex.emplace_at(true, 0, 5);
    testComplex.emplace_at(true, 2, 5);

    // std::vector<double> testFilt{0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
    std::vector<double> testFilt(6, 0.01);
    std::vector<indsz_t> eulerVec{3, 3};

    //    std::cerr << testComplex << std::endl;

    eirene::Complex comp;
    eirene::toComplex(testComplex, testFilt, eulerVec, comp, true);

    return eirene::cellularPersistence(comp);
}

auto testSphere()
{
    spla::SparseBool testComplex(6, 6);

    testComplex.emplace_at(true, 0, 2);
    testComplex.emplace_at(true, 1, 2);

    testComplex.emplace_at(true, 0, 3);
    testComplex.emplace_at(true, 1, 3);

    testComplex.emplace_at(true, 2, 4);
    testComplex.emplace_at(true, 3, 4);

    testComplex.emplace_at(true, 2, 5);
    testComplex.emplace_at(true, 3, 5);

    printSparseBool(testComplex);

    auto testFilt = std::vector<double>{0.01, 0.02, 0.03, 0.04, 0.05, 0.06};
    auto eulerVec = std::vector<indsz_t>{2, 2, 2};

    eirene::Complex comp;
    eirene::toComplex(testComplex, testFilt, eulerVec, comp, true);

    return eirene::cellularPersistence(comp);
}

auto lingareddyComplex()
{
    /*
     * Random complex, no filtration, from https://math.uchicago.edu/~may/REU2018/REUPapers/Lingareddy.pdf.
     */
    spla::SparseBool testComplex(13, 13);

    testComplex.emplace_at(true, 0, 5);
    testComplex.emplace_at(true, 0, 6);
    testComplex.emplace_at(true, 0, 7);

    testComplex.emplace_at(true, 1, 5);
    testComplex.emplace_at(true, 1, 8);

    testComplex.emplace_at(true, 2, 6);
    testComplex.emplace_at(true, 2, 9);

    testComplex.emplace_at(true, 3, 10);
    testComplex.emplace_at(true, 3, 9);
    testComplex.emplace_at(true, 3, 8);
    testComplex.emplace_at(true, 3, 7);
    
    testComplex.emplace_at(true, 4, 10);

    testComplex.emplace_at(true, 5, 11);
    testComplex.emplace_at(true, 7, 11);
    testComplex.emplace_at(true, 8, 11);

    testComplex.emplace_at(true, 6, 12);
    testComplex.emplace_at(true, 7, 12);
    testComplex.emplace_at(true, 9, 12);

    std::cerr << printSparseBool(testComplex) << std::endl;

    auto testFilt = std::vector<double>(13, 0.01);
    testFilt[1] = 0.02;
    testFilt[5] = 0.02;
    testFilt[8] = 0.02;
    testFilt[11] = 0.03;
    auto eulerVec = std::vector<indsz_t>{5, 6, 2};

    eirene::Complex comp;
    eirene::toComplex(testComplex, testFilt, eulerVec, comp, true);
    return comp;
}

auto barbellComplexRelativeHom()
{
    /*
     * Fold da barbell.
     */
    spla::SparseBool testComplex(5, 5);

    testComplex.emplace_at(true, 0, 3);

    testComplex.emplace_at(true, 1, 3);
    testComplex.emplace_at(true, 1, 4);

    testComplex.emplace_at(true, 2, 4);

    std::cerr << printSparseBool(testComplex) << std::endl;

    auto testFilt = std::vector<double>(5, 0.01);
    testFilt[3] = 0.02;
    auto eulerVec = std::vector<indsz_t>{3, 2};

    eirene::Complex comp;
    eirene::toComplex(testComplex, testFilt, eulerVec, comp, true);

    auto subset = std::vector<indsz_t>{1, 3, 4};

    const auto star = eirene::openStar(subset, comp);
    const std::vector<eirene::MorseReducedResult> res = eirene::relativeCellularPersistence(comp, star);

    return res;
}

int main()
{
    auto lingComp = lingareddyComplex();
    auto subset = std::vector<indsz_t>{1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12};

    const auto openComp = eirene::openStar(subset, lingComp);
    const std::vector<eirene::MorseReducedResult> res = eirene::relativeCellularPersistence(lingComp, openComp);

    const std::vector<eirene::MorseReducedResult> res2 = barbellComplexRelativeHom();

    return 0;
}
