A C++ port of https://github.com/Eetion/Eirene.jl. Don't use this for point clouds/vr complexes, since it doesn't actually output anything useful in that case. 
This would take a little more work but there is other software that does that in most languages anyways.

Check out src/eireneTest.cpp for examples on how to use this from C++. The the thing this does that Eirene.jl doesn't is it can compute relative homology, like in the paper "Local homology of abstract simplicial complexes".
I've checked some examples but there are probably still bugs I have yet to hit.

```
auto barbellComplexRelativeHom()
{
    spla::SparseBool testComplex(5, 5);

    testComplex.emplace_at(true, 0, 3);

    testComplex.emplace_at(true, 1, 3);
    testComplex.emplace_at(true, 1, 4);

    testComplex.emplace_at(true, 2, 4);

    std::cerr << printSparseBool(testComplex) << std::endl;

    auto testFilt = std::vector<double>(5, 0.01);
    auto eulerVec = std::vector<indsz_t>{3, 2};

    eirene::Complex comp;
    eirene::toComplex(testComplex, testFilt, eulerVec, comp, true);

    auto subset = std::vector<indsz_t>{1, 3, 4};

    const auto star = eirene::openStar(subset, comp);
    const std::vector<eirene::MorseReducedResult> res = eirene::relativeCellularPersistence(comp, star);

    return res;
}
```

Here we simply take a cellular complex consisting of 3 0-cells and two 1-cells connected in a chain, like an open necklace. The vertices have ('global') indices 0, 1, 2, and the 
1-cells have indices 3, 4. We then take calculate the open star (add in all the cofaces) of the middle vertex and the two 1-cells (which is just that set again in this case), and compute 
the homology of the whole complex relative to the subcomplex (subcomplex <==> closed in Alexandrov <==> contains all it's faces) (X - star). Thus, we'd expect a nontrivial cycle corresponding 
to closing the necklace.

The nontrivial 1-dimensional relative homology can be seen in res:

```
[2] = {
    leftFactor = {
      rowVals = size=0 {}
      colPtrs = size=2 {
        [0] = 0
        [1] = 0
      }
      numRows = 0
    }
    tid = size=1 {
      [0] = 0
    }
    pairVec = size=0 {}
    cycleReps = size=1 {
      [0] = size=2 {
        [0] = 1
        [1] = 0
      }
    }
    birthDeath = size=2 {
      [0] = 0.02
      [1] = +Inf
    }
  }
```

Note the cell indexes in cycleReps are not counting from 0-dimensional cells ('global') but from the first 1-dimensional cells ('local'/'graded'), just like in Eirene.jl.

It does check some aspects of the complex being well defined but not all. For example, make sure the filtration respects the complex structure! 
A face being born after a coface makes no sense.

What python API there is, is copied from ripser and you do currently need to make sure all the numpy dtypes are double. Check out test.py for what that looks like. In the python case it takes the open star
for you so no need to worry about that.

To build the cython shared object do this. And in general you might want to change the compile flags to enable optimization (setup.py, Makefile).

```
python setup.py build_ext --inplace --debug
```







