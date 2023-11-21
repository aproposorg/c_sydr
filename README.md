# c_sydr
C/C++ implementation of the SyDR library.

Written to target C++14-capable compilers.

# TODOs

## Code

The current implementation assumes acquisition and tracking of legacy GPS L1/CA signals and, as such, only generates their codes at compile-time, see the [constants.h](./include/constants.h) header. More codes can be added as needed.

## Acquisition

The library currently implements two versions of the PCPS algorithm: a trivial one that assembles an entire two-dimensional correlation map and searches for its maximum afterwards (see `PCPS`), and a space-reduced version that searches for the maximum in one row of the map at a time (see `NoMapPCPS`). The latter avoids some space requirements of the former but prevents potential benefits of a two-dimensional search. Alternatives could save space by:

1. changing the precision of the correlation map to a much smaller type than `float`,
2. storing only a small window of the correlation map with associated code and frequency offsets, or
3. store a handful of rows of the correlation map with associated frequency offsets.
