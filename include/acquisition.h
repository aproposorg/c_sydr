
#ifndef ACQUISITION_H
#define ACQUISITION_H

#include <iostream>
#include <complex>
#include <cmath>
#include "constants.h"
//#include "Eigen/Dense"
#include "fft.h"
#include "utils.h"

using namespace std;

//void SerialSearch(Eigen::MatrixXcd, size_t, double*, int, int, float, double*);
void TwoCorrelationPeakComparison(float*, size_t, int, int, int*, float*);
int  FindMaxIndex(float*, size_t);
int  FindMaxIndexWithExclude(float*, size_t, int, int);
void PCPS(double*, size_t, int*, size_t, double, int, int, double, double, int, int, float*);

// struct pad {
//   Index size() const { return out_size; }
//   Index operator[] (Index i) const { return std::max<Index>(0,i-(out_size-in_size)); }
//   Index in_size, out_size;
// };
 

#endif


