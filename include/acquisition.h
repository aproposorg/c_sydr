
#ifndef ACQUISITION_H
#define ACQUISITION_H

#include <iostream>
#include <complex>
#include <cmath>
#include "constants.h"
#include "Eigen/Dense"

using namespace std;

void SerialSearch(Eigen::MatrixXcd, size_t, double*, int, int, float, double*);
void TwoCorrelationPeakComparison(double*, size_t, int*, float*);
int  FindMaxIndex(double*, size_t, int);


// struct pad {
//   Index size() const { return out_size; }
//   Index operator[] (Index i) const { return std::max<Index>(0,i-(out_size-in_size)); }
//   Index in_size, out_size;
// };
 

#endif


