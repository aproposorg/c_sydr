
#ifndef ACQUISITION_H
#define ACQUISITION_H

#include <iostream>
#include <complex>
#include <cmath>
#include "constants.h"
#include "fft.h"
#include "utils.h"

using namespace std;

template<typename T>
size_t FindMaxIndex(const T*, size_t);

template<typename T>
size_t FindMaxIndexWithExclude(const T*, size_t, size_t, size_t);

void TwoCorrelationPeakComparison(const float*, size_t, size_t, size_t, size_t*, float*);

void PCPS(const int*, size_t, const char*, size_t, double, size_t, size_t, double, double, size_t, size_t, float*);

void NoMapPCPS(const int*, size_t, const char *, size_t, double, size_t, size_t, double, double, size_t, size_t, size_t*, float*);

// void SerialSearch(Eigen::MatrixXcd, size_t, double*, int, int, float, double*);

#endif
