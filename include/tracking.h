
#ifndef TRACKING_H
#define TRACKING_H

#include <iostream>
#include <complex>
#include <cmath>
#include "constants.h"
#include "fft.h"
#include "utils.h"

using namespace std;

// Correlators
void EPL(const int*, size_t, const char*, size_t, double, double, double, double, double, double, double*);

// Lock loops
double DLL_NNEML(double, double, double, double);
double PLL_costa(double, double);

// Discriminators
double LoopFilterTau1(double, double, double);
double LoopFilterTau2(double, double);
double BorreLoopFilter(double, double, double, double, double);

// Indicators
double PLLIndicator(double, double, double, double);
double CN0_Baulieu(double, double, double, double);

#endif
