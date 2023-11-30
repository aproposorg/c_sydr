
#include "tracking.h"

using namespace std;

// ====================================================================================================================

void EPL(
    const int* rfdata,
    size_t nbSamples,
    const char* code, // Assuming code of 1025 bits, wrapping previous/last bits
    size_t sizeCode,
    double samplingFreq,
    double carrierFreq,
    double remCarrier,
    double remCode,
    double codeStep,
    double corrSpacing,
    double* r_corrResults)
{
    // Init correlators
    double spacings[] = {-corrSpacing, 0.0, corrSpacing};
    for (size_t i {}; i < 3; ++i)
    {
        r_corrResults[i*2]   = 0.0;
        r_corrResults[i*2+1] = 0.0;
    }

    // Perform correlation
    for (size_t i {}; i < nbSamples; ++i)
    {
        // Generate signal replica
        double phase = -(carrierFreq * 2.0 * PI * (i / samplingFreq)) + remCarrier;
        complex<double> replica = exp(1i * phase);

        // Mix carrier
        replica *= complex<double>(rfdata[2*i], rfdata[2*i+1]);

        // Mix with code
        for (size_t j {}; j < 3; ++j)
        {
            size_t idx = ceil(remCode + spacings[j] + i * codeStep);
            r_corrResults[j*2]   += code[idx] * replica.real();
            r_corrResults[j*2+1] += code[idx] * replica.imag();
        }
    }
    
    return;
}

// ====================================================================================================================
// LOCK LOOPS

double DLL_NNEML(double iEarly, double qEarly, double iLate, double qLate)
{
    double earlySqare = sqrt(iEarly * iEarly + qEarly * qEarly);
    double lateSquare = sqrt(iLate * iLate + qLate * qLate);
    double codeError = (earlySqare - lateSquare) / (earlySqare + lateSquare);

    return codeError;
}

// --------------------------------------------------------------------------------------------------------------------

double PLL_costa(double iPrompt, double qPrompt)
{
    double phaseError = atan(qPrompt / iPrompt) / TWO_PI;

    return phaseError;
}

// ====================================================================================================================
// FILTERS

double LoopFilterTau1(double loopNoiseBandwidth, double dampingRatio, double loopGain)
{
    double Wn = loopNoiseBandwidth * 8.0 * dampingRatio / (4.0 * dampingRatio*dampingRatio + 1);
    double tau1 = loopGain / (Wn*Wn);

    return tau1;
}

// --------------------------------------------------------------------------------------------------------------------

double LoopFilterTau2(double loopNoiseBandwidth, double dampingRatio)
{
    double Wn = loopNoiseBandwidth * 8.0 * dampingRatio / (4.0 * dampingRatio*dampingRatio + 1);
    double tau2 = 2.0 * dampingRatio / Wn;

    return tau2;
}

// --------------------------------------------------------------------------------------------------------------------

double BorreLoopFilter(double input, double memory, double tau1, double tau2, double pdi)
{
    double output = tau2 / tau1 * (input - memory) + pdi / tau1 * input;

    return output;
}

// ====================================================================================================================
// INDICATORS

double PLLIndicator(double iprompt, double qprompt, double previous, double alpha)
{
    // Narrow Band Difference
    double nbd = iprompt*iprompt - qprompt*qprompt;

    // Narrow Band Power
    double nbp = iprompt*iprompt + qprompt*qprompt;

    // Pass through low-pass filter
    double pllLock = nbd / nbp;
    pllLock = (1 - alpha) * previous + alpha * pllLock;

    return pllLock;
}

// --------------------------------------------------------------------------------------------------------------------

double CN0_Baulieu(double pdpnRatio, double previous, double nbSamples, double alpha)
{    
    double lambda_c = 1 / (pdpnRatio / nbSamples);

    double cn0 = lambda_c / (nbSamples * 1e-3); // Divide by number of milliseconds
    cn0 = (1 - alpha) * previous + alpha * cn0;

    return cn0;
}
