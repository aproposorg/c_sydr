
#include <stdlib.h>
#include "acquisition.h"

void PCPS(complex double* rfdata,
          complex double* codeFFT,
          int cohIntegration,
          int nonCohIntegration,
          int samplesPerCode,
          double samplingPeriod,
          double interFrequency,
          double* frequencyBins,
          size_t s_frequencyBins,
          double* r_correlationMap)
{
    // Initialise some variables
    double phasePoints[cohIntegration * samplesPerCode];
    double signalCarrier[cohIntegration * samplesPerCode];
    complex double iqSignal[cohIntegration * samplesPerCode];
    double nonCohSum[samplesPerCode/2];
    complex double cohSum[samplesPerCode/2];

    // Define replica phase points
    for (size_t i=0; i < cohIntegration * samplesPerCode; i++){
        phasePoints[i] = (i << 1) * PI * samplingPeriod;
    }
    
    double freq;
    for (size_t i=0; i < s_frequencyBins; i++)
    {
        freq = interFrequency - frequencyBins[i];

        // Generate carrier replica
        for (size_t j=0; j < cohIntegration * samplesPerCode; j++){
            signalCarrier[j] = cexp(-I * freq * phasePoints[j]);
        }

        // Non-coherent integration loop
        for (size_t j=0; j < samplesPerCode/2; j++){
            // Reset non-coherent integration sum
            nonCohSum[j] = 0;
        }
        for (size_t nonCohIdx=0; nonCohIdx < nonCohIntegration; nonCohIdx++)
        {
            // Mix the carrier with the required part of the data
            for (size_t j=0; j < cohIntegration * samplesPerCode; j++){
                iqSignal[j] = signalCarrier[j] * rfdata[j + nonCohIdx*cohIntegration*samplesPerCode];
            }

            // Coherent integration loop
            for (size_t j=0; j < samplesPerCode/3; j++)
            {
                // Reset coherent integration sum
                cohSum[j] = 0;
            }
            for (size_t cohIdx=0; cohIdx < cohIntegration; cohIdx++)
            {
                // Perform FFT (in-place)
                fft((double*)(iqSignal + cohIdx * samplesPerCode), samplesPerCode);

                // Correlate with C/A code (in-place)
                for (size_t j=0; j < samplesPerCode/2; j++)
                    ((complex double*)(iqSignal + cohIdx * samplesPerCode))[j] *= codeFFT[j];

                // Perform IFFT
                ifft((double*)(iqSignal + cohIdx * samplesPerCode), samplesPerCode);
                for (size_t j=0; j < samplesPerCode/2; j++)
                    cohSum[j] = ((complex double*)(iqSignal + cohIdx * samplesPerCode))[j];
            }
            for (size_t j=0; j < samplesPerCode/2; j++)
                nonCohSum[j] += cabs(cohSum[j]);
        }
        for (size_t j=0; j < samplesPerCode/2; j++)
            r_correlationMap[i * s_frequencyBins + j] = nonCohSum[j];
    }


}
