
#include "acquisition.h"

using namespace std;

void SerialSearch(
    Eigen::MatrixXcd rfdata, 
    size_t size,
    double* code, 
    int dopplerRange, 
    int dopplerStep, 
    float samplingFrequency,
    double* r_correlation)
{
    int frequencyBins = dopplerRange * 2 / dopplerStep + 1;
    Eigen::MatrixXcd replica = Eigen::MatrixXcd::Zero(1, size);
    Eigen::MatrixXd code_eigen = Eigen::Map<Eigen::MatrixXd>(code, 1, GPS_L1CA_CODE_SIZE_BITS);
    Eigen::MatrixXd _code = Eigen::Map<Eigen::MatrixXd>(code, 1, GPS_L1CA_CODE_SIZE_BITS);
    Eigen::MatrixXd _upcode = Eigen::MatrixXd::Zero(1, size);
    double iSignal;
    double qSignal; 
    double codeStep = GPS_L1CA_CODE_FREQ / samplingFrequency;

    // Gerenate phase points
    //Eigen::MatrixXd phases = Eigen::MatrixXd::LinSpaced(size, 0, (2 * PI * size / samplingFrequency));

    // Code roll indices
    int codeRollIdx[GPS_L1CA_CODE_SIZE_BITS];
    codeRollIdx[0] = GPS_L1CA_CODE_SIZE_BITS - 1;
    for(int idx=1; idx < GPS_L1CA_CODE_SIZE_BITS; idx++){
        codeRollIdx[idx] = idx - 1;
    }

    double freq = -dopplerRange;
    for(int idxFreq=0; idxFreq <= frequencyBins; idxFreq++){

        // Generate replica
        for(int idxSamples=0; idxSamples < size; idxSamples++){
            replica(0, idxSamples) = exp(1i * (double) (freq * 2.0 * PI * idxSamples/samplingFrequency));
            //cout << replica(0, idxSamples) << endl;
        }
        replica = replica.cwiseProduct(rfdata);

        for(int idxCode=0; idxCode < GPS_L1CA_CODE_SIZE_BITS; idxCode++){

            // Shift the code
            if(idxCode > 0){
                _code = code_eigen(Eigen::all, codeRollIdx);
                code_eigen = _code;
            }

            // Code upsamples
            int codeUpIdx = 0;
            for(int idx=0; idx < size; idx++){
                codeUpIdx = trunc(GPS_L1CA_CODE_FREQ * idx / samplingFrequency);
                _upcode(0, idx) = _code(0, codeUpIdx);
            }

            iSignal = real(replica.real().cwiseProduct(_upcode).sum());
            qSignal = real(replica.imag().cwiseProduct(_upcode).sum());

            // Accumlate
            r_correlation[idxFreq * GPS_L1CA_CODE_SIZE_BITS + idxCode] = iSignal * iSignal + qSignal * qSignal;
        }

        freq += dopplerStep;
    }
    return;
}

// ====================================================================================================================

void TwoCorrelationPeakComparison(double* correlationMap, size_t sizeMap, int* r_idxPeak, float* r_acqMetric) {
    
    // Find largest value
    int idxPeak1 = FindMaxIndex(correlationMap, sizeMap, -1);

    // Find second largest
    int idxPeak2 = FindMaxIndex(correlationMap, sizeMap, idxPeak1);

    *r_acqMetric = correlationMap[idxPeak1] / correlationMap[idxPeak2];
    *r_idxPeak = idxPeak1;

    return;
}

// ====================================================================================================================

int FindMaxIndex(double* array, size_t size, int excludeIdx){
    int idxMax = 0; 
    if (excludeIdx == 0){
        idxMax = 1;
    }
    for(int i=0; i < size; i++){   
        if (excludeIdx == i){
            continue;
        }
        if(array[idxMax] < array[i]){
            idxMax = i;
        }
    }

    return idxMax;
}

// ====================================================================================================================
// /!\ NOT TESTED 
// void PCPS(complex double* rfdata,
//           complex double* codeFFT,
//           int cohIntegration,
//           int nonCohIntegration,
//           int samplesPerCode,
//           double samplingPeriod,
//           double interFrequency,
//           double* frequencyBins,
//           size_t s_frequencyBins,
//           double* r_correlationMap)
// {
//     // Initialise some variables
//     double phasePoints[cohIntegration * samplesPerCode];
//     double signalCarrier[cohIntegration * samplesPerCode];
//     complex double iqSignal[cohIntegration * samplesPerCode];
//     double nonCohSum[samplesPerCode/2];
//     complex double cohSum[samplesPerCode/2];

//     // Define replica phase points
//     for (size_t i=0; i < cohIntegration * samplesPerCode; i++){
//         phasePoints[i] = (i << 1) * PI * samplingPeriod;
//     }
    
//     double freq;
//     for (size_t i=0; i < s_frequencyBins; i++)
//     {
//         freq = interFrequency - frequencyBins[i];

//         // Generate carrier replica
//         for (size_t j=0; j < cohIntegration * samplesPerCode; j++){
//             signalCarrier[j] = cexp(-I * freq * phasePoints[j]);
//         }

//         // Non-coherent integration loop
//         for (size_t j=0; j < samplesPerCode/2; j++){
//             // Reset non-coherent integration sum
//             nonCohSum[j] = 0;
//         }
//         for (size_t nonCohIdx=0; nonCohIdx < nonCohIntegration; nonCohIdx++)
//         {
//             // Mix the carrier with the required part of the data
//             for (size_t j=0; j < cohIntegration * samplesPerCode; j++){
//                 iqSignal[j] = signalCarrier[j] * rfdata[j + nonCohIdx*cohIntegration*samplesPerCode];
//             }

//             // Coherent integration loop
//             for (size_t j=0; j < samplesPerCode/3; j++)
//             {
//                 // Reset coherent integration sum
//                 cohSum[j] = 0;
//             }
//             for (size_t cohIdx=0; cohIdx < cohIntegration; cohIdx++)
//             {
//                 // Perform FFT (in-place)
//                 fft((double*)(iqSignal + cohIdx * samplesPerCode), samplesPerCode);

//                 // Correlate with C/A code (in-place)
//                 for (size_t j=0; j < samplesPerCode/2; j++)
//                     ((complex double*)(iqSignal + cohIdx * samplesPerCode))[j] *= codeFFT[j];

//                 // Perform IFFT
//                 ifft((double*)(iqSignal + cohIdx * samplesPerCode), samplesPerCode);
//                 for (size_t j=0; j < samplesPerCode/2; j++)
//                     cohSum[j] = ((complex double*)(iqSignal + cohIdx * samplesPerCode))[j];
//             }
//             for (size_t j=0; j < samplesPerCode/2; j++)
//                 nonCohSum[j] += cabs(cohSum[j]);
//         }
//         for (size_t j=0; j < samplesPerCode/2; j++)
//             r_correlationMap[i * s_frequencyBins + j] = nonCohSum[j];
//     }
// }
