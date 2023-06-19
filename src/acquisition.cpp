
#include "acquisition.h"

using namespace std;

// ====================================================================================================================

void TwoCorrelationPeakComparison(float* correlationMap, size_t sizeMap, int samplesPerCode, int samplesPerCodeChip, 
                                  int* r_idxPeak, float* r_acqMetric) {
    
    // Find largest value
    int idxPeak1 = FindMaxIndex(correlationMap, sizeMap);
    
    // Find second largest
    int startIdx = (int) idxPeak1 - (idxPeak1 % samplesPerCode);
    int idxPeak2 = FindMaxIndexWithExclude((float*)(correlationMap + startIdx), samplesPerCode, 
                                           idxPeak1 % samplesPerCode, samplesPerCodeChip);
    idxPeak2 += startIdx;
    *r_acqMetric = correlationMap[idxPeak1] / correlationMap[idxPeak2];
    *r_idxPeak = idxPeak1;

    return;
}

// ====================================================================================================================

int FindMaxIndex(float* array, size_t size){
    int idxMax = 0; 
    for(int i=0; i < size; i++){
        if(array[idxMax] < array[i]){
            idxMax = i;
        }
    }

    return idxMax;
}

// ====================================================================================================================

int FindMaxIndexWithExclude(float* array, size_t size, int excludeIdx, int excludeChipSize){
    int idxMax = 0; 
    if (excludeIdx == 0){
        idxMax = 1;
    }
    for(int i=0; i < size; i++){
        if (abs(excludeIdx - i) < excludeChipSize){
            continue;
        }
        if(array[idxMax] < array[i]){
            idxMax = i;
        }
    }

    return idxMax;
}

// ====================================================================================================================

void PCPS(double* rfdata,
          size_t sizeData,
          int* code,
          size_t sizeCode,
          double codeFrequency,
          int cohIntegration,
          int nonCohIntegration,
          double samplingFrequency,
          double interFrequency,
          int dopplerRange,
          int dopplerStep,
          float* r_correlationMap)
{
    // Initialise some variables
    int samplesPerCode = samplingFrequency * sizeCode / codeFrequency;
    double phasePoints[cohIntegration * samplesPerCode];
    complex<double> signalCarrier[cohIntegration * samplesPerCode];
    double iqSignal[cohIntegration * samplesPerCode * 2];

    double nonCohSum[samplesPerCode];
    complex<double> cohSum[samplesPerCode];

    int dopplerBins = dopplerRange * 2 / dopplerStep + 1;

    // Define replica phase points
    for (size_t i=0; i < cohIntegration * samplesPerCode; i++){
        phasePoints[i] = (i << 1) * PI / samplingFrequency;
    }

    // Prepare code FFT
    double codeUpsampled[samplesPerCode];
    complex<double> codeFFT[samplesPerCode];
    upsampleCode(code, sizeCode, GPS_L1CA_CODE_FREQ, samplingFrequency, codeUpsampled); // Upsampling code to sampling frequency
    rfft(codeUpsampled, samplesPerCode); // in-place FFT

    // Need to reconstruct the real-to-complex FFT, as pocket FFT only compute the real part
    // For this, we mirror the FFT coefficients with their conjugate. See FFT theory for more details.
    // We also apply another conjugate here for mixing later.
    // Note: in case the number of samples is odd, there might be an issue with this...
    int j = 1;
    for (size_t i=1; i < samplesPerCode/2; i++){
        complex<double> _complex = codeUpsampled[j] + 1i * codeUpsampled[j+1];
        codeFFT[i] = conj(_complex);
        codeFFT[samplesPerCode-i] = _complex;
        j += 2;
    }
    codeFFT[0] = codeUpsampled[0];
    codeFFT[samplesPerCode/2] = codeUpsampled[samplesPerCode-1];
    
    double freq = -dopplerRange;
    for (size_t idxDoppler=0; idxDoppler < dopplerBins; idxDoppler++)
    {
        freq -= interFrequency;

        // Generate carrier replica
        for (size_t j=0; j < cohIntegration * samplesPerCode; j++){
            signalCarrier[j] = exp(-1i * freq * phasePoints[j]);
        }

        // Non-coherent integration loop
        for (size_t j=0; j < samplesPerCode; j++){
            nonCohSum[j] = 0; // Reset non-coherent integration sum
        }
        for (size_t nonCohIdx=0; nonCohIdx < nonCohIntegration; nonCohIdx++)
        {
            // Mix the carrier with the required part of the data
            for (size_t j=0; j < cohIntegration * samplesPerCode; j++){
                int _idx = j + nonCohIdx * cohIntegration * samplesPerCode;
                complex<double> _complex = signalCarrier[j] * (rfdata[2*_idx] + 1i * rfdata[2*_idx+1]);
                iqSignal[2*j] = real(_complex);
                iqSignal[2*j+1] = imag(_complex);
            }

            // Coherent integration loop
            for (size_t j=0; j < samplesPerCode; j++)
            {
                cohSum[j] = 0;  // Reset coherent integration sum
            }
            for (size_t cohIdx=0; cohIdx < cohIntegration; cohIdx++)
            {
                // Perform FFT (in-place)
                cfft((double*)(iqSignal + cohIdx * samplesPerCode), samplesPerCode);

                // Correlate with C/A code (in-place)
                for (size_t j=0; j < samplesPerCode; j++){
                    int _idx = j + nonCohIdx * cohIntegration * samplesPerCode;
                    complex<double> _complex = iqSignal[2*j] + 1i * iqSignal[2*j+1];
                    _complex *= codeFFT[j];
                    iqSignal[2*j] = real(_complex);
                    iqSignal[2*j+1] = imag(_complex); 
                }

                // Perform IFFT
                cifft((double*)(iqSignal + cohIdx * samplesPerCode), samplesPerCode);
                for (size_t j=0; j < samplesPerCode; j++)
                    cohSum[j] += ((complex<double>*)(iqSignal + cohIdx * samplesPerCode))[j];
            }
            for (size_t j=0; j < samplesPerCode; j++)
                nonCohSum[j] += abs(cohSum[j]);
        }
        for (size_t j=0; j < samplesPerCode; j++)
            r_correlationMap[idxDoppler * samplesPerCode + j] = nonCohSum[j];

        // Increment frequency search step
        freq += dopplerStep;
    }

    return;
}

// ====================================================================================================================

// void SerialSearch(
//     Eigen::MatrixXcd rfdata, 
//     size_t size,
//     double* code, 
//     int dopplerRange, 
//     int dopplerStep, 
//     float samplingFrequency,
//     double* r_correlation)
// {
//     int frequencyBins = dopplerRange * 2 / dopplerStep + 1;
//     Eigen::MatrixXcd replica = Eigen::MatrixXcd::Zero(1, size);
//     Eigen::MatrixXd code_eigen = Eigen::Map<Eigen::MatrixXd>(code, 1, GPS_L1CA_CODE_SIZE_BITS);
//     Eigen::MatrixXd _code = Eigen::Map<Eigen::MatrixXd>(code, 1, GPS_L1CA_CODE_SIZE_BITS);
//     Eigen::MatrixXd _upcode = Eigen::MatrixXd::Zero(1, size);
//     double iSignal;
//     double qSignal; 
//     double codeStep = GPS_L1CA_CODE_FREQ / samplingFrequency;

//     // Gerenate phase points
//     //Eigen::MatrixXd phases = Eigen::MatrixXd::LinSpaced(size, 0, (2 * PI * size / samplingFrequency));

//     // Code roll indices
//     int codeRollIdx[GPS_L1CA_CODE_SIZE_BITS];
//     codeRollIdx[0] = GPS_L1CA_CODE_SIZE_BITS - 1;
//     for(int idx=1; idx < GPS_L1CA_CODE_SIZE_BITS; idx++){
//         codeRollIdx[idx] = idx - 1;
//     }

//     double freq = -dopplerRange;
//     for(int idxFreq=0; idxFreq <= frequencyBins; idxFreq++){

//         // Generate replica
//         for(int idxSamples=0; idxSamples < size; idxSamples++){
//             replica(0, idxSamples) = exp(1i * (double) (freq * 2.0 * PI * idxSamples/samplingFrequency));
//             //cout << replica(0, idxSamples) << endl;
//         }
//         replica = replica.cwiseProduct(rfdata);

//         for(int idxCode=0; idxCode < GPS_L1CA_CODE_SIZE_BITS; idxCode++){

//             // Shift the code
//             if(idxCode > 0){
//                 _code = code_eigen(Eigen::all, codeRollIdx);
//                 code_eigen = _code;
//             }

//             // Code upsamples
//             int codeUpIdx = 0;
//             for(int idx=0; idx < size; idx++){
//                 codeUpIdx = trunc(GPS_L1CA_CODE_FREQ * idx / samplingFrequency);
//                 _upcode(0, idx) = _code(0, codeUpIdx);
//             }

//             iSignal = real(replica.real().cwiseProduct(_upcode).sum());
//             qSignal = real(replica.imag().cwiseProduct(_upcode).sum());

//             // Accumlate
//             r_correlation[idxFreq * GPS_L1CA_CODE_SIZE_BITS + idxCode] = iSignal * iSignal + qSignal * qSignal;
//         }

//         freq += dopplerStep;
//     }
//     return;
// }