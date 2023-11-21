
#include "acquisition.h"

using namespace std;

// ====================================================================================================================

template <typename T>
size_t FindMaxIndex(const T* array, size_t size)
{
    size_t idxMax{};

    for (size_t i{1}; i < size; ++i)
        idxMax = array[idxMax] < array[i] ? i : idxMax;

    return idxMax;
}

// ====================================================================================================================

template <typename T>
size_t FindMaxIndexWithExclude(const T* array, size_t size, size_t excludeIdx, size_t excludeChipSize)
{
    size_t idxMax = excludeIdx == 0 ? 1 : 0;

    for (size_t i{idxMax + 1}; i < size; ++i)
    {
        if ((size_t)abs((int)(excludeIdx - i)) < excludeChipSize)
            continue;

        idxMax = array[idxMax] < array[i] ? i : idxMax;
    }

    return idxMax;
}

// ====================================================================================================================

void TwoCorrelationPeakComparison(const float* correlationMap, size_t sizeMap,
                                  size_t samplesPerCode, size_t samplesPerCodeChip, 
                                  size_t* r_idxPeak, float* r_acqMetric)
{
    // Find largest value
    size_t idxPeak1 = FindMaxIndex<float>(correlationMap, sizeMap);

    // Find second largest in the same frequency
    size_t freqOffset = idxPeak1 - (idxPeak1 % samplesPerCode);
    size_t idxPeak2 = FindMaxIndexWithExclude<float>(correlationMap + freqOffset, samplesPerCode,
                                                     idxPeak1 % samplesPerCode, samplesPerCodeChip);
    idxPeak2 += freqOffset;

    *r_idxPeak = idxPeak1;
    *r_acqMetric = correlationMap[idxPeak1] / correlationMap[idxPeak2];
    return;
}

// ====================================================================================================================

void PCPS(
    const int* rfdata,
    size_t sizeData,
    const char *code,
    size_t sizeCode,
    double codeFrequency,
    size_t cohIntegration,
    size_t nonCohIntegration,
    double samplingFrequency,
    double interFrequency,
    size_t dopplerRange,
    size_t dopplerStep,
    float *r_correlationMap)
{
    // Initialise some variables
    const size_t samplesPerCode = samplingFrequency * sizeCode / codeFrequency;
    const size_t dopplerBins = dopplerRange * 2 / dopplerStep + 1;

    // Allocate buffers before the FFT
    double iqSignal[samplesPerCode * 2];
    double nonCohSum[samplesPerCode];
    complex<double> cohSum[samplesPerCode];

    double codeUpsampled[samplesPerCode];
    complex<double> codeFFT[samplesPerCode];

    upsampleCode(code, GPS_L1CA_CODE_SIZE_BITS, GPS_L1CA_CODE_FREQ, samplingFrequency, codeUpsampled); // Upsampling code to sampling frequency
    rfft(codeUpsampled, samplesPerCode); // in-place FFT

    // Need to reconstruct the real-to-complex FFT, as pocket FFT only compute the real part
    // For this, we mirror the FFT coefficients with their conjugate. See FFT theory for more details.
    // We also apply another conjugate here for mixing later.
    // Note: in case the number of samples is odd, there might be an issue with this...
    for (size_t i {1}; i < samplesPerCode/2; ++i)
    {
        complex<double> _complex(codeUpsampled[2*i-1], codeUpsampled[2*i]);
        codeFFT[i] = conj(_complex);
        codeFFT[samplesPerCode-i] = _complex;
    }
    codeFFT[0] = codeUpsampled[0];
    codeFFT[samplesPerCode/2] = codeUpsampled[samplesPerCode-1];

    double freq = -(int)dopplerRange;
    for (size_t idxDoppler {}; idxDoppler < dopplerBins; ++idxDoppler)
    {
        freq -= interFrequency;

        // Reset non-coherent integration sum
        for (size_t j {}; j < samplesPerCode; ++j)
            nonCohSum[j] = 0;

        // Non-coherent integration loop
        for (size_t nonCohIdx {}; nonCohIdx < nonCohIntegration; ++nonCohIdx)
        {
            // Reset coherent integration sum
            for (size_t j {}; j < samplesPerCode; ++j)
                cohSum[j] = 0;

            // Coherent integration loop
            for (size_t cohIdx {}; cohIdx < cohIntegration; ++cohIdx)
            {
                 // Mix the carrier with the required part of the data
                for (size_t j {}; j < samplesPerCode; ++j)
                {
                    // Generate carrier
                    double phase = (cohIdx * samplesPerCode + j) * 2 * PI / samplingFrequency;
                    complex<double> signalCarrier = exp(-1i * freq * phase);

                    // Mix
                    size_t _idx = j + (nonCohIdx * cohIntegration + cohIdx) * samplesPerCode;
                    complex<double> _complex = signalCarrier * complex<double>((double)rfdata[2*_idx], (double)rfdata[2*_idx+1]);
                    iqSignal[2*j] = _complex.real();
                    iqSignal[2*j+1] = _complex.imag();
                }

                // Perform FFT (in-place)
                cfft(iqSignal, samplesPerCode);

                // Correlate with C/A code (in-place)
                for (size_t j {}; j < samplesPerCode; ++j)
                {
                    complex<double> _complex = codeFFT[j] * complex<double>(iqSignal[2*j], iqSignal[2*j+1]);
                    iqSignal[2*j] = _complex.real();
                    iqSignal[2*j+1] = _complex.imag(); 
                }

                // Perform IFFT
                cifft(iqSignal, samplesPerCode);

                // Sum coherent integration
                for (size_t j {}; j < samplesPerCode; ++j)
                    cohSum[j] += complex<double>(iqSignal[2*j], iqSignal[2*j+1]);
            }

            // Sum non-coherent integration
            for (size_t j {}; j < samplesPerCode; ++j)
                nonCohSum[j] += abs(cohSum[j]);
        }

        // Accumulate correlation map
        for (size_t j {}; j < samplesPerCode; ++j)
            r_correlationMap[idxDoppler * samplesPerCode + j] = nonCohSum[j];

        // Increment frequency search step
        freq += dopplerStep;
    }

    return;
}

// ====================================================================================================================

void NoMapPCPS(
    const int *rfdata,
    size_t sizeData,
    const char *code,
    size_t sizeCode,
    double codeFrequency,
    size_t cohIntegration,
    size_t nonCohIntegration,
    double samplingFrequency,
    double interFrequency,
    size_t dopplerRange,
    size_t dopplerStep,
    size_t* r_idxPeak,
    float* r_acqMetric)
{
    // Initialise some variables
    const size_t samplesPerCode = samplingFrequency * sizeCode / codeFrequency;
    const size_t samplesPerCodeChip = ceil((double)samplesPerCode / GPS_L1CA_CODE_SIZE_BITS);
    const size_t dopplerBins = dopplerRange * 2 / dopplerStep + 1;

    // Allocate buffers before the FFT
    double iqSignal[samplesPerCode * 2];
    double nonCohSum[samplesPerCode];
    complex<double> cohSum[samplesPerCode];

    double codeUpsampled[samplesPerCode];
    complex<double> codeFFT[samplesPerCode];

    upsampleCode(code, GPS_L1CA_CODE_SIZE_BITS, GPS_L1CA_CODE_FREQ, samplingFrequency, codeUpsampled); // Upsampling code to sampling frequency
    rfft(codeUpsampled, samplesPerCode);                                                               // in-place FFT

    // Need to reconstruct the real-to-complex FFT, as pocket FFT only compute the real part
    // For this, we mirror the FFT coefficients with their conjugate. See FFT theory for more details.
    // We also apply another conjugate here for mixing later.
    // Note: in case the number of samples is odd, there might be an issue with this...
    for (size_t i{1}; i < samplesPerCode / 2; ++i)
    {
        complex<double> _complex(codeUpsampled[2 * i - 1], codeUpsampled[2 * i]);
        codeFFT[i] = conj(_complex);
        codeFFT[samplesPerCode - i] = _complex;
    }
    codeFFT[0] = codeUpsampled[0];
    codeFFT[samplesPerCode / 2] = codeUpsampled[samplesPerCode - 1];

    double freq = -(int)dopplerRange;
    for (size_t idxDoppler {}; idxDoppler < dopplerBins; ++idxDoppler)
    {
        freq -= interFrequency;

        // Reset non-coherent integration sum
        for (size_t j{}; j < samplesPerCode; ++j)
            nonCohSum[j] = 0;

        // Non-coherent integration loop
        for (size_t nonCohIdx{}; nonCohIdx < nonCohIntegration; ++nonCohIdx)
        {
            // Reset coherent integration sum
            for (size_t j{}; j < samplesPerCode; ++j)
                cohSum[j] = 0;

            // Coherent integration loop
            for (size_t cohIdx{}; cohIdx < cohIntegration; ++cohIdx)
            {
                // Mix the carrier with the required part of the data
                for (size_t j{}; j < samplesPerCode; ++j)
                {
                    // Generate carrier
                    double phase = (cohIdx * samplesPerCode + j) * 2 * PI / samplingFrequency;
                    complex<double> signalCarrier = exp(-1i * freq * phase);

                    // Mix
                    size_t _idx = j + (nonCohIdx * cohIntegration + cohIdx) * samplesPerCode;
                    complex<double> _complex = signalCarrier * complex<double>(rfdata[2 * _idx], rfdata[2 * _idx + 1]);
                    iqSignal[2 * j] = _complex.real();
                    iqSignal[2 * j + 1] = _complex.imag();
                }

                // Perform FFT (in-place)
                cfft(iqSignal, samplesPerCode);

                // Correlate with C/A code (in-place)
                for (size_t j{}; j < samplesPerCode; ++j)
                {
                    complex<double> _complex = codeFFT[j] * complex<double>(iqSignal[2 * j], iqSignal[2 * j + 1]);
                    iqSignal[2 * j] = _complex.real();
                    iqSignal[2 * j + 1] = _complex.imag();
                }

                // Perform IFFT
                cifft(iqSignal, samplesPerCode);

                // Sum coherent integration
                for (size_t j{}; j < samplesPerCode; ++j)
                    cohSum[j] += complex<double>(iqSignal[2 * j], iqSignal[2 * j + 1]);
            }

            // Sum non-coherent integration
            for (size_t j{}; j < samplesPerCode; ++j)
                nonCohSum[j] += abs(cohSum[j]);
        }

        // Find and update the peak in the non-coherent integration sum as required
        size_t idxPeak1 = FindMaxIndex<double>(nonCohSum, samplesPerCode),
               idxPeak2 = FindMaxIndexWithExclude<double>(nonCohSum, samplesPerCode,
                                                          idxPeak1, samplesPerCodeChip);
        double acqMetric = nonCohSum[idxPeak1] / nonCohSum[idxPeak2];

        if ((idxDoppler == 0) || (*r_acqMetric < acqMetric))
        {
            // On the first iteration, we store results to have something to update
            // later. On other interations, we update if a new peak is found
            *r_idxPeak   = idxDoppler * samplesPerCode + idxPeak1;
            *r_acqMetric = acqMetric;
        }

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