
#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "structures.h"

using namespace std;

#define PI     3.1415926535898 // GPS definition of Pi
#define TWO_PI 6.2831853071796

#define GPS_L1CA_CODE_SIZE_BITS 1023
#define GPS_L1CA_CODE_FREQ      1023e3
#define LNAV_MS_PER_BIT         20

/// ===========================================================================
// Channel-related constants
#define RESULTS_SIZE 32
#define NB_CORRELATORS 3
#define I_EARLY_IDX 0
#define Q_EARLY_IDX 1
#define I_PROMPT_IDX 2
#define Q_PROMPT_IDX 3
#define I_LATE_IDX 4
#define Q_LATE_IDX 5

/// ===========================================================================
// Auto-generate default configurations for the signals and channels
/// Default signal parameters
#define SAMPLING_FREQ     10e6
#define INTERMEDIATE_FREQ 0.0

/// Auto-generate a populated signal configuration
constexpr st_SignalConfig DefaultSignalConfig {
    .samplingFreq = SAMPLING_FREQ,
    .intermediateFreq = INTERMEDIATE_FREQ,

    .codeFreqBasis = GPS_L1CA_CODE_FREQ,
    .codeLengthBits = GPS_L1CA_CODE_SIZE_BITS
};

/// Default channel parameters
//// Acquisition
#define DOPPLER_RANGE      5000
#define DOPPLER_STEP       250
#define ACQ_THRESHOLD      1.0
#define COH_INTEGRATION    1
#define NONCOH_INTEGRATION 1

//// Tracking
#define CORRELATOR_SPACING 0.5
#define DLL_DAMP_RATIO     0.7
#define DLL_NOISE_BW       1.0
#define DLL_LOOP_GAIN      1.0
#define DLL_PDI            1e-3

#define PLL_DAMP_RATIO      0.7
#define PLL_NOISE_BW        15.0
#define PLL_LOOP_GAIN       0.25
#define PLL_PDI             1e-3
#define PLL_INDICATOR_ALPHA 0.01

#define CN0_ALPHA 0.1

/// Auto-generate a populated channel configuration
constexpr st_ChannelConfig DefaultChannelConfig {
    .signalConfig = DefaultSignalConfig,

    .dopplerRange = DOPPLER_RANGE,
    .dopplerStep = DOPPLER_STEP,
    .cohIntegration = COH_INTEGRATION,
    .nonCohIntegration = NONCOH_INTEGRATION,
    .acqThreshold = ACQ_THRESHOLD,

    .correlatorSpacing = CORRELATOR_SPACING,
    .dllNoiseBW = DLL_NOISE_BW,
    .dllDampRatio = DLL_DAMP_RATIO,
    .dllLoopGain = DLL_LOOP_GAIN,
    .dllPDI = DLL_PDI,
    .pllNoiseBW = PLL_NOISE_BW,
    .pllDampRatio = PLL_DAMP_RATIO,
    .pllLoopGain = PLL_LOOP_GAIN,
    .pllPDI = PLL_PDI,

    .pllIndicatorAlpha = PLL_INDICATOR_ALPHA,
    .cn0Alpha = CN0_ALPHA
};

/// ===========================================================================
// Auto-generate the PRN codes for all GPS L1 C/A satellites at compile-time
#define NB_GPS_L1CA_SATS 31

/**
 * The following is adapted from:
 * 
 * Project Title: GNSS-R SDR
 * Author: John Bagshaw
 * Co-author: Surabhi Guruprasad
 * Contact: jotshaw@yorku.ca
 * Supervisor: Prof. Sunil Bisnath
 * Project Manager: Junchan Lee
 * Institution: Lassonde School of Engineering, York University, Canada.
 **/
const short G2Shifts[33] = {0, 5, 6, 7, 8, 17, 18, 139, 140, 141,
                            251, 252, 254, 255, 256, 257, 258, 469,
                            470, 471, 472, 473, 474, 509, 512, 513,
                            514, 515, 516, 859, 860, 861, 862};

template<unsigned N>
struct st_SatConfigs {
private:
    int g1 [GPS_L1CA_CODE_SIZE_BITS];
    int g2 [GPS_L1CA_CODE_SIZE_BITS];
    int g2b[GPS_L1CA_CODE_SIZE_BITS];

    char codes[N][GPS_L1CA_CODE_SIZE_BITS + 2];

public:
    constexpr st_SatConfigs() : g1(), g2(), g2b(), codes() {
        for (unsigned char s {}; s < N; ++s)
        {
            int reg1[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
            int reg2[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
            int saveBit {};

            /* --- Generate G1 code ----------------------------------------------------- */
            /* --- Generate all G1 signal chips based on the G1 feedback polynomial ----- */
            for (size_t nx {}; nx < GPS_L1CA_CODE_SIZE_BITS; ++nx)
            {
                g1[nx] = reg1[9];
                saveBit = reg1[9] ^ reg1[2];
                for (size_t i = 9; i >= 1; i--)
                {
                    reg1[i] = reg1[i - 1];
                }
                reg1[0] = saveBit;
            }

            /* --- Generate G2 code ----------------------------------------------------- */
            /* --- Initialize g2 output to speed up the function --- */
            /* --- Generate all G2 signal chips based on the G2 feedback polynomial ----- */
            for (size_t nx {}; nx < GPS_L1CA_CODE_SIZE_BITS; ++nx)
            {
                g2[nx] = reg2[9];
                saveBit = reg2[9];
                saveBit ^= reg2[8];
                saveBit ^= reg2[7];
                saveBit ^= reg2[5];
                saveBit ^= reg2[2];
                saveBit ^= reg2[1];
                for (size_t i = 9; i >= 1; i--)
                    reg2[i] = reg2[i - 1];
                reg2[0] = saveBit;
            }

            /* --- Shift G2 code -------------------------------------------------------- */
            /* The idea: g2 = concatenate[ g2_right_part, g2_left_part ]; */
            /* --- Form single sample C/A code by multiplying G1 and G2 ----------------- */
            const short g2shift = G2Shifts[s + 1];
            if (g2shift)
            {
                size_t idx {};
                for (size_t nx {(size_t)(GPS_L1CA_CODE_SIZE_BITS - g2shift)}; nx < GPS_L1CA_CODE_SIZE_BITS; ++nx)
                    g2b[idx++] = g2[nx];
                for (size_t nx {}; nx < (size_t)(GPS_L1CA_CODE_SIZE_BITS - g2shift); ++nx)
                    g2b[idx++] = g2[nx];
            }

            for (size_t nx {}; nx < GPS_L1CA_CODE_SIZE_BITS; ++nx)
                codes[s][nx + 1] = g1[nx] == g2b[nx] ? -1 : 1;

            /* --- Wrap code ------------------------------------------------------------ */
            codes[s][0] = codes[s][GPS_L1CA_CODE_SIZE_BITS];
            codes[s][GPS_L1CA_CODE_SIZE_BITS+1] = codes[s][1];
        }
    }

    const char* operator[](unsigned char prn) const {
        return codes[prn-1];
    }
};

static constexpr auto GPS_L1CA_CODES = st_SatConfigs<NB_GPS_L1CA_SATS>();

#endif
