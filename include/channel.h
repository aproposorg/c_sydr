
#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <ccomplex>
#include "acquisition.h"
#include "tracking.h"
#include "constants.h"
#include <cassert>
//#include <Eigen/Dense>

/// ===========================================================================

using namespace std;

class Channel{
private:
    // Configuration flag
    bool m_configured {false};

    // Acquisition
    float m_acqMetric;
    size_t m_indexPeak;

    // Tracking
    bool m_isTrackingInitialised;
    size_t m_trackRequiredSamples;
    double m_iPromptSum;
    double m_qPromptSum;
    int m_nbPromptSum;
    double m_remainingCarrier;
    double m_remainingCode;
    double m_pllDiscrim;
    double m_dllDiscrim;
    double m_codeError;
    double m_phaseError;
    double m_dllTau1;
    double m_dllTau2;
    double m_dllPDI;
    double m_pllTau1;
    double m_pllTau2;
    double m_pllPDI;
    unsigned m_codeCounter;

    // Indicators
    double m_pllLock;
    double m_pdpnRatio; // For C/N0 estimation
    double m_cn0;

    // Channel identification and control
    int m_channelID{-1};
    unsigned m_satelliteID;
    unsigned m_channelState;
    st_ChannelConfig m_config;

    // Intermediate results
    int m_codeOffset;
    double m_carrierFrequency;
    double m_codeFrequency;
    double m_correlatorsResults[NB_CORRELATORS * 2];

public:
    // Constructor
    Channel();
    Channel(char, st_ChannelConfig);

    // Destructor
    ~Channel();

    // Configuration
    void setSatellite(unsigned char);
    void configure(char, st_ChannelConfig);

    // General processing
    void run(const int* _rfdata, size_t size);
    void resetChannel();
    void resetCounters();

    // Acquisition
    void runAcquisition(const int*, size_t);
    void runSignalSearch(const int*, size_t, float*);
    void runPeakFinder(const float*, size_t);
    void postAcquisitionUpdate();
    void runNoMapAcquisition(const int*, size_t);

    // Tracking
    void initTracking(size_t);
    void runTracking(const int*, size_t);
    void runCorrelators(const int*);
    void runDiscriminatorsFilters();
    void runLoopIndicators();
    void postTrackingUpdate();
    const double* getCorrelatorResults() const {
        return m_correlatorsResults;
    }
};

#endif
