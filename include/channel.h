#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <ccomplex>
#include "gen_code.h"
#include "acquisition.h"
#include "tracking.h"
#include "structures.h"
#include <cstring>
//#include <Eigen/Dense>

#define RESULTS_SIZE 32
#define NB_CORRELATORS 3 
#define I_EARLY_IDX  0
#define Q_EARLY_IDX  1 
#define I_PROMPT_IDX 2
#define Q_PROMPT_IDX 3 
#define I_LATE_IDX   4
#define Q_LATE_IDX   5 

/// ===========================================================================

using namespace std;

class Channel{

    private:
        int* m_rfdata;
        size_t m_rfdataSize;

        // Tracking
        bool m_isTrackingInitialised;
        int m_trackRequiredSamples;
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
        int m_codeCounter;
        
        // Indicators
        double m_pllLock;
        double m_pdpnRatio; // For C/N0 estimation
        double m_cn0;
    
    public:
        int m_channelID;
        int m_satelliteID;
        int m_channelState;
        double m_results[RESULTS_SIZE];
        int m_code[GPS_L1CA_CODE_SIZE_BITS+2]; // Code also include the last bit from previous code and first bit from next code for tracking purposes
        st_ChannelConfig* m_config;

        int m_codeOffset;
        double m_carrierFrequency;
        double m_codeFrequency;

        double m_correlatorsResults[NB_CORRELATORS*2];
        
        // Constructor
        Channel(int, st_ChannelConfig*);

        // Destructor
        ~Channel();

        // General processing
        void run(int* _rfdata, size_t size);
        void processHandler();
        void setSatellite(int satelliteID);
        void getTimeSinceTOW();
        void prepareChannelUpdate();
        double* prepareResults();
        void resetChannel();
        void resetCounters();

        // Acquisition
        void runAcquisition();
        void runSignalSearch(float*);
        void runPeakFinder(float*, size_t);
        void postAcquisitionUpdate();
        void prepareResultsAcquisition();

        // Tracking
        void initTracking();
        void runTracking();
        void runCorrelators();
        void runDiscriminatorsFilters();
        void runLoopIndicators();
        void postTrackingUpdate();
        void trackingStateUpdate();
        void runFrequencyDiscriminator();
        void runPhaseDiscriminator();
        void runCodeDiscriminator();
        void prepareResultsTracking();

        // Decoding
        void runDecoding();
        void decodeBit();
        void decodeSubframe();
        void postDecodingUpdate();
        void prepareResultsDecoding();

};

#endif

// ====================================================================================================================

template <typename D, typename S> std::complex<D> cast(const std::complex<S> s)
{
    return std::complex<D>(s.real(), s.imag());
}
