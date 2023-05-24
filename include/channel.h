#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <ccomplex>
#include "gen_code.h"
#include "acquisition.h"
#include "structures.h"
#include <cstring>
#include <Eigen/Dense>

#define RESULTS_SIZE 32
#define PRN_SIZE 1024

/// ===========================================================================

using namespace std;

class Channel{

    private:
        Eigen::MatrixXcd m_rfdata;
        size_t m_rfdataSize;
    
    public:
        int m_channelID;
        int m_satelliteID;
        int m_channelState;
        double m_results[RESULTS_SIZE];
        double m_code[PRN_SIZE];
        st_ChannelConfig* m_config;

        int codeOffset;
        double carrierFrequency;
        
        // Constructor
        Channel(int, st_ChannelConfig*);

        // Destructor
        ~Channel();

        // General processing
        void run(complex<double>* _rfdata, size_t size);
        void processHandler();
        void setSatellite(int satelliteID);
        void getTimeSinceTOW();
        void prepareChannelUpdate();
        double* prepareResults();

        // Acquisition
        void runAcquisition();
        void runSignalSearch(double*);
        void runPeakFinder(double*, size_t);
        void postAcquisitionUpdate();
        void prepareResultsAcquisition();

        // Tracking
        void runTracking();
        void runCorrelators();
        void runDiscriminators();
        void runCarrierFrequencyFilter();
        void runCodeFrequencyFilter();
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
