#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <complex>
#include "gen_code.h"

#define RESULTS_SIZE 32
#define PRN_SIZE 1024

/// ===========================================================================

using namespace std;

class Channel{

    private:
        complex<double>* rfdata;
    
    public:
        int m_channelID;
        int m_satelliteID;
        int m_channelState;
        double m_results[RESULTS_SIZE];
        double m_code[PRN_SIZE];
        st_ChannelConfig* m_config;

        
        // Constructor
        Channel(int, st_ChannelConfig*);

        // Destructor
        ~Channel();

        // General processing
        void run(complex<double>* _rfdata);
        void processHandler();
        void setSatellite(int satelliteID);
        void getTimeSinceTOW();
        void prepareChannelUpdate();
        double* prepareResults();

        // Acquisition
        void runAcquisition();
        void runSignalSearch();
        void runPeakFinder();
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

