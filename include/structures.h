#ifndef STRUCTURES_H
#define STRUCTURES_H

/// ===========================================================================

typedef enum {
    OFF,
    IDLE,
    ACQUIRING,
    TRACKING
} ChannelState;

// ----------------------------------------------------------------------------

typedef enum {
    NONE,
    CHANNEL_UPDATE,
    ACQUISITION_UPDATE,
    TRACKING_UPDATE,
    DECODING_UPDATE
} ChannelMessage;

/// ===========================================================================

struct st_SignalConfig{
    //RF 
    double samplingFreq;
    double intermediateFreq;

    // GNSS
    double codeFreqBasis;
    int codeLengthBits;
};

/// ===========================================================================

struct st_ChannelConfig {

    // General
    st_SignalConfig* signalConfig; // RF and GNSS signal parameters
    int bufferSize; 

    // Acquisition
    int dopplerRange;
    int dopplerStep;
    int cohIntegration;   // Coherent integration
    int nonCohIntegration; // Non coherent integration
    float acqThreshold;

};

/// ===========================================================================

#endif
