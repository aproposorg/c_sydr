#ifndef STRUCTURES_H
#define STRUCTURES_H
/// ===========================================================================

struct st_SignalConfig{
    double samplingFreq;
    double codeFreqBasis;
    double codeLength;
};


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

struct st_ChannelConfig {

    st_SignalConfig* signalConfig;

};

/// ===========================================================================

#endif
