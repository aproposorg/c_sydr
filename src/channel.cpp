
#include "channel.h"

using namespace std;

// ----------------------------------------------------------------------------
// Constructor / Destructor

Channel::Channel(int _channelID, st_ChannelConfig* _config){
    
    m_channelID = _channelID;
    m_config    = _config;
    m_channelState = IDLE;

    cout << "Channel created CID " << m_channelID << endl;
}

Channel::~Channel(){
    // Nothing to do yet
}

// ----------------------------------------------------------------------------
// General processing

void Channel::run(double* _rfdata, size_t size){

    m_rfdataSize = size;
    m_rfdata = _rfdata;

    // Proccess new data
    processHandler();

    // Channel update
    //prepareChannelUpdate();

    return;
}

// ----------------------------------------------------------------------------

void Channel::processHandler(){

    // Select processing based on current state
    switch (m_channelState)
    {
    case OFF:
        // TODO Warning? 
        break;
    case IDLE:
        // TODO Warning
        break;
    case ACQUIRING:
        runAcquisition();
        break;
    case TRACKING:
        //runTracking();
        //runDecoding();
        break;

    default:
        // TODO Error
        break;
    }

    return;
}

// ----------------------------------------------------------------------------

void Channel::setSatellite(int satelliteID){
    m_satelliteID = satelliteID;
    generateCAcode(m_satelliteID, m_code);
    m_channelState = ACQUIRING;
    cout << "Channel initialised with PRN " << m_satelliteID << endl;
}

// ----------------------------------------------------------------------------
// ACQUISITION METHODS

void Channel::runAcquisition(){

    // Check if sufficient data in buffer
    if (m_rfdataSize < m_config->acqRequiredSamples){
        // Not enough samples
        return;
    }

    // Initialise arrays
    int samplesPerCode = m_config->signalConfig->samplingFreq * GPS_L1CA_CODE_SIZE_BITS / GPS_L1CA_CODE_FREQ;
    int sizeMap = (m_config->dopplerRange * 2 / m_config->dopplerStep + 1) * samplesPerCode;
    float acqCorrelationMap[sizeMap];
    memset(acqCorrelationMap, 0, sizeMap); // init to 0 

    // Perform signal search 
    runSignalSearch(acqCorrelationMap);

    // Perform peak finding
    runPeakFinder(acqCorrelationMap, sizeMap);
    
}

// ----------------------------------------------------------------------------

void Channel::runSignalSearch(float* r_correlation){

    // SerialSearch(
    //     m_rfdata, m_rfdataSize, m_code, m_config->dopplerRange, m_config->dopplerStep, 
    //     m_config->signalConfig->samplingFreq, r_correlation);

    PCPS(
        m_rfdata, m_rfdataSize, m_code, GPS_L1CA_CODE_SIZE_BITS, GPS_L1CA_CODE_FREQ,
        m_config->cohIntegration, m_config->nonCohIntegration, m_config->signalConfig->samplingFreq, 
        0.0, m_config->dopplerRange, m_config->dopplerStep, r_correlation);

    return;
}

// ----------------------------------------------------------------------------

void Channel::runPeakFinder(float* acqCorrelationMap, size_t sizeMap){

    int idxPeak = 0.0;
    float acqMetric = 0.0;

    int samplesPerCode = m_config->signalConfig->samplingFreq * GPS_L1CA_CODE_SIZE_BITS / GPS_L1CA_CODE_FREQ;    

    // Find the correlation
    TwoCorrelationPeakComparison(acqCorrelationMap, sizeMap, &idxPeak, &acqMetric);

    // Check if peak is above threshold
    if(acqMetric < m_config->acqThreshold){
        // TODO Should count the number of try an only deactivate later
        m_channelState = IDLE;
    }
    else{
        // Unravel index
        int idxFreq = (int) floor(idxPeak / samplesPerCode);
        float dopplerShift = -((-m_config->dopplerRange) + m_config->dopplerStep * idxFreq);
        
        m_carrierFrequency = m_config->signalConfig->intermediateFreq + dopplerShift;
        m_codeOffset = (int) idxPeak % samplesPerCode;

        cout << "Acquisition PRN " << m_satelliteID << ": " 
            << m_carrierFrequency << ", " << m_codeOffset << ", " << acqMetric << endl;
    }
    return;
}