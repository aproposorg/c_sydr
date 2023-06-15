
#include "channel.h"

using namespace std;

// ----------------------------------------------------------------------------
// Constructor / Destructor

Channel::Channel(int _channelID, st_ChannelConfig* _config){
    
    m_channelID = _channelID;
    m_config    = _config;
    m_channelState = IDLE;

    m_rfdata = Eigen::MatrixXcd::Zero(1, m_config->bufferSize);

    cout << "Channel created CID " << m_channelID << endl;
}

Channel::~Channel(){
    // Nothing to do yet
}

// ----------------------------------------------------------------------------
// General processing

void Channel::run(complex<double>* _rfdata, size_t size){

    m_rfdataSize = size;

    // Copy data for matrix structure
    for(int idx=0; idx < size; idx++){
        m_rfdata(0, idx) = _rfdata[idx];
        //cout << m_rfdata(0, idx) << endl;
    }

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
    int sizeMap = (m_config->dopplerRange * 2 / m_config->dopplerStep + 1) * GPS_L1CA_CODE_SIZE_BITS;
    double acqCorrelationMap[sizeMap];
    memset(acqCorrelationMap, 0, sizeMap); // init to 0 

    // Perform signal search 
    runSignalSearch(acqCorrelationMap);

    // Perform peak finding
    runPeakFinder(acqCorrelationMap, sizeMap);
    
}

// ----------------------------------------------------------------------------

void Channel::runSignalSearch(double* r_correlation){

    SerialSearch(
        m_rfdata, m_rfdataSize, m_code, m_config->dopplerRange, m_config->dopplerStep, 
        m_config->signalConfig->samplingFreq, r_correlation);

    return;
}

// ----------------------------------------------------------------------------

void Channel::runPeakFinder(double* acqCorrelationMap, size_t sizeMap){

    int idxPeak = 0.0;
    float acqMetric = 0.0; 

    // Find the correlation
    cout << acqCorrelationMap[0] << endl;
    TwoCorrelationPeakComparison(acqCorrelationMap, sizeMap, &idxPeak, &acqMetric);

    // Check if peak is above threshold
    if(acqMetric < m_config->acqThreshold){
        // TODO Should count the number of try an only deactivate later
        m_channelState = IDLE;
    }
    else{
        // Unravel index
        int idxFreq = (int) floor(idxPeak / m_config->signalConfig->codeLengthBits);
        float dopplerShift = ((-m_config->dopplerRange) + m_config->dopplerStep * idxFreq);
        
        carrierFrequency = m_config->signalConfig->intermediateFreq + dopplerShift;
        codeOffset = (int) idxPeak % m_config->signalConfig->codeLengthBits;
    }
    return;
}