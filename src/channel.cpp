
#include "channel.h"

using namespace std;

// ----------------------------------------------------------------------------
// Constructor / Destructor

Channel::Channel(int _channelID, st_ChannelConfig* _config){
    
    m_channelID = _channelID;
    m_config    = _config;
    m_channelState = IDLE;

    // Init
    m_dllTau1 = LoopFilterTau1(m_config->dllNoiseBW, m_config->dllDampRatio, m_config->dllLoopGain);
    m_dllTau2 = LoopFilterTau2(m_config->dllNoiseBW, m_config->dllDampRatio, m_config->dllLoopGain);
    m_pllTau1 = LoopFilterTau1(m_config->pllNoiseBW, m_config->pllDampRatio, m_config->pllLoopGain);
    m_pllTau2 = LoopFilterTau2(m_config->pllNoiseBW, m_config->pllDampRatio, m_config->pllLoopGain);

    m_trackRequiredSamples = int(1e-3 * m_config->signalConfig->samplingFreq);

    cout << "Channel created CID " << m_channelID << endl;
}

Channel::~Channel(){
    // Nothing to do yet
}

// ----------------------------------------------------------------------------
// General processing

void Channel::run(int* _rfdata, size_t size){

    m_rfdataSize = size;
    m_rfdata = _rfdata;

    // Proccess new data
    processHandler();

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
        runTracking();
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

    // Generate code
    generateCAcode(m_satelliteID, m_code+1);
    m_code[0] = m_code[GPS_L1CA_CODE_SIZE_BITS];   // Wrapping previous/last bits for tracking
    m_code[GPS_L1CA_CODE_SIZE_BITS+1] = m_code[1];

    // Initialise
    m_channelState = ACQUIRING;
    resetChannel();

    cout << "Channel initialised with PRN " << m_satelliteID << endl;
}

// ----------------------------------------------------------------------------

void Channel::resetChannel(){

    m_remainingCarrier = 0.0;
    m_remainingCode = 0.0;

    m_carrierFrequency = 0.0;
    m_codeOffset = 0;
    m_codeFrequency = GPS_L1CA_CODE_FREQ;

    // Tracking
    m_isTrackingInitialised = false;
    m_dllDiscrim = 0.0;
    m_pllDiscrim = 0.0;

    resetCounters();
}

// ----------------------------------------------------------------------------

void Channel::resetCounters(){

    m_iPromptSum = 0.0;
    m_qPromptSum = 0.0;
    m_nbPromptSum = 0;
    m_codeCounter = 0;
    m_pdpnRatio = 0.0;
}

// ----------------------------------------------------------------------------
// ACQUISITION METHODS

void Channel::runAcquisition(){

    int samplesPerCode = m_config->signalConfig->samplingFreq * GPS_L1CA_CODE_SIZE_BITS / GPS_L1CA_CODE_FREQ;

    // Check if sufficient data in buffer
    int requiredSamples = (int) (samplesPerCode * m_config->cohIntegration * m_config->nonCohIntegration);
    if (m_rfdataSize < requiredSamples){
        // Not enough samples
        return;
    }

    // Initialise arrays
    int sizeMap = (m_config->dopplerRange * 2 / m_config->dopplerStep + 1) * samplesPerCode;
    float acqCorrelationMap[sizeMap];
    memset(acqCorrelationMap, 0, sizeMap); // init to 0 

    // Perform signal search 
    runSignalSearch(acqCorrelationMap);

    // Perform peak finding
    runPeakFinder(acqCorrelationMap, sizeMap);

    // Post acquisition update
    postAcquisitionUpdate();
    
    return;
}

// ----------------------------------------------------------------------------

void Channel::runSignalSearch(float* r_correlation){

    // SerialSearch(
    //     m_rfdata, m_rfdataSize, m_code, m_config->dopplerRange, m_config->dopplerStep, 
    //     m_config->signalConfig->samplingFreq, r_correlation);

    // Recast array
    double _rfdata[m_rfdataSize];
    for(int i=0; i<m_rfdataSize; i++){
        _rfdata[i] = (double) m_rfdata[i];
    }

    PCPS(
        _rfdata, m_rfdataSize, m_code+1, GPS_L1CA_CODE_SIZE_BITS, GPS_L1CA_CODE_FREQ,
        m_config->cohIntegration, m_config->nonCohIntegration, m_config->signalConfig->samplingFreq, 
        0.0, m_config->dopplerRange, m_config->dopplerStep, r_correlation);

    return;
}

// ----------------------------------------------------------------------------

void Channel::runPeakFinder(float* acqCorrelationMap, size_t sizeMap){

    int samplesPerCode = m_config->signalConfig->samplingFreq * GPS_L1CA_CODE_SIZE_BITS / GPS_L1CA_CODE_FREQ;    
    int samplesPerCodeChip = ceil((float) samplesPerCode / GPS_L1CA_CODE_SIZE_BITS); // Code per chip round up to the next integer

    // Find the correlation
    TwoCorrelationPeakComparison(acqCorrelationMap, sizeMap, samplesPerCode, samplesPerCodeChip, &m_indexPeak, &m_acqMetric);

    return;
}

// ----------------------------------------------------------------------------

void Channel::postAcquisitionUpdate(){

    int samplesPerCode = m_config->signalConfig->samplingFreq * GPS_L1CA_CODE_SIZE_BITS / GPS_L1CA_CODE_FREQ;  

    // Check if peak is above threshold
    if(m_acqMetric < m_config->acqThreshold){
        // TODO Should count the number of try an only deactivate later
        m_channelState = IDLE;
    }
    else{
        // Unravel index
        int idxFreq = (int) floor(m_indexPeak / samplesPerCode);
        float dopplerShift = (-m_config->dopplerRange) + m_config->dopplerStep * idxFreq;

        m_carrierFrequency = m_config->signalConfig->intermediateFreq + dopplerShift;
        m_codeOffset = int(m_indexPeak % samplesPerCode) + 1;
        m_channelState = TRACKING;

        cout << "PRN " << m_satelliteID << " successfully acquired: " 
                << m_carrierFrequency << ", " << m_codeOffset << ", " << m_acqMetric << endl;
    }
    return;
}

// ----------------------------------------------------------------------------
// TRACKING METHODS

void Channel::runTracking(){

    // Initialise 
    initTracking();
    
    // Correlators
    runCorrelators();

    // Discriminators
    runDiscriminatorsFilters();

    // Loop indicators
    runLoopIndicators();

    // Post tracking update
    postTrackingUpdate();

    return;
}

void Channel::initTracking(){

    if (!m_isTrackingInitialised){
        m_trackRequiredSamples = m_rfdataSize / 2;
        m_isTrackingInitialised = true;
    }
    
}

// ----------------------------------------------------------------------------

void Channel::runCorrelators(){
    
    // Code step is declared as double as float caused differences with Python code
    // This would need to be evaluated to decrease to float.
    double codeStep = m_codeFrequency / m_config->signalConfig->samplingFreq;
    EPL(
        m_rfdata + (m_codeOffset*2), m_trackRequiredSamples, m_code, GPS_L1CA_CODE_SIZE_BITS+2, m_config->signalConfig->samplingFreq, 
        m_carrierFrequency, m_remainingCarrier, m_remainingCode, codeStep, m_config->correlatorSpacing,
        m_correlatorsResults);

    m_iPromptSum += m_correlatorsResults[I_PROMPT_IDX];
    m_qPromptSum += m_correlatorsResults[Q_PROMPT_IDX];
    m_nbPromptSum += 1;

    return;
}

// ----------------------------------------------------------------------------

void Channel::runDiscriminatorsFilters(){

    // DLL
    double dllDiscrim = DLL_NNEML(
        m_correlatorsResults[I_EARLY_IDX], m_correlatorsResults[Q_EARLY_IDX],
        m_correlatorsResults[I_LATE_IDX], m_correlatorsResults[Q_LATE_IDX]
    );
    m_codeError = BorreLoopFilter(dllDiscrim, m_dllDiscrim, m_dllTau1, m_dllTau2, m_config->dllPDI);
    m_dllDiscrim = dllDiscrim;

    // PLL 
    double pllDiscrim = PLL_costa(
        m_correlatorsResults[I_PROMPT_IDX],
        m_correlatorsResults[Q_PROMPT_IDX]
    );
    m_phaseError = BorreLoopFilter(pllDiscrim, m_pllDiscrim, m_pllTau1, m_pllTau2, m_config->pllPDI);
    m_pllDiscrim = pllDiscrim;

    return;
}

// ----------------------------------------------------------------------------

void Channel::runLoopIndicators(){

    // PLL
    m_pllLock = PLLIndicator(
        m_correlatorsResults[I_PROMPT_IDX], m_correlatorsResults[Q_PROMPT_IDX],
        m_pllLock, m_config->pllIndicatorAlpha
    );

    // CN0
    m_pdpnRatio += (m_correlatorsResults[I_PROMPT_IDX] * m_correlatorsResults[I_PROMPT_IDX] 
                  + m_correlatorsResults[Q_PROMPT_IDX] * m_correlatorsResults[Q_PROMPT_IDX])
                  / pow(abs(m_correlatorsResults[I_PROMPT_IDX]) - abs(m_correlatorsResults[Q_PROMPT_IDX]), 2);
    
    if(m_nbPromptSum == LNAV_MS_PER_BIT){
        m_cn0 = CN0_Baulieu(m_pdpnRatio, m_cn0, m_nbPromptSum, m_config->cn0Alpha);
        m_pdpnRatio = 0.0;
        m_nbPromptSum = 0;
    }

    return;
}


// ----------------------------------------------------------------------------

void Channel::postTrackingUpdate(){

    // Update NCO
    double sampFreq     = m_config->signalConfig->samplingFreq;
    double codeStep     = (m_codeFrequency / sampFreq);
    m_remainingCarrier -= m_carrierFrequency * TWO_PI * m_trackRequiredSamples / sampFreq;
    m_remainingCarrier  = fmod(m_remainingCarrier , TWO_PI);
    m_remainingCode    += m_trackRequiredSamples * codeStep - GPS_L1CA_CODE_SIZE_BITS;
    m_codeFrequency    -= m_codeError;
    m_carrierFrequency += m_phaseError;
    
    codeStep = (m_codeFrequency / sampFreq); // Update value
    m_codeOffset += m_trackRequiredSamples - int(1e-3 * sampFreq); // TODO this is to have the same as python but might be wrong...
    m_trackRequiredSamples = (int) ceil((GPS_L1CA_CODE_SIZE_BITS - m_remainingCode) / codeStep);

    m_codeCounter += 1;
    return;
}