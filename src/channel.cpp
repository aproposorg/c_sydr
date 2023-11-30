
#include "channel.h"

using namespace std;

// ----------------------------------------------------------------------------
// Constructor / Destructor

Channel::Channel()
{
    m_channelState = IDLE;

    cout << "Empty channel created" << endl;
}

Channel::Channel(int _channelID, st_ChannelConfig _config)
{    
    m_channelState = IDLE;

    configure(_channelID, _config);

    cout << "Channel created with CID " << m_channelID << endl;
}

Channel::~Channel()
{
    // Nothing to do yet
}

// ----------------------------------------------------------------------------
// Configuration

void Channel::setSatellite(unsigned satelliteID)
{
    m_satelliteID = satelliteID;

    m_channelState = ACQUIRING;
    resetChannel();

    cout << "Channel initialised with PRN " << m_satelliteID << endl;
}

void Channel::configure(int _channelID, st_ChannelConfig _config)
{
    m_channelID = _channelID;
    m_config    = _config;

    // Init
    m_dllTau1 = LoopFilterTau1(m_config.dllNoiseBW, m_config.dllDampRatio, m_config.dllLoopGain);
    m_dllTau2 = LoopFilterTau2(m_config.dllNoiseBW, m_config.dllDampRatio);
    m_pllTau1 = LoopFilterTau1(m_config.pllNoiseBW, m_config.pllDampRatio, m_config.pllLoopGain);
    m_pllTau2 = LoopFilterTau2(m_config.pllNoiseBW, m_config.pllDampRatio);

    m_trackRequiredSamples = 1e-3 * m_config.signalConfig.samplingFreq;

    // Set configured flag
    m_configured = true;

    cout << "Configured channel CID " << m_channelID << endl;
}

void Channel::resetChannel()
{
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

void Channel::resetCounters()
{
    m_iPromptSum = 0.0;
    m_qPromptSum = 0.0;
    m_nbPromptSum = 0;
    m_pdpnRatio = 0.0;
}

// ----------------------------------------------------------------------------
// General processing

void Channel::run(const int* rfdata, size_t size)
{
    assert(m_configured);

    // Select processing based on current state
    switch (m_channelState)
    {
    case OFF: // TODO Warning?
    case IDLE:
        break;
    case ACQUIRING:
        runAcquisition(rfdata, size);
        break;
    case TRACKING:
        runTracking(rfdata, size);
        break;
    default: // TODO Error?
        break;
    }

    return;
}

// ----------------------------------------------------------------------------
// ACQUISITION METHODS

void Channel::runAcquisition(const int* rfdata, size_t size)
{
    size_t samplesPerCode = m_config.signalConfig.samplingFreq * m_config.signalConfig.codeLengthBits / m_config.signalConfig.codeFreqBasis;

    // Check if sufficient data in buffer
    size_t requiredSamples = samplesPerCode * m_config.cohIntegration * m_config.nonCohIntegration;
    if (size < requiredSamples) return;

    // Initialise arrays
    size_t sizeMap = (m_config.dopplerRange * 2 / m_config.dopplerStep + 1) * samplesPerCode;
    float acqCorrelationMap[sizeMap];

    // Perform signal search 
    runSignalSearch(rfdata, size, acqCorrelationMap);

    // Perform peak finding
    runPeakFinder(acqCorrelationMap, sizeMap);

    // Post acquisition update
    postAcquisitionUpdate();
    
    return;
}

void Channel::runNoMapAcquisition(const int* rfdata, size_t size)
{
    // Fetch the pre-generated code
    const char *code = GPS_L1CA_CODES[m_satelliteID];

    // Perform signal search and peak finding in one
    NoMapPCPS(
        rfdata, size, code + 1, m_config.signalConfig.codeLengthBits, m_config.signalConfig.codeFreqBasis,
        m_config.cohIntegration, m_config.nonCohIntegration, m_config.signalConfig.samplingFreq,
        0.0, m_config.dopplerRange, m_config.dopplerStep, &m_indexPeak, &m_acqMetric);

    // Post acquisition update
    postAcquisitionUpdate();

    return;
}

// ----------------------------------------------------------------------------

void Channel::runSignalSearch(const int* rfdata, size_t size, float* r_correlation)
{
    // SerialSearch(
    //     m_rfdata, m_rfdataSize, m_code, m_config.dopplerRange, m_config.dopplerStep, 
    //     m_config.signalConfig.samplingFreq, r_correlation);

    // Fetch the pre-generated code
    const char* code = GPS_L1CA_CODES[m_satelliteID];

    PCPS(
        rfdata, size, code + 1, m_config.signalConfig.codeLengthBits, m_config.signalConfig.codeFreqBasis,
        m_config.cohIntegration, m_config.nonCohIntegration, m_config.signalConfig.samplingFreq,
        0.0, m_config.dopplerRange, m_config.dopplerStep, r_correlation);

    return;
}

// ----------------------------------------------------------------------------

void Channel::runPeakFinder(const float* acqCorrelationMap, size_t sizeMap)
{
    size_t samplesPerCode = m_config.signalConfig.samplingFreq * m_config.signalConfig.codeLengthBits / m_config.signalConfig.codeFreqBasis;
    size_t samplesPerCodeChip = ceil((float) samplesPerCode / m_config.signalConfig.codeLengthBits); // Code per chip round up to the next integer

    // Find the correlation
    TwoCorrelationPeakComparison(acqCorrelationMap, sizeMap, samplesPerCode, samplesPerCodeChip, &m_indexPeak, &m_acqMetric);

    return;
}

// ----------------------------------------------------------------------------

void Channel::postAcquisitionUpdate()
{
    size_t samplesPerCode = m_config.signalConfig.samplingFreq * m_config.signalConfig.codeLengthBits / m_config.signalConfig.codeFreqBasis;

    // Check if peak is above threshold
    if (m_acqMetric < m_config.acqThreshold) {
        // TODO Should count the number of try an only deactivate later
        m_channelState = IDLE;
    } else {
        // Unravel index
        size_t idxFreq = m_indexPeak / samplesPerCode;
        size_t dopplerShift = (-m_config.dopplerRange) + m_config.dopplerStep * idxFreq;

        m_carrierFrequency = m_config.signalConfig.intermediateFreq + dopplerShift;
        m_codeOffset = (m_indexPeak % samplesPerCode) + 1;
        m_channelState = TRACKING;

        cout << "PRN " << m_satelliteID << " successfully acquired: " 
                << m_carrierFrequency << ", " << m_codeOffset << ", " << m_acqMetric << endl;
    }

    return;
}

// ----------------------------------------------------------------------------
// TRACKING METHODS

void Channel::runTracking(const int* rfdata, size_t size)
{
    // Initialise 
    initTracking(size);

    // Correlators
    runCorrelators(rfdata);

    // Discriminators
    runDiscriminatorsFilters();

    // Loop indicators
    runLoopIndicators();

    // Post tracking update
    postTrackingUpdate();

    return;
}

void Channel::initTracking(size_t size)
{
    if (!m_isTrackingInitialised) {
        m_trackRequiredSamples = size / 2;
        m_isTrackingInitialised = true;
    }
}

// ----------------------------------------------------------------------------

void Channel::runCorrelators(const int* rfdata)
{
    // Fetch the pre-generated code
    const char* code = GPS_L1CA_CODES[m_satelliteID];

    // Code step is declared as double as float caused differences with Python code
    // This would need to be evaluated to decrease to float.
    double codeStep = m_codeFrequency / m_config.signalConfig.samplingFreq;
    EPL(
        rfdata + (m_codeOffset*2), m_trackRequiredSamples, code, m_config.signalConfig.codeLengthBits+2, m_config.signalConfig.samplingFreq, 
        m_carrierFrequency, m_remainingCarrier, m_remainingCode, codeStep, m_config.correlatorSpacing,
        m_correlatorsResults);

    m_iPromptSum += m_correlatorsResults[I_PROMPT_IDX];
    m_qPromptSum += m_correlatorsResults[Q_PROMPT_IDX];
    m_nbPromptSum += 1;

    return;
}

// ----------------------------------------------------------------------------

void Channel::runDiscriminatorsFilters()
{
    // DLL
    double dllDiscrim = DLL_NNEML(
        m_correlatorsResults[I_EARLY_IDX], m_correlatorsResults[Q_EARLY_IDX],
        m_correlatorsResults[I_LATE_IDX], m_correlatorsResults[Q_LATE_IDX]
    );
    m_codeError = BorreLoopFilter(dllDiscrim, m_dllDiscrim, m_dllTau1, m_dllTau2, m_config.dllPDI);
    m_dllDiscrim = dllDiscrim;

    // PLL 
    double pllDiscrim = PLL_costa(
        m_correlatorsResults[I_PROMPT_IDX],
        m_correlatorsResults[Q_PROMPT_IDX]
    );
    m_phaseError = BorreLoopFilter(pllDiscrim, m_pllDiscrim, m_pllTau1, m_pllTau2, m_config.pllPDI);
    m_pllDiscrim = pllDiscrim;

    return;
}

// ----------------------------------------------------------------------------

void Channel::runLoopIndicators()
{
    // PLL
    m_pllLock = PLLIndicator(
        m_correlatorsResults[I_PROMPT_IDX], m_correlatorsResults[Q_PROMPT_IDX],
        m_pllLock, m_config.pllIndicatorAlpha
    );

    // CN0
    m_pdpnRatio += (m_correlatorsResults[I_PROMPT_IDX] * m_correlatorsResults[I_PROMPT_IDX] 
                  + m_correlatorsResults[Q_PROMPT_IDX] * m_correlatorsResults[Q_PROMPT_IDX])
                  / pow(abs(m_correlatorsResults[I_PROMPT_IDX]) - abs(m_correlatorsResults[Q_PROMPT_IDX]), 2);

    if (m_nbPromptSum == LNAV_MS_PER_BIT) {
        m_cn0 = CN0_Baulieu(m_pdpnRatio, m_cn0, m_nbPromptSum, m_config.cn0Alpha);
        m_pdpnRatio = 0.0;
        m_nbPromptSum = 0;
    }

    return;
}

// ----------------------------------------------------------------------------

void Channel::postTrackingUpdate()
{
    // Update NCO
    double sampFreq     = m_config.signalConfig.samplingFreq;
    double codeStep     = (m_codeFrequency / sampFreq);
    m_remainingCarrier -= m_carrierFrequency * TWO_PI * m_trackRequiredSamples / sampFreq;
    m_remainingCarrier  = fmod(m_remainingCarrier , TWO_PI);
    m_remainingCode    += m_trackRequiredSamples * codeStep - m_config.signalConfig.codeLengthBits;
    m_codeFrequency    -= m_codeError;
    m_carrierFrequency += m_phaseError;
    
    codeStep = (m_codeFrequency / sampFreq); // Update value
    m_codeOffset += m_trackRequiredSamples - (int)(1e-3 * sampFreq); // TODO this is to have the same as python but might be wrong...
    m_trackRequiredSamples = (int) ceil((m_config.signalConfig.codeLengthBits - m_remainingCode) / codeStep);

    return;
}
