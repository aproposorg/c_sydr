
#include "channel.h"

using namespace std;

// ----------------------------------------------------------------------------
// Constructor / Destructor

Channel::Channel(int _channelID, st_ChannelConfig* _config){
    
    m_channelID = _channelID;
    m_config    = _config;

    cout << "Channel created CID " << m_channelID << endl;
}

Channel::~Channel(){
    // Nothing to do yet
}

// ----------------------------------------------------------------------------
// General processing

// void Channel::run(complex<double>* _rfdata){

//     rfdata = _rfdata;

//     // Proccess new data
//     processHandler();

//     // Channel update
//     prepareChannelUpdate();

//     return;
// }

// ----------------------------------------------------------------------------

// void Channel::processHandler(){

//     // Select processing based on current state
//     switch (m_channelState)
//     {
//     case IDLE:
//         // TODO Warning
//         break;
//     case ACQUIRING:
//         runAcquisition();
//         break;
//     case TRACKING:
//         runTracking();
//         runDecoding();
//         break;
//     default:
//         // TODO Error
//         break;
//     }

//     return;
// }

// ----------------------------------------------------------------------------

void Channel::setSatellite(int satelliteID){
    m_satelliteID = satelliteID;
    generateCAcode(m_satelliteID, m_code);
    cout << "Channel initialised with PRN " << m_satelliteID << endl;
}

// ----------------------------------------------------------------------------
