#include <iostream>
#include "channel.h"

using namespace std;

int main(){
    cout << "Hello World!" << endl;

    // Define signal parameters
    st_SignalConfig signalConfig;
    signalConfig.codeFreqBasis = 1023e6;
    signalConfig.codeLength = 1023;
    signalConfig.samplingFreq = 10e6;
    
    st_ChannelConfig channelConfig;
    channelConfig.signalConfig = &signalConfig;

    // Create channel
    int prn = 2;
    int channelID = 0;
    Channel channel(channelID, &channelConfig);
    channel.setSatellite(prn);

    return 0;
}