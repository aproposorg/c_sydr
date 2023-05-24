#include <iostream>
#include <fstream>
#include <vector>
#include "channel.h"
#include <Eigen/Dense>

using namespace std;
using Eigen::Matrix;

int main(){
    cout << "Hello World!" << endl;

    // Define signal parameters
    st_SignalConfig signalConfig;
    signalConfig.codeFreqBasis = GPS_L1CA_CODE_FREQ;
    signalConfig.codeLengthBits = GPS_L1CA_CODE_SIZE_BITS;
    signalConfig.samplingFreq = 10e6;
    signalConfig.intermediateFreq = 0.0;
    
    st_ChannelConfig channelConfig;
    channelConfig.bufferSize = 10000;
    channelConfig.dopplerRange = 5000;
    channelConfig.dopplerStep = 100;
    channelConfig.acqThreshold = 1.5;
    channelConfig.acqRequiredSamples = 10000;
    channelConfig.signalConfig = &signalConfig;

    // Create channel
    int prn = 9;
    int channelID = 0;
    Channel channel(channelID, &channelConfig);
    channel.setSatellite(prn);

    // Read 1 millisecond
    int size1ms = 10000;
    complex<double> rfdata[size1ms];
    string filepath = "/mnt/c/Users/vmangr/Documents/Datasets/2021_11_30-TAU_Roof_Antenna_Tallysman/Novatel_20211130_resampled_10MHz_8bit_IQ_gain25.bin";
    ifstream ifs(filepath, ios::binary | ios::in);
    vector<complex<int8_t>> v(size1ms);
    ifs.read(reinterpret_cast<char*>(v.data()), size1ms*2*sizeof(int8_t));
    ifs.close();

    // Recast to array
    for(int idx=0; idx < size1ms; idx++){
        rfdata[idx] = cast<double, int8_t>(v[idx]);
    }

    // Run channel 
    channel.run(rfdata, size1ms);

    return 0;
}