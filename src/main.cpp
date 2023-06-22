#include <iostream>
#include <fstream>
#include <vector>
#include "channel.h"
#include "fft.h"

using namespace std;

int main(){
    cout << "Hello World!" << endl;

    /// =======================================================================
    // CONFIGURATION 
    // Define signal parameters
    st_SignalConfig signalConfig;
    signalConfig.codeFreqBasis = GPS_L1CA_CODE_FREQ;
    signalConfig.codeLengthBits = GPS_L1CA_CODE_SIZE_BITS;
    signalConfig.samplingFreq = 10e6;
    signalConfig.intermediateFreq = 0.0;
    
    st_ChannelConfig channelConfig;
    channelConfig.signalConfig = &signalConfig;

    // Acquisition
    channelConfig.dopplerRange      = 5000;
    channelConfig.dopplerStep       = 250;
    channelConfig.acqThreshold      = 1.0;
    channelConfig.cohIntegration    = 1;
    channelConfig.nonCohIntegration = 1;

    // Tracking
    channelConfig.correlatorSpacing = 0.5;
    channelConfig.dllDampRatio      = 0.7;
    channelConfig.dllNoiseBW        = 1.0;
    channelConfig.dllLoopGain       = 1.0;
    channelConfig.dllPDI            = 1e-3;

    channelConfig.pllDampRatio      = 0.7;
    channelConfig.pllNoiseBW        = 15.0;
    channelConfig.pllLoopGain       = 0.25;
    channelConfig.pllPDI            = 1e-3;
    channelConfig.pllIndicatorAlpha = 0.01;

    channelConfig.cn0Alpha          = 0.1;

    /// =======================================================================
    // MAIN PROGRAM
    
    // Create channel
    int prn = 9; 
    int channelID = 0;
    Channel channel(channelID, &channelConfig);
    channel.setSatellite(prn);
    
    // For example purposes, recording the tracking results
    ofstream myfile;
    myfile.open ("track.txt");

    // Play RF file, 10MHz sampling rate, 8 bit quantization, I/Q interleaved
    // Reading and recasting is very inefficient but this is just for debugging purposes
    string filepath = "/mnt/c/Users/vmangr/Documents/Datasets/2021_11_30-TAU_Roof_Antenna_Tallysman/Novatel_20211130_resampled_10MHz_8bit_IQ_gain25.bin";
    ifstream ifs(filepath, ios::binary | ios::in);
    int size1ms = 10000 * 2;
    int8_t rfdata_int8[size1ms];
    int rfdata[size1ms * 2]; // Storing 2ms
    int bufferIdx = 0;
    int nbMilliseconds = 1000;
    for(int i=0; i < nbMilliseconds; i++){
        // Read 1 millisecond
        ifs.read(reinterpret_cast<char*>(rfdata_int8), size1ms*sizeof(int8_t));

        // Recast to array
        for(int idx=0; idx < size1ms; idx++){
            if(i > 1){
                rfdata[idx] = rfdata[bufferIdx + idx]; // Move previous data
            }
            rfdata[bufferIdx + idx] = (int) rfdata_int8[idx]; // Copy new data
        }
        bufferIdx = size1ms;

        // Run channel 
        channel.run(rfdata, size1ms);

        // Write to file
        myfile  << channel.m_correlatorsResults[0] << ","
                << channel.m_correlatorsResults[1] << ","
                << channel.m_correlatorsResults[2] << ","
                << channel.m_correlatorsResults[3] << ","
                << channel.m_correlatorsResults[4] << ","
                << channel.m_correlatorsResults[5]
                <<"\n";
    }
    ifs.close();
    myfile.close();

    return 0;
}