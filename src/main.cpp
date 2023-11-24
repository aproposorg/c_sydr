#include <iostream>
#include <fstream>
#include "channel.h"

using namespace std;

int main()
{
    // Some constants for the test
    const int channelID {};
    const unsigned prn {9};
    const size_t NbMilliseconds {1000};

    // Create channel
    Channel channel {};
    channel.setSatellite(prn);
    channel.configure(channelID, DefaultChannelConfig);

    // For example purposes, recording the tracking results
    ofstream ofs("track.txt");

    // Play RF file, 10MHz sampling rate, 8 bit quantization, I/Q interleaved
    // Reading and recasting is very inefficient but this is just for debugging purposes
    ifstream ifs("./data/Novatel_20211130_resampled_10MHz_8bit_IQ_gain25.bin", ios::binary);
    if (!ifs) {
        cout << "ERROR: RF data file could not be found. Exiting." << endl;
        return -1;
    }

    // Create data buffers
    size_t size1ms {2 * 10000}, bufferIdx {};
    int rfdata[2 * size1ms];

    for (size_t i {}; i < NbMilliseconds; ++i)
    {
        // Read 1 millisecond
        char rfdata_int8[size1ms];
        ifs.read(rfdata_int8, size1ms * sizeof(char));

        // Recast to array
        for (size_t j {}; j < size1ms; ++j)
        {
            if (i > 1)
                rfdata[j] = rfdata[bufferIdx + j]; // Move previous data
            rfdata[bufferIdx + j] = (int)rfdata_int8[j]; // Copy new data
        }
        bufferIdx = size1ms;

        // Run channel 
        channel.run(rfdata, size1ms);

        // Write to file
        const double* const correlatorResults = channel.getCorrelatorResults();
        ofs << correlatorResults[0];
        for (size_t j {1}; j < 2 * NB_CORRELATORS; ++j)
            ofs << "," << correlatorResults[j];
        ofs << endl;
    }

    ifs.close();
    ofs.close();

    return 0;
}
