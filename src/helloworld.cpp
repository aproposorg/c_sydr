#include <iostream>
#include "channel.h"

using namespace std;

int main(){
    cout << "Hello World!" << endl;

    // Create channel
    int prn = 1;
    int channelID = 0;
    st_ChannelConfig config;
    Channel channel(channelID, &config);
    channel.setSatellite(prn);

    return 0;
}