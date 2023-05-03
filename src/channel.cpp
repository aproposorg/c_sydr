
#include "channel.h"

using namespace std;

Channel::Channel(int _channelID, st_ChannelConfig* _config){
    
    channelID = _channelID;
    config    = _config;

    cout << "Channel created for PRN " << _channelID << endl;
}
