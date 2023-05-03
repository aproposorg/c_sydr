#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>

/// ===================================================================================================================

struct st_ChannelConfig {

};

/// ===================================================================================================================

class Channel{
    
    public:
        int channelID;
        st_ChannelConfig* config;

        Channel(int channelID, st_ChannelConfig* config);

};

#endif

