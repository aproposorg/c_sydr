#include <iostream>
#include "channel.h"

using namespace std;

int main(){
    cout << "Hello World!" << endl;

    // Create channel
    int prn = 1;
    st_ChannelConfig config;
    Channel channel(prn, &config);

    return 0;
}