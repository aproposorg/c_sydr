#include <iostream>
#include "fft.h"

using namespace std;

int main(){

    int size1ms = 16;
    double signal[size1ms] = {2, -1, 31, 2, -15, -32, 6, 11, -7, -2, -10, 9, -23, 7, 23, 3};

    // FFT test
    fft(signal, size1ms);

    for(int idx=0; idx < size1ms; idx++){
        cout << signal[2*idx] << " " << signal[2*idx+1] << endl;
    }

    cout << "-----------" << endl;

    // IFFT test
    ifft(signal, size1ms);

    for(int idx=0; idx < size1ms; idx++){
        cout << signal[2*idx] << " " << signal[2*idx+1] << endl;
    }
}
