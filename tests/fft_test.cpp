#include <iostream>
#include "fft.h"
#include <complex>
#include "gen_code.h"
#include "utils.h"

using namespace std;

int main(){

    int size = 8;
    double signal[] = {2., -1., 31., 2., -15., -32., 6., 11., -7., -2., -10., 9., -23., 7., 23., 3.};

    // COMPLEX
    // FFT test
    cfft(signal, size);

    for(int idx=0; idx < size; idx++){
        cout << signal[2*idx] << " " << signal[2*idx+1] << endl;
    }

    cout << "-----------" << endl;

    // IFFT test
    cifft(signal, size);

    for(int idx=0; idx < size; idx++){
        cout << signal[2*idx] << " " << signal[2*idx+1] << endl;
    }

    cout << "-----------" << endl;

    // REAL
    // size = 32;
    // double signal_real[] = { 
    //     1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
    //     1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,
    //     1.,  1.,  1.,  1., -1., -1.};

    size = 10000;
    // generate code
    int code[1023];
    double upcode[10000];
    generateCAcode(9, code);
    upsampleCode(code, 1023, 1.023e6, 10e6, upcode);

    double* signal_real = upcode;

    rfft(signal_real, size);

    cout << signal_real[0] << endl;
    for(int idx=1; idx < size; idx+=2){
        cout << signal_real[idx] << " " << signal_real[idx+1] << endl;
    }

    for(int idx=1; idx < size; idx+=2){
        cout << conj(signal_real[idx] + 1i * signal_real[idx+1]) << endl;
    }
    
}
