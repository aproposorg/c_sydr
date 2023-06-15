
#include "pocketfft.h"

void fft(double* signal, const size_t size)
{
    cfft_plan plan = make_cfft_plan (size);
    cfft_forward(plan, signal, 1.);
    destroy_cfft_plan (plan);
    return;
}


void ifft(double* signal, const size_t size)
{
    cfft_plan plan = make_cfft_plan (size);
    cfft_backward(plan, signal, 1./size);
    destroy_cfft_plan (plan);
    return;
}