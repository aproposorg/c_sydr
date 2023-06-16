
#include "fft.h"

void cfft(double* signal, const size_t size)
{
    cfft_plan plan = make_cfft_plan (size);
    cfft_forward(plan, signal, 1.);
    destroy_cfft_plan (plan);
    return;
}

void cifft(double* signal, const size_t size)
{
    cfft_plan plan = make_cfft_plan (size);
    cfft_backward(plan, signal, 1./size);
    destroy_cfft_plan (plan);
    return;
}

void rfft(double* signal, const size_t size)
{
    rfft_plan plan = make_rfft_plan (size);
    rfft_forward(plan, signal, 1.);
    destroy_rfft_plan (plan);
    return;
}

void rifft(double* signal, const size_t size)
{
    rfft_plan plan = make_rfft_plan (size);
    rfft_backward(plan, signal, 1./size);
    destroy_rfft_plan (plan);
    return;
}
