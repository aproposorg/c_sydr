
#include "utils.h"

void upsampleCode(const char* code, size_t size, double codeFreq, double samplingFreq, double* r_code)
{
    size_t sizeUpCode = round(samplingFreq * size / codeFreq);

    for (size_t i {}; i < sizeUpCode; ++i)
        r_code[i] = code[(size_t) (i * codeFreq / samplingFreq)];

    return;
}
