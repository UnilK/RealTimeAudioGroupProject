#include "StateVariableFilter.h"

#include <cmath>
#include <algorithm>

namespace dsp
{

StateVariableFilter::StateVariableFilter()
{
}

StateVariableFilter::~StateVariableFilter()
{
}

void StateVariableFilter::prepare(double newSampleRate)
{
    sampleRate = newSampleRate;
    state0 = 0.f;
    state1 = 0.f;
}

void StateVariableFilter::process(float* lpfOut, float* bpfOut, float* hpfOut,
                                  const float* audioIn, const float* freqIn, const float* resoIn,
                                  unsigned int numSamples)
{
    for (unsigned int n = 0; n < numSamples; ++n)
    {
        float twoR = 1.f / std::clamp(resoIn[n], 0.1f, 10.f);
        float g    = std::tan(static_cast<float>(M_PI / sampleRate) * std::clamp(freqIn[n], 20.f, 20000.f));
        float g0   = twoR + g;
        float d    = 1.f / (1.f + twoR * g + g * g);

        float hp = (audioIn[n] - state1 - g0 * state0) * d;
        float v0 = g * hp;
        float bp = v0 + state0;
        state0   = bp + v0;
        float v1 = g * bp;
        float lp = v1 + state1;
        state1   = lp + v1;

        lpfOut[n] = lp;
        bpfOut[n] = bp;
        hpfOut[n] = hp;
    }
}

}
