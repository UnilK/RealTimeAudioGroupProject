#include "dsp/psola.h"
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static constexpr int MIN_TRUSTED_PERIOD = 80;

namespace dsp {

void PSolaShifter::prepare(double /*sampleRate*/, int maxGrainSize)
{
    maxGrain = maxGrainSize;
    bufLen   = maxGrain * 8;

    inBuf .assign(bufLen, 0.0f);
    outBuf.assign(bufLen, 0.0f);
    grain .assign(maxGrain, 0.0f);

    inPos            = 0;
    outPos           = 0;
    sampleCount      = 0.0;
    synPos           = 0.0;
    synInit          = false;
    markCount        = 0;
    for (int i = 0; i < PERIOD_HIST_N; ++i) periodHist[i] = 0;
    periodHistIdx    = 0;
    periodLastSeen   = 0;
    voiceStrength    = 0.0f;
    smoothPitchShift = 1.0f;
}

static inline float hannWindow(int i, int grainSize)
{
    if (grainSize < 2) return 0.0f;
    return 0.5f * (1.0f - std::cos(2.0f * (float)M_PI * (float)i
                                       / (float)(grainSize - 1)));
}

void PSolaShifter::extractGrain(int grainSize)
{
    int start = inPos - grainSize;
    for (int i = 0; i < grainSize; ++i)
    {
        float w   = hannWindow(i, grainSize);
        int   idx = ((start + i) % bufLen + bufLen) % bufLen;
        grain[i]  = inBuf[idx] * w;
    }
}

void PSolaShifter::placeGrain(int center, int grainSize)
{
    int start = center - grainSize / 2;
    for (int i = 0; i < grainSize; ++i)
    {
        int idx = ((start + i) % bufLen + bufLen) % bufLen;
        outBuf[idx] += grain[i];
    }
}

float PSolaShifter::process(float input, float pitchShift, int period, bool voiced)
{
    const int rawPeriod = period;

    pitchShift = std::max(0.25f, std::min(pitchShift, 4.0f)); // Limit to +-2 octaves.
    smoothPitchShift += (pitchShift - smoothPitchShift) * 0.005f;
    pitchShift = smoothPitchShift;

    period = std::max(MIN_TRUSTED_PERIOD, std::min(period, maxGrain / 4)); // Limit period to a reasonable range.

    if (rawPeriod < MIN_TRUSTED_PERIOD) voiced = false;

    // Median filter on the detected period rejects octave errors.
    if (period != periodLastSeen && period > 0)
    {
        periodHist[periodHistIdx] = period;
        periodHistIdx = (periodHistIdx + 1) % PERIOD_HIST_N;
        periodLastSeen = period;
    }
    {
        int sorted[PERIOD_HIST_N];
        for (int i = 0; i < PERIOD_HIST_N; ++i)
            sorted[i] = (periodHist[i] > 0) ? periodHist[i] : period;
        std::sort(sorted, sorted + PERIOD_HIST_N);
        period = sorted[PERIOD_HIST_N / 2];
    }

    // grainSize = 2*hop  ->  50% Hann overlap
    const double hop       = (double)period / pitchShift; // hop is fractional 
    const int    grainSize = std::min(maxGrain, std::max(8, 2 * (int)std::ceil(hop)));

    inBuf[inPos] = input;
    inPos = (inPos + 1) % bufLen;

    // 
    if (!voiced) 
    {
        markCount = 0;
        synInit   = false;
        for (int i = 0; i < PERIOD_HIST_N; ++i) periodHist[i] = 0;
        periodLastSeen = 0;
    }
    else if (++markCount >= period)
    {
        markCount -= period;

        const double writeAhead = grainSize * 0.5;
        if (!synInit || synPos < sampleCount + writeAhead)
        {
            synPos  = sampleCount + writeAhead;
            synInit = true;
        }

        extractGrain(grainSize);

        while (synPos < sampleCount + writeAhead + (double)period)
        {
            int center = (int)(std::llround(synPos) % (long long)bufLen); 
            if (center < 0) center += bufLen;
            placeGrain(center, grainSize);
            synPos += hop;
        }
    }

    float out = outBuf[outPos];
    outBuf[outPos] = 0.0f;
    outPos = (outPos + 1) % bufLen;
    sampleCount += 1.0;

    // Smooth crossfade between dry (unvoiced) and harmony (voiced).
    const float targetVoice = voiced ? 1.0f : 0.0f;
    voiceStrength += (targetVoice - voiceStrength) * 0.005f;
    return /* voiceStrength * */ out; // + (1.0f - voiceStrength) * input;
}

} // namespace dsp
