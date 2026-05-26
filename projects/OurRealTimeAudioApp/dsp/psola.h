#pragma once
#include <vector>

namespace dsp {

class PSolaShifter {
public:
    void prepare(double sampleRate, int maxGrain = 4096);

    // pitchShift is a frequency ratio (1.0 = no change, 2.0 = +1 octave).
    float process(float input, float pitchShift, int period, bool voiced);

private:
    int bufLen   = 0;
    int maxGrain = 0;

    std::vector<float> inBuf;
    std::vector<float> outBuf;
    std::vector<float> grain;

    int    inPos       = 0;
    int    outPos      = 0;
    double sampleCount = 0.0;   // absolute read position, never wraps
    double synPos      = 0.0;   // absolute position of next grain centre
    bool   synInit     = false;
    int    markCount   = 0;

    void extractGrain(int grainSize);
    void placeGrain(int center, int grainSize);

    static constexpr int PERIOD_HIST_N = 5;
    int periodHist[PERIOD_HIST_N] = {0,0,0,0,0};
    int periodHistIdx  = 0;
    int periodLastSeen = 0;

    float voiceStrength    = 0.0f;
    float smoothPitchShift = 1.0f;
};

} // namespace dsp
