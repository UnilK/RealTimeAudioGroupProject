/*
 * PSOLA (Pitch Synchronous Overlap-Add) Pitch Shifter
 * Copyright (C) 2024
 */

#include "dsp/psola.h"
#include <cmath>
#include <algorithm>
#include <cassert>

namespace dsp {

PSolaShifter::PSolaShifter()
    : analysisBuffer(8192, 0.0f),
      synthesisBuffer(8192, 0.0f),
      windowedFrame(4096, 0.0f),
      resampledFrame(4096, 0.0f),
      hannWindow(4096, 0.0f)
{
}

void PSolaShifter::prepare(double sampleRate, int maxFrameSize)
{
    this->sampleRate = sampleRate;
    this->maxFrameSize = maxFrameSize;
    this->anaBufferSize = maxFrameSize * 2;

    analysisBuffer.assign(anaBufferSize, 0.0f);
    synthesisBuffer.assign(anaBufferSize, 0.0f);
    windowedFrame.assign(maxFrameSize, 0.0f);
    resampledFrame.assign(maxFrameSize, 0.0f);

    analysisBufferIndex = 0;
    synthesisBufferIndex = 0;
    synthesisPhase = 0.0;
    frameCounter = 0;  
}

void PSolaShifter::fillAnalysisBuffer(float sample)
{
    analysisBuffer[analysisBufferIndex] = sample;
    analysisBufferIndex = (analysisBufferIndex + 1) % anaBufferSize;
}

void PSolaShifter::extractAnalysisFrame(int pitchPeriod)
{
    int frameSize = std::min(pitchPeriod * 2, maxFrameSize);
    frameSize = std::max(frameSize, 64);

    int startIdx = (analysisBufferIndex - frameSize / 2 + anaBufferSize) % anaBufferSize;

    for (int i = 0; i < frameSize; ++i)
    {
        int idx = (startIdx + i) % anaBufferSize;
        float w = 0.5f * (1.0f - std::cos(2.0f * M_PI * i / (frameSize - 1)));
        windowedFrame[i] = analysisBuffer[idx] * w;
    }

    // zero out the rest
    for (int i = frameSize; i < maxFrameSize; ++i)
        windowedFrame[i] = 0.0f;
}

void PSolaShifter::resampleFrame(float pitchShift, int pitchPeriod)
{
    // Resample the windowed frame based on pitch shift
    int inputFrameSize = std::min(static_cast<int>(windowedFrame.size()), pitchPeriod * 2);
    inputFrameSize = std::max(inputFrameSize, 64);

    int outputFrameSize = std::max(1, static_cast<int>(inputFrameSize / pitchShift));
    outputFrameSize = std::min(outputFrameSize, maxFrameSize);

    // Linear interpolation resampling
    for (int i = 0; i < outputFrameSize; ++i) {
        float srcIdx = i * pitchShift;
        int idx = static_cast<int>(srcIdx);
        float frac = srcIdx - idx;

        if (idx + 1 < inputFrameSize) {
            resampledFrame[i] = windowedFrame[idx] * (1.0f - frac) +
                                windowedFrame[idx + 1] * frac;
        } else if (idx < inputFrameSize) {
            resampledFrame[i] = windowedFrame[idx];
        } else {
            resampledFrame[i] = 0.0f;
        }
    }

    // Zero out rest of buffer
    for (int i = outputFrameSize; i < static_cast<int>(resampledFrame.size()); ++i) {
        resampledFrame[i] = 0.0f;
    }

    lastPitchShift = pitchShift;
}

void PSolaShifter::overlapAdd()
{
    int inputFrameSize = std::min(lastPitchPeriod * 2, maxFrameSize);
    inputFrameSize = std::max(inputFrameSize, 64);

    int synthFrameSize = std::max(1, static_cast<int>(inputFrameSize / lastPitchShift));
    synthFrameSize = std::min(synthFrameSize, maxFrameSize);

    for (int i = 0; i < synthFrameSize; ++i)
        synthesisBuffer[(synthesisBufferIndex + i) % anaBufferSize] += resampledFrame[i];
}

float PSolaShifter::process(float inputSample, float pitchShift, int pitchPeriod, bool isVoiced)
{
    if (!std::isfinite(pitchShift)) pitchShift = 1.0f;
    pitchShift = std::max(0.5f, std::min(2.0f, pitchShift));
    pitchPeriod = std::max(20, std::min(pitchPeriod, maxFrameSize / 2)); // ← cap to half buffer

    fillAnalysisBuffer(inputSample);
    lastPitchPeriod = pitchPeriod;

    // removed local declaration — frameCounter is now a member
    frameCounter++;
    if (frameCounter >= pitchPeriod)
    {
        frameCounter = 0;
        if (pitchShift > 0.0f)
        {
            extractAnalysisFrame(pitchPeriod);
            resampleFrame(pitchShift, pitchPeriod);
            overlapAdd();
        }
    }

    float output = synthesisBuffer[synthesisBufferIndex];
    synthesisBuffer[synthesisBufferIndex] = 0.0f;
    synthesisBufferIndex = (synthesisBufferIndex + 1) % anaBufferSize;

    return output;
}
} // namespace dsp
