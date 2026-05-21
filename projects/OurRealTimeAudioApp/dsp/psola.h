/*
 * PSOLA (Pitch Synchronous Overlap-Add) Pitch Shifter
 * Copyright (C) 2024
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 */

#pragma once

#include <vector>
#include <cstddef>

namespace dsp {

class PSolaShifter {
public:
    PSolaShifter();
    ~PSolaShifter() = default;

    // Prepare the shifter with sample rate and max frame size
    void prepare(double sampleRate, int maxFrameSize = 4096);

    // Process a single sample with pitch shifting
    // pitchShift: factor to shift pitch (1.0 = no shift, 2.0 = up octave, 0.5 = down octave)
    // pitchPeriod: detected fundamental period in samples
    // isVoiced: whether the current frame is voiced
    float process(float inputSample, float pitchShift, int pitchPeriod, bool isVoiced);

    // Get the current analysis buffer center position
    const float* getAnalysisBuffer() const { return analysisBuffer.data(); }
    int getAnalysisBufferSize() const { return analysisBuffer.size(); }

private:
    // Internal buffers
    std::vector<float> analysisBuffer;      // Input analysis buffer
    std::vector<float> synthesisBuffer;     // Output synthesis buffer
    std::vector<float> windowedFrame;       // Windowed analysis frame
    std::vector<float> resampledFrame;      // Resampled frame for synthesis
    std::vector<float> hannWindow;          // Hann window for PSOLA

    // State variables
    int analysisBufferIndex = 0;
    int synthesisBufferIndex = 0;
    double synthesisPhase = 0.0;            // Phase accumulator for synthesis
    int lastPitchPeriod = 80;               // Default pitch period (~500Hz at 44.1kHz)
    float lastPitchShift = 1.0f;

    double sampleRate = 44100.0;
    int maxFrameSize = 4096;
    int anaBufferSize = 8192;

    int frameCounter { 0 };

    // Helper functions
    void fillAnalysisBuffer(float sample);
    void extractAnalysisFrame(int pitchPeriod);
    void resampleFrame(float pitchShift, int pitchPeriod);
    void applyHannWindow(std::vector<float>& frame, int size);
    void overlapAdd();
};

} // namespace dsp
