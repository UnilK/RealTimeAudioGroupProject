/*
 * Pitch detection algorithm
 * Copyright (C) 2024 Unto Karila
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
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
 * SOFTWARE.
*/

#pragma once

#include <vector>

#include "dsp/rbuffer.h"

namespace dsp {

struct Ticker {
    int period = 1, state = 0;
    bool tick();
};

struct PitchDetectorConstructor {
    float framerate = -1.0f;
    float minPitchHz = 80.0f;
    float maxPitchHz = 1200.0f;
    float popWidth = 200.0f;
    float periodCalcFrequency = 400.0f;
    float defaultPitch = 200.0f;
    float voicedLimit = 0.9f;
    float halfTime = 0.0225f;
};

class PitchDetector {

public:
    int n, m;
    std::vector<float> ix, iy, x, y, mse;

    int min, max, pop;
    float voicedLimit, framerate;

    Ticker clock;

    PitchDetector(const PitchDetectorConstructor& ctr);

    void prepare(const PitchDetectorConstructor& ctr);

    int get_reguired_buffer_radius();

    // run this for each sample of the input signal
    // The buffer is indexed from bufferCenter[-radius] to bufferCenter[radius-1]
    void update_period(const float* bufferCenter);

    int period, top;
    float similarity, pitch, stablePeriod, halfRate;
    bool isVoiced = false;
};

}