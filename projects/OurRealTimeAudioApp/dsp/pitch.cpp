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

#include <cmath>
#include <cassert>

#include "dsp/pitch.h"
#include "math/commonMath.h"
#include "math/fft.h"

#include <iostream>

namespace dsp {

bool Ticker::tick(){
    state++;
    if(state >= period) state = 0;
    return state == 0;
}

PitchDetector::PitchDetector(const PitchDetectorConstructor& ctr){
    prepare(ctr);
}

void PitchDetector::prepare(const PitchDetectorConstructor& ctr){
    assert(ctr.framerate > 0.0f && ctr.framerate < 1e8f);
    assert(ctr.minPitchHz < ctr.maxPitchHz);

    framerate = ctr.framerate;
    pop = std::max<int>(8, std::ceil(ctr.framerate / (ctr.popWidth * 2)));
    min = std::floor(ctr.framerate / ctr.maxPitchHz);
    max = std::ceil(ctr.framerate / ctr.minPitchHz);
    min = std::max(8, std::min(min, pop));
    max = std::max(2 * pop, max);

    voicedLimit = ctr.voicedLimit;
    period = pop;
    stablePeriod = ctr.framerate / ctr.defaultPitch;
    pitch = ctr.framerate / pop;
    similarity = 0.0;
    clock = Ticker{(int)std::ceil(ctr.framerate / ctr.periodCalcFrequency), -1};
    
    n = max + 2 * pop;
    m = 1;
    while(m < n) m *= 2;
    x.resize(2*m);
    y.resize(2*m);
    ix.resize(n);
    iy.resize(n);
    mse = &x[2*pop];
}

int PitchDetector::get_reguired_buffer_radius(){
    return max + pop;
}

void PitchDetector::update_period(const float* bufferCenter){
    if(!clock.tick()) return;

    int l = -max-pop;
    int r = -pop;

    float avg = 0.0f;
    for(int i=-pop; i<pop; i++) avg += bufferCenter[i];
    avg /= 2*pop;

    for(int i=0; i<n; i++){
        const float noise = math::rnd(1e-4f); 
        ix[i] = bufferCenter[l+i] - avg + noise;
        iy[i] = bufferCenter[r+i] - avg + noise;
    }

    math::energy_mse(
            n, m,
            &ix[0], &iy[0],
            &x[0], &y[0]);

    float best = 1.0f;
    int top = pop;

    for(int i=min+1; i+1<=max; i++){
        if(mse[i] < best && mse[i] < mse[i-1] && mse[i] < mse[i+1]){
            top = i;
            best = mse[i];
        }
    }

    best = 1.0f - best;

    isVoiced = best > voicedLimit;
    if(isVoiced){
        float decay = std::pow(0.5f, (float)clock.period * halfRate);
        stablePeriod = stablePeriod * decay + top * (1.0f - decay);
    }

    int jumps = std::ceil(stablePeriod / (float)top);
    float a = (jumps - 1) * top;
    float b = jumps * top;
    float da = a / stablePeriod;
    float db = stablePeriod / b;

    if(b >= max) period = a;
    else period = da > db ? a : b;

    best = 1.0f;
    int jmin = std::max(min+1, period-5);
    int jmax = std::min(max-1, period+5);
    for(int i=jmin; i<=jmax; i++){
        if(mse[i] < best){
            period = i;
            best = mse[i];
        }
    }

    best = 1.0f - best;
    
    pitch = framerate / period;
    similarity = std::max(0.0f, std::min(best, 1.0f));
}

}