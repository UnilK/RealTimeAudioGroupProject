/*
 * Math and dsp utility functions
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

#include "commonMath.h"

#include "constants.h"

#include <cassert>
#include <cmath>
#include <random>

namespace math {

using std::vector;
using std::complex;

float rnd(const float &d){
    static std::mt19937 rng32;
    return std::uniform_real_distribution<float>(-d, d)(rng32);
}

int rnd_int(int min, int max){
    static std::mt19937 rng32;
    return std::uniform_int_distribution<int>(min, max)(rng32);
}

float sinc(float x){
    if(x*x > 1e-10f) return sin(x * PIF) / (x * PIF);
    return 1.0f;
}

complex<float> unit_complex(const complex<float> &c){
    return c / (std::abs(c) + 1e-18f);
}

vector<float> cos_window(int length){

    double a = 2 * PIF / length;

    vector<float> window(length);
    for(int i=0; i<length; i++) window[i] = 1 - std::cos(a * i);

    return window;
}

vector<complex<float> > cos_wavelet(int length, double spins){
    
    const double a = 2 * PIF / length;
    const double b = 2 * PIF * spins / length;
    vector<complex<float> > w(length);
    
    for(int i=0; i<length; i++) w[i] = std::polar<double>(1 - std::cos(a * i), b * i);

    return w;
}

vector<complex<float> > extract_phase(
        const vector<complex<float> > &frequency)
{
    vector<complex<float> > phase(frequency.size());
    for(unsigned i=0; i<phase.size(); i++) phase[i] = unit_complex(frequency[i]);
    return phase;
}

vector<float> extract_energy(
        const vector<complex<float> > &frequency)
{    
    vector<float> energy(frequency.size());
    for(unsigned i=0; i<energy.size(); i++) energy[i] = std::norm(frequency[i]);
    return energy;
}

vector<complex<float> > create_frequency(
        const vector<float> &energy,
        vector<complex<float> > phase)
{
    assert(energy.size() == phase.size());

    for(unsigned i=0; i<phase.size(); i++) phase[i] *= std::sqrt(energy[i]);
    return phase;
}

}

