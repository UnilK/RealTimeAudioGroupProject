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

#pragma once

#include <complex>
#include <vector>

namespace math {

float rnd(const float &d);

int rnd_int(int min, int max);

float sinc(float x);

std::complex<float> unit_complex(const std::complex<float> &c);

std::vector<float> cos_window(int length);

std::vector<std::complex<float> > cos_wavelet(int length, double spins);

std::vector<std::complex<float> > extract_phase(
        const std::vector<std::complex<float> > &frequency);

std::vector<float> extract_energy(
        const std::vector<std::complex<float> > &frequency);

std::vector<std::complex<float> > create_frequency(
        const std::vector<float> &energy,
        std::vector<std::complex<float> > phase);

}

