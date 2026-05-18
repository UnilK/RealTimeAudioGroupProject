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

#include "fft-impl.h"

#include <array>
#include <complex>
#include <vector>

namespace math {

void in_place_fft(std::vector<std::complex<float> > &v, bool inv = 0);

std::vector<std::complex<float> > fft(
        const std::vector<float> &v);

std::array<std::vector<std::complex<float> >, 2> fft(
        const std::vector<float> &a,
        const std::vector<float> &b);

std::vector<float> inverse_fft(
        const std::vector<std::complex<float> > &v);

std::array<std::vector<float>, 2> inverse_fft(
        std::vector<std::complex<float> > a,
        std::vector<std::complex<float> > b);

void fft2(int n, float *x1, float *y1, float *x2, float *y2);

void energy_mse(int n, int m, float *ix, float *iy, float *x, float *y);

void energy_mse(int n, int m, float *ix, float *iy, std::complex<float> *c, float *x);

std::vector<float> energy_mse(
        const std::vector<float> &a,
        const std::vector<float> &b);

std::vector<float> padded_energy_mse(
        const std::vector<float> &a,
        const std::vector<float> &b,
        int border);

}

