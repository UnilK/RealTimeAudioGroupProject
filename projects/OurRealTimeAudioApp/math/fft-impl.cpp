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

#include "fft-impl.h"

#include "constants.h"

#include <cassert>
#include <vector>

#if defined(_WIN32)

#include <intrin.h>

#else

#if defined(__SSE__)

#include <immintrin.h>

#elif (defined(__ARM_NEON))

#include <arm_neon.h>

#endif

#endif

namespace math {

using std::vector;
using std::complex;

#if (defined(__SSE__) || defined(_M_AMD64) || defined(_M_X64))

// x86 SSE ////////////////////////////////////////////////////////////////////////////////////////

static __m128 fmadd(const __m128 &a, const __m128 &b, const __m128 &c){
    return _mm_add_ps(_mm_mul_ps(a, b), c);
}

static __m128 fmsub(const __m128 &a, const __m128 &b, const __m128 &c){
    return _mm_sub_ps(_mm_mul_ps(a, b), c);
}

static void cmul(
        const __m128 &Xr, const __m128 &Xi,
        const __m128 &Yr, const __m128 &Yi,
        __m128 &Zr, __m128 &Zi)
{
    Zr = fmsub(Xr, Yr, _mm_mul_ps(Xi, Yi));
    Zi = fmadd(Xr, Yi, _mm_mul_ps(Xi, Yr));
}

namespace fft_precalc {
    static int N = -1;
    static vector<int> invbit;
    static vector<vector<complex<float> > > w;
    static vector<complex<float> > vc;

    static __m128 *wr4[64] = {0};
    static __m128 *wi4[64] = {0};
}

void init_fft(int exp2){

    using namespace fft_precalc;

    assert(exp2 < 30 && exp2 > 1);
    
    N = exp2;
    int M = exp2-1;


    w.resize(N);
    invbit = vector<int>(1<<N, 0);

    vc.resize(1<<N);

    w[M].resize(1<<M);
    for(int i=0; i<(1<<M); i++) w[M][i] = std::polar<double>(1.0, -PI * i / (1<<M));
    
    invbit[1] = 1;
    for(int b=M-1; b>=0; b--){

        int x = 1<<b, y = 1<<(M-b);
        w[b].resize(x);
        
        for(int i=0; i<y; i++){
            invbit[i] <<= 1;
            invbit[i+y] = invbit[i] | 1;
        }

        for(int i=0; i<x; i++) w[b][i] = w[b+1][2*i];
    }

    for(int b=2; b<N; b++){

        if(!wr4[b]){
            wr4[b] = new __m128[1<<(b-2)];
            wi4[b] = new __m128[1<<(b-2)];
        }

        for(int i=0; i<(1<<b); i+=4){

            wr4[b][i>>2] = _mm_set_ps(
                    w[b][i+3].real(),
                    w[b][i+2].real(),
                    w[b][i+1].real(),
                    w[b][i].real());
            
            wi4[b][i>>2] = _mm_set_ps(
                    w[b][i+3].imag(),
                    w[b][i+2].imag(),
                    w[b][i+1].imag(),
                    w[b][i].imag());
        }
    }
}

static void small_fft(int n, float *x, float *y, bool inv){

    using namespace fft_precalc;

    int bits = 0;
    while(1ll<<bits < n) bits++;
    assert(1ll<<bits == n && bits <= N);
    if(bits == 0) return;

    int shift = N - bits;

    for(int i=0; i<n; i++){
        if(i < invbit[i]>>shift){
            std::swap(x[i], x[invbit[i]>>shift]);
            std::swap(y[i], y[invbit[i]>>shift]);
        }
    }

    for(int r=0, rd=1; r<bits; r++, rd*=2){
        for(int i=0; i<n; i+=2*rd){
            for(int d=0; d<rd; d++){
                int j = i+d;
                const int jj = j+rd;
                const auto &c = w[r][d];
                float wr = c.real() * x[jj] - c.imag() * y[jj];
                float wi = c.real() * y[jj] + c.imag() * x[jj];
                x[jj] = x[j] - wr;
                y[jj] = y[j] - wi;
                x[j] = x[j] + wr;
                y[j] = y[j] + wi;
            }
        }
    }
    
    if(!inv) return;

    float invn = 1.0 / n;
    x[0] *= invn;
    y[0] *= invn;
    for(int i=1; i<=n/2; i++){
        float tmpx = x[i];
        x[i] = x[n-i] * invn;
        x[n-i] = tmpx * invn;
        float tmpy = y[i];
        y[i] = y[n-i] * invn;
        y[n-i] = tmpy * invn;
    }
}

void in_place_fft(int n, float *x, float *y, bool inv){

    using namespace fft_precalc;

    int bits = 0;
    while(1ll<<bits < n) bits++;
    assert(1ll<<bits == n && bits <= N);
    
    if(bits < 4){
        small_fft(n, x, y, inv);
        return;
    }

    int shift = N - bits;

    for(int i=0; i<n; i++){
        if(i < invbit[i]>>shift){
            std::swap(x[i], x[invbit[i]>>shift]);
            std::swap(y[i], y[invbit[i]>>shift]);
        }
    }

    {
        __m128 m0 = _mm_set_ps(-1.0f, 1.0, -1.0f, 1.0f);
        __m128 m1 = _mm_set_ps(-1.0f, -1.0, 1.0f, 1.0f);
        __m128 m2 = _mm_set_ps(1.0f, -1.0, -1.0f, 1.0f);

        for(int i=0; i<n; i+=4){
            
            __m128 x0 = _mm_loadu_ps(&x[i]);
            __m128 y0 = _mm_loadu_ps(&y[i]);
            
            __m128 x1 = _mm_shuffle_ps(x0, x0, 0b11110101);
            __m128 y1 = _mm_shuffle_ps(y0, y0, 0b11110101);
            x0 = _mm_shuffle_ps(x0, x0, 0b10100000);
            y0 = _mm_shuffle_ps(y0, y0, 0b10100000);

            x0 = fmadd(x1, m0, x0);
            y0 = fmadd(y1, m0, y0);

            __m128 xy23 = _mm_shuffle_ps(x0, y0, 0b11101110);
            x0 = _mm_shuffle_ps(x0, x0, 0b01000100);
            y0 = _mm_shuffle_ps(y0, y0, 0b01000100);
            x1 = _mm_shuffle_ps(xy23, xy23, 0b11001100);
            y1 = _mm_shuffle_ps(xy23, xy23, 0b01100110);

            x0 = fmadd(x1, m1, x0);
            y0 = fmadd(y1, m2, y0);

            _mm_storeu_ps(&x[i], x0);
            _mm_storeu_ps(&y[i], y0);
        }
    }
    
    for(int r=2, rd=4; r<bits; r++, rd*=2){
        for(int i=0; i<n; i+=2*rd){
            for(int d=0; 4*d<rd; d++){
                
                int j = i+4*d;
                int jj = j+rd;
                
                __m128 x0 = _mm_loadu_ps(&x[j]);
                __m128 y0 = _mm_loadu_ps(&y[j]);
                __m128 x1 = _mm_loadu_ps(&x[jj]);
                __m128 y1 = _mm_loadu_ps(&y[jj]);
                __m128 wr = wr4[r][d];
                __m128 wi = wi4[r][d];
                
                __m128 wvr, wvi;
                cmul(x1, y1, wr, wi, wvr, wvi);

                x1 = _mm_sub_ps(x0, wvr);
                y1 = _mm_sub_ps(y0, wvi);
                x0 = _mm_add_ps(x0, wvr);
                y0 = _mm_add_ps(y0, wvi);

                _mm_storeu_ps(&x[j], x0);
                _mm_storeu_ps(&y[j], y0);
                _mm_storeu_ps(&x[jj], x1);
                _mm_storeu_ps(&y[jj], y1);

            }
        }
    }
    
    if(!inv) return;

    float invn = 1.0 / n;
    __m128 invn4 = _mm_set1_ps(invn);
    
    x[0] *= invn;
    y[0] *= invn;

    int i = 1;
    for(; i+4<n/2; i+=4){
        
        __m128 xl = _mm_loadu_ps(&x[i]); 
        __m128 xr = _mm_loadu_ps(&x[n-i-3]); 
        __m128 yl = _mm_loadu_ps(&y[i]); 
        __m128 yr = _mm_loadu_ps(&y[n-i-3]); 
        
        xl = _mm_shuffle_ps(xl, xl, 0b00011011);
        xr = _mm_shuffle_ps(xr, xr, 0b00011011);
        yl = _mm_shuffle_ps(yl, yl, 0b00011011);
        yr = _mm_shuffle_ps(yr, yr, 0b00011011);
        
        xl = _mm_mul_ps(xl, invn4);
        xr = _mm_mul_ps(xr, invn4);
        yl = _mm_mul_ps(yl, invn4);
        yr = _mm_mul_ps(yr, invn4);

        _mm_storeu_ps(&x[i], xr);
        _mm_storeu_ps(&x[n-i-3], xl);
        _mm_storeu_ps(&y[i], yr);
        _mm_storeu_ps(&y[n-i-3], yl);
    }
    
    for(; i<=n/2; i++){
        float tmpx = x[i];
        x[i] = x[n-i] * invn;
        x[n-i] = tmpx * invn;
        float tmpy = y[i];
        y[i] = y[n-i] * invn;
        y[n-i] = tmpy * invn;
    }
}

#elif defined(__ARM_NEON)

static float32x4_t fmadd(const float32x4_t &a, const float32x4_t &b, const float32x4_t &c){
    return vfmaq_f32(c, a, b);
}

static float32x4_t fmsub(const float32x4_t &a, const float32x4_t &b, const float32x4_t &c){
    return vfmsq_f32(c, a, b);
}

static float32x4_t set(float x3, float x2, float x1, float x0){
    return {x0, x1, x2, x3};
}

static float32x4_t loadu(const float *x){
    return vld1q_f32(x);
}

static void storeu(float *y, float32x4_t val){
    vst1q_f32(y, val);
}

static void cmul(
        const float32x4_t &Xr, const float32x4_t &Xi,
        const float32x4_t &Yr, const float32x4_t &Yi,
        float32x4_t &Zr, float32x4_t &Zi)
{
    Zr = fmsub(Xi, Yi, vmulq_f32(Xr, Yr));
    Zi = fmadd(Xr, Yi, vmulq_f32(Xi, Yr));
}

namespace fft_precalc {
    static int N = -1;
    static vector<int> invbit;
    static vector<vector<complex<float> > > w;
    static vector<complex<float> > vc;

    static float32x4_t *wr4[64] = {0};
    static float32x4_t *wi4[64] = {0};
}

void init_fft(int exp2){

    using namespace fft_precalc;

    assert(exp2 < 30 && exp2 > 1);
    
    N = exp2;
    int M = exp2-1;


    w.resize(N);
    invbit = vector<int>(1<<N, 0);

    vc.resize(1<<N);

    w[M].resize(1<<M);
    for(int i=0; i<(1<<M); i++) w[M][i] = std::polar<double>(1.0, -PI * i / (1<<M));
    
    invbit[1] = 1;
    for(int b=M-1; b>=0; b--){

        int x = 1<<b, y = 1<<(M-b);
        w[b].resize(x);
        
        for(int i=0; i<y; i++){
            invbit[i] <<= 1;
            invbit[i+y] = invbit[i] | 1;
        }

        for(int i=0; i<x; i++) w[b][i] = w[b+1][2*i];
    }

    for(int b=2; b<N; b++){

        if(!wr4[b]){
            wr4[b] = new float32x4_t[1<<(b-2)];
            wi4[b] = new float32x4_t[1<<(b-2)];
        }

        for(int i=0; i<(1<<b); i+=4){

            wr4[b][i>>2] = set(
                    w[b][i+3].real(),
                    w[b][i+2].real(),
                    w[b][i+1].real(),
                    w[b][i].real());
            
            wi4[b][i>>2] = set(
                    w[b][i+3].imag(),
                    w[b][i+2].imag(),
                    w[b][i+1].imag(),
                    w[b][i].imag());
        }
    }
}

void in_place_fft(int n, float *x, float *y, bool inv){

    using namespace fft_precalc;

    int bits = 0;
    while(1ll<<bits < n) bits++;
    assert(1ll<<bits == n && bits <= N);

    int shift = N - bits;

    for(int i=0; i<n; i++){
        if(i < invbit[i]>>shift){
            std::swap(x[i], x[invbit[i]>>shift]);
            std::swap(y[i], y[invbit[i]>>shift]);
        }
    }

    for(int r=0, rd=1; r<bits && r<2; r++, rd*=2){
        for(int i=0; i<n; i+=2*rd){
            for(int d=0; d<rd; d++){
                int j = i+d;
                const int jj = j+rd;
                const auto &c = w[r][d];
                float wr = c.real() * x[jj] - c.imag() * y[jj];
                float wi = c.real() * y[jj] + c.imag() * x[jj];
                x[jj] = x[j] - wr;
                y[jj] = y[j] - wi;
                x[j] = x[j] + wr;
                y[j] = y[j] + wi;
            }
        }
    }
      
    for(int r=2, rd=4; r<bits; r++, rd*=2){
        for(int i=0; i<n; i+=2*rd){
            for(int d=0; 4*d<rd; d++){
                
                int j = i+4*d;
                int jj = j+rd;
                
                float32x4_t x0 = loadu(&x[j]);
                float32x4_t y0 = loadu(&y[j]);
                float32x4_t x1 = loadu(&x[jj]);
                float32x4_t y1 = loadu(&y[jj]);
                float32x4_t wr = wr4[r][d];
                float32x4_t wi = wi4[r][d];
                
                float32x4_t wvr, wvi;
                cmul(x1, y1, wr, wi, wvr, wvi);

                x1 = vsubq_f32(x0, wvr);
                y1 = vsubq_f32(y0, wvi);
                x0 = vaddq_f32(x0, wvr);
                y0 = vaddq_f32(y0, wvi);

                storeu(&x[j], x0);
                storeu(&y[j], y0);
                storeu(&x[jj], x1);
                storeu(&y[jj], y1);

            }
        }
    }
    
    if(!inv) return;
  
    float invn = 1.0 / n;
    x[0] *= invn;
    y[0] *= invn;
    for(int i=1; i<=n/2; i++){
        float tmpx = x[i];
        x[i] = x[n-i] * invn;
        x[n-i] = tmpx * invn;
        float tmpy = y[i];
        y[i] = y[n-i] * invn;
        y[n-i] = tmpy * invn;
    }
}

#else

// DEFAULT ////////////////////////////////////////////////////////////////////////////////////////

namespace fft_precalc {
    static int N = -1;
    static vector<int> invbit;
    static vector<vector<complex<float> > > w;
}

void init_fft(int exp2){

    using namespace fft_precalc;

    assert(exp2 < 30 && exp2 > 1);
    
    N = exp2;
    int M = exp2-1;

    w.resize(N);
    invbit = vector<int>(1<<N, 0);

    w[M].resize(2<<M);
    for(int i=0; i<(2<<M); i++) w[M][i] = std::polar<double>(1.0, -PI * i / (1<<M));
    
    invbit[1] = 1;
    for(int b=M-1; b>=0; b--){

        int x = 2<<b, y = 1<<(M-b);
        w[b].resize(x);
        
        for(int i=0; i<y; i++){
            invbit[i] <<= 1;
            invbit[i+y] = invbit[i] | 1;
        }

        for(int i=0; i<x; i++) w[b][i] = w[b+1][2*i];
    }
}

void in_place_fft(int n, float *x, float *y, bool inv){

    using namespace fft_precalc;

    int bits = 0;
    while(1ll<<bits < n) bits++;
    assert(1ll<<bits == n && bits <= N);
    if(bits == 0) return;

    int shift = N - bits;

    for(int i=0; i<n; i++){
        if(i < invbit[i]>>shift){
            std::swap(x[i], x[invbit[i]>>shift]);
            std::swap(y[i], y[invbit[i]>>shift]);
        }
    }

    for(int r=0, rd=1; r<bits; r++, rd*=2){
        for(int i=0; i<n; i+=2*rd){
            for(int d=0; d<rd; d++){
                int j = i+d;
                const int jj = j+rd;
                const auto &c = w[r][d];
                float wr = c.real() * x[jj] - c.imag() * y[jj];
                float wi = c.real() * y[jj] + c.imag() * x[jj];
                x[jj] = x[j] - wr;
                y[jj] = y[j] - wi;
                x[j] = x[j] + wr;
                y[j] = y[j] + wi;
            }
        }
    }
    
    if(!inv) return;

    float invn = 1.0 / n;
    x[0] *= invn;
    y[0] *= invn;
    for(int i=1; i<=n/2; i++){
        float tmpx = x[i];
        x[i] = x[n-i] * invn;
        x[n-i] = tmpx * invn;
        float tmpy = y[i];
        y[i] = y[n-i] * invn;
        y[n-i] = tmpy * invn;
    }
}

#endif

void in_place_fft(int n, complex<float> *v, bool inv){

    float *x = new float[n];
    float *y = new float[n];

    for(int i=0; i<n; i++){
        x[i] = v[i].real();
        y[i] = v[i].imag();
    }

    in_place_fft(n, x, y, inv);

    for(int i=0; i<n; i++) v[i] = {x[i], y[i]};

    delete[] x;
    delete[] y;
}

}
