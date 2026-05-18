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

#include "fft.h"

#include <algorithm>
#include <cassert>

namespace math {

using std::vector;
using std::complex;
using std::array;

void in_place_fft(vector<complex<float> > &v, bool inv){

    int n = v.size();

    in_place_fft(n, &v[0], inv);
}

vector<complex<float> > fft(const vector<float> &v){
    
    int n = v.size();
    int m = n/2 + 1;

    vector<complex<float> > w(n);
    for(int i=0; i<n; i++) w[i] = v[i];
    
    in_place_fft(w, 0);
    w.resize(m);
    
    return w;
}

array<vector<complex<float> >, 2> fft(const vector<float> &a, const vector<float> &b){
    
    assert(a.size() == b.size());

    int n = a.size();
    int m = n/2 + 1;

    vector<complex<float> > c(n), fa(m), fb(m);
    for(int i=0; i<n; i++) c[i] = {a[i], b[i]};

    in_place_fft(c, 0);

    for(int i=0; i<m; i++){
        int j = (n-i)&(n-1);
        fa[i] = 0.5f * (c[i] + conj(c[j]));
        fb[i] = 0.5f * (c[i] - conj(c[j]));
        fb[i] = {fb[i].imag(), -fb[i].real()};
    }

    return {fa, fb};
}

vector<float> inverse_fft(const vector<complex<float> > &v){
    
    int m = v.size();
    int n = 2 * m - 2;

    auto w = v;
    w.resize(n);
    
    for(int i=1; i<m; i++) w[n-i] = std::conj(w[i]);

    in_place_fft(w, 1);
    
    vector<float> r(n);
    for(int i=0; i<n; i++) r[i] = w[i].real();
    
    return r;
}

array<vector<float>, 2> inverse_fft(vector<complex<float> > a, vector<complex<float> > b){

    assert(a.size() == b.size());

    int m = a.size();
    int n = 2 * m - 2;
    
    a.resize(n);
    b.resize(n);
    
    for(int i=1; i<m; i++){
        a[n-i] = std::conj(a[i]);  
        b[n-i] = std::conj(b[i]);
    }

    for(int i=0; i<n; i++){
        a[i].real(a[i].real() - b[i].imag());
        a[i].imag(a[i].imag() + b[i].real());
    }

    in_place_fft(a, 1);

    vector<float> ta(n), tb(n);
    for(int i=0; i<n; i++){
        ta[i] = a[i].real();
        tb[i] = a[i].imag();
    }
    
    return {ta, tb};
}

void fft2(int n, float *x1, float *y1, float *x2, float *y2){
    in_place_fft(n, x1, x2, 0);

    for(int i=1, j=n-1; i<j; i++, j--){
        float ar = 0.5f * (x1[i] + x1[j]);
        float ai = 0.5f * (x2[i] - x2[j]);
        float br = 0.5f * (x1[i] - x1[j]);
        float bi = 0.5f * (x2[i] + x2[j]);
        x1[i] = x1[j] = ar;
        y1[i] = ai;
        y1[j] = -ai;
        x2[i] = x2[j] = bi;
        y2[i] = -br;
        y2[j] = br;
    }
}

void energy_mse(int n, int m, float *ix, float *iy, float *x, float *y){
    
    x[0] = ix[0];
    y[0] = iy[0];
    for(int i=1; i<n; i++){
        x[i] = ix[i];
        y[2*m-i] = iy[i];
    }

    for(int i=n; i<2*m; i++) x[i] = y[2*m-i] = 0.0f;

    in_place_fft(2*m, x, y, 0);
    
    x[0] = x[0] * y[0];
    x[m] = x[m] * y[m];
    y[0] = y[m] = 0.0f;

    for(int i=1, j=2*m-1; i<j; i++, j--){
        float ar = x[i] - x[j];
        float ai = y[i] + y[j];
        float br = x[i] + x[j];
        float bi = y[i] - y[j];
        x[j] = x[i] = 0.25f * (ar * bi + ai * br);
        y[i] = 0.25f * (ai * bi - ar * br);
        y[j] = -y[i];
    }

    in_place_fft(2*m, x, y, 1);

    double a2 = 0, b2 = 0;
    for(int i=n-1; i>=0; i--){
        a2 += ix[i]*ix[i];
        b2 += iy[n-1-i]*iy[n-1-i];
        x[i] = (a2 + b2 - 2.0 * x[i]) / (a2 + b2 + 1e-12f);
    } x[n] = 0.0f;

    std::reverse(x, x+n+1);
}

void energy_mse(int n, int m, float *ix, float *iy, complex<float> *c, float *x){

    c[0] = {ix[0], iy[0]};
    for(int i=1; i<n; i++){
        c[i].real(ix[i]);
        c[2*m-i].imag(iy[i]);
    }

    for(int i=n; i<2*m; i++){
        c[i].real(0.0f);
        c[2*m-i].imag(0.0f);
    }

    in_place_fft(2*m, c, 0);
    
    c[0] = c[0].real() * c[0].imag();
    c[m] = c[m].real() * c[m].imag();

    for(int i=1, j=2*m-1; i<j; i++, j--){
        float ar = c[i].real() - c[j].real();
        float ai = c[i].imag() + c[j].imag();
        float br = c[i].real() + c[j].real();
        float bi = c[i].imag() - c[j].imag();
        float abr = 0.25f * (ar * bi + ai * br);
        float abi = 0.25f * (ai * bi - ar * br);
        c[i] = {abr, abi};
        c[j] = {abr, -abi};
    }

    in_place_fft(2*m, c, 1);
    for(int i=0; i<=n; i++) x[i] = c[i].real();

    double a2 = 0, b2 = 0;
    for(int i=n-1; i>=0; i--){
        a2 += ix[i]*ix[i];
        b2 += iy[n-1-i]*iy[n-1-i];
        x[i] = (a2 + b2 - 2.0 * x[i]) / (a2 + b2 + 1e-12f);
    } x[n] = 0.0f;

    std::reverse(x, x+n+1);
}

vector<float> energy_mse(const vector<float> &a, const vector<float> &b){
    
    int n = a.size();
    
    assert(a.size() == b.size());
    
    int m = 1;
    while(m < n) m *= 2;

    vector<complex<float> > cc(2*m, 0.0f);
    for(int i=0; i<n; i++) cc[i].real(a[i]);
    cc[0].imag(b[0]);
    for(int i=1; i<n; i++) cc[2*m-i].imag(b[i]);

    in_place_fft(cc, 0);
    
    cc[0] = cc[0].real()*cc[0].imag();
    cc[m] = cc[m].real()*cc[m].imag();
    for(int i=1, j=2*m-1; i<j; i++, j--){
        cc[i] = (cc[i]-conj(cc[j]))*(cc[i]+conj(cc[j]));
        cc[i] = {cc[i].imag()*0.25f, -0.25f*cc[i].real()};
        cc[j] = conj(cc[i]);
    }

    in_place_fft(cc, 1);

    vector<float> c(n+1);
    for(int i=0; i<=n; i++) c[i] = cc[i].real();

    double a2 = 0, b2 = 0;
    for(int i=n-1; i>=0; i--){
        a2 += a[i]*a[i];
        b2 += b[n-1-i]*b[n-1-i];
        c[i] = std::max(1e-7, a2 + b2 - 2.0 * c[i]);
        c[i] /= std::max(1e-7, a2 + b2);
    } c[n] = 0.0f;
    
    reverse(c.begin(), c.end());

    return c;
}

vector<float> padded_energy_mse(const vector<float> &a, const vector<float> &b, int border){

    auto mse = energy_mse(a, b);
    
    for(int i=0; i+2*border < (int)mse.size(); i++) mse[i] = mse[i+2*border];
    mse.resize(mse.size() - 2 * border);

    return mse;
}

}
