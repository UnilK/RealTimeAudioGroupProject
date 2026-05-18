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

#include <vector>
#include <cstring>
#include <cassert>

namespace dsp {

template<typename T> class rbuffer {

public:

    rbuffer(){
        pointer = offset = bsize = 0;
    }

    rbuffer(int size){
        pointer = offset = bsize = 0;
        resize(size, (T)0);
    }

    rbuffer(int size, T val){
        pointer = offset = bsize = 0;
        resize(size, val);
    }
    
    void push(T val){
        if(pointer + bsize >= (int)buffer.size()){
            std::memcpy(buffer.data(), buffer.data() + pointer, sizeof(T) * bsize);
            pointer = 0;
        }
        buffer[pointer + bsize] = val;
        pointer++;
    }

    T& operator[](int i){
        return buffer[pointer + offset + i];
    }
    
    std::vector<T> slice(int left, int right) const {
        std::vector<T> s(right-left);
        int po = pointer + offset;
        for(int i=left; i<right; i++) s[i - left] = buffer[po + i];
        return s;
    }

    void set_offset(int o){ offset = o; }
    
    void resize(int size, T val){
        if(4 * size > (int)buffer.size()) buffer.resize(4 * size);
        if(pointer + bsize >= (int)buffer.size()){
            std::memcpy(buffer.data(), buffer.data() + pointer, sizeof(T) * bsize);
            pointer = 0;
        }
        for(int i=bsize; i<size; i++) buffer[pointer + i] = val;
        bsize = size;
    }

    int size() const { return bsize; }
    int left() const { return offset; }
    int right() const { return bsize - offset; }

private:

    int bsize, pointer, offset;
    std::vector<T> buffer;

};

}

