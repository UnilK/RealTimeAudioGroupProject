#pragma once

#include <cmath>

namespace dsp
{

template<typename F>
class Ramp
{
public:
    static constexpr F DefaultRampTime { static_cast<F>(0.05) };

    static_assert(std::is_floating_point<F>::value, "Only supports floating point!");

    Ramp(F rampTimeSec) :
        rampTime { std::fmax(rampTimeSec, minRampTime) }
    { }

    ~Ramp() { }

    Ramp() :
        rampTime { DefaultRampTime }
    { }

    Ramp(const Ramp&) = delete;
    const Ramp& operator=(const Ramp&) = delete;
    Ramp(Ramp&&) = delete;
    const Ramp& operator=(Ramp&&) = delete;

    void prepare(double newSampleRate, bool skipRamp = false, F skipRampToValue = static_cast<F>(0))
    {
        sampleRate = newSampleRate;
        if (skipRamp)
            setTarget(skipRampToValue, true);
        else
            setTarget(targetValue);
    }

    void setTarget(F newTargetValue, bool skipRamp = false)
    {
        if (std::abs(newTargetValue - currentValue) > minDelta)
        {
            targetValue = newTargetValue;
            rampStep = (targetValue - currentValue) / static_cast<F>(sampleRate * rampTime);
        }

        if (skipRamp)
            currentValue = targetValue = newTargetValue;
    }

    void setRampTime(F newRampTimeSec)
    {
        rampTime = std::fmax(newRampTimeSec, 0.f);
    }

    void applySum(F* buffers, unsigned int numChannels)
    {
        const F targetDelta { std::fabs(targetValue - currentValue) };
        if ((targetDelta > std::fabs(static_cast<F>(2) * rampStep)) && (std::fabs(rampStep) > minDelta))
            currentValue += rampStep;
        else
            currentValue = targetValue;

        for (unsigned int ch = 0; ch < numChannels; ++ch)
            buffers[ch] += currentValue;
    }

    void applySum(F* const* buffers, unsigned int numChannels, unsigned int numSamples)
    {
        for (unsigned int n = 0; n < numSamples; ++n)
        {
            const F targetDelta { std::fabs(targetValue - currentValue) };
            if ((targetDelta > std::fabs(static_cast<F>(2) * rampStep)) && (std::fabs(rampStep) > minDelta))
                currentValue += rampStep;
            else
                currentValue = targetValue;

            for (unsigned int ch = 0; ch < numChannels; ++ch)
                buffers[ch][n] += currentValue;
        }
    }

    void applySum(F* const* output, const F* const* input, unsigned int numChannels, unsigned int numSamples)
    {
        for (unsigned int n = 0; n < numSamples; ++n)
        {
            const F targetDelta { std::fabs(targetValue - currentValue) };
            if ((targetDelta > std::fabs(static_cast<F>(2) * rampStep)) && (std::fabs(rampStep) > minDelta))
                currentValue += rampStep;
            else
                currentValue = targetValue;

            for (unsigned int ch = 0; ch < numChannels; ++ch)
                output[ch][n] = currentValue + input[ch][n];
        }
    }

    void applyGain(F* buffers, unsigned int numChannels)
    {
        const F targetDelta { std::fabs(targetValue - currentValue) };
        if ((targetDelta > std::fabs(static_cast<F>(2) * rampStep)) && (std::fabs(rampStep) > minDelta))
            currentValue += rampStep;
        else
            currentValue = targetValue;

        for (unsigned int ch = 0; ch < numChannels; ++ch)
            buffers[ch] *= currentValue;
    }

    void applyGain(F* const* buffers, unsigned int numChannels, unsigned int numSamples)
    {
        for (unsigned int n = 0; n < numSamples; ++n)
        {
            const F targetDelta { std::fabs(targetValue - currentValue) };
            if ((targetDelta > std::fabs(static_cast<F>(2) * rampStep)) && (std::fabs(rampStep) > minDelta))
                currentValue += rampStep;
            else
                currentValue = targetValue;

            for (unsigned int ch = 0; ch < numChannels; ++ch)
                buffers[ch][n] *= currentValue;
        }
    }

    void applyGain(F* const* output, const F* const* input, unsigned int numChannels, unsigned int numSamples)
    {
        for (unsigned int n = 0; n < numSamples; ++n)
        {
            const F targetDelta{ std::fabs(targetValue - currentValue) };
            if ((targetDelta > std::fabs(static_cast<F>(2) * rampStep)) && (std::fabs(rampStep) > minDelta))
                currentValue += rampStep;
            else
                currentValue = targetValue;

            for (unsigned int ch = 0; ch < numChannels; ++ch)
                output[ch][n] = currentValue * input[ch][n];
        }
    }

    float getNext()
    {
        const F targetDelta { std::fabs(targetValue - currentValue) };
        if ((targetDelta > std::fabs(static_cast<F>(2) * rampStep)) && (std::fabs(rampStep) > minDelta))
            currentValue += rampStep;
        else
            currentValue = targetValue;

        return currentValue;
    }

    static constexpr F minRampTime { static_cast<F>(1e-3) };
    static constexpr F minDelta    { static_cast<F>(1e-9) };

private:
    double sampleRate { 48000.0 };
    F rampTime;
    F rampStep    { static_cast<F>(0) };
    F targetValue { static_cast<F>(0) };
    F currentValue{ static_cast<F>(0) };
};

}
