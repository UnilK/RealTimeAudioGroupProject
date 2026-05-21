#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <algorithm>
#include <cmath>
#include "math/constants.h"
#include "math/fft.h"

static const std::vector<mrta::ParameterInfo> ParameterInfos
{
    { Param::ID::PostGain,   Param::Name::PostGain,   "dB",   0.0f, -60.f,  12.f,  0.1f, 3.8018f },
    { Param::ID::PitchShift, Param::Name::PitchShift, "semi", 0.0f, -24.0f, 24.0f, 1.0f, 1.0001f },
};

MainProcessor::MainProcessor() :
    mrta::BaseProcessor(ParameterInfos),
    pitchDetector({ .framerate = 44100.0f }),
    ibuff(1)
{
    math::init_fft(18);

    registerParameterCallback(Param::ID::PostGain,
        [this] (float value, bool forced)
        {
            float dbValue { 0.f };
            if (value > -60.f) dbValue = std::pow(10.f, value * 0.05f);
            gain = dbValue;
        });

    registerParameterCallback(Param::ID::PitchShift,
        [this] (float value, bool forced)
        {
            DBG("raw callback value = " << value);
            if (!std::isfinite(value)) value = 0.0f;
            value = std::max(-24.0f, std::min(24.0f, value));
            pitchShift = std::pow(2.0f, value / 12.0f);
        });
}

MainProcessor::~MainProcessor()
{
}

void MainProcessor::prepare(double sampleRate, int samplesPerBlock)
{
    juce::uint32 numChannels { static_cast<juce::uint32>(std::max(getMainBusNumInputChannels(), getMainBusNumOutputChannels())) };
    pitchDetector.prepare({ .framerate = (float)sampleRate });
    // psolaShifter.prepare(sampleRate, samplesPerBlock);
    psolaShifter.prepare(sampleRate, 4096);
    int radius = pitchDetector.get_reguired_buffer_radius() + 1;
    ibuff.resize(radius * 2, 0.0f);
    ibuff.set_offset(radius);
}

void MainProcessor::process(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& /*midiMessages*/)
{
    juce::ScopedNoDenormals noDenormals;
    int n = buffer.getNumSamples();
    int m = std::min(buffer.getNumChannels(), 2);
    const float* x = buffer.getReadPointer(0);

    for (int i = 0; i < n; i++)
    {
        ibuff.push(x[i]);
        pitchDetector.update_period(&ibuff[0]);

        float shiftedSample = psolaShifter.process(
            x[i],
            pitchShift,
            pitchDetector.period,
            pitchDetector.isVoiced
        );

        float sample = std::max(-1.0f, std::min<float>(1.0f, shiftedSample * gain));
        for (int j = 0; j < m; j++) buffer.getWritePointer(j)[i] = sample;
    }
}

juce::AudioProcessorEditor* MainProcessor::createEditor()
{
    return new MainProcessorEditor(*this);
}

CREATE_PLUGIN(MainProcessor)