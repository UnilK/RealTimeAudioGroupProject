#include "PluginProcessor.h"
#include "PluginEditor.h"

#include <algorithm>
#include <cmath>

#include "math/constants.h"
#include "math/fft.h"


void modeMix(float mode, float& lpf, float& bpf, float& hpf)
{
    mode = std::clamp(mode, -1.f, 1.f);
    lpf = std::fmax(-mode, 0.f);
    bpf = std::fmax(1.f - std::fabs(mode), 0.f);
    hpf = std::fmax(mode, 0.f);
}

static const std::vector<mrta::ParameterInfo> ParameterInfos
{
    { Param::ID::PostGain,  Param::Name::PostGain,  "dB", 0.0f, -60.f, 12.f, 0.1f, 3.8018f },
    { Param::ID::Freq, Param::Name::Freq, Param::Unit::Hz, 500.0f, Param::Ranges::FreqMin, Param::Ranges::FreqMax, Param::Ranges::FreqInc, Param::Ranges::FreqSkw },
    { Param::ID::Reso, Param::Name::Reso, "", 1.0f, Param::Ranges::ResoMin, Param::Ranges::ResoMax, Param::Ranges::ResoInc, Param::Ranges::ResoSkw },
    { Param::ID::Mode, Param::Name::Mode, "", 0.0f, Param::Ranges::ModeMin, Param::Ranges::ModeMax, Param::Ranges::ModeInc, Param::Ranges::ModeSkw}
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
            DBG(Param::Name::PostGain + ": " + juce::String { value });
            float dbValue { 0.f };
            if (value > -60.f) dbValue = std::pow(10.f, value * 0.05f);
            gain = dbValue;

        });

    registerParameterCallback(Param::ID::Freq,
        [this] (float value, bool force)
        {
            freqHz = value;
            freqRamp.setTarget(value, force);
        });

    registerParameterCallback(Param::ID::Reso,
        [this] (float value, bool force)
        {
            reso = value;
            resoRamp.setTarget(value, force);
        });

    registerParameterCallback(Param::ID::Mode,
        [this] (float value, bool force)
        {
            mode = value;
            float lpf(0.f), bpf(0.f), hpf(0.f);
            modeMix(mode, lpf, bpf, hpf);
            lpfRamp.setTarget(lpf, force);
            bpfRamp.setTarget(bpf, force);
            hpfRamp.setTarget(hpf, force);
        });   
}

MainProcessor::~MainProcessor()
{
}

void MainProcessor::prepare(double sampleRate, int samplesPerBlock)
{   
    juce::uint32 numChannels { static_cast<juce::uint32>(std::max(getMainBusNumInputChannels(), getMainBusNumOutputChannels())) };
    
    pitchDetector.prepare({ .framerate = (float)sampleRate });


    int radius = pitchDetector.get_reguired_buffer_radius()+ 1;
    
    ibuff.resize(radius * 2, 0.0f);
    ibuff.set_offset(radius);


    svf.prepare(sampleRate);

    float lpf, bpf, hpf;
    modeMix(mode, lpf, bpf, hpf);
    lpfRamp.prepare(sampleRate, true, lpf);
    bpfRamp.prepare(sampleRate, true, bpf);
    hpfRamp.prepare(sampleRate, true, hpf);

    freqInBuffer.setSize(1, samplesPerBlock);
    resoInBuffer.setSize(1, samplesPerBlock);
    lpfOutBuffer.setSize(2, samplesPerBlock);
    bpfOutBuffer.setSize(2, samplesPerBlock);
    hpfOutBuffer.setSize(2, samplesPerBlock);
}

void MainProcessor::process(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& /*midiMessages*/)
{
    juce::ScopedNoDenormals noDenormals;

    int n = buffer.getNumSamples();
    int m = std::min(buffer.getNumChannels(), 2);
    const float *x = buffer.getReadPointer(0);

    float ifs = 1.0f / (float)getSampleRate();

    for(int i=0; i<n; i++){
        ibuff.push(x[i]);
        pitchDetector.update_period(&ibuff[0]);

        if(pitchDetector.isVoiced){
            phaseState = std::fmod(phaseState + pitchDetector.pitch * ifs * 2 * PI, 2 * PI);
        }

        float sample = std::max(-1.0f, std::min<float>(1.0f, std::sin(phaseState) * gain));

        for(int j=0; j<m; j++) buffer.getWritePointer(j)[i] = sample;
    }

    freqInBuffer.clear();
    resoInBuffer.clear();
    lpfOutBuffer.clear();
    bpfOutBuffer.clear();
    hpfOutBuffer.clear();

    freqRamp.applySum(freqInBuffer.getWritePointer(0), n);
    resoRamp.applySum(resoInBuffer.getWritePointer(0), n);

    svf.process(lpfOutBuffer.getWritePointer(0),
                    bpfOutBuffer.getWritePointer(0),
                    hpfOutBuffer.getWritePointer(0),
                    buffer.getReadPointer(0),
                    freqInBuffer.getReadPointer(0),
                    resoInBuffer.getReadPointer(0),
                    n);

    lpfRamp.applyGain(lpfOutBuffer.getArrayOfWritePointers(), m, n);
    bpfRamp.applyGain(bpfOutBuffer.getArrayOfWritePointers(), m, n);
    hpfRamp.applyGain(hpfOutBuffer.getArrayOfWritePointers(), m, n);
    buffer.clear();

    for (int ch = 0; ch < m; ++ch)
    {
        buffer.addFrom(ch, 0, lpfOutBuffer, ch, 0, n);
        buffer.addFrom(ch, 0, bpfOutBuffer, ch, 0, n);
        buffer.addFrom(ch, 0, hpfOutBuffer, ch, 0, n);
    }

}

juce::AudioProcessorEditor* MainProcessor::createEditor()
{
    return new MainProcessorEditor(*this);
}

CREATE_PLUGIN(MainProcessor)
