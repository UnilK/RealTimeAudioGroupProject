#include "PluginProcessor.h"
#include "PluginEditor.h"
#include <algorithm>
#include <cmath>

#include "math/constants.h"
#include "math/fft.h"

static void modeMix(float mode, float& lpf, float& bpf, float& hpf)
{
    mode = std::clamp(mode, -1.f, 1.f);
    lpf  = std::fmax(-mode, 0.f);
    bpf  = std::fmax(1.f - std::fabs(mode), 0.f);
    hpf  = std::fmax(mode, 0.f);
}

static const std::vector<mrta::ParameterInfo> ParameterInfos
{
    { Param::ID::PostGain,   Param::Name::PostGain,   "dB",    -6.0f, -60.f,  12.f,  0.1f, 3.8018f },
    { Param::ID::PitchShift, Param::Name::PitchShift, "semi",   7.0f, -24.0f, 24.0f, 1.0f, 1.0001f },
    { Param::ID::HarmonyMix, Param::Name::HarmonyMix, "%",     75.0f,   0.0f, 100.0f, 1.0f, 1.0f   },
    { Param::ID::Freq, Param::Name::Freq, Param::Unit::Hz, 500.0f, Param::Ranges::FreqMin, Param::Ranges::FreqMax, Param::Ranges::FreqInc, Param::Ranges::FreqSkw },
    { Param::ID::Reso, Param::Name::Reso, "",           1.0f, Param::Ranges::ResoMin, Param::Ranges::ResoMax, Param::Ranges::ResoInc, Param::Ranges::ResoSkw },
    { Param::ID::Mode, Param::Name::Mode, "",           0.0f, Param::Ranges::ModeMin, Param::Ranges::ModeMax, Param::Ranges::ModeInc, Param::Ranges::ModeSkw },
};

MainProcessor::MainProcessor() :
    mrta::BaseProcessor(ParameterInfos),
    pitchDetector({ .framerate = 44100.0f }),
    ibuff(1)
{
    math::init_fft(18);

    registerParameterCallback(Param::ID::PostGain,
        [this] (float value, bool /*forced*/)
        {
            gain = (value > -60.f) ? std::pow(10.f, value * 0.05f) : 0.f;
        });

    registerParameterCallback(Param::ID::PitchShift,
        [this] (float value, bool /*forced*/)
        {
            pitchShift = std::pow(2.0f, std::clamp(value, -24.0f, 24.0f) / 12.0f);
        });

    registerParameterCallback(Param::ID::HarmonyMix,
        [this] (float value, bool /*forced*/)
        {
            harmonyMix = value / 100.0f;
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

MainProcessor::~MainProcessor() {}

void MainProcessor::prepare(double sampleRate, int samplesPerBlock)
{
    pitchDetector.prepare({ .framerate = (float)sampleRate, .halfTime = 0.05f });
    psolaShifter.prepare(sampleRate, 4096);

    int radius = pitchDetector.get_reguired_buffer_radius() + 1;
    ibuff.resize(radius * 2, 0.0f);
    ibuff.set_offset(radius);

    svfLeft.prepare(sampleRate);
    svfRight.prepare(sampleRate);

    float lpf, bpf, hpf;
    modeMix(mode, lpf, bpf, hpf);
    lpfRamp.prepare(sampleRate, true, lpf);
    bpfRamp.prepare(sampleRate, true, bpf);
    hpfRamp.prepare(sampleRate, true, hpf);
    freqRamp.prepare(sampleRate, true, freqHz);
    resoRamp.prepare(sampleRate, true, reso);

    psolaBuffer.setSize(2, samplesPerBlock);

    freqInBuffer.setSize(1, samplesPerBlock);
    resoInBuffer.setSize(1, samplesPerBlock);
    lpfOutBuffer.setSize(2, samplesPerBlock);
    bpfOutBuffer.setSize(2, samplesPerBlock);
    hpfOutBuffer.setSize(2, samplesPerBlock);
}

void MainProcessor::process(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& /*midiMessages*/)
{

    /*
    juce::ScopedNoDenormals noDenormals;

    int n = buffer.getNumSamples();
    int m = std::min(buffer.getNumChannels(), 2);
    const float *x = buffer.getReadPointer(0);

    float ifs = 1.0f / (float)getSampleRate();
    static double phaseState = 0.0;

    for(int i=0; i<n; i++){
        ibuff.push(x[i]);
        pitchDetector.update_period(&ibuff[0]);

        if(pitchDetector.isVoiced){
            phaseState = std::fmod(phaseState + pitchDetector.pitch * ifs * 2 * PI, 2 * PI);
        }

        float sample = std::max(-1.0f, std::min<float>(1.0f, std::sin(phaseState) * gain));

        for(int j=0; j<m; j++) buffer.getWritePointer(j)[i] = sample;
    }
    */

    juce::ScopedNoDenormals noDenormals;
    const int    n = buffer.getNumSamples();
    const int    m = std::min(buffer.getNumChannels(), 2);
    const float* x = buffer.getReadPointer(0);

    psolaBuffer.clear();

    // Harmonizer
    for (int i = 0; i < n; ++i)
    {
        const float dry = x[i];

        ibuff.push(dry);
        pitchDetector.update_period(&ibuff[0]);

        const float harmony = psolaShifter.process(
            dry, pitchShift, pitchDetector.period, pitchDetector.isVoiced);

        // const float mixed = (dry + harmonyMix * harmony) / (1.0f + harmonyMix);
        const float out   = std::clamp(gain, -1.0f, 1.0f);

        for (int j = 0; j < m; ++j)
            psolaBuffer.getWritePointer(j)[i] = out;
    }

    // Filter
    freqInBuffer.clear();
    resoInBuffer.clear();
    lpfOutBuffer.clear();
    bpfOutBuffer.clear();
    hpfOutBuffer.clear();

    freqRamp.applySum(freqInBuffer.getWritePointer(0), n);
    resoRamp.applySum(resoInBuffer.getWritePointer(0), n);

    svfLeft.process(lpfOutBuffer.getWritePointer(0),
                bpfOutBuffer.getWritePointer(0),
                hpfOutBuffer.getWritePointer(0),
                psolaBuffer.getReadPointer(0),
                freqInBuffer.getReadPointer(0),
                resoInBuffer.getReadPointer(0),
                n);

    if (m > 1)
    {
        svfRight.process(lpfOutBuffer.getWritePointer(1),
                         bpfOutBuffer.getWritePointer(1),
                         hpfOutBuffer.getWritePointer(1),
                         psolaBuffer.getReadPointer(1),
                         freqInBuffer.getReadPointer(0),
                         resoInBuffer.getReadPointer(0),
                         n);
    }

    lpfRamp.applyGain(lpfOutBuffer.getArrayOfWritePointers(), m, n);
    bpfRamp.applyGain(bpfOutBuffer.getArrayOfWritePointers(), m, n);
    hpfRamp.applyGain(hpfOutBuffer.getArrayOfWritePointers(), m, n);
    psolaBuffer.clear();
    for (int ch = 0; ch < m; ++ch)
    {
        psolaBuffer.addFrom(ch, 0, lpfOutBuffer, ch, 0, n);
        psolaBuffer.addFrom(ch, 0, bpfOutBuffer, ch, 0, n);
        psolaBuffer.addFrom(ch, 0, hpfOutBuffer, ch, 0, n);
    }

    for (int i = 0; i < n; ++i)
    {
        const float dry = x[i];

        const float harmony = psolaBuffer.getWritePointer(0)[i];

        const float mixed = (dry + harmonyMix * harmony) / (1.0f + harmonyMix);
        const float out   = std::clamp(mixed * gain, -1.0f, 1.0f);

        for (int j = 0; j < m; ++j)
            buffer.getWritePointer(j)[i] = out;
    }
}

juce::AudioProcessorEditor* MainProcessor::createEditor()
{
    return new MainProcessorEditor(*this);
}

CREATE_PLUGIN(MainProcessor)
