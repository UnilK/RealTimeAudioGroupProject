#pragma once

#include <juce_dsp/juce_dsp.h>

#include <BaseProcessor.h>
#include "dsp/pitch.h"
#include "dsp/psola.h"
#include "dsp/rbuffer.h"
#include "dsp/StateVariableFilter.h"
#include "dsp/Ramp.h"

namespace Param
{
    namespace ID
    {
        static const juce::String PostGain   { "post_gain" };
        static const juce::String PitchShift { "pitch_shift" };
        static const juce::String HarmonyMix { "harmony_mix" };
        static const juce::String Freq       { "freq" };
        static const juce::String Reso       { "reso" };
        static const juce::String Mode       { "mode" };
    }

    namespace Name
    {
        static const juce::String PostGain   { "Post-Gain" };
        static const juce::String PitchShift { "Pitch Shift" };
        static const juce::String HarmonyMix { "Harmony Mix" };
        static const juce::String Freq       { "Frequency" };
        static const juce::String Reso       { "Resonance" };
        static const juce::String Mode       { "Mode" };
    }

    namespace Unit
    {
        static const juce::String Hz { "Hz" };
    }

    namespace Ranges
    {
        static constexpr float FreqMin { 100.f };
        static constexpr float FreqMax { 10000.f };
        static constexpr float FreqInc { 1.f };
        static constexpr float FreqSkw { 0.4f };

        static constexpr float ResoMin { 0.5f };
        static constexpr float ResoMax { 5.f };
        static constexpr float ResoInc { 0.01f };
        static constexpr float ResoSkw { 0.4f };

        static constexpr float ModeMin { -1.f };
        static constexpr float ModeMax { 1.f };
        static constexpr float ModeInc { 0.01f };
        static constexpr float ModeSkw { 1.f };
    }
}

class MainProcessor final : public mrta::BaseProcessor
{
public:
    MainProcessor();
    ~MainProcessor() override;

    void prepare(double sampleRate, int samplesPerBlock) override;
    void process(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    juce::AudioProcessorEditor* createEditor() override;

private:
    // Harmonizer
    dsp::PitchDetector  pitchDetector;
    dsp::PSolaShifter   psolaShifter;
    dsp::rbuffer<float> ibuff;

    float gain       = 1.0f;
    float pitchShift = 1.0f;
    float harmonyMix = 0.75f;

    // Filter
    dsp::StateVariableFilter svf;

    float mode   { 0.f };
    float reso   { 1.0f };
    float freqHz { 500.f };

    dsp::Ramp<float> freqRamp;
    dsp::Ramp<float> resoRamp;
    dsp::Ramp<float> lpfRamp;
    dsp::Ramp<float> bpfRamp;
    dsp::Ramp<float> hpfRamp;

    juce::AudioBuffer<float> freqInBuffer;
    juce::AudioBuffer<float> resoInBuffer;
    juce::AudioBuffer<float> lpfOutBuffer;
    juce::AudioBuffer<float> bpfOutBuffer;
    juce::AudioBuffer<float> hpfOutBuffer;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainProcessor)
};
