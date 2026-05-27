#include "PluginEditor.h"
#include <cassert>

// Width of the whole GUI
static constexpr int WIDTH { 250 };
static constexpr int PITCH_WIDTH { 750 };

// Height of each paramter knob on the paramEditor
static const int PARAM_HEIGHT { 100 };

MainProcessorEditor::MainProcessorEditor(MainProcessor& p) :
    juce::AudioProcessorEditor(p),
    pluginProcessor { p },
    paramEditor(pluginProcessor.getParameterManager(), PARAM_HEIGHT)
{
    addAndMakeVisible(paramEditor);

    // Calculate window height based on number of parameters
    const auto height { pluginProcessor.getParameterManager().getParameters().size() * PARAM_HEIGHT };
    setSize(WIDTH + PITCH_WIDTH, height);
}

MainProcessorEditor::~MainProcessorEditor()
{
}

void MainProcessorEditor::paint(juce::Graphics &g)
{
    auto &pd = pluginProcessor.pitchDetector;

    // static std::vector<float> ix, iy, x;
    // static const float* mse = nullptr;

    // {
    //     int n = pd.get_reguired_buffer_radius();
    //     float mi = 0.0f;
    //     for(int i=0; i<n; i++) mi = std::min(mi, pd.mse[i]);
    //     if(mi > -0.1f){
    //         repaint();
    //         return;
    //     }
    // }
    
    g.fillAll(juce::Colour(0xFF201918));

    {
        int n = pd.ix.size();
        const auto height { pluginProcessor.getParameterManager().getParameters().size() * PARAM_HEIGHT };
        float wi = (float)PITCH_WIDTH / n;
        float hi = height * 0.5;

        juce::Path path;
        for(int i=1; i<n; i++){
            path.addLineSegment(juce::Line<float>(
                (i-1) * wi + WIDTH,
                pd.ix[i-1] * hi + height/2,
                i * wi + WIDTH, 
                pd.ix[i] * hi + height/2), 2);
        }

        g.setColour(juce::Colour(0xFF3B7ABD));
        g.fillPath(path);
    }

    {
        int n = pd.iy.size();
        const auto height { pluginProcessor.getParameterManager().getParameters().size() * PARAM_HEIGHT };
        float wi = (float)PITCH_WIDTH / n;
        float hi = height * 0.5;

        juce::Path path;
        for(int i=1; i<n; i++){
            path.addLineSegment(juce::Line<float>(
                (i-1) * wi + WIDTH,
                pd.iy[i-1] * hi + height/4*3,
                i * wi + WIDTH, 
                pd.iy[i] * hi + height/4*3), 1);
        }

        g.setColour(juce::Colour(0xFFCC321F));
        g.fillPath(path);
    }


    {
        int n = pd.get_reguired_buffer_radius();
        const auto height { pluginProcessor.getParameterManager().getParameters().size() * PARAM_HEIGHT };
        float wi = (float)PITCH_WIDTH / n;
        float hi = height * 0.9;

        juce::Path path;
        for(int i=1; i<n; i++){
            path.addLineSegment(juce::Line<float>(
                (i-1) * wi + WIDTH,
                height * 0.95 - (1.0f - pd.mse[i-1]) * hi,
                i * wi + WIDTH, 
                height * 0.95 - (1.0f - pd.mse[i]) * hi), 2);
        }

        g.setColour(juce::Colour(0xFFF7E13A));
        g.fillPath(path);
    }

    {
        int n = pd.get_reguired_buffer_radius();
        const auto height { pluginProcessor.getParameterManager().getParameters().size() * PARAM_HEIGHT };
        float wi = (float)PITCH_WIDTH / n;
        float r = 7;
        g.setColour(juce::Colour(0xFF7DD625));
        g.fillEllipse(WIDTH + pd.period * wi-r, height*0.05-r, 2*r, 2*r);
    }

    {
        int n = pd.get_reguired_buffer_radius();
        const auto height { pluginProcessor.getParameterManager().getParameters().size() * PARAM_HEIGHT };
        float wi = (float)PITCH_WIDTH / n;
        float r = 7;
        g.setColour(juce::Colour(0xFFE14444));
        g.fillEllipse(WIDTH + pd.stablePeriod * wi - r, height*0.05-r + 2*r, 2*r, 2*r);
    }
    
    {
        int n = pd.get_reguired_buffer_radius();
        const auto height { pluginProcessor.getParameterManager().getParameters().size() * PARAM_HEIGHT };
        float wi = (float)PITCH_WIDTH / n;
        float r = 7;
        g.setColour(juce::Colour(0xFF3B7ABD));
        g.fillEllipse(WIDTH + pd.top * wi - r, height*0.05-r + 4*r, 2*r, 2*r);
    }

    repaint();
}

void MainProcessorEditor::resized()
{
    auto b = getLocalBounds();
    b.setWidth(b.getWidth() - PITCH_WIDTH);
    paramEditor.setBounds(b);
}
