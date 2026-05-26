cmake_minimum_required(VERSION 3.25)

# Windows static runtime
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

# MacOS unversal binary
set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64")

# MacOS minimum deployment target
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")

# Main project
project(MRTA VERSION 0.1.0 LANGUAGES C CXX)

# Get 'add_plugin' function and add JUCE to project
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(add_plugin)


# General variables
# You should change this to something more interesting
set(company_name "Modern Real-Time Audio")
set(company_copyright "Modern Real-Time Audio Oy")
set(company_code "Mrta")
set(company_bundle_id "com.ModernRealTimeAudio")


## Plugin projects

# Example project
set(ortaa_source ${CMAKE_CURRENT_SOURCE_DIR}/projects/OurRealTimeAudioApp)

add_plugin(ortaa
    VERSION 0.1.0
    PLUGIN_NAME "Harmonizer"
    PROD_NAME Harmonizer
    PROD_CODE Harm
    SYNTH FALSE
    SOURCES
        ${ortaa_source}/PluginEditor.cpp
        ${ortaa_source}/PluginProcessor.cpp
        ${ortaa_source}/math/commonMath.cpp
        ${ortaa_source}/math/fft-impl.cpp
        ${ortaa_source}/math/fft.cpp
        ${ortaa_source}/dsp/pitch.cpp
        ${ortaa_source}/dsp/psola.cpp
        ${ortaa_source}/dsp/StateVariableFilter.cpp
    INCLUDE_DIRS
        ${ortaa_source}
        ${ortaa_source}/math
        ${ortaa_source}/dsp)
