/*
  ==============================================================================

    This code is based on the contents of the book: "Audio Effects: Theory,
    Implementation and Application" by Joshua D. Reiss and Andrew P. McPherson.

    Code by Juan Gil <https://juangil.com/>.
    Copyright (C) 2017-2020 Juan Gil.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "PluginParameter.h"

//==============================================================================

TemplateFrequencyDomainAudioProcessor::TemplateFrequencyDomainAudioProcessor():
#ifndef JucePlugin_PreferredChannelConfigurations
    AudioProcessor (BusesProperties()
                    #if ! JucePlugin_IsMidiEffect
                     #if ! JucePlugin_IsSynth
                      .withInput  ("Input",  AudioChannelSet::stereo(), true)
                     #endif
                      .withOutput ("Output", AudioChannelSet::stereo(), true)
                    #endif
                   ),
#endif
    parameters (*this)
    , paramFftSize (parameters, "FFT size", fftSizeItemsUI, fftSize512,
                    [this](float value){
                        const ScopedLock sl (lock);
                        value = (float)(1 << ((int)value + 5));
                        paramFftSize.setCurrentAndTargetValue (value);
                        stft.updateParameters((int)paramFftSize.getTargetValue(),
                                              (int)paramHopSize.getTargetValue());
                        return value;
                    })
    , paramHopSize (parameters, "Hop size", hopSizeItemsUI, hopSize8,
                    [this](float value){
                        const ScopedLock sl (lock);
                        value = (float)(1 << ((int)value + 1));
                        paramHopSize.setCurrentAndTargetValue (value);
                        stft.updateParameters((int)paramFftSize.getTargetValue(),
                                              (int)paramHopSize.getTargetValue());
                        return value;
                    })
    , paramWindowType (parameters, "Window type", windowTypeItemsUI, STFT::windowTypeHann,
                       [this](float value){
                           const ScopedLock sl (lock);
                           paramWindowType.setCurrentAndTargetValue (value);
                           stft.updateWindowShape((float)paramWindowType.getTargetValue());
                           return value;
                       })
    , paramcontinuousfftsize (parameters, "hopp size", "samples", 0.f, 1.f, 0.f, 
                        [this](float value){
                            const ScopedLock sl (lock);
                            paramcontinuousfftsize.setCurrentAndTargetValue(value);
                           stft.updateBlur(value);
                           return value;
                        })
    , paramThreshold (parameters, "Threshold", "units", 0.f, 1.f, 0.f, 
                        [this](float value){
                            const ScopedLock sl(lock);
                            paramThreshold.setCurrentAndTargetValue(value);
                            stft.updateThreshold(value);
                            return value;
                        })
    , paramRandomizePhase (parameters, "random phase", "uni", 0.f, 10.f, 0.f, 
                        [this](float value){
                            const ScopedLock sl(lock);
                            paramRandomizePhase.setCurrentAndTargetValue(value);
                            stft.updateRandomPhase(value);
                            return value;
                        })
    , paramNoiseFiltering (parameters, "noise filtering", " units",0.f, 1.f, 0.f, 
                        [this](float value){
                            const ScopedLock sl(lock);
                            paramNoiseFiltering.setCurrentAndTargetValue(value);
                            stft.updateNoiseFiltering(value);
                            return value;
                        })                        
    , paramPhaseAmplitude (parameters, "phase amplitude", " units",0.f, 1.f, 1.f, 
                        [this](float value){
                            const ScopedLock sl(lock);
                            paramPhaseAmplitude.setCurrentAndTargetValue(value);
                            stft.updatePhaseAmplitude(value);
                            return value;
                        })
    , paramRandomizeAmplitude (parameters, "random amlitude", "uni", 0.f, 2.f, 0.f, 
                        [this](float value){
                            const ScopedLock sl(lock);
                            paramRandomizeAmplitude.setCurrentAndTargetValue(value);
                            stft.updateRandomAmplitude(value);
                            return value;
                        })
    , paramWindowShape (parameters, "win shape", "t", 0.f, 7.f, 0.f, 
                        [this](float value){
                            const ScopedLock sl(lock);
                            paramWindowShape.setCurrentAndTargetValue(value);
                            stft.updateWindowShape(value);
                            return value;
                        })
    , paramFrequencyFolding (parameters, "folding", "t", 0.f, 100.f, 0.f, 
                        [this](float value){
                            const ScopedLock sl(lock);
                            paramFrequencyFolding.setCurrentAndTargetValue(value);
                            stft.updateFrequencyFolding((int)value);
                            return value;
                        })
{
    parameters.apvts.state = ValueTree (Identifier (getName().removeCharacters ("- ")));
}

TemplateFrequencyDomainAudioProcessor::~TemplateFrequencyDomainAudioProcessor()
{
}

//==============================================================================

void TemplateFrequencyDomainAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    const double smoothTime = 1e-3;
    paramFftSize.reset (sampleRate, smoothTime);
    paramHopSize.reset (sampleRate, smoothTime);
    paramWindowType.reset (sampleRate, smoothTime);
    paramThreshold.reset (sampleRate, smoothTime);
    paramcontinuousfftsize.reset (sampleRate, smoothTime);
    paramRandomizePhase.reset(sampleRate, smoothTime);
    paramRandomizeAmplitude.reset(sampleRate, smoothTime);
    paramPhaseAmplitude.reset(sampleRate, smoothTime);
    paramWindowShape.reset(sampleRate, smoothTime);
    paramFrequencyFolding.reset(sampleRate, smoothTime);

    //======================================

    stft.updateRandomAmplitude(paramPhaseAmplitude.getTargetValue());
    stft.updateWindowShape(paramWindowShape.getTargetValue());
    stft.setup(getTotalNumInputChannels());
    stft.updateParameters((int)paramFftSize.getTargetValue(),
                          (int)paramHopSize.getTargetValue());
    stft.updateThreshold(paramThreshold.getTargetValue());
    stft.updateFrequencyFolding(paramFrequencyFolding.getTargetValue());
}

void TemplateFrequencyDomainAudioProcessor::releaseResources()
{
}

void TemplateFrequencyDomainAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
    const ScopedLock sl (lock);

    ScopedNoDenormals noDenormals;

    const int numInputChannels = getTotalNumInputChannels();
    const int numOutputChannels = getTotalNumOutputChannels();
    const int numSamples = buffer.getNumSamples();

    //======================================

    stft.processBlock (buffer);

    //======================================

    for (int channel = numInputChannels; channel < numOutputChannels; ++channel)
        buffer.clear (channel, 0, numSamples);
}

//==============================================================================






//==============================================================================

void TemplateFrequencyDomainAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    auto state = parameters.apvts.copyState();
    std::unique_ptr<XmlElement> xml (state.createXml());
    copyXmlToBinary (*xml, destData);
}

void TemplateFrequencyDomainAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    std::unique_ptr<XmlElement> xmlState (getXmlFromBinary (data, sizeInBytes));

    if (xmlState.get() != nullptr)
        if (xmlState->hasTagName (parameters.apvts.state.getType()))
            parameters.apvts.replaceState (ValueTree::fromXml (*xmlState));
}

//==============================================================================

AudioProcessorEditor* TemplateFrequencyDomainAudioProcessor::createEditor()
{
    return new TemplateFrequencyDomainAudioProcessorEditor (*this);
}

bool TemplateFrequencyDomainAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

//==============================================================================

#ifndef JucePlugin_PreferredChannelConfigurations
bool TemplateFrequencyDomainAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

//==============================================================================

const String TemplateFrequencyDomainAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool TemplateFrequencyDomainAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool TemplateFrequencyDomainAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool TemplateFrequencyDomainAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double TemplateFrequencyDomainAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

//==============================================================================

int TemplateFrequencyDomainAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int TemplateFrequencyDomainAudioProcessor::getCurrentProgram()
{
    return 0;
}

void TemplateFrequencyDomainAudioProcessor::setCurrentProgram (int index)
{
}

const String TemplateFrequencyDomainAudioProcessor::getProgramName (int index)
{
    return {};
}

void TemplateFrequencyDomainAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================

// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new TemplateFrequencyDomainAudioProcessor();
}

//==============================================================================
