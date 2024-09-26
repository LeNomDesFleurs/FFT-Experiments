/*
 ==============================================================================

 Code by Juan Gil <https://juangil.com/>.
 Copyright (C) 2017-2020 Juan Gil.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.

 ==============================================================================
 */

#pragma once

#include <cmath>

#include "../JuceLibraryCode/JuceHeader.h"
#include "Tools.hpp"

class STFT {
 public:
  enum windowTypeIndex {
    windowTypeRectangular = 0,
    windowTypeBartlett,
    windowTypeHann,
    windowTypeHamming,
    windowTypeFalling,
    windowTypeRaising,
    windowTypeHammingReversed,
    windowTypeNoise,
  };

  //======================================

  STFT() : numChannels(1) {}

  virtual ~STFT() {}

  //======================================

  void setup(const int numInputChannels) {
    numChannels = (numInputChannels > 0) ? numInputChannels : 1;
  }

  void updateParameters(const int newFftSize, const int newOverlap) {
    updateFftSize(newFftSize);
    updateHopSize(newOverlap);
  }

  void updateThreshold(const float new_value) { threshold = new_value; }

  void updateBlur(const float new_value) { blur = new_value; }

  void updateRandomPhase(const float new_value) { randomize_phase = new_value; }

  void updatePhaseAmplitude(const float new_value) {
    phase_amplitude = new_value;
  }

  void updateNoiseFiltering(const float new_value) {
    noise_filtering = new_value;
  }

  void updateRandomAmplitude(const float new_value) {
    randomize_amplitude = new_value;
  }

  void updateFrequencyFolding(const int new_value){
    frequency_folding = new_value;
  }

  //======================================

  void processBlock(AudioSampleBuffer& block) {
    numSamples = block.getNumSamples();
    for (int channel = 0; channel < numChannels; ++channel) {
      float* channelData = block.getWritePointer(channel);

      currentInputBufferWritePosition = inputBufferWritePosition;
      currentOutputBufferWritePosition = outputBufferWritePosition;
      currentOutputBufferReadPosition = outputBufferReadPosition;
      currentSamplesSinceLastFFT = samplesSinceLastFFT;

      for (int sample = 0; sample < numSamples; ++sample) {
        const float inputSample = channelData[sample];
        inputBuffer.setSample(channel, currentInputBufferWritePosition,
                              inputSample);
        if (++currentInputBufferWritePosition >= inputBufferLength)
          currentInputBufferWritePosition = 0;

        channelData[sample] =
            outputBuffer.getSample(channel, currentOutputBufferReadPosition);
        outputBuffer.setSample(channel, currentOutputBufferReadPosition, 0.0f);
        if (++currentOutputBufferReadPosition >= outputBufferLength)
          currentOutputBufferReadPosition = 0;

        if (++currentSamplesSinceLastFFT >= hopSize) {
          currentSamplesSinceLastFFT = 0;

          random_value = noi::Tools::linearCrossfade(
              random.getNormalizedValue(), random_value, noise_filtering);
          analysis(channel);
          modification();
          synthesis(channel);
        }
      }
    }

    inputBufferWritePosition = currentInputBufferWritePosition;
    outputBufferWritePosition = currentOutputBufferWritePosition;
    outputBufferReadPosition = currentOutputBufferReadPosition;
    samplesSinceLastFFT = currentSamplesSinceLastFFT;
  }

  void updateWindowShape(const float window_shape) {
    windowTypeIndex newWindowType = (windowTypeIndex)window_shape;
    float entier;
    float fraction = modf(window_shape, &entier);

    // for (int sample = 0; sample < fftSize; ++sample) fftWindow[sample]
    // = 1.0f;

    if (fraction == 0) {
      fillWindow(fftWindow, newWindowType);
    } else {
      fillWindow(fftWindowA, newWindowType);
      windowTypeIndex next_window = (windowTypeIndex)(((int)newWindowType) + 1);
      fillWindow(fftWindowB, next_window);
      for (int sample = 0; sample < fftSize; sample++)
        fftWindow[sample] = (fftWindowA[sample] * fraction) +
                            (fftWindowB[sample] * (1.f - fraction));
    }
    float windowSum = 0.0f;
    for (int sample = 0; sample < fftSize; ++sample)
      windowSum += fftWindow[sample];

    windowScaleFactor = 0.0f;
    if (overlap != 0 && windowSum != 0.0f)
      windowScaleFactor = 1.0f / (float)overlap / windowSum * (float)fftSize;
  }

 private:
  //======================================

  void updateFftSize(const int newFftSize) {
    fftSize = newFftSize;
    fft = std::make_unique<dsp::FFT>(log2(fftSize));

    inputBufferLength = fftSize;
    inputBuffer.clear();
    inputBuffer.setSize(numChannels, inputBufferLength);

    outputBufferLength = fftSize;
    outputBuffer.clear();
    outputBuffer.setSize(numChannels, outputBufferLength);

    fftWindow.realloc(fftSize);
    fftWindow.clear(fftSize);
    fftWindowA.realloc(fftSize);
    fftWindowA.clear(fftSize);
    fftWindowB.realloc(fftSize);
    fftWindowB.clear(fftSize);

    bin_phases.realloc(fftSize);
    bin_phases.clear(fftSize);
    bin_magnitudes.realloc(fftSize);
    bin_magnitudes.clear(fftSize);

    timeDomainBuffer.realloc(fftSize);
    timeDomainBuffer.clear(fftSize);

    frequencyDomainBuffer.realloc(fftSize);
    frequencyDomainBuffer.clear(fftSize);

    inputBufferWritePosition = 0;
    outputBufferWritePosition = 0;
    outputBufferReadPosition = 0;
    samplesSinceLastFFT = 0;
  }

  void updateHopSize(const int newOverlap) {
    overlap = newOverlap;
    if (overlap != 0) {
      hopSize = fftSize / overlap;
      outputBufferWritePosition = hopSize % outputBufferLength;
    }
  }

  void fillWindow(HeapBlock<float>& window,
                  const windowTypeIndex newWindowType) {
    switch (newWindowType) {
      case windowTypeRectangular: {
        for (int sample = 0; sample < fftSize; ++sample) window[sample] = 1.0f;
        break;
      }
      case windowTypeBartlett: {
        for (int sample = 0; sample < fftSize; ++sample)
          window[sample] =
              1.0f - fabs(2.0f * (float)sample / (float)(fftSize - 1) - 1.0f);
        break;
      }
      case windowTypeHann: {
        for (int sample = 0; sample < fftSize; ++sample)
          window[sample] = 0.5f - 0.5f * cosf(2.0f * M_PI * (float)sample /
                                              (float)(fftSize - 1));
        break;
      }
      case windowTypeHamming: {
        for (int sample = 0; sample < fftSize; ++sample)
          window[sample] = 0.54f - 0.46f * cosf(2.0f * M_PI * (float)sample /
                                                (float)(fftSize - 1));
        break;
      }

      case windowTypeFalling: {
        for (int sample = 0; sample < fftSize; ++sample)
          window[sample] = 1.f - ((float)sample / (float)(fftSize - 1));
        break;
      }
      case windowTypeRaising: {
        for (int sample = 0; sample < fftSize; ++sample)
          window[sample] = (float)sample / (float)(fftSize - 1);
        break;
      }
      case windowTypeHammingReversed: {
        for (int sample = 0; sample < fftSize; ++sample)
          window[sample] =
              1.f - (0.54f - 0.46f * cosf(2.0f * M_PI * (float)sample /
                                          (float)(fftSize - 1)));
        break;
      }
      case windowTypeNoise: {
        for (int sample = 0; sample < fftSize; ++sample)
          window[sample] = random.getNormalizedValue();
        break;
      }
    }
  }

  //======================================

  void analysis(const int channel) {
    int inputBufferIndex = currentInputBufferWritePosition;
    for (int index = 0; index < fftSize; ++index) {
      timeDomainBuffer[index].real(
          fftWindow[index] * inputBuffer.getSample(channel, inputBufferIndex));
      timeDomainBuffer[index].imag(0.0f);

      if (++inputBufferIndex >= inputBufferLength) inputBufferIndex = 0;
    }
  }

  virtual void modification() {
    fft->perform(timeDomainBuffer, frequencyDomainBuffer, false);

    frequency_folding_counter = (++frequency_folding_counter) % frequency_folding;
    if (frequency_folding_counter == 0 || frequency_folding == 0){
      // collect all bin values
      for (int index = 0; index < fftSize / 2 + 1; ++index) {
        bin_magnitudes[index] = abs(frequencyDomainBuffer[index]);
        bin_phases[index] = arg(frequencyDomainBuffer[index]);
      }
    }
    for (int index = 0; index < fftSize / 2 + 1; ++index) {
      float phase = bin_phases[index];
      float magnitude = bin_magnitudes[index];
      phase += (random_value * randomize_phase);
      phase *= phase_amplitude;
      magnitude += random.getNormalizedValue() * randomize_amplitude;
      float prev_magnitude =
          bin_magnitudes[(int)noi::Tools::clipValue(index - 1, 0, fftSize - 1)];
      magnitude = noi::Tools::equalPowerCrossfade(magnitude, prev_magnitude, blur);
      // threshold
      if (magnitude < threshold) {
        magnitude = 0;
      }
      // lock magnitude amplitude
      else {
        //   magnitude = 0.5;
      }
      // safety net
      magnitude = noi::Tools::clipValue(magnitude, 0.F, 1.F);

      float real = magnitude * cosf(phase);
      frequencyDomainBuffer[index].real(real);
      // frequencyDomainBuffer[index].real (magnitude * cosf (phase));

      frequencyDomainBuffer[index].imag(magnitude * sinf(phase));

      if (index > 0 && index < fftSize / 2) {
        frequencyDomainBuffer[fftSize - index].real(magnitude * cosf(phase));
        frequencyDomainBuffer[fftSize - index].imag(magnitude * sinf(-phase));
      }
    }

    fft->perform(frequencyDomainBuffer, timeDomainBuffer, true);
  }

  void synthesis(const int channel) {
    int outputBufferIndex = currentOutputBufferWritePosition;
    for (int index = 0; index < fftSize; ++index) {
      float outputSample = outputBuffer.getSample(channel, outputBufferIndex);
      outputSample += timeDomainBuffer[index].real() * windowScaleFactor;
      outputBuffer.setSample(channel, outputBufferIndex, outputSample);

      if (++outputBufferIndex >= outputBufferLength) outputBufferIndex = 0;
    }

    currentOutputBufferWritePosition += hopSize;
    if (currentOutputBufferWritePosition >= outputBufferLength)
      currentOutputBufferWritePosition = 0;
  }

 protected:
  //======================================
  int numChannels;
  int numSamples;

  int fftSize;
  std::unique_ptr<dsp::FFT> fft;

  int inputBufferLength;
  AudioSampleBuffer inputBuffer;

  int outputBufferLength;
  AudioSampleBuffer outputBuffer;

  HeapBlock<float> fftWindow;
  HeapBlock<float> fftWindowA;
  HeapBlock<float> fftWindowB;
  HeapBlock<float> bin_phases;
  HeapBlock<float> bin_magnitudes;
  HeapBlock<dsp::Complex<float>> timeDomainBuffer;
  HeapBlock<dsp::Complex<float>> frequencyDomainBuffer;

  int overlap;
  int hopSize;
  float windowScaleFactor;

  int inputBufferWritePosition;
  int outputBufferWritePosition;
  int outputBufferReadPosition;
  int samplesSinceLastFFT;

  int currentInputBufferWritePosition;
  int currentOutputBufferWritePosition;
  int currentOutputBufferReadPosition;
  int currentSamplesSinceLastFFT;

  float threshold{};
  float blur{};
  float phase_amplitude{};
  float noise_filtering{};
  float randomize_phase{};
  float randomize_amplitude{};
  float random_value{};
  int frequency_folding_counter{};
  int frequency_folding{};

  noi::Tools::RandomGenerator random{};
};

//==============================================================================
