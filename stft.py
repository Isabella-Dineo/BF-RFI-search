#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt 

#-----------------
# define function:
#-----------------


def stft(signal, samplingFreq, spectraNum, fftSize, windowSize=None):
    """Function to perform short-time fourier transform of the signal
    
    Args:
    -----
    signal       : ndarray containing the signal to process
    samplingFreq : sampling frequency (Hz)
    spectraNum   : number of spectra in each channel
    fftSize      : number of samples in FT
    windowSize   : size of the segment of the signal (optimal windowSize = 256).
    
    """
    # 1. Pick a segment of the signal
    N = len(signal)
    if windowSize == None:
        print ("Window size not definedi for STFT. Using the whole spectrum...")
        windowSize = N # use the whole spectra if segment size is not defined

    segments = np.int(np.ceil(signal) / np.float(windowSize)) # Total number of segments
    
    # 2. Multipy by half-cosine function
    windowFunc = np.hanning(fftSize)
    # 3. Zeros to double the size of each segment.
    innerPad = np.zeros(fftSize)
    padSignal = np.concatenate(signal, innerPad) # padded signal 

    result = np.empty((segments, fftSize), dtype=float)
    for sid in xrange(segments):
        hop = fftSize * sid                              # current segment ofset
        segment = padSignal[hop:hop+fftsize]             # current segment
        windowed = segment * windowFunc                  # multiply by tapering function
        padded = np.append(windowed, innerPad)           # pad signal to double the data
        spectrum = np.fft.fft(padded) / fft_size         # FFT and scale by number of samples
        autopower = np.abs(spectrum * np.conj(spectrum)) # Find the autopower spectrum
        result[i, :] = autopower[:fft_size]

    result = 20*np.log10(result)          # scale to db
    result = np.clip(result, -40, 200)    # clip values

    return result

