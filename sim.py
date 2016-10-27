#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Simulating a sinusoidal wave observed at 800 MHz and its FFT

# Define Functions.

def stft(signal, samplingFreq, fftSize, windowSize=None):
    """Function to perform short-time fourier transform of the signal
    
    Args:
    -----
    signal     : ndarray containing the signal to process
    spectraNum : number of spectra in each channel
    fftSize    : number of samples in FT
    windowSize : size of the segment of the signal (optimal windowSize = 256).
    
    """
    # 1. Pick a segment of the signal
    N = len(signal)
    if windowSize == None:
        print ("Window size not defined. Using the whole spectrum...")
        windowSize = N # use the whole spectra if segment size is not defined

    #segments = np.int(np.ceil(dataH5["Data/bf_raw"].shape[1]) / np.float(256))
    segments = np.int64(np.ceil(signal) / np.float(windowSize)) # Total number of segments
    # 2. Multipy by half-cosine function
    windowFunc = np.hanning(fftSize)
    # 3. Zeros to double the size of each segment.
    innerPad = np.zeros((fftSize, signal.shape[0]))
    data = np.concatenate([signal, innerPad]) # padded signal 

    #for i in range(len(se)):
    result = np.empty((np.int64(segments), np.int64(fftSize)), dtype=float)
    for sid in xrange(segments):
        hop = fftSize * sid                              # current segment ofset
        segment = data[hop:hop+fftsize]                  # current segment
        windowed = segment * windowFunc                  # multiply by tapering function
        padded = np.append(windowed, innerPad)           # pad signal to double the data
        spectrum = np.fft.fft(padded) / fft_size         # FFT and scale by number of samples
        autopower = np.abs(spectrum * np.conj(spectrum)) # Find the autopower spectrum
        result[i, :] = autopower[:fft_size]

    result = 20*np.log10(result)          # scale to db
    result = np.clip(result, -40, 200)    # clip values

    return result

#----------------------------
# MAIN BODY:
#----------------------------

f = 100 # frequency of the wave in Hz
samplingRate = 300 # sampling rate
fs = samplingRate * f # sampling frequency
phase = 0 # phase shift
ncycle = 5 # Cycles of the wave
time = np.arange(0, ncycle * (1/float(f)), 1/float(fs))
signal = np.sin(2 * np.pi * f * time + phase)
signal = np.reshape(signal, np.product(signal.shape))
plt.plot(time, signal)
plt.title("Signal in time domain")
plt.xlabel('time (s)')
plt.ylabel('amplitude')

#--------------------------------------
# Representation in frequency domain.
#--------------------------------------
#dft = stft(signal, fs, fftSize=1024, windowSize=int(ncycle/float(1/float(fs))))
dft = np.fft.fft(signal)
samplep = np.arange(0,len(dft))
plt.figure()
plt.plot(samplep, abs(dft.real))
plt.title('Raw DFT values')
plt.xlabel('Sample points')
plt.ylabel('DFT values')

# Normalised frequency
normFreq = np.linspace(0, (len(dft) - 1)/len(dft), num=len(dft))
plt.figure()
plt.plot(normFreq, dft)

plt.show()
