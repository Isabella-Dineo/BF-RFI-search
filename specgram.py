#!/usr/bin/env python

import scipy
import matplotlib.pyplot as plt

#===============================
# DEFINE FUNCTIONS:
#===============================

def stft(x, fs, feamesz, hop):
    """Function to perform a short time fourier transform

       Args:
       -----
       x       : signal to proces (ndarray)
       fs      : sampling frequency (Hz)
       framesz : frame size (Hz)
       hop     : hop size (seconds)

       Return:
       -------
       X       : short-time fourier transform of x
    """
    framesamp = int(framesz*fs) # sample frame size
    hopsamp = int(hop*fs) 
    w = scipy.hanning(framesamp) # Hann wondow (a tapering window)
    X = scipy.array([scipy.fft(w*x[i:i+framesamp]) 
                     for i in range(0, len(x)-framesamp, hopsamp)])
    return X

# TEST

if __name__ == '__main__':
    f0 = 440         # Compute the STFT of a 440 Hz sinusoid
    fs = 8000        # sampled at 8 kHz
    T = 5            # lasting 5 seconds (length of the signal)
    framesz = 0.050  # with a frame size of 50 milliseconds
    hop = 0.025      # and hop size of 25 milliseconds.

    # Create test signal and STFT.
    t = scipy.linspace(0, T, T*fs, endpoint=False)
    x = scipy.sin(2*scipy.pi*f0*t)
    X = stft(x, fs, framesz, hop)

    # Plot the magnitude spectrogram.
    plt.figure()
    plt.imshow(scipy.absolute(X.T), origin='lower', aspect='auto',
                 interpolation='nearest')
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.show()

    # Compute the ISTFT.
#    xhat = istft(X, fs, T, hop)

    # Plot the input and output signals over 0.1 seconds.
#    T1 = int(0.1*fs)

#    pylab.figure()
#    pylab.plot(t[:T1], x[:T1], t[:T1], xhat[:T1])
#    pylab.xlabel('Time (seconds)')

#    pylab.figure()
#    pylab.plot(t[-T1:], x[-T1:], t[-T1:], xhat[-T1:])
#    pylab.xlabel('Time (seconds)')
