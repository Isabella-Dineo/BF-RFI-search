#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt 

#=======================================
# FUNCTION:

def chunkit(seq, num):
    ''' Devide an array sequence into chunks of given size'''
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

#=======================================
# TIME SERIES:
#-------------
N = 46022
x = np.random.random_sample((N,))  # time series samples for the duration of the obs.
x = x - x.max()/2.                 
dt = 8.69140625 * 1e-5                # sampling interval in seconds
T = N * dt                         # duration of the observation in seconds
P = 0.089                          # rotation period 
Npulses = int(T/float(P))          # Number of pulses observed

# ADDING PERIODIC SIGNAL:
#------------------------
fp = 1/P            # frequency of the periodic signal
fs = 1/dt           # Sampling frequency in Hz (number of samples per second)
step = N/(T*fp)     # step of the position of the periodic signal
x[::int(step)] =  5 # Adding a periodic signal

''' Devide an array sequence into chunks of given size'''
print '--------------------------------------------------------------'
print 'I                      TIME SERIES INFO                      I'
print '--------------------------------------------------------------\n'
print 'Sampling inteval (s):', dt
print 'Duration of the observation (s):', T
print 'Period of the pulsar (s):', P
print 'Expected number of pulses observed:', Npulses
print '\n--------------------------------------------------------------' 

#    DFT:
#------------
X = np.fft.fft(x)    # dft sample (signal in the frequency domain)
df = 1/float(P)      # Fundamental frequency in Hz (first harmonic)
dw = 2 * np.pi * df  # sampling rate of the dft (rad/sec) [of the dft elememts on the frequency axis]
ny = (dw * N) / 2.   # nyquist frequency (2*max freq in the signal observed = Top frequency)


'''
# FOLDING:
#-------------
n = np.ceil(N/(P * dt * 1e6))        # Number of samples/bins in one pulse period
segmentNumber = int(N/float(n))      # Total number of samples to devide N
folded = np.zeros(int(n))            # An array to store the folded signal
for i in range(int(n)):
    offset = n * i                            # current offset
    pulse = X[offset : offset + segmentSize]  # current segment (first period profile)
    phase = pulse * ((dt * 1e6)/P)            # phase in fractions of a second
    
'''
'''segments = N/1000.
segmentSize = int(N/fs)
full = []

    offset = segmentSize * i
    segment = X[offset : offset + segmentSize]
    full.append(segment)
stacked = np.asarray(full).reshape((segments, segmentSize))
print stacked'''


#FREQUENCY ASSOCIATED WITH A PARTICULAR ELEMENT IN THE DFT (IN X):
#----------------------------------------------------------------

f = np.fft.fftfreq(N, d=dt)
#f = np.fft.fftfreq(N.d)          # dimentionless frequency associated with N samples fft
#f = f * N * df                   # dimensional frequencies in Hz
w = np.fft.fftfreq(N) * N * dw   # dimensional frequencies in rad/s
print '---------------------------------------------------------------'
print 'I                 FREQUENCY DOMAIN INFO                       I'
print '---------------------------------------------------------------\n'
print 'Fundamental frequency (Hz):', df
print 'DFT sampling rate: (rad/s)', dw   
print 'Top frequency (Nyquist frequency Hz)', ny

#------------------
#   POWER SPECTRUM:
#------------------
power = X * np.conj(X)/N

# VELA EVERY 11 Hz:
#onesided = f[0:int(N/2.)]
    
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(x)
plt.title('Time series samples')
plt.subplot(2, 1, 2)
plt.plot(f, abs(X.real))
plt.title('Frequency domain signal')
plt.ylabel('DFT value')
plt.xlabel('Frequency bin centers (cycle per second)')

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(f[0:int(N/2.)], abs(X[0:int(N/2.)]))
plt.ylabel('DFT magnitude')
plt.title('Raw DFT valuesi[one sided]')
plt.subplot(2, 1, 2)
plt.plot(f[0:int(N/2.)], abs(power[0:int(N/2.)]))
plt.ylabel('Power')
plt.xlabel('Frequency bin centers (cycle per second)')
plt.show()
