#!/usr/bin/env python 

__author__ = 'IsabellaRammala'

import numpy as np 
import glob, sys, os, time 
import argparse
import h5py
import struct
#import ephem
#import katpoint
import matplotlib.pyplot as plt
import scipy 
from scipy.fftpack import fft

#==============================================
# Define Functions
#==============================================
def stft(signal, samplingFreq, spectraNum, fftSize, windowSize=None):
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
    innerPad = np.zeros((fftSize, 2))
    padSignal = np.concatenate([signal, innerPad]) # padded signal 
    
    for i in range(len(se)): 
    result = np.empty((np.int64(segments), np.int64(fftSize)), dtype=float)
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
    
#==============================================
# Initialize parameters from the command line:
#==============================================
parser = argparse.ArgumentParser(description="AR1 BF RFI characterization script.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--freq", metavar="<freqCent>", type=float, default="1391.0", help="Center frequency in MHz." )
parser.add_argument("--direc", metavar="<path/to/file>", type=str, default=None, help="Directory to look for the .h5 file.")
parser.add_argument("--h5file", metavar="<.h5File>", type=str, default=None, help="Input .h5 file.")
parser.add_argument("--chunk", metavar="<chunkSize>", type=int, default="256", help="Number of samples to process.")
parser.add_argument("--out", metavar="<outFileName>", type=str, default="outFile", help="Output file name.")
args = parser.parse_args()
directory = args.direc
h5file = args.h5file
chunkSize = args.chunk
freqCent = args.freq

#===============================================
# Get data from h5 input file:
#===============================================
t = time.time() # script start time
# Checking if h5 file or location is specified.
if h5file != None: 
    if directory != None:
        file = directory + h5file
    else:
        file = h5file
else:
    if directory != None:
        files = np.array(glob.glob(directory + "/*.h5"))
        if files.size > 1:
            print '\n List of files to choose from:\n', files
            print '\n Select a file to search...'
            file = raw_input()
        else:
            print '\n File selected = ', files[0]
            file = files[0]
    else:
        print '\n Give .h5 file. Exiting...\n'
        sys.exit(1) 
# Global variables.
# Sampling time in microseconds .
samplingTime = 4.78504672897196 * 1e-6

# Sampling frequency in Hz.
samplingFreq = 1712.0e6
print ("samplingFreq : 1712.0 MHz")

# Channel bandwidth.
channelBW =  1712.0/(1712.0*4.78504672897196)
#channelBW = 1712.0/8192.0
print ("channelBW: %f") % channelBW

# Load the h5 file.
dataH5 = h5py.File(h5file, "r")
print ("h5 file: %s") % h5file
print ("h5 Data: %s") % dataH5

# Getting number of channels.
channelNumber = dataH5["Data/bf_raw"].shape[0]
print ("Channel Number: %d") % channelNumber

# Determine maximum and minimum frequency.
freqMax = (((channelNumber/2) - 1) * channelBW) + freqCent
freqMin = freqCent - ((channelNumber/2) * channelBW)
print ("Bottom freq: %.2f") % freqMin
print ("Top freq: %.2f") % freqMax
# Getting number of dump products.
spectraNumber = dataH5["Data/bf_raw"].shape[1]
print ("Spectra Number: %d") % spectraNumber
print ("SHAPE = ") + str(np.shape(dataH5["Data/bf_raw"]))
# Find observation start time.
startTime = h5file.split(".")[0]
#print ("File start time: ") 
#os.system("date --utc -d @'startTime'")

# Getting the ADC counts/ timestamps.
countsADC = dataH5["Data/timestamps"][:]
print ("contADC: %d") % countsADC[0]

# Extracting data from h5 file.
#spectraChunk = dataH5["Data/bf_raw"][:].size
#print ("SpectraChunk: %f") % spectraChunk

#======================================================
#  	 	TRIAL ON 1 CHANNEL
#======================================================
print ("I==========================================I")
print ("I	 Searching for RFI in channel 0    I")
print ("I==========================================I")
# Starting with channel 0
print ("Analysing channel #0")

# Extract spectrum in channel 0.
spectrum = dataH5["Data/bf_raw"][0]
print ("Spectrum shape:" + str(np.shape(spectrum)))

# STFT algorithm (using the concept of 'sliding window'):
fftSize = 500 # testvalue
spectrumFT = stft(spectrum, samplingFreq/1e6, spectraNumber, fftSize, windowSize=chunkSize)

#plot the spectrum in that channel.
#print ("Plotting spectrum in channel 0")
plt.plot(countsADC[0:chunkSize], spectrum[0:chunkSize])
plt.title("real spectrum value spectrum[0]")

plt.figure()
plt.imshow(spectrumFT, origin='lower', cmap='jet', interpolation='nearest', aspect='auto')

# Fourier transform the spectra.
print ("\n Fourier Transforming the time-series...")
fourierSpectrum = scipy.fft(spectrum)
#freqs = scipy.fftpack.fftfreq(chunkSize)
print ("Plotting the FFT spectra")
plt.figure()
plt.plot(range(chunkSize), fourierSpectrum[0:chunkSize])
#plt.grid()
plt.show()
