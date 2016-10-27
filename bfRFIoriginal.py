#!/usr/bin/env python

import numpy as np
import glob, sys, os, time
import h5py
import katdal
import matplotlib.pyplot as plt
import argparse

#==============================================
# Initialize parameters from the command line:
#==============================================
parser = argparse.ArgumentParser(description="AR1 BF RFI characterization script. Uses the corresponding CBF .h5 file to determine which channels are RFI flagged.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--freq", metavar="<freqCent>", type=float, default="1391.0", help="Center frequency in MHz." )
parser.add_argument("--bf0_direc", metavar="<path/to/bf_file>", type=str, default=None, help="Directory to look for the BF .h5 file (h pol).")
parser.add_argument("--bf1_direc", metavar="<path/to/bf_file>", type=str, default=None, help="Directory to look for the BF .h5 file (v pol).")
parser.add_argument("--cbf_direc", metavar="<path/to/cbf_file>", type=str, default=None, help="Directory to look for the CBF .h5 file.")
parser.add_argument("--bf0", metavar="<bf_h5File>", type=str, default=None, help="Input BF .h5 file (h pol).")
parser.add_argument("--bf1", metavar="<bf_h5File>", type=str, default=None, help="Input BF .h5 file (v pol).")
parser.add_argument("--cbf_h5", metavar="<cbf_h5File>", type=str, default=None, help="Input CBF .h5 file.")
parser.add_argument("--out", metavar="<outfile name>", type=str, default=None, help="Output file names.")
args = parser.parse_args()
bf0_direc = args.bf0_direc
bf1_direc = args.bf1_direc
cbf_direc = args.cbf_direc 
bf0 = args.bf0
bf1 = args.bf1
cbf_h5 = args.cbf_h5
freqCent = args.freq
outfile = args.out

#================================================
# 		MAIN BODY:
#================================================
#-----------------------------------------------------------------------------------------------------
# Load the cbf file to search for flagged channels.
#-----------------------------------------------------------------------------------------------------

# Locate the cbf file.
if cbf_h5 != None:
    if cbf_direc != None:
        file = cbf_direc + cbf_h5
    else:
        file = cbf_h5
else:
    if cbf_direc != None:
        files = np.array(glob.glob(cbf_direc + "/*.h5"))
        if files.size > 1:
            print '\n List of files to choose from:\n', files
            print '\n Select a file to search...'
            file = raw_input()
        else:
            print '\n File selected = ', files[0]
            file = files[0]
    else:
        print '\n Give the corresponding visibility .h5 file to check get flagged channels. Exiting...\n'
        sys.exit(1)

# Load the cbf file with katdal and get flagged channels.
cbfFile = katdal.open(cbf_h5)
threshold = 0.8
flags = cbfFile.flags()[:,:,0]
flagChans = np.where(np.mean(flags, axis=0) > threshold)
chansRFI = flagChans[0]
print "Channels with RFI: \n", chansRFI

#-----------------------------------------------------------------------------------------------------
# Load the bf file.
#-----------------------------------------------------------------------------------------------------

# Locate the bf file.
if bf0 != None:
    if bf0_direc != None:
        file = direc + bf0
    else:
        file = bf0
else:
    if direc != None:
        files = np.array(glob.glob(direc + "/*.h5"))
        if files.size > 1:
            print '\n List of files to choose from:\n', files
            print '\n Select a file to search...'
            file = raw_input()
        else:
            print '\n File selected = ', files[0]
            file = files[0]
    else:
        print '\n Give the bf .h5 file. Exiting...\n'
        sys.exit(1)

# Load the bf file with h5py.
bf0 = h5py.File(bf0, "r") # h pol
bf1 = h5py.File(bf1, "r") # v pol

# Getting number of channels.
channels0 = bf0["Data/bf_raw"].shape[0]
channels1 = bf1["Data/bf_raw"].shape[0]

if channels0 != channels1:
    print "Channels in the BF files are not the same..."
    sys.exit(1)

# Total number of samples.
samples = bf0["Data/bf_raw"].shape[1]
 
# Time interval between samples (seconds).
samplingTime = 4.78504672897196 *1e-6 # in microseconds

obsTime = samples*samplingTime
print 'obsTime:', obsTime
# Frequency spacing (frequency resolution in Hz).
#freqRes = 1/(samples * samplingTime)

# Sampling rate (sampling frequency Hz).
samplingRate = 1 / samplingTime 
print samplingRate
print samplingTime

# Channel bandwidth in MHz.
channelBW =  1712.0/(1712.0 * 4.78504672897196)

# Maximum channel frequency (MHz).
freqMax = ((channels0/2 - 1) * channelBW) + freqCent

# Minimum channel frequency (MHz).
freqMin = freqCent - ((channels0/2 - 1) * channelBW)

#Frequency range (MHz).
freqRange = np.linspace(freqMin, freqMax, num=channelBW, endpoint=True)


#-------------------------------------------------------------------------------------------------------
# Fourier transforms.
#-------------------------------------------------------------------------------------------------------
# Number of segments to divide the samples equvalent to samples every 4 seconds.
segments = np.ceil(samples/(4*float(samplingRate))) # Total number of segments

# Length of the Fourier transform:
segmentSize = int(samples/float(segments)) # samples every 4 seconds

# Test on a 'clean channel'
cleanChannel = 2048/2
fullChannelPol0 = []
for i in range(int(segments)):
#for i in range(2):
    offset = segmentSize * i                                                          # current offset
   # segmentPol0 = bf0['Data/bf_raw'][cleanChannel, offset : offset + segmentSize, :] * np.conj(bf0['Data/bf_raw'][cleanChannel, offset : offset + segmentSize, :])  # current segment spectrum
    segmentPol0 = bf0['Data/bf_raw'][cleanChannel, offset : offset + segmentSize, 0]
    #dftPol0 = np.fft.fft(segmentPol0)
    dftPol0 = np.abs(np.fft.fft(segmentPol0))
    #fullChannelPol0 = np.hstack((fullChannelPol0, dftPol0))
    fullChannelPol0.append(dftPol0)
stacked = np.asarray(fullChannelPol0).reshape((segments, segmentSize))
avg = np.average(stacked, axis=0)
plt.plot(avg)
plt.figure()
plt.plot(avg[len(avg)/2:])
#freqRange = (1/(samplingTime * 1e6)) #* np.arange(len(amplitudes))/len(amplitudes) 
#freqRange = 11*np.arange(len(amplitudes))/len(amplitudes) 
#freqRange = np.arange(len(amplitudes))/len(amplitudes) 
#power = amplitudes * np.conj(amplitudes)/(len(amplitudes)**2)
#plt.figure()
#plt.plot(freqRange, power)
#plt.title('Power spectral density in "clean channel", channel: %d' %cleanChannel)
#plt.xlabel('normalized frequency')
#plt.ylabel('Power')


# Test on 1 channel with RFI:


channelRFI = 384
rfiChannelPol0 = []
for i in range(int(segments)):
    offset = segmentSize * i                                                      
    segmentPol0 = bf0['Data/bf_raw'][channelRFI, offset : offset + segmentSize, 0]
    dftPol0 = np.abs(np.fft.fft(segmentPol0))                              
    rfiChannelPol0.append(dftPol0)
stacked0 = np.asarray(rfiChannelPol0).reshape((segments, segmentSize))
avg0 = np.average(stacked0, axis=0)
amp = rfiChannelPol0[0:len(rfiChannelPol0)/2]   # half the dft size
freqR = (1/(samplingTime * 1e6)) * np.arange(len(amp))/len(amp)   # Normalized frequency axis
power0 = amp * np.conj(amp)/(len(amp)**2)
plt.figure()
plt.plot(freqR, power0) 
plt.title('Power spectral denstiy plot of a channel with RFI, channel: %d' %channelRFI)
plt.xlabel('Normalized frequency')
plt.ylabel('power')
plt.figure(avg0)
# Test on flagged channels
# 1. Perform FFTs on each channelised time-series (channels with RFI)
#fig = plt.figure() 
#allChannelspol0 = []
#allChannelspol1 = []
#for cid in chansRFI:
#    if cid > channels0:
#        break
#    print 'Channel:', cid
#    fullChannelpol0 = []
#    fullChannelpol1 = []
#    for sid in range(int(segments)):
        # Determine the next segment offset.
#        offset = segmentSize * sid         
        # current segment of samples.        
#        segmentpol0 = bf0['Data/bf_raw'][cid, offset : offset + segmentSize, 0]
#        segmentpol1 = bf1['Data/bf_raw'][cid, offset : offset + segmentSize, 0]
        # Perform FFT on each segment.
#        dftpol0 = np.abs(np.fft.fft(segmentpol0)[0:len(segmentpol0)/2]) # amplitudes (half of the dft results)
#        dftpol1 = np.abs(np.fft.fft(segmentpol1)[0:len(segmentpol1)/2]) 
#        halfdft = dftpol0[0:len(dftpol0)/2]
#        fullChannelpol0 = np.hstack((fullChannelpol0, dftpol0))
#    allChannelspol0.append(fullChannelpol0)
#    print np.shape(fullChannelpol0)
#    plt.figure()
#    plt.plot(fullChannelpol0) # amplitude plots
#    plt.title('Amplitude plots, channel: %d' %cid) 
#    plt.xlabel('samples')
#    plt.ylabel('dft values')
#    plt.plot(np.arange(len(fullChannelpol0)), fullChannelpol0)
#    plt.plot(np.arange(len(fullChannelpol0))/len(fullChannelpol0), fullChannelpol0) # normalized frequency range
#    plt.show()
#    allChannelspol0.append(fullChannelpol0)
#print len(chansRFI)
#print np.shape(allChannelspol0)
#        full_channel.append(np.abs(np.fft.fft(segment1))) 
#    np.concatenate((full_channel))
#print np.shape(full_channel)

#plt.imshow(full_channel)
plt.show()
