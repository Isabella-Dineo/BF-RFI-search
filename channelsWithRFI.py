#!/usr/bin/env python

import katdal
import numpy as np
import argparse 

# Command line parser:
parser = argparse.ArgumentParser(description="Check the correlator .h5 files for channels that are RFI flagged .", prog='Command', usage='%(prog)s [option]', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--h5file", metavar="<.h5File>", type=str, default=None, help="Input .h5 file.")

# initialize parameters.
args = parser.parse_args()
threshold = 0.8

# Load the .h5 file and get flagged channels.
f = katdal.open(args.h5file)
flags = f.flags()[:,:,0]
flag_chans = np.where(np.mean(flags, axis=0) > threshold)[0]
for chanID, chan in enumerate(flag_chans):
    print "Channel flagged: ", chan

#for indices in zip(*np.where(np.mean(f.flags()[:], axis=0) > threshold)[0]):
#    print "Channel: ", indices[0], "Frequency: ", f.channel_freqs[indices[0]]/1e6, "Baseline: ", f.corr_products[indices[1]]


