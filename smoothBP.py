import sys
import numpy as np
from interferometry import splineComplex
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-b', dest='smoothWidth', metavar='smoothWidth',
    help='Smoothing width [MHz]', default='')
parser.add_option('-c', dest='scan', metavar='scan',
    help='Scan ID  e.g. 3', default='')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA42', default='')
parser.add_option('-s', dest='spw', metavar='spw',
    help='Spectral Window', default='')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
refant  = options.refant
smoothWidth=  float(options.smoothWidth)
spw = int(options.spw)
scan= int(options.scan)
#-------- Procedures
Freq = np.load('%s-SPW%d-Freq.npy' % (prefix, spw))
Fwidth = np.median(np.diff(Freq))
chBin = int(smoothWidth / abs(Fwidth) * 1.0e6)
BPant   = np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refant, scan, spw))
XYspec  = np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, refant, scan, spw))
antNum, parapolNum, chNum  = BPant.shape[0], BPant.shape[1], BPant.shape[2]
if chBin > 1:
    smXY, smBP = np.ones(chNum, dtype=complex), np.ones([antNum, parapolNum, chNum], dtype=complex)
    smXY = splineComplex(np.arange(chNum), XYspec, chBin)
    for ant_index in list(range(antNum)):
        for pol_index in list(range(parapolNum)):
            smBP[ant_index, pol_index] = splineComplex(np.arange(chNum), BPant[ant_index, pol_index], chBin)
else:
    smXY, smBP = XYspec, BPant
np.save('%s-REF%s-SC-1-SPW%d-BPant.npy'  % (prefix, refant, spw), smBP) 
np.save('%s-REF%s-SC-1-SPW%d-XYspec.npy' % (prefix, refant, spw), smXY) 
