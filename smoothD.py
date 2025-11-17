import os
import sys
import numpy as np
from interferometry import GetAntName, splineComplex
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antList', metavar='antList',
    help='Antenna List (default = all)', default='')
parser.add_option('-b', dest='smoothWidth', metavar='smoothWidth',
    help='Smoothing width [MHz]', default='')
parser.add_option('-s', dest='spw', metavar='spw',
    help='Spectral Window', default='')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
msfile = prefix + '.ms'
if options.antList == '':
    antList = GetAntName(msfile).tolist()
else: 
    antList  = [ant for ant in options.antList.split(',')]
smoothWidth=  float(options.smoothWidth)
spw = int(options.spw)
#-------- Procedures
for ant_index, ant in enumerate(antList):
    DfileName = '%s-SPW%d-%s.DSpec.npy' % (prefix, spw, ant)
    if not os.path.isfile(DfileName): continue
    Dspec = np.load(DfileName)
    Freq = Dspec[0]     # Frequency [GHz]
    Fwidth = np.median(np.diff(Freq))
    chBin = int(smoothWidth / abs(1.0e3* Fwidth))
    chNum = Dspec.shape[1]
    if chBin > 1:
        Dx, Dy = Dspec[1] + (0.0 + 1.0j)* Dspec[2], Dspec[3] + (0.0 + 1.0j)* Dspec[4]
        smDx, smDy = splineComplex(np.arange(chNum), Dx, chBin), splineComplex(np.arange(chNum), Dx, chBin)
        Dspec[1], Dspec[2], Dspec[3], Dspec[4]  = smDx.real, smDx.imag, smDy.real, smDy.imag
        np.save(DfileName, Dspec) 
#
