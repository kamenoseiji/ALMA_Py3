import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from Plotters import plotDSpec
#exec(open(SCR_DIR + 'Plotters.py').read())
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA42', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
spwList = [int(spw) for spw in  options.spwList.split(',')]
refantName = options.refant
#
FreqList, DxList, DyList = [], [], []
for spw_index, spw in enumerate(spwList):
    antList = np.load('%s-SPW%d-%s.Ant.npy' % (prefix, spw, refantName))
    for ant_index, antName in enumerate(antList):
        Dterm = np.load('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, antName))
        DxList   = DxList   + [Dterm[1] + (0.0 + 1.0j)* Dterm[2]]
        DyList   = DyList   + [Dterm[3] + (0.0 + 1.0j)* Dterm[4]]
    #
    FreqList = FreqList + [Dterm[0]]
#
pp = PdfPages('D_%s-REF%s-Dspec.pdf' % (prefix, refantName))
plotDSpec(pp, prefix, antList, spwList, FreqList, DxList, DyList)
