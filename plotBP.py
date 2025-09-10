import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import indexList
from Plotters import plotSP
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan ID  e.g. 3,5,7', default='')
parser.add_option('-p', dest='plotAnt', metavar='plotAnt',
    help='Antennas to plot', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPWs to plot e.g. 17,19,21', default='')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA41', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
scanList= [int(scan) for scan in options.scanList.split(',')]
spwList = [int(spw) for spw in options.spwList.split(',')]
plotAntList = [] if options.plotAnt == '' else [ant for ant in options.plotAnt.split(',')]
refantName = options.refant
'''
prefix = 'uid___A002_X12dca4b_Xcf44'
scanList = [3,5,7]
spwList = [5,17,21,23]
plotAntList = ['DA41','DA42']
refantName = 'DA51'
'''
antList = np.load(prefix + '.Ant.npy')
if plotAntList == []: plotAntList = antList
antMap = indexList(plotAntList, antList)
plotAntList = antList[antMap].tolist()
plotAntNum, columnNum = len(plotAntList), len(spwList)
#-------- Bandpass Table
antList = np.load('%s-REF%s.Ant.npy' % (prefix, refantName))
BPList, FreqList = [], []
for spw_index, spw in enumerate(spwList):
    #-- First Target scan
    for scan_index, scan in enumerate(scanList):
        print('--Loading Bandpass Table for spw%d scan%d' % (spw, scan))
        BPList = BPList + [np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refantName, scan, spw))]
        FreqList = FreqList + [np.load('%s-SPW%d-Freq.npy' % (prefix, spw))]
    #
#
ppolNum = BPList[0].shape[1]
PolList = ['X', 'Y']
#-------- Plots
if 'plotMin' not in locals(): plotMin = 0.0
if 'plotMax' not in locals(): plotMax = 1.2
for scan_index, scan in enumerate(scanList):
    print('--Plotting scan%d' % (scan))
    pp = PdfPages('BP_%s_REF%s_Scan%d.pdf' % (prefix, refantName, scan))
    plotSP(pp, prefix, antList, spwList, FreqList, BPList, plotMin, plotMax) 
#
