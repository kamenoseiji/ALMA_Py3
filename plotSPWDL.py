import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from casatools import quanta as qatool
import analysisUtils as au
import datetime
from interferometry import indexList, RADDEG
from Plotters import polColor, polName
qa = qatool()
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-p', dest='plotAnt', metavar='plotAnt',
    help='Antennas to plot', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPWs to plot e.g. 17,19,21', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scans to plot e.g. 3,5,7', default='')
parser.add_option('-R', dest='refant', metavar='refant',
    help='reference antenna', default='')
parser.add_option('-M', dest='plotRange', metavar='plotRange',
    help='Max y-axis range in [ns] e.g. 0.01', default='0.01')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
spwList = [int(spw) for spw in options.spwList.split(',')]
scanList = [int(scan) for scan in options.scanList.split(',')]
refant = options.refant
plotMax = float(options.plotRange)
plotAntList = [] if options.plotAnt == '' else [ant for ant in options.plotAnt.split(',')]
antList = np.load(prefix + '.Ant.npy')
if plotAntList == []: plotAntList = antList
antMap = indexList(plotAntList, antList)
plotAntList = antList[antMap].tolist()
plotAntNum, columnNum = len(plotAntList), len(spwList)
#-------- Load tables
polLabel = ['X', 'Y']
polSymbol = ['x', '1']
cmap = plt.get_cmap("tab10")
pp = PdfPages('DL_%s.pdf' %  (prefix))
figDL = plt.figure(figsize = (8, 11))
figDL.suptitle(prefix + ' Single-Band Delay')
for row_index, ant in enumerate(plotAntList):
    ant_index = np.where(antList == ant)[0][0]
    DLPL = figDL.add_subplot( plotAntNum, 1, row_index + 1 )
    DLPL.yaxis.offsetText.set_fontsize(3)
    DLPL.tick_params(labelsize=4); DLPL.tick_params(axis='x')
    DLPL.set_ylabel('Delay [ns] w.r.t. Scan %d' % (scanList[0]), fontsize=6)
    DLPL.axis([np.min(scanList), np.max(scanList), -plotMax, plotMax])
    for pol_index, polchar in enumerate(polSymbol):
        for spw_index, spw in enumerate(spwList):
            plotDL = np.zeros([len(spwList), 2, len(scanList)])
            for scan_index, scan in enumerate(scanList):
                DL =  np.load('%s-REF%s-SC%d-SPW%d-DL.npy' % (prefix, refant, scan, spw))
                if scan_index == 0:
                    DL0 = DL[ant_index, pol_index]
                else:
                    plotDL[spw_index, pol_index, scan_index] = DL[ant_index, pol_index] - DL0
    DLPL.plot( np.array(scanList), np.mean(plotDL, axis=(0,1)), '-', color=cmap(spw_index))
    DLPL.plot( np.array(scanList), np.mean(plotDL, axis=(0,1)), 'o', markersize=2, color=cmap(spw_index), label='SPW%d' % (spw))
    if row_index == 0: DLPL.legend(loc = 'best', prop={'size' :4}, numpoints = 1)
    DLPL.text(scanList[0], 0.9*plotMax, ant, fontsize=6)
figDL.text(0.45, 0.05, 'Scan')
figDL.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/plotSPWDL.py to generate this plot', fontsize=4)
figDL.savefig(pp, format='pdf')
plt.close('all')
pp.close()
