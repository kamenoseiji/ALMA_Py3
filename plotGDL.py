#-------- Plot Multi-BB Group Delay
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
parser.add_option('-M', dest='plotRange', metavar='plotRange',
    help='Max y-axis range in [ps] e.g. 10.0', default='10.0')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
plotMax = float(options.plotRange)
plotAntList = [] if options.plotAnt == '' else [ant for ant in options.plotAnt.split(',')]
'''
plotMax = 10.0
plotAntList = [] 
prefix = 'uid___A002_X12e95ca_X7383'
'''
DL = np.load(prefix + '-GDL.npy')
antList = np.load(prefix + '.Ant.npy')
antNum = len(antList)
if plotAntList == []: plotAntList = antList
antMap = indexList(plotAntList, antList)
plotAntList = antList[antMap].tolist()
plotAntNum = len(plotAntList)
rowNum = int(min(4, np.floor(np.sqrt(plotAntNum))))
#-------- Load tables
DT, TS = [], DL[0]
for mjdSec in TS.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
polLabel = ['X', 'Y']
polSymbol = ['x', '1']
cmap = plt.get_cmap("tab10")
pp = PdfPages('DL_%s.pdf' %  (prefix))
figDL = plt.figure(figsize = (8, 11))
figDL.suptitle(prefix + ' Multi-BB Group Residual Delay')
for row_index, ant in enumerate(plotAntList):
    ant_index = np.where(antList == ant)[0][0]
    DLPL = figDL.add_subplot( int(np.ceil(plotAntNum/rowNum)), rowNum, row_index + 1 )
    DLPL.yaxis.offsetText.set_fontsize(3)
    DLPL.tick_params(labelsize=4); DLPL.tick_params(axis='x')
    #DLPL.set_ylabel('Delay [ps]', fontsize=6)
    DLPL.axis([np.min(DT), np.max(DT), -plotMax, plotMax])
    DLPL.grid(axis='y', linestyle='--', color='gray')
    text_sd = ant
    for pol_index, polchar in enumerate(polSymbol):
        plotDL = 1.0e3* DL[1+ (pol_index + 2)*antNum + ant_index]   # Delay in [ps] # Delay in [ps]
        DLPL.plot( DT, plotDL, '-', color=cmap(pol_index), linewidth=0.25)
        DLPL.plot( DT, plotDL, polSymbol[pol_index], markersize=4, color=cmap(pol_index), label='%s-%s' % (ant, polLabel[pol_index]))
        text_sd = text_sd + ' %.1f ' % (np.std(plotDL))
    DLPL.text(DT[0], 0.8*plotMax, text_sd, fontsize=5)
    #if row_index == 0: DLPL.legend(loc = 'best', prop={'size' :4}, numpoints = 1)
figDL.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
figDL.text(0.03, 0.45, 'Multi-BB Group Delay [ps]', rotation=90)
figDL.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/plotGDL.py to generate this plot', fontsize=4)
figDL.savefig(pp, format='pdf')
plt.close('all')
pp.close()
