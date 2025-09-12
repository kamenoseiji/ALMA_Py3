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
parser.add_option('-M', dest='plotRange', metavar='plotRange',
    help='Max y-axis range e.g. 60', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
spwList = [int(spw) for spw in options.spwList.split(',')]
plotMax =  180 if options.plotRange == '' else float(options.plotRange)
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
pp = PdfPages('PH_%s.pdf' %  (prefix))
figPhs = plt.figure(figsize = (8, 11))
figPhs.suptitle(prefix + ' Gain Phase')
for row_index, ant in enumerate(plotAntList):
    ant_index = np.where(antList == ant)[0][0]
    timeFile = '%s-SPW%d.TS.npy' % (prefix, spwList[0])
    timeStamp = np.load(timeFile)
    DT = [datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f') for mjdSec in timeStamp.tolist()]
    for pol_index, polchar in enumerate(polSymbol):
        PhsPL = figPhs.add_subplot( plotAntNum, 2, 2*row_index + pol_index + 1 )
        PhsPL.yaxis.offsetText.set_fontsize(3)
        PhsPL.tick_params(labelsize=4); PhsPL.tick_params(axis='x', rotation=45)
        PhsPL.axis([np.min(DT), np.max(DT), -plotMax, plotMax])
        PhsPL.set_ylabel('Phase [deg] w.r.t. SPW%d' % (spwList[0]), fontsize=6)
        for spw_index, spw in enumerate(spwList):
            Gain, FG = np.load('%s-SPW%d.GA.npy' % (prefix, spw)), np.load('%s-SPW%d.FG.npy' % (prefix, spw))
            flag_index = np.where(FG[ant_index] > 0.01)[0].tolist()
            if spw_index == 0: GA0 = Gain[ant_index, pol_index, flag_index]*Gain[ant_index, pol_index, 0].conjugate()
            plotPhase = np.angle(Gain[ant_index, pol_index, flag_index]*Gain[ant_index, pol_index, 0].conjugate() / GA0)*RADDEG
            PhsPL.plot( np.array(DT)[flag_index], plotPhase, '-', linewidth=0.5, color=cmap(spw_index))
            PhsPL.plot( np.array(DT)[flag_index], plotPhase, polchar, markersize=2, color=cmap(spw_index), label='SPW%d Pol%s' % (spw, polLabel[pol_index]))
        if row_index == 0: PhsPL.legend(loc = 'best', prop={'size' :4}, numpoints = 1)
    #
    PhsPL.text(np.min(DT), 0.9*plotMax, ant, fontsize=6)
figPhs.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
figPhs.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/plotSPWphase.py to generate this plot', fontsize=4)
figPhs.savefig(pp, format='pdf')
plt.close('all')
pp.close()
