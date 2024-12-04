import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
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
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
spwList = [int(spw) for spw in options.spwList.split(',')]
plotAntList = [] if options.plotAnt == '' else [ant for ant in options.plotAnt.split(',')]
antList = np.load(prefix + '.Ant.npy')
if plotAntList == []: plotAntList = antList
antMap = indexList(plotAntList, antList)
plotAntList = antList[antMap].tolist()
plotAntNum, columnNum = len(plotAntList), len(spwList)
#-------- Load tables
pp = PdfPages('GA_%s.pdf' %  (prefix))
figAmp, figPhs = plt.figure(figsize = (8, 11)), plt.figure(figsize = (8, 11))
figAmp.suptitle(prefix + ' Gain Amplitude'); figPhs.suptitle(prefix + ' Gain Phase')
for col_index, spw in enumerate(spwList):
    timeFile = '%s-SPW%d.TS.npy' % (prefix, spw)
    timeStamp = np.load(timeFile)
    DT = [datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f') for mjdSec in timeStamp.tolist()]
    Gain, FG = np.load('%s-SPW%d.GA.npy' % (prefix, spw)), np.load('%s-SPW%d.FG.npy' % (prefix, spw))
    polList = polName[:Gain.shape[1]]
    plotMin, plotMax = 0.0, 1.1* np.percentile(np.max(abs(Gain[antMap]), axis=1)* FG[antMap], 90)
    for row_index, ant in enumerate(plotAntList):
        #-------- Plot Gain
        plot_index = row_index* columnNum + col_index
        ant_index = np.where(antList == ant)[0][0]
        flag_index = np.where(FG[ant_index] > 0.01)[0].tolist()
        amp_text, phs_text = 'SPW%d Gain(median) = (' % (spw), 'SPW%d Phase(rms) = (' % (spw)
        for pol_index, pol in enumerate(polList):
            amp_text = amp_text + '%.2f%% ' % (100.0* np.median(abs(Gain[ant_index, pol_index, flag_index])) )
            phs_text = phs_text + '%.1f deg ' % (RADDEG* np.std(np.angle(Gain[ant_index, pol_index, flag_index])) )
        amp_text, phs_text = amp_text[:-1] + ')', phs_text[:-1] + ')'
        print(amp_text)
        AmpPL = figAmp.add_subplot( plotAntNum, columnNum, plot_index + 1 )
        PhsPL = figPhs.add_subplot( plotAntNum, columnNum, plot_index + 1 )
        for pol_index, pol in enumerate(polList): AmpPL.plot( np.array(DT)[flag_index], abs(Gain[ant_index, pol_index, flag_index]), '.', markersize=3, color=polColor[pol_index], label=pol)
        for pol_index, pol in enumerate(polList): PhsPL.plot( np.array(DT)[flag_index], np.angle(Gain[ant_index, pol_index, flag_index]*Gain[ant_index, pol_index, 0].conjugate())*RADDEG, '.', markersize=3, color=polColor[pol_index], label=pol)
        AmpPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
        AmpPL.yaxis.offsetText.set_fontsize(3)
        PhsPL.yaxis.offsetText.set_fontsize(3)
        AmpPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
        AmpPL.tick_params(labelsize=4); AmpPL.tick_params(axis='x', rotation=45)
        PhsPL.tick_params(labelsize=4); PhsPL.tick_params(axis='x', rotation=45)
        AmpPL.axis([np.min(DT), np.max(DT), plotMin, plotMax])
        PhsPL.axis([np.min(DT), np.max(DT), -180.0, 180.0])
        AmpPL.text( 0.05, 1.02, amp_text, transform=AmpPL.transAxes, fontsize=5)
        PhsPL.text( 0.05, 1.02, phs_text, transform=PhsPL.transAxes, fontsize=5)
        if col_index == 0: AmpPL.set_ylabel(ant, fontsize=6); PhsPL.set_ylabel(ant, fontsize=6)
        if row_index < plotAntNum - 1: AmpPL.set_xticks([]); PhsPL.set_xticks([])
        if plot_index == 0:
            AmpPL.legend(loc = 'best', prop={'size' :4}, numpoints = 1)
            PhsPL.legend(loc = 'best', prop={'size' :4}, numpoints = 1)
figAmp.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d'))); figPhs.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
figAmp.text(0.03, 0.45, 'Gain Amplitude = sqrt(correlated flux / SEFD)', rotation=90); figPhs.text(0.03, 0.45, 'Gain Phase [deg]', rotation=90)
figAmp.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/plotGain.py to generate this plot', fontsize=4)
figPhs.text(0.1, 0.03, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/plotGain.py to generate this plot', fontsize=4)
figAmp.savefig(pp, format='pdf'); figPhs.savefig(pp, format='pdf')
plt.close('all')
pp.close()
