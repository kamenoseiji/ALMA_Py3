#---- Script for Band-3 Astroholograpy Data
import math
import analysisUtils as au
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from interferometry import GetBaselineIndex, CrossCorrAntList, GetAntName, GetTimerecord, GetChNum, GetVisAllBL, Bl2Ant, bunchVec
from casatools import quanta as qatool
qa = qatool()
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-b', dest='chBunch', metavar='chBunch',
    help='Channel binning  e.g. 2', default='1')
parser.add_option('-c', dest='scanID', metavar='scanID',
    help='Scan ID  e.g. 2', default='2')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
parser.add_option('-S', dest='startTime', metavar='startTime',
    help='Start time e.g. 2020-03-03T14:11:25', default='')
parser.add_option('-i', dest='integTime', metavar='integTime',
    help='Integration [s] e.g. 60', default='')
parser.add_option('-M', dest='plotMax', metavar='plotMax',
    help='Max amplitude to plot', default='')
(options, args) = parser.parse_args()
#-------- Definitions
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
BPscan  =  int(options.scanID)
spwList = [int(spw) for spw in options.spwList.split(',')]
chBunch =  int(options.chBunch)
if options.plotMax != ''  : plotMax = float(options.plotMax)
if options.startTime != '': startMJD = qa.convert(options.startTime, 's')['value']
if options.integTime != '': integTime = float(options.integTime)
msfile = prefix + '.ms'
Antenna1, Antenna2 = GetBaselineIndex(msfile, spwList[0], BPscan)
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)[UseAntList]
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
spwNum  = len(spwList)
#-------- Procedures
blMap = list(range(blNum))
#-------- Procedures
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], BPscan)
integDuration = np.median(interval)
timeNum = len(timeStamp)
if 'integTime' in locals(): timeNum = int(np.ceil(integTime / integDuration))
timeNum = min(timeNum, len(timeStamp))
#
#-------- Prepare Plots
for spw_index in range(spwNum):
    figInch = max(16,antNum)
    fontSize = min(32, figInch)
    #w, h = plt.figaspect(1)
    figSPW = plt.figure(figsize=(figInch,figInch))
    figSPW.text(0.475, 0.05, 'Frequency [GHz]', fontsize=fontSize)
    figSPW.text(0.05, 0.5, 'Phase [rad]', rotation=90, fontsize=fontSize)
    figSPW.text(0.95, 0.5, 'Amplitude', rotation=-90, fontsize=fontSize)
    #-------- Plot BP
    print(' Loading SPW = %d' % (spwList[spw_index]))
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPscan)
    if chBunch > 1:
        chNum, chWid, Freq = int(chNum / chBunch), chWid* chBunch, bunchVec(Freq, chBunch)
        for ch_index in list(range(chNum)):
            sourceRange = list(range(ch_index*chBunch, (ch_index+1)*chBunch))
            Pspec[:,ch_index] = np.mean(Pspec[:,sourceRange], axis=1)
            Xspec[:,ch_index] = np.mean(Xspec[:,sourceRange], axis=1)
        Pspec = Pspec[:,0:chNum]
        Xspec = Xspec[:,0:chNum]
    #---- integration timerange
    startMJD = min(max(startMJD, timeStamp[0]), timeStamp[-1]) if 'startMJD' in locals() else timeStamp[0]
    endMJD = min(timeStamp[-1], startMJD + timeNum* integDuration)
    st_index, timeNum = np.argmin(abs(timeStamp - startMJD)), int((endMJD - startMJD + 0.1*integDuration) / integDuration)
    timeRange = list(range(st_index, st_index + timeNum))
    text_timerange = au.call_qa_time('%fs' % (startMJD), form='fits', prec=6) + ' - ' + au.call_qa_time('%fs' % (endMJD), form='fits', prec=6)
    print('Integration in %s (%.1f sec)' % (text_timerange, timeNum* integDuration))
    figSPW.suptitle('%s SPW=%d Scan=%d Integration in %s (%.1f sec)' % (prefix, spwList[spw_index], BPscan, text_timerange, timeNum* integDuration), fontsize=fontSize)
    #---- polarization format
    polNum = Xspec.shape[0]
    if polNum == 4: pol = [0,3]; polName = ['X', 'Y']         # parallel pol in full-pol correlations
    if polNum == 2: pol = [0,1]; polName = ['X', 'Y']         # XX and YY
    if polNum == 1: pol = [0]; polName = ['X']           # Only XX
    polColor = ['b','g']
    polNum = len(pol)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap][:,:,:,timeRange]
    Pspec = Pspec[pol]; Pspec = Pspec[:,:,:,timeRange]
    tempVis = np.mean(Xspec, axis=3)    # tempVis[pol, ch, bl]
    tempAC  = np.mean(Pspec, axis=3)    # tempVis[pol, ch, ant]
    pMax = np.percentile(abs(tempVis), 98) if 'plotMax' not in locals() else plotMax
    aMax = np.percentile(abs(tempAC), 98)
    polColor = ['b', 'g']
    for bl_index in list(range(blNum)):
        ants = Bl2Ant(bl_index)
        BLamp = figSPW.add_subplot(antNum, antNum, ants[1]*antNum + ants[0] + 1)
        BLphs = figSPW.add_subplot(antNum, antNum, ants[0]*antNum + ants[1] + 1)
        for pol_index in list(range(polNum)):
            plotVis = tempVis[pol_index, :, bl_index]
            BLamp.step(Freq, abs(plotVis), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
            BLphs.plot( Freq, np.angle(plotVis), '.', color=polColor[pol_index], label = 'Pol=' + polName[pol_index])
        #
        if np.max(abs(plotVis)) > 0.1:
            print('Out of Range : %s - %s' % (antList[ants[0]], antList[ants[1]]))
        BLamp.axis([np.min(Freq), np.max(Freq), 0.0, 1.25*pMax])
        BLphs.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
        BLamp.xaxis.set_major_locator(plt.NullLocator())
        if ants[1] == 0:    # Antenna label in the top and leftside
            BLamp.set_title( antList[ants[0]] )
            BLphs.set_ylabel( antList[ants[0]] )
        else:
            BLphs.yaxis.set_major_locator(plt.NullLocator())
        #
        if ants[0] == antNum - 1:    # Antenna at rightside
            BLamp.yaxis.tick_right()
        else:
            BLamp.yaxis.set_major_locator(plt.NullLocator())
        if ants[0] < antNum - 1:    # except bottom panel : skip drawing X-axis
            BLphs.xaxis.set_major_locator(plt.NullLocator())
        #
    #-------- Plot autocorrelations
    for ant_index in list(range(antNum)):
        BLamp = figSPW.add_subplot(antNum, antNum, ant_index*antNum + ant_index + 1)
        BLamp.patch.set_facecolor('pink')
        for pol_index in list(range(polNum)):
            BLamp.step(Freq, abs(tempAC[pol_index, :, ant_index]), color=polColor[pol_index], where='mid', label = polName[pol_index])
        #
        BLamp.axis([np.min(Freq), np.max(Freq), 0.0, 1.25*aMax])
        if ant_index < antNum-1:
            BLamp.xaxis.set_major_locator(plt.NullLocator())
        else:
            BLamp.yaxis.tick_right()
        if ant_index > 0:
            BLamp.yaxis.set_major_locator(plt.NullLocator())
        else: 
            BLamp.set_title( antList[0])
            BLamp.set_ylabel( antList[0])
            BLamp.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
        #
    #
    figSPW.text(0.1, 0.05, 'See https://github.com/kamenoseiji/ALMA_Py3/wiki/checkXspec.py to generate this plot', fontsize=fontSize)
    plt.show()
    pngFile = 'PS_%s_Scan%d_SPW%d' % (prefix, BPscan, spwList[spw_index])
    #pdfFile = pngFile + '.pdf'
    figSPW.savefig(pngFile + '.png', format='png', dpi=72)
    #figSPW.savefig(pdfFile, format='pdf', dpi=144)
    plt.close('all')
    #os.system('pdftoppm -png %s %s' % (pdfFile, pngFile))
#
