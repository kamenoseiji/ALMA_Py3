#---- Multi-band bandpass
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import GetBaselineIndex, CrossCorrAntList, Ant2Bl, Ant2BlD, Bl2Ant, indexList, ANT0, ANT1, bestRefant, bunchVec, GetAntName, GetUVW, GetChNum, BPtable, GetVisAllBL, gainComplexVec, ParaPolBL, RADDEG
from ASDM_XML import SPW_FULL_RES, spwMS, spwIDMS
from Plotters import plotXYP, plotBP, plotSP, polColor, polColor
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
parser.add_option('-B', dest='BPprefix', metavar='BPprefix',
    help='Bandpass table to apply', default='')
parser.add_option('-b', dest='chBin', metavar='chBin',
    help='Channel binning', default='1')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List  e.g. 3,6', default='')
parser.add_option('-f', dest='FG', metavar='FG',
    help='Apply flagging', action="store_true")
parser.add_option('-g', dest='GAspw', metavar='GAspw',
    help='SPW for gain table to apply', default='')
parser.add_option('-m', dest='plotMin', metavar='plotMin', type="float",
    help='Plot range minimum', default=0.0)
parser.add_option('-M', dest='plotMax', metavar='plotMax', type="float",
    help='Plot range minimum', default=1.2)
parser.add_option('-r', dest='refscan', metavar='refscan',
    help='Reference (bandpass) scan', default='3')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
parser.add_option('-P', dest='BPPLOT', metavar='BPPLOT',
    help='Plot PDF', action="store_true")
parser.add_option('-R', dest='refant', metavar='refant',
    help='reference antenna', default='')
parser.add_option('-D', dest='delayMessage', metavar='delayMessage',
    help='Print residual delay)', action="store_true")
parser.add_option('-X', dest='XYLog', metavar='XYLog',
    help='Record XY delay in a logfile (*.XYdelay.log)', action="store_true")
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
BPprefix= options.BPprefix
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
if not os.path.isdir(prefix + '.ms'):
    if not os.path.isdir(prefix):
        os.system('rm -rf %s.ms' % (prefix))
        os.system('asdmExport %s' % (prefix))
    importasdm(prefix)
antFlag = [ant for ant in options.antFlag.split(',')]
scanList  = [] if options.scanList == '' else [int(scan) for scan in options.scanList.split(',')]
chBin=  int(options.chBin)
spwList = [] if options.spwList == '' else [int(spw) for spw in options.spwList.split(',')]
plotMin = options.plotMin
plotMax = options.plotMax
BPPLOT  = options.BPPLOT
delayMessage  = options.delayMessage
FG      = options.FG
XYLog   = options.XYLog
refant  = options.refant
GAspw   = int(options.GAspw)
GAprefix = BPprefix
refScan = int(options.refscan)
'''
BPprefix = 'uid___A002_X12e95ca_X7383.WVR'
GAprefix = 'uid___A002_X12e95ca_X7383.WVR'
prefix = 'uid___A002_X12e95ca_X7383.WVR'
antFlag = []
refScan = 3
scanList = [26]
#spwList = [15,17,19,21]
spwList = [0,2,4,6]
GAspw   = 0
chBin = 1
plotMin = 0.0
plotMax = 1.5
XYLog = False
BPPLOT = True
FG     = False
refant = 'DA45'
'''
#--------
def delayFit(freq, spec):
    chNum = len(freq)
    normspec = np.array((spec.real/abs(spec)).tolist() +  (spec.imag/abs(spec)).tolist())
    refFreq  = np.array((freq - np.median(freq)).tolist() + (freq - np.median(freq)).tolist())                     # Frequency relative to the center
    param = np.array([np.median(np.angle(spec)),0.0 ])  # initial phase and delay
    P = np.array([(-normspec[chNum:2*chNum]).tolist() + normspec[0:chNum].tolist(), (-normspec[chNum:2*chNum]).tolist() + normspec[0:chNum].tolist()])
    P[1] = P[1]* 2.0* np.pi* refFreq
    for iteration in list(range(3)):
        resid = np.array((normspec[0:chNum] - np.cos( param[0] + 2.0* np.pi* param[1]* refFreq[0:chNum])).tolist() + (normspec[chNum:2*chNum] - np.sin( param[0] + 2.0* np.pi* param[1]* refFreq[chNum:2*chNum])).tolist())
        correction = np.linalg.inv(P.dot(P.T)).dot(P.dot(resid))
        param = param + correction
    return param
#
#-------- bandpass and gain tables to apply
msfile = prefix + '.ms'
BPantList = np.load('%s-REF%s.Ant.npy' % (BPprefix, refant))
BPfileList, BPfreqList = [], []
for spw in spwList:
    BPfileList = BPfileList + ['%s-REF%s-SC%d-SPW%d-BPant.npy' % (BPprefix, refant, refScan, spw)]
    BPfreqList = BPfreqList + ['%s-SPW%d-Freq.npy' % (BPprefix, spw)]
GAfileName = '%s-SPW%d.GA.npy' % (GAprefix, GAspw)
TSfileName = '%s-SPW%d.TS.npy' % (GAprefix, GAspw)
#-------- Procedures
if os.path.isdir(prefix + '/SpectralWindow.xml'):
    SPWdic = SPW_FULL_RES(prefix)
    SPWdic = spwIDMS(SPWdic, msfile)
else:
    SPWdic = spwMS(msfile)
if len(spwList) == 0:   # Use all available spws
    spwList = [spw for spw in SPWdic.keys() if SPWdic[spw]['chNum'] > 4 and len(SPWdic[spw]['scanList']) > 0]
#
if len(scanList) == 0:  # Use all available scans
    for spw in spwList: scanList = scanList + SPWdic[spw]['scanList']
    scanList = list(set(scanList)); scanList.sort()
#
antList = GetAntName(msfile)
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
#-------- Configure Array
print('---Checking array configulation in scan %d' % (refScan))
BPantMap = indexList(np.delete(BPantList, indexList(antFlag, BPantList)), BPantList)
antMap = indexList( BPantList[BPantMap], antList)
flagAnt = indexList(antFlag, antList)
UseAnt = list(set(range(antNum)) - set(flagAnt)); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
#-------- Baseline Mapping
print('---Baseline Mapping')
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print( '  %d baselines are inverted.' % (len(np.where( blInv )[0])))
#-------- Gain Table
GA = np.load(GAfileName)
TS = np.load(TSfileName)
#-------- Bandpass Table
print('---Loading antenna-based bandpass table')
SideBand = ['LSB', 'USB']
FreqList, BPList, Freq1D = [], [], []
for spw_index, spw in enumerate(spwList):
    BPList   = BPList + [np.load(BPfileList[spw_index])]
    FreqList = FreqList + [np.load(BPfreqList[spw_index])]
    SPWdic[spw]['chRange'] = list(range(int(0.05*SPWdic[spw]['chNum']), int(0.95*SPWdic[spw]['chNum'])))
    Freq1D = Freq1D + np.load(BPfreqList[spw_index])[SPWdic[spw]['chRange']].tolist()
Freq1D = 1.0e-9* np.array(Freq1D)
#-------- Load Visiblities
GDL = np.zeros([4* antNum + 1,len(scanList)])
for scan_index, scan in enumerate(scanList):
    caledBPList = []
    pp = PdfPages('PHS_%s_REF%s_Scan%d.pdf' % (prefix, BPantList[0], scan))
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle('%s Scan %d Antenna-based phase' % (prefix, scan))
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Phase [deg]', rotation=90)
    for spw_index, spw in enumerate(spwList):
        print('---Loading visibilities Scan %d SPW %d' % (scan, spw))
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)    # Xspec[pol, ch, bl, time]
        GAscan = GA[:,:,indexList(timeStamp, TS)]                   # GAscan[ant, pol, time]
        Xspec  = ParaPolBL(Xspec[:,:,blMap], blInv)
        BPant  = BPList[spw_index][BPantMap]                        # BPant[ant, pol, ch]
        CaledXspec = np.mean(GAscan[ant0].conjugate()* GAscan[ant1]* Xspec.transpose(1,2,0,3), axis=3).transpose(1,2,0) # CaledXspec[bl, pol, ch]
        CaledXspec = CaledXspec / (BPant[ant0]* BPant[ant1].conjugate())
        caledBPList = caledBPList + [np.array([gainComplexVec(CaledXspec[:,0]), gainComplexVec(CaledXspec[:,1])]).transpose(1,0,2)]
    GDL[0, scan_index] = np.median(timeStamp)
    for ant_index, ant in enumerate(antList[antMap]):
        PhsPL = figAnt.add_subplot(1, 1, 1 )
        PhsPL.set_title(ant)
        BPX, BPY = [], []
        for spw_index, spw in enumerate(spwList):
            BPX = BPX + caledBPList[spw_index][ant_index,0][SPWdic[spw]['chRange']].tolist()
            BPY = BPY + caledBPList[spw_index][ant_index,1][SPWdic[spw]['chRange']].tolist()
        BPX, BPY = np.array(BPX), np.array(BPY)
        paramX, paramY = delayFit(Freq1D, BPX), delayFit(Freq1D, BPY)
        GDL[1 +            ant_index, scan_index] = paramX[0]
        GDL[1 +   antNum + ant_index, scan_index] = paramY[0]
        GDL[1 + 2*antNum + ant_index, scan_index] = paramX[1]
        GDL[1 + 3*antNum + ant_index, scan_index] = paramY[1]
        PhsPL.plot(Freq1D, RADDEG* np.angle(BPX), '.', color=polColor[0], label = 'Pol-%s DL=%+.3f [ps] phs=%.1f [deg]' % ('X', 1e3*paramX[1], RADDEG*paramX[0]))
        PhsPL.plot(Freq1D, RADDEG* np.angle(BPY), '.', color=polColor[1], label = 'Pol-%s DL=%+.3f [ps] phs=%.1f [deg]' % ('Y', 1e3*paramY[1], RADDEG*paramY[0]))
        text_sd = 'Scan%d %s : multi-delay = %6.3f %6.3f [ps]  %6.1f %6.1f [deg]' % (scan, ant, 1e3*paramX[1], 1e3*paramY[1], RADDEG*paramX[0], RADDEG*paramY[0])
        print(text_sd)
        PhsPL.axis([np.min(Freq1D), np.max(Freq1D), -30.0, 30.0])
        PhsPL.tick_params(axis='both', labelsize=9)
        PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
        figAnt.savefig(pp, format='pdf')
        figAnt.delaxes(PhsPL)
    #
    plt.close('all')
    pp.close()
np.save('%s-GDL.npy' % (prefix), GDL)
