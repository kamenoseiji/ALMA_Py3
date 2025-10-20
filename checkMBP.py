#---- Multi-band bandpass
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import GetBaselineIndex, CrossCorrAntList, Ant2Bl, Ant2BlD, Bl2Ant, indexList, ANT0, ANT1, bestRefant, bunchVec, GetAntName, GetUVW, GetChNum, BPtable, GetVisAllBL, gainComplexVec, ParaPolBL
from ASDM_XML import SPW_FULL_RES, spwMS, spwIDMS
from Plotters import plotXYP, plotBP, plotSP
from optparse import OptionParser
'''
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
parser.add_option('-m', dest='plotMin', metavar='plotMin', type="float",
    help='Plot range minimum', default=0.0)
parser.add_option('-M', dest='plotMax', metavar='plotMax', type="float",
    help='Plot range minimum', default=1.2)
parser.add_option('-n', dest='NPY', metavar='NPY',
    help='Generate numpy files (.npy)', action="store_true")
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
NPY     = options.NPY
XYLog   = options.XYLog
refant  = options.refant
'''
BPprefix = 'uid___A002_X12e95ca_X7383'
GAprefix = 'uid___A002_X12e95ca_X7383'
prefix = 'uid___A002_X12e95ca_X7383'
antFlag = ['DA50','DV24']
refScan = 3
scanList = [28]
spwList = [15,17,19,21]
GAspw   = 15
chBin = 1
plotMin = 0.0
plotMax = 1.5
XYLog = False
BPPLOT = True
FG     = False
refant = 'DA45'
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
GA = np.load(GAfileName)[antMap]
TS = np.load(TSfileName)
#-------- Bandpass Table
print('---Loading antenna-based bandpass table')
SideBand = ['LSB', 'USB']
FreqList, BPList = [], []
for spw_index, spw in enumerate(spwList):
    BPList   = BPList + [np.load(BPfileList[spw_index])]
    FreqList = FreqList + [np.load(BPfreqList[spw_index])]
#-------- Load Visiblities
print('---Loading visibilities')
for scan_index, scan in scanList:
    for spw_index, spw in enumerate(spwList):
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)    # Xspec[pol, ch, bl, time]
        Xspec  = ParaPolBL(Xspec[:,:,blMap], blInv)
        GAscan = GA[:,:,indexList(timeStamp, TS)]                   # GAscan[ant, pol, time]
        BPant  = BPList[spw_index][BPantMap]                        # BPant[ant, pol, ch]
        CaledXspec = np.mean(GAscan[ant0].conjugate()* GAscan[ant1]* Xspec.transpose(1,2,0,3), axis=3).transpose(1,2,0) # CaledXspec[bl, pol, ch]
        CaledXspec = CaledXspec / (BPant[ant0]* BPant[ant1].conjugate())
'''



        BP_ant[:,0], BP_ant[:,1] = (gainComplexVec(CaledXspec[:,:,0].T), gainComplexVec(CaledXspec[:,:,1].T))
    #
    ppolNum = BPList[0].shape[1]
    PolList = ['X', 'Y']
    #
    #-------- Save CalTables
    if NPY:
        np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '.Ant.npy', antList[antMap]) 
        for spw_index, spw in enumerate(spws):
            np.save('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spw), BPList[spw_index]) 
            np.save('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spw), XYList[spw_index]) 
            np.save('%s-REF%s-SC%d-SPW%d-XYdelay.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spw), XYdelayList[spw_index]) 
    #
    #-------- Plots
    if BPPLOT:
        if XYsnr > 0.0:
            pp = PdfPages('XYP_%s_REF%s_Scan%d.pdf' % (prefix, antList[UseAnt[refantID]], BPscan))
            plotXYP(pp, prefix, spws, XYList, XYdelayList, chBin) 
        #
        pp = PdfPages('BP_%s_REF%s_Scan%d.pdf' % (prefix, antList[UseAnt[refantID]], BPscan))
        if 'spurRFLists' in locals():
            delay_ant = plotBP(pp, prefix, antList[antMap], spws, BPscan, BPList, chBin, 1.2, spurRFLists) 
        else:
            delay_ant = plotSP(pp, prefix, antList[antMap], spws, FreqList, BPList, plotMin, plotMax, delayMessage)
        if NPY:
            for spw_index, spw in enumerate(spws): np.save('%s-REF%s-SC%d-SPW%d-DL.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spw), delay_ant[:,spw_index])
    #
if XYLog: xyLog.close()
'''
