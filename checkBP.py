#---- Script for Band-3 Astroholograpy Data
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import GetBaselineIndex, CrossCorrAntList, Ant2Bl, Ant2BlD, Bl2Ant, indexList, ANT0, ANT1, bestRefant, bunchVec, GetAntName, GetUVW, GetChNum, BPtable
from ASDM_XML import SPW_FULL_RES, spwMS, spwIDMS
from Plotters import plotXYP, plotBP, plotSP
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
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
parser.add_option('-X', dest='XYLog', metavar='XYLog',
    help='Record XY delay in a logfile (*.XYdelay.log)', action="store_true")
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
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
FG      = options.FG
NPY     = options.NPY
XYLog   = options.XYLog
refant  = options.refant
#-------- Procedures
if XYLog: xyLog = open(prefix + '.XYdelay.log', 'w')
msfile = prefix + '.ms'
SPWdic = SPW_FULL_RES(prefix) if os.path.isdir(prefix) else spwMS(msfile)
SPWdic = spwIDMS(SPWdic, msfile)
if len(spwList) == 0:   # Use all available spws
    spwList = [spw for spw in SPWdic.keys() if SPWdic[spw]['chNum'] > 4 and len(SPWdic[spw]['scanList']) > 0]
#
if len(scanList) == 0:  # Use all available scans
    for spw in spwList: scanList = scanList + SPWdic[spw]['scanList']
    scanList = list(set(scanList)); scanList.sort()
#
refSPW = spwList[0]
Antenna1, Antenna2 = GetBaselineIndex(msfile, refSPW, SPWdic[refSPW]['scanList'][0])
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)[UseAntList]
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
#-------- Configure Array
print('---Checking array configulation in scan %d' % (scanList[0]))
flagAnt = indexList(antFlag, antList)
UseAnt = list(set(range(antNum)) - set(flagAnt)); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
if refant not in antList[UseAnt]: refant = ''
if refant == '':
    timeStamp, UVW = GetUVW(msfile, refSPW, SPWdic[refSPW]['scanList'][0])
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
else:
    refantID = np.where(antList[UseAnt] == refant)[0][0]
print( '  Use %s as the refant.' % (antList[UseAnt[refantID]]))
#-------- Baseline Mapping
print('---Baseline Mapping')
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
for bl_index in range(UseBlNum):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
print( '  %d baselines are inverted.' % (len(np.where( blInv )[0])))
#-------- Bandpass Table
print('---Generating antenna-based bandpass table')
SideBand = ['LSB', 'USB']
for BPscan in scanList:
    spws = [spw for spw in spwList if BPscan in SPWdic[spw]['scanList']]
    FreqList, BPList, XYList, XYdelayList = [], [], [], []
    for spw_index, spw in enumerate(spws):
        if FG:  # Flag table
            FG = np.load('%s-SPW%d.FG.npy' % (prefix, spw))
            TS = np.load('%s-SPW%d.TS.npy' % (prefix, spw))
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, chBin, FG, TS)
        else:
            if spw_index == 0:
                BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, chBin)
            else :
                BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, chBin, np.array([]), np.array([]), Gain)
        #
        BPList = BPList + [BP_ant]
        XYList = XYList + [XY_BP]
        chNum, chWid, Freq = GetChNum(msfile, spw)
        chNum, chWid, Freq = int(chNum / chBin), chWid* chBin, bunchVec(Freq, chBin)
        if NPY: np.save('%s-SPW%d-Freq.npy' % (prefix, spw), Freq) 
        FreqList = FreqList + [Freq]
        BW = chNum* np.median(chWid)    # Bandwidth
        XYdelayList = XYdelayList + [0.5e9* XYD / BW]
        text_sd = 'Scan%2d SPW%2d BB%d: [%s] XY delay = %+.3f [ns] : SNR = %.1f' % (BPscan, spw, SPWdic[spw]['BB']+1, SideBand[int((np.sign(np.median(chWid))+1)/2)], -0.5* XYD / (BW * 1.0e-9), XYsnr)
        print(text_sd)
        if XYLog: xyLog.write(text_sd + '\n')
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
            plotBP(pp, prefix, antList[antMap], spws, BPscan, BPList, chBin, 1.2, spurRFLists) 
        else:
            plotSP(pp, prefix, antList[antMap], spws, FreqList, BPList, plotMin, plotMax)
    #
if XYLog: xyLog.close()
