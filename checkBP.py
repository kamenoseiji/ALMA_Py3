#---- Script for Band-3 Astroholograpy Data
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import GetBaselineIndex, CrossCorrAntList, Ant2Bl, Ant2BlD, Bl2Ant, indexList, ANT0, ANT1, bestRefant, bunchVec, GetAntName, GetUVW, GetChNum, BPtable
from Plotters import plotXYP, plotBP, plotSP
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
parser.add_option('-b', dest='bunchNum', metavar='bunchNum',
    help='Channel binning', default='1')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan List  e.g. 3,6', default='3')
parser.add_option('-f', dest='FG', metavar='FG',
    help='Apply flagging', action="store_true")
parser.add_option('-m', dest='plotMin', metavar='plotMin', type="float",
    help='Plot range minimum', default=0.0)
parser.add_option('-M', dest='plotMax', metavar='plotMax', type="float",
    help='Plot range minimum', default=1.2)
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
parser.add_option('-P', dest='BPPLOT', metavar='BPPLOT',
    help='Plot PDF', action="store_true")
parser.add_option('-R', dest='refant', metavar='refant',
    help='reference antenna', default='')
#
(options, args) = parser.parse_args()
#-------- BB_spw : BB power measurements for list of spws, single antenna, single scan
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
antFlag = [ant for ant in options.antFlag.split(',')]
scanList  = [int(scan) for scan in options.scanList.split(',')]
bunchNum=  int(options.bunchNum)
spwList = [] if options.spwList == '' else [int(spw) for spw in options.spwList.split(',')]
plotMin = options.plotMin
plotMax = options.plotMax
BPPLOT  = options.BPPLOT
FG      = options.FG
refant  = options.refant
#-------- spwList
def spwListForBandpass( msfile, scan ):
    msmd.open(msfile)
    scanSPWs    = msmd.spwsforscan(scan)
    fullResSPWs = msmd.tdmspws().tolist() + msmd.fdmspws().tolist()
    msmd.close()
    spwList = list(set(scanSPWs) & set(fullResSPWs))
    spwList.sort()
    return spwList
#-------- Procedures
msfile = prefix + '.ms'
refSPW = spwListForBandpass(msfile, scanList[0])[0] if spwList == [] else spwList[0]
Antenna1, Antenna2 = GetBaselineIndex(msfile, refSPW, scanList[0])
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
    #ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
    #for bl_index in list(range(UseBlNum)): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
    timeStamp, UVW = GetUVW(msfile, refSPW, scanList[0])
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
    spws = spwListForBandpass(msfile, BPscan) if spwList == [] else spwList
    FreqList, BPList, XYList, XYdelayList = [], [], [], []
    for spw_index, spw in enumerate(spws):
        if FG:  # Flag table
            FG = np.load('%s-SPW%d.FG.npy' % (prefix, spw))
            TS = np.load('%s-SPW%d.TS.npy' % (prefix, spw))
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, bunchNum, FG, TS)
        else:
            if spw_index == 0:
                BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, bunchNum)
            else :
                BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, bunchNum, np.array([]), np.array([]), Gain)
        #
        BPList = BPList + [BP_ant]
        XYList = XYList + [XY_BP]
        XYdelayList = XYdelayList + [XYD]
        chNum, chWid, Freq = GetChNum(msfile, spw)
        chNum, chWid, Freq = int(chNum / bunchNum), chWid* bunchNum, bunchVec(Freq, bunchNum)
        np.save('%s-SPW%d-Freq.npy' % (prefix, spw), Freq) 
        FreqList = FreqList + [Freq]
        BW = chNum* np.median(chWid)    # Bandwidth
        print('SPW%2d: [%s] XY delay = %+f [ns] : SNR = %f' % (spw, SideBand[int((np.sign(np.median(chWid))+1)/2)], 0.5* XYD / (BW * 1.0e-9), XYsnr))
    #
    ppolNum = BPList[0].shape[1]
    PolList = ['X', 'Y']
    #
    #-------- Save CalTables
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
            plotXYP(pp, prefix, spws, XYList, bunchNum) 
        #
        pp = PdfPages('BP_%s_REF%s_Scan%d.pdf' % (prefix, antList[UseAnt[refantID]], BPscan))
        if 'spurRFLists' in locals():
            plotBP(pp, prefix, antList[antMap], spws, BPscan, BPList, bunchNum, 1.2, spurRFLists) 
        else:
            plotSP(pp, prefix, antList[antMap], spws, FreqList, BPList, plotMin, plotMax)
    #
#
