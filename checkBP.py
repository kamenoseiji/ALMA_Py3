#---- Script for Band-3 Astroholograpy Data
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import GetBaselineIndex, CrossCorrAntList, Ant2Bl, Ant2BlD, Bl2Ant, indexList, ANT0, ANT1, bestRefant, bunchVec, GetAntName, GetUVW, GetChNum, BPtable
from Plotters import plotXYP, plotBP, plotSP
#-------- Procedures
msfile = wd + prefix + '.ms'
Antenna1, Antenna2 = GetBaselineIndex(msfile, spwList[0], BPscan)
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)[UseAntList]
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
#-------- Configure Array
print('---Checking array configulation in scan %d' % (BPscan))
if 'BPPLOT' not in locals(): BPPLOT = False
if 'antFlag' not in locals(): antFlag = []
flagAnt = indexList(antFlag, antList)
UseAnt = list(set(range(antNum)) - set(flagAnt)); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
try:
    refantID = np.where(antList[UseAnt] == refant )[0][0]
except:
    ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
    for bl_index in list(range(UseBlNum)): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
    timeStamp, UVW = GetUVW(msfile, spwList[0], BPscan)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
#
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
FreqList, BPList, XYList, XYdelayList = [], [], [], []
if 'bunchNum' not in locals(): bunchNum = 1
for spw_index, spw in enumerate(spwList):
    if 'FGprefix' in locals():  # Flag table
        try:
            FG = np.load('%s-SPW%d.FG.npy' % (FGprefix, spw))
            TS = np.load('%s-SPW%d.TS.npy' % (FGprefix, spw))
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, bunchNum, FG, TS)
        except:
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, bunchNum )
        #
    else:
        if spw_index == 0:
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, bunchNum)
        else :
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spw, BPscan, blMap, blInv, bunchNum, np.array([]), np.array([]), Gain)
        #
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
for spw_index, spw in enumerate(spwList):
    np.save('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spw), BPList[spw_index]) 
    np.save('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spw), XYList[spw_index]) 
    np.save('%s-REF%s-SC%d-SPW%d-XYdelay.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spw), XYdelayList[spw_index]) 
#
#-------- Plots
if BPPLOT:
    if XYsnr > 0.0:
        pp = PdfPages('XYP_%s_REF%s_Scan%d.pdf' % (prefix, antList[UseAnt[refantID]], BPscan))
        plotXYP(pp, prefix, spwList, XYList, bunchNum) 
    #
    pp = PdfPages('BP_%s_REF%s_Scan%d.pdf' % (prefix, antList[UseAnt[refantID]], BPscan))
    if 'spurRFLists' in locals():
        plotBP(pp, prefix, antList[antMap], spwList, BPscan, BPList, bunchNum, 1.2, spurRFLists) 
    else:
        if 'plotMin' not in locals(): plotMin = 0.0
        if 'plotMax' not in locals(): plotMax = 1.2
        plotSP(pp, prefix, antList[antMap], spwList, FreqList, BPList, plotMin, plotMax)
#
