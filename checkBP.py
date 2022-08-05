#---- Script for Band-3 Astroholograpy Data
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
#
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
spwNum = len(spwList)
if 'bunchNum' not in locals(): bunchNum = 1
for spw_index in list(range(spwNum)):
    if 'FGprefix' in locals():  # Flag table
        try:
            FG = np.load('%s-SPW%d.FG.npy' % (FGprefix, spwList[spw_index]))
            TS = np.load('%s-SPW%d.TS.npy' % (FGprefix, spwList[spw_index]))
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPscan, blMap, blInv, bunchNum, FG, TS)
            #BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPscan, blMap, blInv, bunchNum, FG)
        except:
            BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPscan, blMap, blInv, bunchNum )
        #
    else:
        BP_ant, XY_BP, XYD, Gain, XYsnr = BPtable(msfile, spwList[spw_index], BPscan, blMap, blInv, bunchNum)
    #
    BPList = BPList + [BP_ant]
    XYList = XYList + [XY_BP]
    XYdelayList = XYdelayList + [XYD]
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index])
    chNum, chWid, Freq = int(chNum / bunchNum), chWid* bunchNum, bunchVec(Freq, bunchNum)
    np.save('%s-SPW%d-Freq.npy' % (prefix, spwList[spw_index]), Freq) 
    FreqList = FreqList + [Freq]
    BW = chNum* np.median(chWid)    # Bandwidth
    print('SPW%2d: [%s] XY delay = %+f [ns] : SNR = %f' % (spwList[spw_index], SideBand[int((np.sign(np.median(chWid))+1)/2)], 0.5* XYD / (BW * 1.0e-9), XYsnr))
#
ppolNum = BPList[0].shape[1]
PolList = ['X', 'Y']
#
#-------- Save CalTables
np.save(prefix + '-REF' + antList[UseAnt[refantID]] + '.Ant.npy', antList[antMap]) 
for spw_index in list(range(spwNum)):
    np.save('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spwList[spw_index]), BPList[spw_index]) 
    np.save('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spwList[spw_index]), XYList[spw_index]) 
    np.save('%s-REF%s-SC%d-SPW%d-XYdelay.npy' % (prefix, antList[UseAnt[refantID]], BPscan, spwList[spw_index]), XYdelayList[spw_index]) 
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
