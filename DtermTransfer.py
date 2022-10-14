#exec(open(SCR_DIR + 'interferometry.py').read())
#exec(open(SCR_DIR + 'Grid.py').read())
#exec(open(SCR_DIR + 'Plotters.py').read())
import os
import re
import numpy as np
import scipy
from interferometry import ANT0, ANT1, Ant2BlD, ALMA_lat, BANDPA, indexList, GetAzEl, AzElMatch, GetBaselineIndex, CrossCorrAntList, CrossPolBL, GetAntName, GetSourceList, GetVisAllBL, AzEl2PA, polariGain, PAVector, TransferD
from Grid import sourceRename
from matplotlib.backends.backend_pdf import PdfPages
import pickle
import analysisUtils as au
pattern = r'RB_..'
#----------------------------------------- Procedures
if 'antFlag' not in locals():   antFlag = []
spwNum = len(spwList)
polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
#
msfile = wd + prefix + '.ms'
sourceList = []
sources, posList = GetSourceList(msfile); sourceList = sourceList + sourceRename(sources)
sourceList = np.unique(sourceList).tolist()
azelTime, AntID, AZ, EL = GetAzEl(msfile)
Antenna1, Antenna2 = GetBaselineIndex(msfile, spwList[0], scanList[0])
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)
refAntID = np.where(antList == refantName)[0][0]
antList = antList[UseAntList]
trkAntList = [refantName] + [ant for ant in antList if ant not in scanAntList + [refantName]]
antNum  = len(antList); blNum = int(antNum * (antNum - 1)/2)
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
trkAntMap = indexList(trkAntList, antList)
scnAntMap = indexList(scanAntList, antList)
msmd.open(msfile)
spwName = msmd.namesforspws(spwList)[0]; BandName = re.findall(pattern, spwName)[0]; bandID = int(BandName[3:5])
BandPA = (BANDPA[bandID] + 90.0)*np.pi/180.0
for spw_index, spw in enumerate(spwList):
    #-------- XY phase
    XYphase = np.load('%s-SPW%d-%s.XYPH.npy' % (prefix, spw, refantName))
    #-------- Stokes parameters
    StokesFP = open('Stokes.%s-SPW%d.dic' % (prefix, spw), 'rb')
    StokesDic = pickle.load(StokesFP)
    StokesFP.close()
    sourceList = list(StokesDic.keys())
    #-------- Load D-term for trakking antennas
    DxSpec, DySpec = [], []
    for antName in trkAntList:
        DtermData = np.load('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, antName))
        DxSpec = DxSpec + [DtermData[1] + (0.0 + 1.0j)* DtermData[2]]
        DySpec = DySpec + [DtermData[3] + (0.0 + 1.0j)* DtermData[4]]
    DxSpec, DySpec = np.array(DxSpec), np.array(DySpec)
    #-------- Load BP table 
    BPantList, BP_ant = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (BPprefix, refantName, BPscan, spw))
    XYspec = np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (XYprefix, refantName, BPscan, spw))
    GA = np.load('%s-SPW%d-%s.GA.npy' % (XYprefix, spw, refantName))
    TS = np.load('%s-SPW%d-%s.TS.npy' % (XYprefix, spw, refantName))
    print('Apply XY phase into Y-pol Bandpass.'); BP_ant[:,1] *= XYspec  # XY phase cal
    BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()
    BPantMap = indexList(BPantList,antList)
    BPblMap, BPblInv= list(range(blNum)), [False]* blNum
    for bl_index in list(range(blNum)):
        BPblMap[bl_index], BPblInv[bl_index]  = Ant2BlD(BPantMap[ant0[bl_index]], BPantMap[ant1[bl_index]])
    #-------- Load visibilities
    mjdSec, Az, El, PA, timeNum, scanST, XspecList = [], [], [], [], [], [], []
    scanST = scanST + [0]
    for scan_index, scan in enumerate(scanList):
        sourceName = sourceRename(sources)[msmd.sourceidforfield(msmd.fieldsforscan(scan)[0])]
        if(len(StokesDic[sourceName]) < 4): continue
        print('Scan%d %s : %.2f %.2f %.2f %.2f' % (scan, sourceName, StokesDic[sourceName][0],StokesDic[sourceName][1],StokesDic[sourceName][2],StokesDic[sourceName][3]))
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)  # Xspec[POL, CH, BL, TIME]
        del Pspec
        scanAz, scanEl = AzElMatch(timeStamp, azelTime, AntID, refAntID, AZ, EL)
        scanPA = AzEl2PA(scanAz, scanEl, ALMA_lat) + BandPA
        mjdSec, Az, El, PA = mjdSec + timeStamp.tolist(), Az + scanAz.tolist(), El + scanEl.tolist(), PA + scanPA.tolist()
        XspecList = XspecList + [Xspec]
        del Xspec
        timeNum = timeNum + [len(timeStamp)]
        if scan_index != 0: scanST = scanST + [scanST[scan_index - 1] + timeNum[scan_index - 1]]
    #
    #-------- Stack visibilities
    chNum = XspecList[0].shape[1]; chRange = list(range(int(0.05*chNum), int(0.95*chNum)))
    VisSpec = np.zeros([4, chNum, blNum, len(mjdSec)], dtype=complex)
    timeIndex = 0
    for scan_index, Xspec in enumerate(XspecList):
        VisSpec[:,:,:,timeIndex:timeIndex + timeNum[scan_index]] = CrossPolBL(Xspec[:,:,BPblMap], BPblInv)
        timeIndex = timeIndex + timeNum[scan_index]
    #-------- Apply bandpass and gain cal
    CS, SN = np.cos(2.0* np.array(PA)), np.sin(2.0* np.array(PA))
    QCpUS = CS* StokesDic[sourceName][1]/StokesDic[sourceName][0] - SN* StokesDic[sourceName][2]/StokesDic[sourceName][0]
    PS = PAVector(np.array(PA), np.ones(len(mjdSec))).transpose(2,0,1).dot(StokesDic[sourceName]).T/StokesDic[sourceName][0]
    VisSpec = (VisSpec.transpose(3, 2, 0, 1) / BP_bl).transpose(3,2,1,0)
    GainCaledVisSpec = VisSpec / (GA[polYindex][:,ant0][:,:,indexList( np.array(mjdSec), TS )]* GA[polXindex][:,ant1][:,:,indexList( np.array(mjdSec), TS )].conjugate()) # GainCaledVisSpec[ch, pol, bl, time]
    #-------- For each scan antenna
    for ant_index, scanAnt in enumerate(scanAntList):
        Dx, Dy = np.zeros(chNum, dtype=complex), np.zeros(chNum, dtype=complex)
        print('---- Dterm transfer for %s' % (scanAnt))
        scanAntIndex = np.where(antList[BPantMap] == scanAnt)[0][0]
        trkScnAntMap = indexList(np.array(trkAntList), antList[BPantMap])
        trkScnBLMap, trkScnInv = [], []
        for trkAntIndex in trkScnAntMap:
            blIndex, blInv = Ant2BlD(scanAntIndex, trkAntIndex)
            trkScnBLMap, trkScnInv = trkScnBLMap + [blIndex], trkScnInv + [blInv]
        #
        trcScanVisSpec = CrossPolBL( GainCaledVisSpec[:,:,trkScnBLMap].transpose(1,0,2,3), trkScnInv)
        for ch_index in list(range(chNum)):
            Dx[ch_index], Dy[ch_index] = TransferD(np.mean(trcScanVisSpec[:,ch_index], axis=1), np.mean(DxSpec[:, ch_index]), np.mean(DySpec[:, ch_index]), PS)
        #
        DtermData[1:5] = np.array([Dx.real, Dx.imag, Dy.real, Dy.imag])
        np.save('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, scanAnt), DtermData)
    #
#
msmd.done()
msmd.close()
