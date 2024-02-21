#import sys
#import pickle
import math
import numpy as np
import analysisUtils as au
import xml.etree.ElementTree as ET
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import Tcmb, kb, BANDPA, BANDFQ, indexList, subArrayIndex, GetAntName, GetAntD, GetSourceList, GetBandNames, GetAeff, GetDterm, quadratic_interpol, GetAtmSPWs, GetBPcalSPWs, GetSPWFreq, GetOnSource, GetUVW, loadScanSPW, AzElMatch, gainComplexErr, bestRefant, ANT0, ANT1, Ant2Bl, Ant2BlD, Bl2Ant, gainComplexVec, CrossPolBL, CrossPolBP, SPWalign, delay_search, linearRegression, VisMuiti_solveD, AllanVarPhase, specBunch
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from Plotters import plotSP, plotXYP
from Grid import *
from ASDM_XML import CheckCorr, BandList
from PolCal import GetAMAPOLAStokes, GetSSOFlux, PolResponse
msfile = wd + prefix + '.ms'
#-------- Check Correlator Type
BLCORR = True
if 'ACA' in CheckCorr(prefix): BLCORR = False
#-------- Check Receivers
RXList = BandList(prefix)
antList = GetAntName(msfile)
#-------- Check SPWs of atmCal and bandpass
print('---Checking SPWs and Scan information')
atmSPWs, bpSPWs = GetAtmSPWs(msfile), GetBPcalSPWs(msfile)
if 'spwFlag' in locals():
    atmSPWs = list(set(atmSPWs)- set(spwFlag));atmSPWs.sort()
    bpSPWs  = list(set(bpSPWs) - set(spwFlag)); bpSPWs.sort()
bandNameList = GetBandNames(msfile, atmSPWs)
for BandName in RXList:
    if BandName not in bandNameList: RXList.remove(BandName)
#
BandPA  = dict(zip(RXList, [[]]*len(RXList)))    # Band PA
BandbpSPW  = dict(zip(RXList, [[]]*len(RXList))) # Band SPW for visibilitiies
BandatmSPW = dict(zip(RXList, [[]]*len(RXList))) # Band SPW for atmCal
BandScanList = dict(zip(RXList, [[]]*len(RXList))) # Band scan list
#
OnScanList = GetOnSource(msfile)
if len(OnScanList) == 0: RXList = []
if 'antFlag' not in locals(): antFlag = []
msmd.open(msfile)
for BandName in RXList:
    BandPA[BandName] = (BANDPA[int(BandName[3:5])] + 90.0)*math.pi/180.0
    BandbpSPW[BandName]  = {'spw': np.array(bpSPWs)[np.where(np.array(bandNameList) == BandName)[0].tolist()].tolist()}
    BandatmSPW[BandName] = {'spw': np.array(atmSPWs)[np.where(np.array(bandNameList) == BandName)[0].tolist()].tolist()}
    BandScanList[BandName] = list(set(msmd.scansforspw(BandbpSPW[BandName]['spw'][0])) & set(OnScanList))
    BandScanList[BandName].sort()
    #---- Bandpass scan to check Allan Variance
    def AV(vis): return AllanVarPhase(np.angle(vis), 1)
    checkScan = msmd.scansforintent('*BANDPASS*')
    if len(checkScan) == 0: checkScan = msmd.scansforintent('*POINTING*')
    if len(checkScan) == 0: checkScan = [msmd.scansforintent('*PHASE*')[0]]
    checkScan = checkScan[-1]
    print('---Checking usable antennas by ASD in Scan %d' % (checkScan))
    chavSPWs = list((set(msmd.chanavgspws()) - set(msmd.almaspws(sqld=True)) - set(msmd.almaspws(wvr=True))) & set(msmd.spwsforscan(checkScan)))
    timeStampList, XspecList = loadScanSPW(msfile, chavSPWs, [checkScan])  # XspecList[spw][scan] [corr, ch, bl, time]
    parapolIndex = [0,3] if XspecList[0][0].shape[0] == 4 else [0,1]
    bunchNum = 8 if XspecList[0][0].shape[3] > 30 else int(XspecList[0][0].shape[3]/3)
    for spw_index, spw in enumerate(chavSPWs):
        checkVis = XspecList[spw_index][0][parapolIndex][:,0]
        timeRange = list(range(checkVis.shape[2] % bunchNum, checkVis.shape[2]))
        checkVis = [specBunch(checkVis[0][:,timeRange], 1, bunchNum), specBunch(checkVis[1][:,timeRange], 1, bunchNum)]
        AV_bl = np.apply_along_axis(AV, 1, checkVis[0]) + np.apply_along_axis(AV, 1, checkVis[1])
        errBL = list(set(np.where(AV_bl > 2.0)[0]) | set(np.where(np.median(abs(checkVis[0]), axis=1) > 15.0*np.median(abs(checkVis[0])))[0]) | set(np.where(np.median(abs(checkVis[1]), axis=1) > 15.0*np.median(abs(checkVis[1])))[0]))
        errCount = np.zeros(Bl2Ant(len(AV_bl))[0])
        for bl in errBL: errCount[list(Bl2Ant(bl))] += 1
        antFlag = list(set(antFlag + antList[np.where(errCount > len(antFlag)+2 )[0].tolist()].tolist()))
    #
msmd.close()
BandbpSPW = GetSPWFreq(msfile, BandbpSPW)   # BandbpSPW[BandName] : [[SPW List][freqArray][chNum][BW]]
BandatmSPW = GetSPWFreq(msfile, BandatmSPW)
#
#-------- Tsys measurement
if len(antFlag) < len(antList) - 3: exec(open(SCR_DIR + 'TsysCal.py').read())
else: RXList = []
if 'Tau0med' in locals():
    if any(Tau0med < 0.0): RXList = []
#-------- Check Antenna List
antDia = GetAntD(antList)
antNum = len(antList)
blNum = int(antNum* (antNum - 1) / 2)
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0; useAntMap = np.where(flagAnt > 0.1)[0].tolist()
flagBL = flagAnt[ANT0[0:blNum]]* flagAnt[ANT1[0:blNum]]; useBlMap = np.where(flagBL > 0.1)[0].tolist()
polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
#-------- Check source list
print('---Checking source list')
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
FscaleDic = dict(zip(np.array(sourceList)[SSOList].tolist(), [[]]* len(SSOList)))
#-------- Loop for Bands
for BandName in RXList:
    if BandName not in ['RB_03', 'RB_04', 'RB_06', 'RB_07']: continue
    #-------- Prepare log files
    logfile = open(prefix + '-' + BandName + '-Flux.log', 'w')
    ingestFile = open(prefix + '-' + BandName + '-Ingest.log', 'w')
    text_corr = ' %s BLCORR' % (prefix) if BLCORR else ' %s ACACORR' % (prefix)
    logfile.write(text_corr + '\n')
    print('-----%s----' % (BandName))
    #-------- D-term from history
    Dcat = GetDterm(TBL_DIR, antList, int(BandName[3:5]), np.mean(azelTime))
    #-------- Load Aeff file
    etaA = 0.01* GetAeff(TBL_DIR, antList, int(UniqBands[band_index][3:5]), np.mean(azelTime)).T
    print('-----Estimation from AMAPOLA and Butler-JPL-Horizons')
    #-------- Load Visibilities into memory
    timeStampList, XspecList = loadScanSPW(msfile, BandbpSPW[BandName]['spw'], BandScanList[BandName])  # XspecList[spw][scan] [corr, ch, bl, time]
    StokesDic = GetAMAPOLAStokes(R_DIR, SCR_DIR, sourceList, qa.time('%fs' % (timeStampList[0][0]), form='ymd')[0], BANDFQ[int(BandName[3:5])])
    StokesDic, SSODic = GetSSOFlux(StokesDic, qa.time('%fs' % (timeStampList[0][0]), form='ymd')[0], [1.0e-9* np.median(BandbpSPW[BandName]['freq'][spw_index]) for spw_index, spw in enumerate(BandbpSPW[BandName]['spw'])])
    #-------- Polarization responses per scan
    scanDic = PolResponse(msfile, StokesDic, BandPA[BandName], BandScanList[BandName], timeStampList)
    QSOscanList = [scan for scan in scanDic.keys() if scanDic[scan]['source'][0] == 'J' and str.isdigit(scanDic[scan]['source'][1])]
    AzScanList, ElScanList = [], []
    #-------- Check AZEL
    for scan_index, scan in enumerate(BandScanList[BandName]):
        AzScan, ElScan = AzElMatch(timeStampList[scan_index], azelTime, AntID, 0, AZ, EL)
        AzScanList, ElScanList = AzScanList + [AzScan], ElScanList + [ElScan]
        scanDic[scan]['SA'] = SunAngleSourceList[sourceList.index(scanDic[scan]['source'])]
    #-------- Apply Tsys calibration
    scanDic, XspecList = applyTsysCal(prefix, BandName, BandbpSPW[BandName], scanDic, SSODic, XspecList)
    #-------- Check usable antennas and refant
    print('-----Filter usable antennas')
    chRange = BandbpSPW[BandName]['chRange'][0]
    checkScan   = QSOscanList[np.argmax(np.array([np.median(abs(scanDic[scan]['UCmQS'])* np.sin(scanDic[scan]['EL'] - ELshadow))* np.exp(-np.median(scanDic[scan]['Tau'][0])) for scan in QSOscanList]))]
    if not np.mean(np.array(scanDic[checkScan]['Tau'])) > -0.5 : continue
    checkSource = scanDic[checkScan]['source']
    print('-----Check Scan %d : %s' % (checkScan, checkSource))
    Xspec       = XspecList[spw_index][BandScanList[BandName].index(checkScan)][:,:,useBlMap]
    checkVis    = np.mean(Xspec[[0,3]][:,chRange], axis=1) / scanDic[checkScan]['I']
    Gain =  np.array([gainComplexVec(checkVis[0]), gainComplexVec(checkVis[1])])
    antCoh = np.array([abs(Gain[1,ant_index].dot(Gain[0,ant_index].conjugate())) for ant_index, ant in enumerate(useAntMap)]) / Gain.shape[2]
    Aeff = 8.0* kb* antCoh / (np.pi* antDia[useAntMap]**2)
    for ant_index, ant in enumerate(useAntMap):
        if abs(Aeff[ant_index] - np.median(Aeff)) > 0.25: flagAnt[ant] *= 0.0
    useAntMap = np.where(flagAnt > 0.0)[0].tolist(); useAntNum = len(useAntMap)
    useBlNum  = int(useAntNum* (useAntNum - 1) / 2); flagBL = flagAnt[ANT0[0:blNum]]* flagAnt[ANT1[0:blNum]]; useBlMap = np.where(flagBL > 0.1)[0].tolist()
    text_sd = '  Usable antennas (%d) : ' % (len(useAntMap))
    for ants in antList[useAntMap].tolist(): text_sd = text_sd + ants + ' '
    print(text_sd)
    if useAntNum < 4: continue
    print('-----Select reference antenna')
    timeStamp, UVW = GetUVW(msfile, BandbpSPW[BandName]['spw'][0], checkScan)
    uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2); uvDist[useBlMap] *= 0.001
    refantID = bestRefant(uvDist)
    print('Use %s as refant' % (antList[refantID]))
    antMap = [refantID] + list(set(useAntMap) - set([refantID]))
    useAntNum = len(antMap); useBlNum  = int(useAntNum* (useAntNum - 1) / 2)
    blMap, blInv = list(range(useBlNum)), np.ones(useBlNum)
    for bl_index, bl in enumerate(useBlMap): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ANT0[bl_index]], antMap[ANT1[bl_index]])
    ant0, ant1 = ANT0[0:useBlNum], ANT1[0:useBlNum]
    #-------- Remap baseline ordering
    for scan_index, scan in enumerate(BandScanList[BandName]):
        timeStamp, UVW = GetUVW(msfile, BandbpSPW[BandName]['spw'][0], scan)
        scanDic[scan]['UVW'] = UVW[:,blMap]
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            Xspec = CrossPolBL(XspecList[spw_index][scan_index][:,:,blMap], blInv)
            XspecList[spw_index][scan_index] = Xspec
    print('-----Bandpass to align SPWs and polarization')
    #-------- Bandpass using checkScan
    FreqList, BPList, spwGainList = [], [], []
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        Xspec = XspecList[spw_index][BandScanList[BandName].index(checkScan)]
        BP_ant, BPCaledXYSpec, XYdelay, Gain, XYsnr = CrossPolBP(Xspec)
        BPList = BPList + [BP_ant]
        spwGainList = spwGainList + [Gain]
    #-------- SPW phase offsets
    spwTwiddle = SPWalign(np.array(spwGainList))
    #-------- Phase-aligned bandpass table
    print('-----SPW-aligned bandpass in scan %d : %s' % (checkScan, scanDic[checkScan]['source']))
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        BPList[spw_index] = (BPList[spw_index].transpose(2,0,1) * spwTwiddle[:,:,spw_index]).transpose(1,2,0)
    pp = PdfPages('BP-%s-%s.pdf' % (prefix,BandName))
    plotSP(pp, prefix, antList[antMap], BandbpSPW[BandName]['spw'], BandbpSPW[BandName]['freq'], BPList)
    os.system('rm -rf B0')
    bandpass(msfile, caltable='B0', spw=','.join(map(str, BandbpSPW[BandName]['spw'])), scan=str(checkScan), refant=antList[refantID], solnorm=True, minblperant=3)
    #-------- SPW-combined phase calibration
    print('-----Antenna-based gain correction')
    text_sd = '        coherence loss % :'
    for antName in antList[antMap]: text_sd = text_sd + ' ' +  antName
    print(text_sd)
    for scan_index, scan in enumerate(BandScanList[BandName]):
        text_sd = '----- Scan%3d %10s :' % (scan, scanDic[scan]['source'])
        chAvgList = []
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            chRange = BandbpSPW[BandName]['chRange'][spw_index]
            BP_ant = BPList[spw_index]
            BPcaledSpec = XspecList[spw_index][scan_index][[0,3]].transpose(3,2,0,1) / (BP_ant[ant0]* BP_ant[ant1].conjugate())
            chAvgVis    = np.mean(BPcaledSpec[:,:,:,chRange], axis=(2,3))
            chAvgList = chAvgList + [chAvgVis]
        #
        scanGain = gainComplexVec(np.mean(np.array(chAvgList), axis=0).T)
        scanFlag = np.unique(np.where(np.max(abs(scanGain), axis=0)  <  min([1.0, 5.0*np.median(abs(scanGain))]))).tolist()
        scanDic[scan]['Flag'] = scanFlag
        scanDic[scan]['Gain'] = scanGain[:,scanFlag]
        scanDic[scan]['mjdSec'] = scanDic[scan]['mjdSec'][scanFlag]
        scanDic[scan]['EL']     = scanDic[scan]['EL'][scanFlag]
        scanDic[scan]['PA']     = scanDic[scan]['PA'][scanFlag]
        scanDic[scan]['QCpUS']  = scanDic[scan]['QCpUS'][scanFlag]
        scanDic[scan]['UCmQS']  = scanDic[scan]['UCmQS'][scanFlag]
        scanDic[scan]['UVW']    = scanDic[scan]['UVW'][:,:,scanFlag]
        coh = abs(np.mean(scanGain, axis=1)) / np.mean(abs(scanGain), axis=1)
        for ant_index, ant in enumerate(antList[antMap]): text_sd = text_sd + ' %.2f' % (100.0* (1.0 - coh[ant_index]))
        print(text_sd)
    #
    #-------- Scan-by-scan bandpass
    BPavgScanList, BPList, XYList, XYWList = [], [], [], []
    text_sd = '-----Scan-by-scan BP     :'
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']): text_sd = text_sd + '     SPW%2d     ' % (spw)
    print(text_sd)
    text_sd = '----- XY delay           :'
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']): text_sd = text_sd + '  [ns]  ( SNR )'
    print(text_sd)
    for scan_index, scan in enumerate(BandScanList[BandName]):
        if scan not in QSOscanList : continue              # filter by QSO
        text_sd = '----- Scan%3d %10s :' % (scan, scanDic[scan]['source'])
        scanFlag  = scanDic[scan]['Flag']
        if len(scanFlag) == 0: continue
        scanPhase = scanDic[scan]['Gain'] / abs(scanDic[scan]['Gain'])
        BPSPWList, XYSPWList, XYsnrList = [], [], []
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            chRange = BandbpSPW[BandName]['chRange'][spw_index]
            Xspec = XspecList[spw_index][scan_index][:,:,:,scanFlag]         # Xspec[pol, ch, bl, time]
            XPspec = np.mean(Xspec/ (scanPhase[ant0]* scanPhase[ant1].conjugate()), axis=3)   # antenna-based phase correction
            BP_ant = np.array([gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[3].T)])
            #---- Amplitude normalization
            for pol_index in [0,1]:
                BP_eq_gain = np.mean(abs(BP_ant[pol_index][:,chRange]), axis=1)
                BP_ant[pol_index] = (BP_ant[pol_index].T / BP_eq_gain).T
            BPCaledXspec = XPspec.transpose(0, 2, 1)/(BP_ant[polYindex][:,ant0]* BP_ant[polXindex][:,ant1].conjugate())
            BPCaledXYSpec = np.mean(BPCaledXspec[1], axis=0) +  np.mean(BPCaledXspec[2], axis=0).conjugate()
            XYdelay, XYsnr = delay_search( BPCaledXYSpec[chRange] )
            XYdelay = (float(chNum) / float(len(chRange)))* XYdelay
            text_sd = text_sd + ' %6.3f (%5.1f)' % (0.5* XYdelay/(BandbpSPW[BandName]['BW'][spw_index]*1.0e-9), XYsnr)
            XYsnrList = XYsnrList + [XYsnr if XYsnr > 7.0 else 0.001]
            BPSPWList = BPSPWList + [BP_ant.transpose(1,0,2)]
            XYSPWList = XYSPWList + [BPCaledXYSpec]
        #
        print(text_sd)
        BPavgScanList = BPavgScanList + [scan]
        BPList = BPList + [BPSPWList]
        XYList = XYList + [XYSPWList]
        XYWList=XYWList + [XYsnrList]
        #pp = PdfPages('BP-%s-%s-%d.pdf' % (prefix, BandName, scan))
        #plotSP(pp, prefix, antList[antMap], BandbpSPW[BandName]['spw'], BandbpSPW[BandName]['freq'], BPSPWList)
        #pp = PdfPages('XYP_%s_REF%s_Scan%d.pdf' % (prefix, antList[antMap[0]], scan))
        #plotXYP(pp, prefix, BandbpSPW[BandName]['spw'], XYSPWList)
    #-------- Average bandpass
    XYW = np.array(XYWList)**2
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        chRange = BandbpSPW[BandName]['chRange'][spw_index]
        refXY = XYList[np.argmax( XYW[:,spw_index])][spw_index]
        XY = 0.0* refXY
        BP = 0.0* BPList[0][spw_index]
        for scan_index, scan in enumerate(BPavgScanList):
            if len(scanDic[scan]['Flag']) < 3: continue
            #-------- average BP
            BPW = abs(np.mean(scanDic[scan]['Gain'], axis=1)) / np.median(np.std(np.angle(scanDic[scan]['Gain']), axis=1))
            BPW[ np.where(BPW < 0.2)[0].tolist() ] *= 0.1 
            BP  = BP + (BPList[scan_index][spw_index].transpose(1,2,0)* BPW**2).transpose(2,0,1)
            #-------- average XY
            XYspec = XYList[scan_index][spw_index]
            Weight = XYW[scan_index][spw_index] * np.sign( XYspec[chRange].dot(refXY[chRange].conjugate()))
            XY     = XY + Weight* XYspec
        #---- Save averaged BP
        BPSPWList[spw_index] = (BP.transpose(2,0,1) / np.mean(abs(BP[:,:,chRange]), axis=2)).transpose(1,2,0)
        BPSPWList[spw_index][:,1] *= (XY / abs(XY))
        #---- Save into CASA caltable
        if os.path.isfile('B0'):
            tb.open('B0', nomodify=False)
            SPWQ = tb.query('SPECTRAL_WINDOW_ID == %d'%(spw))
            BPtable = SPWQ.getcol('CPARAM')
            for ant_index, ant in enumerate(antMap): BPtable[:,:,ant] = BPSPWList[spw_index][ant_index]
            SPWQ.putcol('CPARAM', BPtable)
            SPWQ.close()
            tb.close()
        #
    #
    pp = PdfPages('BP-%s-%s-%d.pdf' % (prefix, BandName, 0))
    plotSP(pp, prefix, antList[antMap], BandbpSPW[BandName]['spw'], BandbpSPW[BandName]['freq'], BPSPWList, 0.0, 1.2, True)
    del BPavgScanList, BPList, XYList, XYWList, XYW, XY, refXY, BP, XYSPWList, XYsnrList
    #-------- XY sign in checkScan
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        BP_ant = BPSPWList[spw_index].transpose(1,2,0)
        scanFlag  = scanDic[checkScan]['Flag']
        scanPhase = scanDic[checkScan]['Gain'] / abs(scanDic[checkScan]['Gain'])
        Xspec = XspecList[spw_index][list(scanDic.keys()).index(checkScan)][:,:,:,scanFlag]* (scanPhase[ant1]* scanPhase[ant0].conjugate())
        BPCaledXspec = (Xspec.transpose(3, 0, 1, 2) / (BP_ant[polYindex][:,:,ant0]* BP_ant[polXindex][:,:,ant1].conjugate())).transpose(1,2,3,0)
        BPCaledXY    = np.mean(BPCaledXspec[1][chRange], axis=(0,1)) +  np.mean(BPCaledXspec[2][chRange], axis=(0,1)).conjugate()
        XYphase = np.angle(scanDic[checkScan]['UCmQS'][spw_index]*np.mean(BPCaledXY.conjugate()))
        XYsign = np.sign(np.cos(XYphase))
        print('SPW[%d] : XY phase = %6.1f [deg] sign = %3.0f' % (spw, 180.0*XYphase/np.pi, XYsign))
        BPSPWList[spw_index][:,1] *= XYsign
    #-------- Apply Bandpass and Phase Correction
    for scan_index, scan in enumerate(BandScanList[BandName]):
        scanFlag  = scanDic[scan]['Flag']
        if len(scanFlag) == 0: continue
        scanPhase = scanDic[scan]['Gain'] / abs(scanDic[scan]['Gain'])
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            # print('---Applying Gain and Bandpass Correction for scan %d, spw %d' % (scan, spw))
            BP_ant = BPSPWList[spw_index].transpose(1,2,0)
            Xspec = XspecList[spw_index][scan_index][:,:,:,scanFlag]* (scanPhase[ant1]* scanPhase[ant0].conjugate())
            XspecList[spw_index][scan_index] = (Xspec.transpose(3, 0, 1, 2) / (BP_ant[polYindex][:,:,ant0]* BP_ant[polXindex][:,:,ant1].conjugate())).transpose(1,2,3,0)
        #
    #-------- Aperture Efficiencies Determination using Solar System Objects
    for scan_index, scan in enumerate(BandScanList[BandName]):
        if len(scanDic[scan]['Flag']) == 0: continue
        if scan in QSOscanList : continue              # filter QSO out
        uvw = np.mean(scanDic[scan]['UVW'], axis=2) #; uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        FscaleDic[scanDic[scan]['source']] = SSOAe(antList[antMap], BandbpSPW[BandName], uvw, scanDic[scan], SSODic, [XspecList[spw_index][scan_index][0::3] for spw_index in list(range(spwNum))])
    SSOList = list(FscaleDic.keys())
    for SSO in SSOList:
        if len(FscaleDic[SSO]) == 0: del FscaleDic[SSO]
    #-------- A priori Aperture Efficiencies (if no SSO) 
    WgSum = 0.0
    for SSO in FscaleDic.keys():
        SSOWG = 1.0
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']): SSOWG *= np.median(FscaleDic[SSO]['Wg'][spw_index])
        WgSum += SSOWG
    if WgSum > 1.0e-8: Aeff = averageAe(FscaleDic, BandbpSPW[BandName]['spw'])   # Aeff[ant, pol, spw]
    else: 
        Aeff = np.ones([useAntNum, 2, spwNum])
        for spw_index in list(range(spwNum)): Aeff[:,:,spw_index] = etaA[:,antMap].T
    text_sd = ' Aeff: '
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']): text_sd = text_sd + 'SPW%02d-X SPW%02d-Y ' % (spw, spw)
    fluxCalText = ''
    for source in FscaleDic.keys() : fluxCalText = fluxCalText + ' %s,' % (source)
    text_sd = text_sd + fluxCalText
    logfile.write(text_sd +'\n'); print(text_sd)
    logjy = open(prefix + '-' + BandName + '-JyK.log', 'w')
    timeLabel = qa.time('%fs' % np.median(scanDic[checkScan]['mjdSec']), form='ymd')[0]
    for ant_index, ant in enumerate(antList[antMap]):
        text_sd = '%s : ' % (ant)
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            for pol_index, polName in enumerate(['X','Y']):
                text_sd = text_sd + '  ----- ' if Aeff[ant_index][pol_index][spw_index] < 0.1 else text_sd + '  %4.1f%% ' % (100.0* Aeff[ant_index][pol_index][spw_index])
                text_jy = '%s %s %d %s %f %s %s %e %e %5.1f' % (prefix, ant, spw, polName, 2.0* kb / (0.25* Aeff[ant_index][pol_index][spw_index]* np.pi*antDia[ant_index]**2), timeLabel, BandName, np.median(BandbpSPW[BandName]['freq'][spw_index]), abs(BandbpSPW[BandName]['BW'][spw_index]), tempAtm); logjy.write(text_jy + '\n')
        logfile.write(text_sd +'\n'); print(text_sd)
    logjy.write('\n'); logjy.close()
    #-------- Gain transfer and equalization
    QSONonShadowScanList = [scan for scan in QSOscanList if np.median(scanDic[scan]['EL']) > ELshadow]
    if len(QSONonShadowScanList) > 0: checkScan = QSONonShadowScanList[np.argmax([scanDic[scan]['I'] for scan in QSONonShadowScanList])]
    scan_index = list(scanDic.keys()).index(checkScan)
    newAeff = np.ones([useAntNum, 2, spwNum])
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        newAeff[:,:,spw_index] = AeTransfer(np.mean(XspecList[spw_index][scan_index][0::3][:,chRange], axis=(1,3)), Aeff[:,:,spw_index], antDia[antMap] )
    for ant_index, ant in enumerate(antMap): Aeff[ant_index] = newAeff[ant_index]
    #-------- Stokes visibilities for each scan
    ScanFlux, ScanSlope, ErrFlux = np.zeros([len(BandScanList[BandName]), spwNum, 4]), np.zeros([len(BandScanList[BandName]), spwNum, 4]), np.zeros([len(BandScanList[BandName]), spwNum, 4])
    text_Stokes = np.repeat('',spwNum).tolist()
    #-------- Prepare plot files
    pp = PdfPages('FL_' + prefix + '_' + BandName + '.pdf')
    polLabel = ['I', 'Q', 'U', 'V']
    Pcolor   = ['black', 'blue', 'red', 'green']
    for scan_index, scan in enumerate(BandScanList[BandName]):
        if len(scanDic[scan]['Flag']) == 0: continue
        figFL, axes = plt.subplots(3, 4, figsize = (11, 8))
        figFL.suptitle(prefix + ' ' + BandName)
        figFL.text(0.45, 0.05, 'Projected baseline [m]')
        figFL.text(0.03, 0.85, 'Phase [deg]', rotation=90)
        figFL.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
        UVW = scanDic[scan]['UVW']; uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        #timeStamp, UVW = GetUVW(msfile, BandbpSPW[BandName]['spw'][1], scan)
        #uvw = np.mean(UVW, axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)[blMap]
        sourceName = scanDic[scan]['source']
        text_src  = ' %02d %010s EL=%4.1f deg' % (scan, sourceName, 180.0* np.median(scanDic[scan]['EL'])/np.pi)
        timeLabel = qa.time('%fs' % np.median(scanDic[scan]['mjdSec']), form='ymd')[0] + ' SA=%.1f' % (scanDic[scan]['SA']) + ' deg.'
        visChavList = []
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            text_Stokes[spw_index] = ' SPW%02d %5.1f GHz ' % (spw, 1.0e-9* np.median(BandbpSPW[BandName]['freq'][spw_index]))
            visChav = GainScale(newAeff[:,:,spw_index], antDia[antMap], np.mean(XspecList[spw_index][scan_index][:,chRange], axis=1))
            visChavList = visChavList + [visChav]
            StokesVis = Vis2Stokes(visChav, Dcat[antMap][:,:,spw_index], scanDic[scan]['PA'])
            #-------- SSO visibility to correct by model
            if sourceName in FscaleDic.keys(): StokesVis *= (SSODic[sourceName][1][spw_index] / FscaleDic[sourceName]['model'][spw_index])
            #-------- Linear regression to determine zero-spacing visibilities
            ScanFlux[scan_index, spw_index], ScanSlope[scan_index, spw_index], ErrFlux[scan_index, spw_index] = lmStokes(StokesVis, uvDist)
            for pol_index in list(range(4)):
                text_Stokes[spw_index] = text_Stokes[spw_index] + ' %7.4f (%.4f) ' % (ScanFlux[scan_index, spw_index, pol_index], ErrFlux[scan_index, spw_index, pol_index])
            text_Stokes[spw_index] = text_Stokes[spw_index] + '%6.3f   %6.1f ' % (100.0* np.sqrt(ScanFlux[scan_index, spw_index, 1]**2 + ScanFlux[scan_index, spw_index, 2]**2)/ScanFlux[scan_index, spw_index, 0], np.arctan2(ScanFlux[scan_index, spw_index, 2],ScanFlux[scan_index, spw_index, 1])*90.0/np.pi)
            #-------- Plot Stokes visibilities
            axes[0,spw_index].plot( uvDist, 180.0* np.angle(StokesVis[0])/ np.pi, '.', label=polLabel[0], color=Pcolor[0])
            axes[1,spw_index].plot( uvDist, StokesVis[0].real, '.', label=polLabel[0], color=Pcolor[0])
            axes[2,spw_index].plot( uvDist, StokesVis[1].real, '.', label=polLabel[1], color=Pcolor[1])
            axes[2,spw_index].plot( uvDist, StokesVis[2].real, '.', label=polLabel[2], color=Pcolor[2])
            axes[2,spw_index].plot( uvDist, StokesVis[3].real, '.', label=polLabel[3], color=Pcolor[3])
        #---- Update scanDic record entry
        CS, SN = np.cos(2.0* scanDic[scan]['PA']), np.sin(2.0* scanDic[scan]['PA'])
        scanDic[scan]['I'] = [ScanFlux[scan_index, spw_index, 0]* np.ones(len(CS)) for spw_index, spw in enumerate(BandbpSPW[BandName]['spw'])]
        scanDic[scan]['QCpUS'] = [ScanFlux[scan_index, spw_index, 1]* CS + ScanFlux[scan_index, spw_index, 2]* SN for spw_index, spw in enumerate(BandbpSPW[BandName]['spw'])]
        scanDic[scan]['UCmQS'] = [ScanFlux[scan_index, spw_index, 2]* CS - ScanFlux[scan_index, spw_index, 1]* SN for spw_index, spw in enumerate(BandbpSPW[BandName]['spw'])]
        scanDic[scan]['visChav'] = visChavList
        #---- Display results
        uvMin, uvMax, IMax = min(uvDist), max(uvDist), max(ScanFlux[scan_index,:,0])
        axes[0,0].text(0.0, 1.15*180, text_src)
        axes[0,3].text(0.0, 1.15*180, timeLabel)
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            axes[0,spw_index].plot( np.array([0.0, uvMax]), np.array([0.0, 0.0]), '-', color='grey')
            axes[1,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 0], ScanFlux[scan_index, spw_index, 0]+ uvMax* ScanSlope[scan_index, spw_index, 0]]), '-', color=Pcolor[0])
            axes[2,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 1], ScanFlux[scan_index, spw_index, 1]+ uvMax* ScanSlope[scan_index, spw_index, 1]]), '-', color=Pcolor[1])
            axes[2,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 2], ScanFlux[scan_index, spw_index, 2]+ uvMax* ScanSlope[scan_index, spw_index, 2]]), '-', color=Pcolor[2])
            axes[2,spw_index].plot( np.array([0.0, uvMax]), np.array([ScanFlux[scan_index, spw_index, 3], ScanFlux[scan_index, spw_index, 3]+ uvMax* ScanSlope[scan_index, spw_index, 3]]), '-', color=Pcolor[3])
            axes[0,spw_index].axis([0.0, uvMax, -180, 180])
            axes[1,spw_index].axis([0.0, uvMax, 0.0, 1.25*IMax])
            axes[2,spw_index].axis([0.0, uvMax, -0.25*IMax, 0.25*IMax])
            axes[2,spw_index].text(0.0, 1.02*180, 'SPW%2d %5.1f GHz' % (spw, 1.0e-9* np.median(BandbpSPW[BandName]['freq'][spw_index])))
        #
        axes[1,spw_index].legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        axes[2,spw_index].legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        plt.show()
        figFL.savefig(pp, format='pdf')
        #
        refFreq = 1.0e-9* np.mean(BandbpSPW[BandName]['freq'])
        relFreq = 1.0e-9* np.array([np.median( BandbpSPW[BandName]['freq'][spw_index]) for spw_index, spw in enumerate(BandbpSPW[BandName]['spw'])]) - refFreq
        text_mean = ' mean  %5.1f GHz ' % (refFreq)
        pflux, pfluxerr = np.zeros(4), np.zeros(4)
        for pol_index in list(range(4)):
            sol, solerr = linearRegression(relFreq, ScanFlux[scan_index, :, pol_index], ErrFlux[scan_index, :, pol_index] )
            pflux[pol_index], pfluxerr[pol_index] = sol[0], solerr[0]
            text_mean = text_mean + ' %7.4f (%.4f) ' % (pflux[pol_index], pfluxerr[pol_index])
        text_mean = text_mean + '%6.3f   %6.1f' % (100.0* np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/np.pi) 
        text_sd = text_src + ' ' + timeLabel
        logfile.write(text_sd +'\n'); print(text_sd)
        text_sd = ' SPW  Frequency     I                 Q                 U                 V             |  %Pol     EVPA '; logfile.write(text_sd + '\n'); print(text_sd)
        text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print(text_sd)
        for spw_index in list(range(spwNum)): logfile.write(text_Stokes[spw_index] + '\n'); print(text_Stokes[spw_index])
        text_sd = ' --------------------------------------------------------------------------------------------------------'; logfile.write(text_sd + '\n'); print(text_sd)
        if spwNum > 1: logfile.write(text_mean + '\n'); print(text_mean)
        text_sd = ' UV_min_max  %6.1f  %6.1f ' % (uvMin, uvMax); logfile.write(text_sd + '\n'); print(text_sd)
        text_ingest = '%10s, NE, NE, NE, NE, %e, %6.3f, %6.4f, %6.3f, %6.4f, %5.1f, %5.1f, NE, NE, %s, %s, %s\n' % (scanDic[scan]['source'], np.median(BandbpSPW[BandName]['freq'][spw_index]), pflux[0], np.sqrt(0.0004*pflux[0]**2 + pfluxerr[0]**2), np.sqrt(pflux[1]**2 + pflux[2]**2)/pflux[0], np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/pflux[0], np.arctan2(pflux[2],pflux[1])*90.0/np.pi, np.sqrt(pfluxerr[1]**2 + pfluxerr[2]**2)/np.sqrt(pflux[1]**2 + pflux[2]**2)*90.0/np.pi, timeLabel.replace('/','-'), fluxCalText, prefix); ingestFile.write(text_ingest)
    #-------- Review D-term
    Dterm = np.zeros([useAntNum, spwNum, 2], dtype=complex)
    DtermDic = {'mjdSec': np.median(timeStamp)}
    DtermDic['Freq'] = [np.median(BandbpSPW[BandName]['freq'][spw_index]) for spw_index, spw in enumerate(BandbpSPW[BandName]['spw'])]
    text_sd = 'D-term        '; logfile.write(text_sd); print(text_sd, end='')
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        IList, QCpUSList, UCmQSList, visChavList = [], [], [], []
        for scan_index, scan in enumerate(BandScanList[BandName]):
            if len(scanDic[scan]['Flag']) == 0: continue
            if min(scanDic[scan]['EL']) < ELshadow : continue
            IList       = IList + scanDic[scan]['I'][spw_index].tolist()
            QCpUSList   = QCpUSList + scanDic[scan]['QCpUS'][spw_index].tolist()
            UCmQSList   = UCmQSList + scanDic[scan]['UCmQS'][spw_index].tolist()
            visChavList = visChavList + scanDic[scan]['visChav'][spw_index].transpose(2, 0, 1).tolist()
        #
        Dterm[:,spw_index, 0], Dterm[:,spw_index,1] = VisMuiti_solveD(np.array(visChavList).transpose(1, 2, 0), np.array(QCpUSList), np.array(UCmQSList), Dcat[antMap,0,spw_index], Dcat[antMap,1,spw_index], np.array(IList))
    #---- Display D-terms
    for ant_index, ant in enumerate(antList[antMap]): DtermDic[ant] = Dterm[ant_index]
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        text_sd = '    SPW%02d                    ' % (spw);  logfile.write(text_sd); print(text_sd, end=' ')
    logfile.write('\n'); print('')
    text_sd = 'Ant  '; logfile.write(text_sd); print(text_sd, end='')
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        text_sd = '       Dx             Dy     ';  logfile.write(text_sd); print(text_sd, end='')
    logfile.write('\n'); print('')
    for ant_index, ant in enumerate(antList[antMap]):
        text_sd = '%s  ' % (ant)
        text_fd = text_sd
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            for pol_index in [0,1]:
                if abs(Dterm[ant_index,spw_index,pol_index]) > 0.1:
                    text_sd = text_sd + '\033[91m %+.3f%+.3fi \033[0m' % (Dterm[ant_index,spw_index,pol_index].real, Dterm[ant_index,spw_index,pol_index].imag)
                else:
                    text_sd = text_sd + ' %+.3f%+.3fi ' % (Dterm[ant_index,spw_index,pol_index].real, Dterm[ant_index,spw_index,pol_index].imag)
                #
                text_fd = text_fd + ' %+.3f%+.3fi ' % (Dterm[ant_index,spw_index,pol_index].real, Dterm[ant_index,spw_index,pol_index].imag)
            #
        #
        logfile.write(text_fd + '\n'); print(text_sd)
    #----
    logfile.close()
    ingestFile.close()
    plt.close('all')
    pp.close()
    del text_fd,text_sd,text_ingest,UCmQSList,QCpUSList,IList,DtermDic,Dterm,sol,solerr,pflux,pfluxerr,refFreq,relFreq,uvMin,uvMax,IMax,CS,SN,StokesVis,visChav,XspecList,scanDic,SSODic,FscaleDic,BandbpSPW,visChavList,ScanFlux,timeStamp,Xspec,BPCaledXspec,BPCaledXY,XPspec,BP_eq_gain,BPW,XYspec,Weight,pp,scanPhase,XYphase,XYsign,Aeff,newAeff,ScanSlope,ErrFlux,BPSPWList,scanGain,QSONonShadowScanList,BPcaledSpec,chAvgList,RXList,OnScanList,antList
#
