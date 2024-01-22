#import sys
#import pickle
import math
import numpy as np
import analysisUtils as au
import xml.etree.ElementTree as ET
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import Tcmb, kb, BANDPA, BANDFQ, indexList, subArrayIndex, GetAntName, GetAntD, GetSourceList, GetBandNames, GetAeff, quadratic_interpol, GetAtmSPWs, GetBPcalSPWs, GetSPWFreq, GetOnSource, GetUVW, loadScanSPW, AzElMatch, gainComplexErr, bestRefant, ANT0, ANT1, Ant2Bl, Ant2BlD, Bl2Ant, gainComplexVec, CrossPolBL, CrossPolBP, SPWalign, BPaverage, delay_search
import matplotlib.pyplot as plt
from Plotters import plotSP, plotXYP
from Grid import *
from ASDM_XML import CheckCorr, BandList
from PolCal import GetAMAPOLAStokes, GetSSOFlux, PolResponse
msfile = wd + prefix + '.ms'
#
#-------- Check Correlator Type
BLCORR = True
if 'ACA' in CheckCorr(prefix): BLCORR = False
#-------- Check Receivers
RXList = BandList(prefix)
BandPA  = dict(zip(RXList, [[]]*len(RXList)))    # Band PA
BandbpSPW  = dict(zip(RXList, [[]]*len(RXList))) # Band SPW for visibilitiies
BandatmSPW = dict(zip(RXList, [[]]*len(RXList))) # Band SPW for atmCal
BandScanList = dict(zip(RXList, [[]]*len(RXList))) # Band scan list
#-------- Check SPWs of atmCal and bandpass
print('---Checking SPWs and Scan information')
atmSPWs, bpSPWs = GetAtmSPWs(msfile), GetBPcalSPWs(msfile)
if 'spwFlag' in locals():
    atmSPWs = list(set(atmSPWs)- set(spwFlag));atmSPWs.sort()
    bpSPWs  = list(set(bpSPWs) - set(spwFlag)); bpSPWs.sort()
bandNameList = GetBandNames(msfile, atmSPWs)
OnScanList = GetOnSource(msfile)
msmd.open(msfile)
for BandName in RXList:
    BandPA[BandName] = (BANDPA[int(BandName[3:5])] + 90.0)*math.pi/180.0
    BandbpSPW[BandName]  = {'spw': np.array(bpSPWs)[np.where(np.array(bandNameList) == BandName)[0].tolist()].tolist()}
    BandatmSPW[BandName] = {'spw': np.array(atmSPWs)[np.where(np.array(bandNameList) == BandName)[0].tolist()].tolist()}
    BandScanList[BandName] = list(set(msmd.scansforspw(BandbpSPW[BandName]['spw'][0])) & set(OnScanList))
    BandScanList[BandName].sort()
#
msmd.close()
BandbpSPW = GetSPWFreq(msfile, BandbpSPW)   # BandbpSPW[BandName] : [[SPW List][freqArray][chNum][BW]]
BandatmSPW = GetSPWFreq(msfile, BandatmSPW)
#-------- Tsys measurement
exec(open(SCR_DIR + 'TsysCal.py').read())
#-------- Check Antenna List
if 'SNR_THRESH' not in locals(): SNR_THRESH = 0.0
if 'antFlag' not in locals(): antFlag = []
antList = GetAntName(msfile)
antDia = GetAntD(antList)
antNum = len(antList)
antDic = dict(zip(antList, [[]]* len(antList)))
blNum = int(antNum* (antNum - 1) / 2)
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
#-------- Check source list
print('---Checking source list')
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
FscaleDic = dict(zip(np.array(sourceList)[SSOList].tolist(), [[]]* len(SSOList)))
#-------- Check AZEL
#azelTime, AntID, AZ, EL = GetAzEl(msfile)
#azelTime_index = np.where( AntID == 0 )[0].tolist() 
#-------- Loop for Bands
for BandName in RXList:
    print('-----%s----' % (BandName))
    #-------- Load Aeff file
    etaA = 0.01* GetAeff(TBL_DIR, antList, int(UniqBands[band_index][3:5]), np.mean(azelTime)).T
    #Ae = 0.0025* np.pi* etaA* antDia[antMap]**2
    print('-----Estimation from AMAPOLA and Butler-JPL-Horizons')
    #-------- Load Visibilities into memory
    timeStampList, XspecList = loadScanSPW(msfile, BandbpSPW[BandName]['spw'], BandScanList[BandName])
    StokesDic = GetAMAPOLAStokes(R_DIR, SCR_DIR, sourceList, qa.time('%fs' % (timeStampList[0][0]), form='ymd')[0], BANDFQ[int(BandName[3:5])])
    if len(SSOList) > 0:
        StokesDic, SSODic = GetSSOFlux(StokesDic, qa.time('%fs' % (timeStampList[0][0]), form='ymd')[0], [1.0e-9* np.median(BandbpSPW[BandName]['freq'][spw_index]) for spw_index, spw in enumerate(BandbpSPW[BandName]['spw'])])
    #-------- Polarization responses per scan
    scanDic = PolResponse(msfile, StokesDic, BandPA[BandName], BandScanList[BandName], timeStampList)
    QSOscanList = [scan for scan in scanDic.keys() if scanDic[scan]['source'][0] == 'J' and str.isdigit(scanDic[scan]['source'][1])]
    AzScanList, ElScanList = [], []
    #-------- Check AZEL
    for scan_index, scan in enumerate(BandScanList[BandName]):
        AzScan, ElScan = AzElMatch(timeStampList[scan_index], azelTime, AntID, 0, AZ, EL)
        AzScanList, ElScanList = AzScanList + [AzScan], ElScanList + [ElScan]
    #-------- Trx and Zenith optical depth
    TauE = np.load('%s-%s.TauE.npy' % (prefix, BandName))   #  TauE[spw,scan]: time-variable excexs of zenith optical depth
    atmTime = np.load('%s-%s.atmTime.npy' % (prefix, BandName))#  atmTime[scan] : mjdSed at TauE measurements
    Tau0List, TrxList, TaNList, TrxFreq = [], [], [], []
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        Tau0List = Tau0List + [np.load('%s-%s-SPW%d.Tau0.npy' % (prefix, BandName, spw))]   # Tau0List[spw] [ch]
        TrxList  = TrxList  + [np.load('%s-%s-SPW%d.Trx.npy'  % (prefix, BandName, spw))]   # TrxList[spw] [pol, ch, ant, scan]
        TaNList  = TaNList  + [np.load('%s-%s-SPW%d.TantN.npy'% (prefix, BandName, spw))]   # TaNList[spw] [ant, ch]
        TrxFreq  = TrxFreq  + [np.load('%s-%s-SPW%d.TrxFreq.npy'% (prefix, BandName, spw))] # TrxFreq[spw] [ch]
    TrxAntList = np.load('%s-%s.TrxAnt.npy' % (prefix, BandName)) 
    #-------- Put Tsys into scanDic and Tsys correction
    atmReltime = atmTime - atmTime[0]
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    for scan_index, scan in enumerate(scanDic.keys()):
        scanTau = []
        TsysScanDic = dict(zip(TrxAntList, [[]]* len(TrxAntList)))
        source = scanDic[scan]['source']
        for spw_index in range(spwNum):
            TrxAnt = (np.median(TrxList[spw_index], axis=3) + TaNList[spw_index].T).transpose(1, 2, 0)  # [ch, ant, pol]
            Tant = np.zeros([chNum, antNum, 2])
            if source in SSOCatalog:
                Tant = Tant + (SSODic[source][1][spw_index]* etaA* np.pi* antDia**2 / (2.0* kb)).T    # Tant[ch, ant, pol] Antenna temperature of SSO
            SP = tauSMTH(atmReltime, TauE[spw_index] )
            Tau0SP = Tau0[spw_index] + np.median(scipy.interpolate.splev(scanDic[scan]['mjdSec'] - atmTime[0], SP))
            secZ = np.mean(1.0 / np.sin(scanDic[scan]['EL']))
            zenithTau = Tau0SP + Tau0Coef[spw_index][0] + Tau0Coef[spw_index][1]*secZ
            scanTau = scanTau + [zenithTau * secZ]
            exp_Tau = np.exp(-zenithTau * secZ )
            atmCorrect = 1.0 / exp_Tau
            TsysScan = atmCorrect* TrxAnt.transpose(1,2,0) + (Tcmb + Tant.transpose(1, 2, 0))*exp_Tau + tempAtm* (1.0 - exp_Tau)
            #-------- Tsys correction
            Xspec = XspecList[spw_index][scan_index].transpose(3, 2, 0, 1)
            XspecList[spw_index][scan_index] = (Xspec * np.sqrt(TsysScan[ant0][:,polXindex] * TsysScan[ant1][:,polYindex])).transpose(2,3,1,0)
            for ant_index, ant in enumerate(TrxAntList):
                TsysScanDic[ant] = TsysScanDic[ant] + [TsysScan[ant_index]]
        #
        scanDic[scan]['Tau']  = scanTau
        scanDic[scan]['Tsys'] = TsysScanDic
    #
    #-------- Check usable antennas and refant
    print('-----Filter usable antennas and determine reference antenna')
    #chRange = list(range(int(0.1*chNum), int(0.95*chNum)))
    chRange = BandbpSPW[BandName]['chRange'][0]
    checkScan   = QSOscanList[np.argmax(np.array([scanDic[scan]['I'] for scan in QSOscanList]))]
    checkSource = scanDic[checkScan]['source']
    Xspec       = XspecList[0][BandScanList[BandName].index(checkScan)]
    checkVis    = np.mean(Xspec[[0,3]][:,chRange], axis=1)
    GainX, tempErrX = np.apply_along_axis(gainComplexErr, 0, checkVis[0])
    GainY, tempErrY = np.apply_along_axis(gainComplexErr, 0, checkVis[1])
    SNRX = np.median(abs(GainX), axis=1) / np.sqrt(np.mean(abs(tempErrX)**2, axis=1))*np.sqrt(tempErrX.shape[1]) 
    SNRY = np.median(abs(GainY), axis=1) / np.sqrt(np.mean(abs(tempErrY)**2, axis=1))*np.sqrt(tempErrY.shape[1]) 
    flagAnt[np.where(SNRX < SNR_THRESH)[0].tolist()] = 0.0
    flagAnt[np.where(SNRY < SNR_THRESH)[0].tolist()] = 0.0
    UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
    text_sd = '  Usable antennas (%d) : ' % (len(UseAnt))
    for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
    print(text_sd)
    blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
    ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
    for bl_index in list(range(UseBlNum)): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
    timeStamp, UVW = GetUVW(msfile, BandbpSPW[BandName]['spw'][0], checkScan)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
    refantID = 12
    print('Use %s as refant' % (antList[UseAnt[refantID]]))
    antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
    for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
    print('-----Bandpass to align SPWs and polarization')
    #-------- Bandpass using checkScan
    FreqList, BPList, spwGainList = [], [], []
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        Xspec = CrossPolBL(XspecList[spw_index][BandScanList[BandName].index(checkScan)][:,:,blMap], blInv)
        BP_ant, BPCaledXYSpec, XYdelay, Gain, XYsnr = CrossPolBP(Xspec)
        BPList = BPList + [BP_ant]
        #-------- Bandpass-corrected cross-power spectrum
        spwGainList = spwGainList + [Gain]
    #-------- SPW phase offsets
    spwTwiddle = SPWalign(np.array(spwGainList))
    #-------- Phase-aligned bandpass table
    print('-----SPW-aligned bandpass')
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        BPList[spw_index] = (BPList[spw_index].transpose(2,0,1) * spwTwiddle[:,:,spw_index]).transpose(1,2,0)
    pp = PdfPages('BP-%s-%s.pdf' % (prefix,BandName))
    plotSP(pp, prefix, antList[antMap], BandbpSPW[BandName]['spw'], BandbpSPW[BandName]['freq'], BPList)
    os.system('rm -rf B0')
    bandpass(msfile, caltable='B0', spw=','.join(map(str, BandbpSPW[BandName]['spw'])), scan=str(checkScan), refant=antList[antMap[0]], solnorm=True)
    #-------- SPW-combined phase calibration
    print('-----Antenna-based gain correction')
    text_sd = '        coherence loss % :'
    for antName in antList[antMap]: text_sd = text_sd + ' ' +  antName
    print(text_sd)
    for scan_index, scan in enumerate(BandScanList[BandName]):
        #if scan not in QSOscanList : continue              # filter by QSO
        text_sd = '----- Scan%3d %10s :' % (scan, scanDic[scan]['source'])
        chAvgList = []
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            chRange = BandbpSPW[BandName]['chRange'][spw_index]
            BP_ant = BPList[spw_index]
            #medBP = np.median(abs(BP_ant), axis=(0,1))
            BPcaledSpec = CrossPolBL(XspecList[spw_index][scan_index][:,:,blMap], blInv)[[0,3]].transpose(3,2,0,1) / (BP_ant[ant0]* BP_ant[ant1].conjugate())
            chAvgList = chAvgList + [np.mean( BPcaledSpec[:,:,:,chRange], axis=(2,3))]
        #
        scanGain = gainComplexVec(np.mean(np.array(chAvgList), axis=0).T)
        scanDic[scan]['Gain'] = scanGain
        coh = abs(np.mean(scanGain, axis=1)) / np.mean(abs(scanGain), axis=1)
        for ant_index in list(range(antNum)): text_sd = text_sd + ' %.2f' % (100.0* (1.0 - coh[ant_index]))
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
        scanGain = scanDic[scan]['Gain'] / abs(scanDic[scan]['Gain'])
        BPSPWList, XYSPWList, XYsnrList = [], [], []
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            chRange = BandbpSPW[BandName]['chRange'][spw_index]
            Xspec = CrossPolBL(XspecList[spw_index][scan_index][:,:,blMap], blInv)          # Xspec[pol, ch, bl, time]
            XPspec = np.mean(Xspec/ (scanGain[ant0]* scanGain[ant1].conjugate()), axis=3)   # antenna-based phase correction
            BP_ant = np.array([gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[3].T)])
            #---- Amplitude normalization
            for pol_index in [0,1]:
                BP_eq_gain = np.mean(abs(BP_ant[pol_index][:,chRange]), axis=1)
                BP_ant[pol_index] = (BP_ant[pol_index].T / BP_eq_gain).T
            #
            BPCaledXspec = XPspec.transpose(0, 2, 1)/(BP_ant[polYindex][:,ant0]* BP_ant[polXindex][:,ant1].conjugate())
            BPCaledXYSpec = np.mean(BPCaledXspec[1], axis=0) +  np.mean(BPCaledXspec[2], axis=0).conjugate()
            XYdelay, XYsnr = delay_search( BPCaledXYSpec[chRange] )
            XYdelay = (float(chNum) / float(len(chRange)))* XYdelay
            text_sd = text_sd + ' %6.3f (%5.1f)' % (0.5* XYdelay/(BandbpSPW[BandName]['BW'][spw_index]*1.0e-9), XYsnr)
            XYsnrList = XYsnrList + [XYsnr**2]
            BPSPWList = BPSPWList + [BP_ant.transpose(1,0,2)]
            XYSPWList = XYSPWList + [BPCaledXYSpec]
        #
        print(text_sd)
        BPavgScanList = BPavgScanList + [scan]
        BPList = BPList + [BPSPWList]
        XYList = XYList + [XYSPWList]
        XYWList=XYWList + [XYsnrList]
        pp = PdfPages('BP-%s-%s-%d.pdf' % (prefix, BandName, scan))
        plotSP(pp, prefix, antList[antMap], BandbpSPW[BandName]['spw'], BandbpSPW[BandName]['freq'], BPSPWList)
        pp = PdfPages('XYP_%s_REF%s_Scan%d.pdf' % (prefix, antList[antMap[0]], scan))
        plotXYP(pp, prefix, BandbpSPW[BandName]['spw'], XYSPWList)
    #
    #-------- Average bandpass
    XYW = np.array(XYWList)**2
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        chRange = BandbpSPW[BandName]['chRange'][spw_index]
        refXY = XYList[np.argmax( XYW[:,spw_index])][spw_index]
        XY = 0.0* refXY
        BP = 0.0* BPList[0][spw_index]
        for scan_index, scan in enumerate(BPavgScanList):
            #-------- average BP
            BPW = abs(np.mean(scanDic[scan]['Gain'], axis=1)) / np.median(np.std(np.angle(scanDic[scan]['Gain']), axis=1))
            BP  = BP + (BPList[scan_index][spw_index].transpose(1,2,0)* BPW**2).transpose(2,0,1)
            #-------- average XY
            XYspec = XYList[scan_index][spw_index]
            Weight = XYW[scan_index][spw_index] * np.sign( XYspec.dot(refXY.conjugate()))
            XY     = XY + Weight* XYspec
        #---- Save averaged BP
        BPSPWList[spw_index] = (BP.transpose(2,0,1) / np.mean(abs(BP[:,:,chRange]), axis=2)).transpose(1,2,0)
        BPSPWList[spw_index][:,1] *= (XY / abs(XY))
        #---- Save into CASA caltable
        tb.open('B0', nomodify=False)
        SPWQ = tb.query('SPECTRAL_WINDOW_ID == %d'%(spw))
        SPWQ.putcol('CPARAM', BPSPWList[spw_index][indexList(np.array(range(antNum)), np.array(antMap))].transpose(1,2,0))
        SPWQ.close()
    #
    tb.close()
    pp = PdfPages('BP-%s-%s-%d.pdf' % (prefix, BandName, 0))
    plotSP(pp, prefix, antList[antMap], BandbpSPW[BandName]['spw'], BandbpSPW[BandName]['freq'], BPSPWList, 0.0, 1.2, True)
    del BPavgScanList, BPList, XYList, XYWList, XYW, XY, refXY, BP, XYSPWList, XYsnrList
    #---- 
    # Now we have
    #  Visibilities : XspecList [spw][scan][pol, ch, bl, time]
    #  Bandpass     : BPSPWList [spw][ant, pol, ch]
    #  scanDic      : 
    #-------- Apply Bandpass and Gain Correction
    for scan_index, scan in enumerate(BandScanList[BandName]):
        for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
            print('---Applying Gain and Bandpass Correction for scan %d, spw %d' % (scan, spw))
            BP_ant = BPSPWList[spw_index].transpose(1,2,0)
            scanPhase = scanDic[scan]['Gain']/abs(scanDic[scan]['Gain'])
            Xspec = CrossPolBL(XspecList[spw_index][scan_index][:,:,blMap], blInv)* (scanPhase[ant1]* scanPhase[ant0].conjugate())
            XspecList[spw_index][scan_index] = (Xspec.transpose(3, 0, 1, 2) / (BP_ant[polYindex][:,:,ant0]* BP_ant[polXindex][:,:,ant1].conjugate())).transpose(1,2,3,0)
        #
    #
    #-------- Aperture Efficiencies Determination using Solar System Objects
    AeList, WgList = [], []
    for scan_index, scan in enumerate(BandScanList[BandName]):
        if scan in QSOscanList : continue              # filter QSO out
        timeStamp, UVW = GetUVW(msfile, BandbpSPW[BandName]['spw'][1], scan)
        uvw = np.mean(UVW, axis=2) #; uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
        FscaleDic[scanDic[scan]['source']] = SSOAe(antList, antMap, BandbpSPW[BandName], uvw[:,blMap], scanDic[scan], SSODic, [XspecList[spw_index][scan_index][0::3] for spw_index in list(range(spwNum))])
    #-------- A priori Aperture Efficiencies (if no SSO) 
    if len(FscaleDic.keys()) == 0:
        Aeff = np.ones([antNum, 2, spwNum])
        for spw_index in list(range(spwNum)): Aeff[:,:,spw_index] = etaA.T
    else:
        Aeff = averageAe(FscaleDic, antList, BandbpSPW[BandName]['spw'])   # Aeff[ant, pol, spw]
    #
    #-------- Gain transfer and equalization
    QSONonShadowScanList = [scan for scan in QSOscanList if np.median(scanDic[scan]['EL']) > ELshadow]
    checkScan = QSONonShadowScanList[np.argmax([scanDic[scan]['I'] for scan in QSONonShadowScanList])]
    scan_index = list(scanDic.keys()).index(checkScan)
    newAeff = np.ones([antNum, 2, spwNum])
    for spw_index, spw in enumerate(BandbpSPW[BandName]['spw']):
        newAeff[antMap][:,:,spw_index] = AeTransfer(np.mean(XspecList[spw_index][scan_index][0::3][:,chRange], axis=(1,3)), Aeff[antMap][:,:,spw_index], antDia[antMap] )
    #
#
