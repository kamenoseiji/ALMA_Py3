#import sys
#import pickle
import math
import numpy as np
import analysisUtils as au
import xml.etree.ElementTree as ET
from interferometry import BANDPA, indexList, GetAntName, GetSourceList, GetBandNames, GetAtmSPWs, GetBPcalSPWs, GetOnSource, GetVisAllBL
from Grid import *
from ASDM_XML import CheckCorr, BandList
#exec(open(SCR_DIR + 'interferometry.py').read())
#exec(open(SCR_DIR + 'Grid.py').read())
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
atmSPWs, bpSPWs = GetAtmSPWs(msfile), GetBPcalSPWs(msfile)
bandNameList = GetBandNames(msfile, atmSPWs)
OnScanList = GetOnSource(msfile)
msmd.open(msfile)
for BandName in RXList:
    BandPA[BandName] = (BANDPA[int(BandName[3:5])] + 90.0)*math.pi/180.0
    BandbpSPW[BandName]= np.array(bpSPWs)[np.where(np.array(bandNameList) == BandName)[0].tolist()].tolist()
    BandatmSPW[BandName]= np.array(atmSPWs)[np.where(np.array(bandNameList) == BandName)[0].tolist()].tolist()
    BandScanList[BandName] = list(set(msmd.scansforspw(BandbpSPW[BandName][0])) & set(OnScanList))
    BandScanList[BandName].sort()
#
msmd.close()
#-------- Tsys measurement
#exec(open(SCR_DIR + 'TsysCal.py').read())
#execfile(SCR_DIR + 'TsysCal.py')
#-------- Check Antenna List
antList = GetAntName(msfile)
antNum = len(antList)
blNum = int(antNum* (antNum - 1) / 2)
#-------- Check source list
print('---Checking source list')
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
#-------- Loop for Bands
for BandName in RXList:
    scanList = BandScanList[BandName]
    spwList  = BandbpSPW[BandName]
    atmspwList  = BandatmSPW[BandName]
    #-------- Load Visibilities into memory
    timeStampList, XspecList = [], []
    for spw in spwList:
        XscanList = []
        for scan in scanList:
            timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
            XscanList = XscanList + [Xspec]
            if spw == spwList[0]: timeStampList = timeStampList + [timeStamp]
        #
        XspecList = XspecList + [XscanList]
    #
#
a = 1.0/0.0
msmd.open(msfile)
ONScans = sort(np.array(list(set(msmd.scansforintent("*CALIBRATE_POLARIZATION*")) | set(msmd.scansforintent("*CALIBRATE_AMPLI*")) | set(msmd.scansforintent("*CALIBRATE_BANDPASS*")) | set(msmd.scansforintent("*CALIBRATE_FLUX*")) | set(msmd.scansforintent("*CALIBRATE_PHASE*")) | set(msmd.scansforintent("*OBSERVE_TARGET*")) | set(msmd.scansforintent("*CALIBRATE_DELAY*")) )))
msmd.close(); msmd.done()
#-------- Loop for Bands
for band_index in list(range(NumBands)):
    msmd.open(msfile)
    timeSum = 0
    bandName = UniqBands[band_index]; bandID = int(UniqBands[band_index][3:5])
    #if Tau0Max[band_index] < 0.01: print('Failed in Trx/Tsys measurements...'); os.system('touch ' + prefix + '-' + UniqBands[band_index] + '-Flux.log'); continue
    if Tau0Max[band_index] > TauLimit[bandID]: print('Too high optical depth...'); os.system('touch ' + prefix + '-' + UniqBands[band_index] + '-Flux.log'); continue
    ONScan = np.array(bpscanLists[band_index])[indexList( ONScans, np.array(bpscanLists[band_index]))]
    ATMScan= np.array(atmscanLists[band_index])[indexList( ONScans, np.array(atmscanLists[band_index]))]
    onsourceScans, atmScans = ONScan.tolist(), ATMScan.tolist()
    scanNum, atmscanNum = len(onsourceScans), len(atmScans)
    SSOScanIndex = []
    StokesDic = dict(zip(sourceList, [[]]*numSource))   # Stokes parameters for each source
    #-------- Check AZEL
    azelTime, AntID, AZ, EL = GetAzEl(msfile)
    azelTime_index = np.where( AntID == 0 )[0].tolist() 
    azel = np.r_[AZ[azelTime_index], EL[azelTime_index]].reshape(2, len(azelTime_index))
    OnAZ, OnEL, OnPA, BPquality, EQquality, PQquality, sourceIDscan, FLscore, refTime = [], [], [], [], [], [], [], np.zeros(scanNum), []
    #-------- Check QU catalog
    if QUMODEL: # A priori QU model from Rdata
        os.system('rm -rf CalQU.data')
        text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, qa.time('%fs' % (azelTime[0]), form='ymd')[0], BANDFQ[bandID])
        for source in sourceList: text_sd = text_sd + ' ' + source
        print(text_sd)
        os.system(text_sd)
        fp = open('CalQU.data')
        lines = fp.readlines()
        fp.close()
        for eachLine in lines:
            sourceName = eachLine.split()[0]
            StokesDic[sourceName] = [float(eachLine.split()[1]), float(eachLine.split()[2]), float(eachLine.split()[3]), 0.0]
        #
    #
    scanList = []
    for scan_index in range(scanNum):
        sourceID = msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0])
        #sourceIDscan.append( msmd.sourceidforfield(msmd.fieldsforscan(onsourceScans[scan_index])[0]))
        sourceIDscan = sourceIDscan + [sourceID]
        SunAngle = SunAngleSourceList[sourceID]
        interval, timeStamp = GetTimerecord(msfile, 0, 0, bpspwLists[band_index][0], onsourceScans[scan_index])
        timeSum += len(timeStamp)
        try:
            AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 0, AZ, EL)
        except:
            AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, 1, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]; dPA = np.std(np.sin(PA)) #dPA = abs(np.sin(max(PA) - min(PA)))
        OnAZ.append(np.median(AzScan)); OnEL.append(np.median(ElScan)); OnPA.append(np.median(PA))
        refTime = refTime + [np.median(timeStamp)]
        sourceName = sourceList[sourceID]
        if len(StokesDic[sourceName]) == 4:
            CS, SN = np.cos(2.0* OnPA[scan_index]), np.sin(2.0* OnPA[scan_index])
            QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
            UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
            if QUMODEL:
                BPquality = BPquality + [1000.0* StokesDic[sourceName][0]* abs(UCmQS)* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(StokesDic[sourceName][0] + 0.5)]
            else:
                #BPquality = BPquality + [1000.0* abs(UCmQS)* dPA* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(StokesDic[sourceName][0])]
                BPquality = BPquality + [1000.0* abs(UCmQS)* np.sin(OnEL[scan_index] - 0.5*ELshadow) / np.sqrt(StokesDic[sourceName][0])]
            #
            EQquality = EQquality + [(StokesDic[sourceName][0] - abs(QCpUS))**2 * np.sin(OnEL[scan_index] - ELshadow) ]
        else:
            QCpUS, UCmQS = 0.0, 0.0
            BPquality = BPquality + [-1]
            EQquality = EQquality + [-1]
        #
        if SunAngle < SunAngleThresh: BPquality[scan_index], EQquality[scan_index] = -1, -1
        print('Scan%02d : %10s AZ=%6.1f EL=%4.1f SA=%5.1f PA=%6.1f dPA=%5.2f pRes=%5.2f BPquality=%7.4f EQquality=%6.0f' % (onsourceScans[scan_index], sourceList[sourceIDscan[scan_index]], 180.0*OnAZ[scan_index]/np.pi, 180.0*OnEL[scan_index]/np.pi, SunAngle, 180.0*OnPA[scan_index]/np.pi, 180.0*dPA/np.pi, UCmQS, BPquality[scan_index], EQquality[scan_index]))
        if sourceIDscan[scan_index] in SSOList: FLscore[scan_index] = np.exp(np.log(math.sin(OnEL[scan_index])-0.34))* SSOscore[bandID-1][SSOCatalog.index(sourceList[sourceIDscan[scan_index]])]
    #
    #-------- Select Bandpass Calibrator
    if 'BPScans' not in locals():
        BPscanIndex = np.argmax(BPquality)
    else:
        if len(BPScans) < NumBands:
            BPscanIndex = np.argmax(BPquality)
        else:
            BPscanIndex = onsourceScans.index(BPScans[band_index])
        #
    #
    BPScan = onsourceScans[BPscanIndex]; BPcal = sourceList[sourceIDscan[BPscanIndex]]; timeLabelBP = qa.time('%fs' % (refTime[BPscanIndex]), form='ymd')[0]
    BPEL = OnEL[onsourceScans.index(BPScan)]
    #---- Bandpass table
    scanList = np.array(onsourceScans)[np.where(np.array(BPquality) > 10.0)[0].tolist()].tolist()
    if len(scanList) > 1:
        BPPLOT = True
        spwList = bpspwLists[band_index]
        XYwgt = []
        '''
        SNR_THRESH = 3
        for spw in spwList:
            exec(open(SCR_DIR + 'checkGain.py').read())
        '''
        for BPscan in scanList:
            exec(open(SCR_DIR + 'checkBP.py').read())
            refant = antList[UseAnt[refantID]]
            XYwgt = XYwgt + [BPquality[onsourceScans.index(BPscan)]]
        #
        BPscan = onsourceScans[np.argmax(BPquality)]
        exec(open(SCR_DIR + 'averageBP.py').read())
        BPscan = 0
    else:
        exec(open(SCR_DIR + 'checkBP.py').read())
    #
    #-------- Select Equalization Calibrator
    if 'EQScans' not in locals():
        EQscanIndex = np.argmax(EQquality)
    else:
        if len(EQScans) < NumBands:
            EQscanIndex = np.argmax(EQquality)
        else:
            EQscanIndex = onsourceScans.index(EQScans[band_index])
        #
    #
    EQScan = onsourceScans[EQscanIndex]; EQcal = sourceList[sourceIDscan[EQscanIndex]]; timeLabelEQ = qa.time('%fs' % (refTime[EQscanIndex]), form='ymd')[0]
    BPcalText = 'Use %s [Scan%d EL=%4.1f deg] %s as Bandpass Calibrator' % (BPcal, BPScan, 180.0* OnEL[onsourceScans.index(BPScan)]/np.pi, timeLabelBP); print(BPcalText)
    EQcalText = 'Use %s [Scan%d EL=%4.1f deg] %s as Gain Equalizer' % (EQcal, EQScan, 180.0* OnEL[onsourceScans.index(EQScan)]/np.pi, timeLabelEQ); print(EQcalText)
    #-------- Check Baseline Limit
    if 'uvLimit' in locals(): antFlag = angFlagBL(msfile, uvLimit, bpspwLists[band_index][0], BPScan, antFlag)
    #-------- SSO in observed source list
    BandSSOList = list( set(SSOList) & set(sourceIDscan) )
    if len(BandSSOList) == 0: Apriori = True
    #-------- Polarization setup
    atmspw = atmspwLists[band_index]; spwNum = len(atmspw)
    scnspw = bpspwLists[band_index]; scnspwNum = len(scnspw)
    #polNum = msmd.ncorrforpol(msmd.polidfordatadesc(scnspw[0]))
    polNum = 4
    PolList = ['X', 'Y']
    msmd.done()
    if polNum == 4:
        pPol, cPol = [0,3], [1,2];  ppolNum, cpolNum = len(pPol), len(cPol)
        exec(open(SCR_DIR + 'SSO_Stokes.py').read()) # Flux calibration using SSO
        #exec(open(SCR_DIR + 'aprioriStokes.py').read())
        '''
        if Apriori:
            try:
                exec(open(SCR_DIR + 'aprioriStokes.py').read())
            except:
                print('  --A priori flux calibration falied.')
        else:
            try:
                exec(open(SCR_DIR + 'SSO_Stokes.py').read())
            except:
                print('  --SSO-based flux calibration falied. Switch to a priori (SEFD) calibration.')
                try:
                    exec(open(SCR_DIR + 'aprioriStokes.py').read())
                except:
                    print('  --A priori flux calibration falied.')
            #
        #
        '''
    #
    if polNum == 2:
        cPol = [0,1], []; ppolNum, cpolNum = len(pPol), len(cPol)
        if not Apriori:
            try:
                exec(open(SCR_DIR + 'checkSEFD.py').read())
            except:
                print(' --SSO-based flux calibration falied. Switch to a priori (SEFD) calibration.')
                exec(open(SCR_DIR + 'aprioriFlux.py').read())
        else:
            exec(open(SCR_DIR + 'aprioriFlux.py').read())
    #
#---- end of band loop
del msfile, UniqBands, useAnt, atmspwLists, atmSPWs
if 'spwFlag' in locals(): del spwFlag
if 'BPScans' in locals(): del BPScans
if 'EQScans' in locals(): del EQScans
if 'antFlag' in locals(): del antFlag
if 'flagAnt' in locals(): del flagAnt
