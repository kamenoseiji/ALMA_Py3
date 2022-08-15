#import sys
#import pickle
import math
import numpy as np
import analysisUtils as au
import xml.etree.ElementTree as ET
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import BANDPA, BANDFQ, indexList, GetAntName, GetSourceList, GetBandNames, GetAtmSPWs, GetBPcalSPWs, GetSPWFreq, GetOnSource, GetAzEl, GetUVW, loadScanSPW, AzElMatch, gainComplexErr, bestRefant, ANT0, ANT1, Ant2Bl, Ant2BlD, CrossPolBL, gainComplexVec, CrossPolBL, CrossPolBP, SPWalign
import matplotlib.pyplot as plt
from Plotters import plotBP
from Grid import *
from ASDM_XML import CheckCorr, BandList
from PolCal import GetAMAPOLAStokes, GetSSOFlux, PolResponse
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
BandbpSPW = GetSPWFreq(msfile, BandbpSPW)   # BandbpSPW[BandName] : [[SPW List][freqArray][chNum][BW]]
BandatmSPW = GetSPWFreq(msfile, BandatmSPW)
#-------- Tsys measurement
#exec(open(SCR_DIR + 'TsysCal.py').read())
#execfile(SCR_DIR + 'TsysCal.py')
#-------- Check Antenna List
if 'SNR_THRESH' not in locals(): SNR_THRESH = 0.0
if 'antFlag' not in locals(): antFlag = []
antList = GetAntName(msfile)
antNum = len(antList)
blNum = int(antNum* (antNum - 1) / 2)
chRange = list(range(3, 60))
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
#-------- Check source list
print('---Checking source list')
sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList); numSource = len(sourceList)
SSOList   = indexList( np.array(SSOCatalog), np.array(sourceList))
#-------- Check AZEL
azelTime, AntID, AZ, EL = GetAzEl(msfile)
azelTime_index = np.where( AntID == 0 )[0].tolist() 
#-------- Loop for Bands
for BandName in RXList:
    print('-----%s----' % (BandName))
    print('-----Estimation from AMAPOLA and Butler-JPL-Horizons')
    #-------- Load Visibilities into memory
    timeStampList, XspecList = loadScanSPW(msfile, BandbpSPW[BandName][0], BandScanList[BandName])
    StokesDic = GetAMAPOLAStokes(R_DIR, SCR_DIR, sourceList, qa.time('%fs' % (timeStampList[0][0]), form='ymd')[0], BANDFQ[int(BandName[3:5])])
    StokesDic, SSODic = GetSSOFlux(StokesDic, qa.time('%fs' % (timeStampList[0][0]), form='ymd')[0], [np.median(BandbpSPW[BandName][1][spw_index]) for spw_index, spw in enumerate(BandbpSPW[BandName][0])])
    PAList, CSList, SNList, QCpUSList, UCmQSList = [], [], [], [], []
    #-------- Check AZEL
    print('-----AZ, EL, PA')
    AzScanList, ElScanList = [], []
    for scan_index, scan in enumerate(BandScanList[BandName]):
        AzScan, ElScan = AzElMatch(timeStampList[scan_index], azelTime, AntID, 0, AZ, EL)
        AzScanList, ElScanList = AzScanList + [AzScan], ElScanList + [ElScan]
    #-------- Polarization responses
    print('-----Estimated polarization responses')
    PAList, CSList, SNList, QCpUSList, UCmQSLis, scanDic = PolResponse(msfile, StokesDic, BandPA[BandName], BandScanList[BandName], AzScanList, ElScanList)
    #-------- Check usable antennas and refant
    print('-----Filter usable antennas and determine reference antenna')
    checkScan   = BandScanList[BandName][np.argmax(np.array([scanDic[scan][3] for scan in BandScanList[BandName]]))]
    checkSource = scanDic[checkScan][0]
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
    timeStamp, UVW = GetUVW(msfile, BandbpSPW[BandName][0][0], checkScan)
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
    print('Use %s as refant' % (antList[UseAnt[refantID]]))
    antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
    for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
    print('-----Bandpass to align SPWs and polarization')
    #-------- Bandpass using checkScan
    FreqList, BPList, spwGainList = [], [], []
    for spw_index, spw in enumerate(BandbpSPW[BandName][0]):
        Xspec = CrossPolBL(XspecList[spw_index][BandScanList[BandName].index(checkScan)][:,:,blMap], blInv)
        BP_ant, BPCaledXYSpec, XYdelay, Gain, XYsnr = CrossPolBP(Xspec)
        BPList = BPList + [BP_ant]
        #-------- Bandpass-corrected cross-power spectrum
        BPcaled = (Xspec.transpose(3, 2, 0, 1) / (BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate())).transpose(3,2,1,0)
        checkVis= np.mean(BPcaled[chRange], axis=0)[[0,3]]
        spwGainList = spwGainList + [np.array([gainComplexVec(checkVis[0]), gainComplexVec(checkVis[1])])]
    #-------- SPW phase offsets
    spwTwiddle = SPWalign(np.array(spwGainList))
    #-------- Phase-aligned bandpass table
    for spw_index, spw in enumerate(BandbpSPW[BandName][0]):
        BPList[spw_index] = (BPList[spw_index].transpose(2,0,1)* spwTwiddle[:,:,spw_index]).transpose(1,2,0)
    pp = PdfPages('BP-%s-%s.pdf' % (prefix,BandName))
    plotBP(pp, prefix, antList[antMap], BandbpSPW[BandName][0], checkScan, BPList)
    #-------- SPW-combined phase calibration
    GainList = []
    for scan_index, scan in enumerate(BandScanList[BandName]):
        chAvgList = []
        for spw_index, spw in enumerate(BandbpSPW[BandName][0]):
            BP_ant = BPList[spw_index][:,:,chRange]
            Xspec = CrossPolBL(XspecList[spw_index][scan_index][:,:,blMap], blInv)
            chAvgList = chAvgList + [np.mean(Xspec[:,chRange].transpose(3,2,0,1) / (BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()), axis=3).transpose(2,1,0)[[0,3]]]
        #
        GainList = GainList + [gainComplexVec(np.mean(np.array(chAvgList), axis=(0, 1)))]
    #
    #-------- Scan-by-scan bandpass
    BPList, XYList = [], []
    for scan_index, scan in enumerate(BandScanList[BandName]):
        BPSPWList, XYSPWList = [], []
        for spw_index, spw in enumerate(BandbpSPW[BandName][0]):
            Xspec = CrossPolBL(XspecList[spw_index][scan_index][:,:,blMap], blInv)  # Xspec[pol, ch, bl, time]
            XPspec = np.mean(Xspec* GainList[scan_index][ant1]* GainList[scan_index][ant0].conjugate(), axis=3)
            BP_ant = np.array([gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[3].T)])
            BP_ant = (BP_ant.transpose(2,0,1) / abs(np.mean(BP_ant[:,:,chRange], axis=2))).transpose(2,1,0)
            BPCaledXspec = XPspec.transpose(2, 0, 1)/(BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate())
            BPCaledXYSpec = np.mean(BPCaledXspec[:,1], axis=0) +  np.mean(BPCaledXspec[:,2], axis=0).conjugate()
            BPCaledXYSpec = BPCaledXYSpec / abs(BPCaledXYSpec)
            BPSPWList = BPSPWList + [BP_ant]
            XYSPWList = XYSPWList + [BPCaledXYSpec]
        #
        BPList = BPList + [BPSPWList]
        XYList = XYList + [XYSPWList]
        #pp = PdfPages('BP-%s-%s-%d.pdf' % (prefix, BandName, scan))
        #plotBP(pp, prefix, antList[antMap], BandbpSPW[BandName][0], scan, BPList)
    #
    # XY reference scan
    BPscanIndex = np.argmax(np.array([scanDic[scan][3] for scan in BandScanList[BandName]]))
    XYscanIndex = np.argmax(np.array([scanDic[scan][4] for scan in BandScanList[BandName]]))
    for spw_index, spw in enumerate(BandbpSPW[BandName][0]):
        BPScanList, XYScanList = [], []
        for scan_index, scan in enumerate(BandScanList[BandName]):
            BPScanList = BPScanList + [BPList[scan_index][spw_index]]
            XYScanList = XYScanList + [XYList[scan_index][spw_index]]
        BPant  = np.array(BPScanList)
        XYspec = np.array(XYScanList)
        BPweight = np.zeros([XYspec.shape[0], antNum, 2], dtype=complex)
        BPmean = BPant[BPscanIndex]
        XYmean = XYspec[BPscanIndex]
        for iter in list(range(10)):
            #-------- BP table 
            for ant_index in list(range(antNum)):
                for pol_index in list(range(2)):
                    BPpower = np.sum(BPant[:, ant_index, pol_index]* BPant[:, ant_index, pol_index].conjugate(), axis=1).real
                    BPcorr = BPant[:, ant_index, pol_index].dot(BPmean[ant_index, pol_index].conjugate()) / np.sqrt(BPpower* (BPmean[ant_index, pol_index].dot(BPmean[ant_index, pol_index].conjugate()).real))
                    BPvar  = -np.log(abs(BPcorr))
                    BPweight[:, ant_index, pol_index]  = BPcorr.conjugate() / (BPvar + np.percentile(BPvar, 100/len(BandScanList[BandName])))
                    BPweight[:, ant_index, pol_index] = BPweight[:, ant_index, pol_index] / np.sum(abs(BPweight[:, ant_index, pol_index]))
                #
            #
            #
            BPmean = np.sum(BPant.transpose(3,0,1,2)* BPweight, axis=1).transpose(1,2,0)
            BPscale = np.mean(abs(BPmean[:,:,chRange]), axis=2)
            BPmean = (BPmean.transpose(2,0,1) / BPscale).transpose(1,2,0)
            #-------- XY phase 
            XYcorr = XYspec.dot(XYmean.conjugate()) / len(XYmean)
            XYsign = np.sign(XYcorr.real)
            if 'XYwgt' in locals():
                XYweight = XYsign* np.array(XYwgt)
            else:
                XYvar  = -np.log(abs(XYcorr))
                XYweight =  XYsign / (XYvar + np.percentile(XYvar, 100/len(BandScanList[BandName])))
            #
            XYmean   = (XYspec.T).dot(XYweight); XYmean = XYmean / abs(XYmean)
        #
        text_BPwgt, text_XYwgt, text_scan = 'BP wgt:', 'XY wgt:', 'Scan  :'
        for scan_index, scan in enumerate(BandScanList[BandName]):
            text_scan   = text_scan   + '    %3d ' % (scan)
            text_BPwgt  = text_BPwgt + '%7.3f ' % (np.median(abs(BPweight), axis=(1,2))[scan_index])
            text_XYwgt  = text_XYwgt + '%7.1f ' % (XYweight[scan_index])
        #
        print(text_scan)
        print(text_BPwgt)
        print(text_XYwgt)


    '''
        




    cmap = plt.get_cmap("tab10")
    for scan_index, scan in enumerate(BandScanList[BandName]):
        for ant_index in UseAnt:
            plt.plot(timeStampList[scan_index], np.angle(GainList[scan_index][ant_index]), '.', color=cmap(ant_index))
            #plt.plot(timeStampList[scan_index], abs(GainList[scan_index][ant_index]), '.', color=cmap(ant_index))
    '''
    #



    '''

    #-------- Gain table for all scans
    for scan_index, scan in enumerate(BandScanList[BandName]):
        for spw_index, spw in enumerate(BandbpSPW[BandName]):
            BP_ant = BPList[spw_index][:,:,chRange]
            chAvgVis = np.mean(XspecList[spw_index][scan_index][[0,3]][:,chRange].transpose(3,2,0,1) / (BP_ant[ant0]* BP_ant[ant1].conjugate()), axis=3).transpose(2,1,0)



            XspecList[spw_index][scan_index][[0,3]].transpose(3,2,0,1) / (BP_ant[ant0]* BP_ant[ant1].conjugate())
            , axis=1)
            Gain, tempErr = np.apply_along_axis(gainComplexErr, 0, chAvgVis[0])

        Xspec  = XspecList[spw_index][list(StokesDic.keys()).index(checkSource)]
        chNum = Xspec.shape[1]
        BP_ant  = np.ones([antNum, 2, chNum], dtype=complex)
        Xspec  = CrossPolBL(Xspec[:,:,blMap], blInv)
        Gain = np.array([gainComplexVec(np.mean(Xspec[0,chRange], axis=0)), gainComplexVec(np.mean(Xspec[3,chRange], axis=0))])
        CaledXspec = (abs(Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1])* Xspec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())).transpose(1,0,2,3)
        XPspec = np.mean(CaledXspec, axis=3)
        BP_ant[:,0], BP_ant[:,1] = gainComplexVec(XPspec[0].T), gainComplexVec(XPspec[3].T)
        #---- Amplitude normalization
        for pol_index in [0,1]:
            ant_index = np.where( abs(np.mean(BP_ant[:,pol_index], axis=1)) > 0.1* np.median( abs(np.mean(BP_ant[:,pol_index], axis=1)) ))[0].tolist()
            BP_ant[ant_index, pol_index] = (BP_ant[ant_index, pol_index].T / np.mean( abs(BP_ant[ant_index, pol_index]), axis=1)).T
        #---- XY delay
        BPCaledXspec = XPspec.transpose(2, 0, 1)* abs(BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex])**2 /(BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate())
        BPCaledXYSpec = np.mean(BPCaledXspec[:,1], axis=0) +  np.mean(BPCaledXspec[:,2], axis=0).conjugate()
        XYdelay, XYsnr = delay_search( BPCaledXYSpec[chRange] )
        XYdelay = (float(chNum) / float(len(chRange)))* XYdelay
        BPCaledXYSpec = BPCaledXYSpec / abs(BPCaledXYSpec)
        '''



    #-------- Gain table
    '''


    #-------- Polarization 
    print('        Source     :    I     p%     EVPA  QCpUS  UCmQS')
    print('-------+-----------+-------+------+------+------+------')
    for scan_index, scan in enumerate(BandScanList[BandName]):
        AzScan, ElScan = AzElMatch(timeStampList[scan_index], azelTime, AntID, 0, AZ, EL)
        PA = AzEl2PA(AzScan, ElScan) + BandPA[BandName]
        CS, SN, QCpUS, UCmQS = np.cos(2.0* PA), np.sin(2.0* PA), np.zeros(len(PA)), np.zeros(len(PA))
        sourceName = sourceList[msmd.sourceidforfield(msmd.fieldsforscan(scan)[0])]
        if StokesDic[sourceName] != []:
            QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
            UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
            print('Scan%3d %s : %6.2f %6.2f %6.2f %6.2f %6.2f' % (scan, sourceName, StokesDic[sourceName][0], 100.0*np.sqrt(StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2)/StokesDic[sourceName][0], 90.0* np.arctan2(StokesDic[sourceName][2], StokesDic[sourceName][1])/np.pi, np.median(QCpUS), np.median(UCmQS)))
        #
        PAList, CSList, SNList, QCpUSList, UCmQSList = PAList + [PA], CSList + [CS], SNList + [SN], QCpUSList + [QCpUS], UCmQSList + [UCmQS]
    msmd.close(); msmd.done()
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
    '''
#
'''
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
'''
