# TsGysCal.py 
# inputs
#    SCR_DIR    : (character) path to the script (e.g. '/home/skameno/ALMA_SV/' )
#    PLOTTAU    : (boolean)   plot optical depths (True or False)
#    PLOTTSYS   : (boolean)   plot Trx and Tsys spectra (True or False)
#    prefix     : (character) UID name (e.g. 'uid___A002_Xc02418_X3f67' )
#
# outputs
#  chAvgTsys[band* scan* ant* spw* pol] : List of channel-averaged system noise temperature
#  TsysSpec[band* scan* ant* spw* pol][ch] : List of Tsys spectrum
#  TsysFlag[ant, spw, pol, scan] : 0 = invalid, 1=valid
#  Tau0med[spw] : median-value of the zenith optical depth
#  onTau[spw, scan] : on-source optical depth
#
#  They include all of antennas (even if flagged) in MS order
#
import analysisUtils as au
import sys
import scipy
import numpy as np
from interferometry import indexList, AzElMatch, GetTemp, GetAntName, GetAtmSPWs, GetBPcalSPWs, GetBandNames, GetAzEl, GetLoadTemp, GetPSpec, GetPSpecScan, GetSourceDic, GetChNum, get_progressbar_str
from atmCal import scanAtmSpec, residTskyTransfer, residTskyTransfer0, residTskyTransfer2, tau0SpecFit, TrxTskySpec, LogTrx, concatScans, ATTatm
from Plotters import plotTauSpec, plotTauFit, plotTau0E, plotTsys, plotTauEOn
from ASDM_XML import BandList
'''
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
parser.add_option('-T', dest='PLOTTAU', metavar='PLOTTAU',
    help='Plot opacity spectra', action="store_true")
parser.add_option('-t', dest='PLOTTSYS', metavar='PLOTTSYS',
    help='Plot Tsys and Trx spectra', action="store_true")
parser.add_option('-o', dest='ONTAU', metavar='ONTAU',
    help='Online Tsys correction', action="store_true")
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
antFlag = [ant for ant in options.antFlag.split(',')]
PLOTTAU = options.PLOTTAU
PLOTTSYS= options.PLOTTSYS
ONTAU   = options.ONTAU
'''
prefix = 'uid___A002_X122dee6_X129a0'
antFlag = []
PLOTTAU = True
PLOTTSYS = True
ONTAU = True
SunAngleTsysLimit = 5.0 # [deg] 
if 'PLOTTAU'  not in locals(): PLOTTAU  = False
if 'PLOTTSYS' not in locals(): PLOTTSYS = False
if 'ONTAU' not in locals(): ONTAU = False   # on-source real-time optical depth correction using channel-averaged autocorr power
#-------- Check MS file
msfile = prefix + '.ms'
tempAtm = GetTemp(msfile)
if tempAtm != tempAtm: tempAtm = 270.0; print('Cannot get ambient-load temperature ... employ 270.0 K, instead.')
antList = GetAntName(msfile)
antNum = len(antList)
if 'flagAnt' not in locals(): flagAnt = np.ones(antNum)
if 'antFlag' in locals():
    index =  indexList(np.array(antFlag), antList)
    if len(index) > 0: flagAnt[index] = 0.0
useAnt = np.where(flagAnt == 1.0)[0].tolist(); useAntNum = len(useAnt)
#-------- Check SPWs
print('---Checking spectral windows and scans with atmCal for ' + prefix)
if 'atmSPWs' not in locals():
    bpSPWs  = GetBPcalSPWs(msfile)
    atmSPWs = GetAtmSPWs(msfile)
    atmSPWs = list(set(bpSPWs) & set(atmSPWs)) if len(set(bpSPWs) & set(atmSPWs)) > 3 else atmSPWs
    atmSPWs.sort()
atmBandNames = GetBandNames(msfile, atmSPWs); UniqBands = list(set(atmBandNames))
if UniqBands == []: UniqBands = BandList(prefix)
NumBands = len(UniqBands)
msmd.open(msfile)
spwNameList = msmd.namesforspws()
atmspwLists, atmscanLists, SQLDspwLists, CHAVspwLists, OnScanLists = [], [], [], [], []
for band_index, bandName in enumerate(UniqBands):
    spwBandNameList = [spwName for spwName in spwNameList if bandName in spwName]
    #-------- BB List
    SQLDspwList = msmd.almaspws(sqld=True)
    SQLDspwNameList = msmd.namesforspws(SQLDspwList)
    BBList = [int(re.split(r'BB_',spwName)[1].split('#')[0]) for spwName in SQLDspwNameList if bandName in spwName]
    if 'BBFlag' in locals(): [BB for BB in BBList if BB not in BBFlag]
    BBList.sort()
    #-------- atm scan list
    atmscans = np.sort(msmd.scansforintent("CALIBRATE_ATMOSPHERE*"))
    atmscanList = [scan for scan in atmscans if bandName in msmd.namesforspws(msmd.spwsforscan(scan))[0]]
    if 'scanFlag' in locals(): atmscanList = [scan for scan in atmscanList if scan not in scanFlag]
    atmscanLists = atmscanLists + [atmscanList]
    #-------- SQLD List
    atmspwList, SQLDspwList, CHAVspwList = [], [], []
    for BB_index, BB in enumerate(BBList):
        atmspwList = atmspwList + [spw_index for spw_index, spwName in enumerate(spwNameList) if spw_index in msmd.spwsforscan(atmscanLists[band_index][0]) and 'FULL_RES' in spwName and 'BB_%d' % BB in spwName]
        SQLDspwList = SQLDspwList + [spw_index for spw_index, spwName in enumerate(spwNameList) if spw_index in msmd.spwsforscan(atmscanLists[band_index][0]) and 'SQLD' in spwName and 'BB_%d' % BB in spwName]
        CHAVspwList = CHAVspwList + [spw_index for spw_index, spwName in enumerate(spwNameList) if spw_index in msmd.spwsforscan(atmscanLists[band_index][0]) and 'CH_AVG' in spwName and 'BB_%d' % BB in spwName]
    atmspwLists  = atmspwLists  + [atmspwList]
    SQLDspwLists = SQLDspwLists + [SQLDspwList]
    CHAVspwLists = CHAVspwLists + [CHAVspwList]
    OnScanLists  = OnScanLists  + [list( (set(msmd.scansforintent('*ON_SOURCE')) - set(msmd.scansforintent('*ATMOSPHERE*'))) & set(msmd.scansforspw(SQLDspwLists[band_index][0])))]
    print(' %s: atmSPW=' % (UniqBands[band_index]), end=''); print(atmspwLists[band_index])
    TsysDigitalCorrection = True
#
# atmSPWs[band] : SPWs used in atmCal scans
# bpSPWs[band]  : SPWs used in bandpass scan (i.e. SPWs for OBS_TARGET)
# atmscanLists[band] : atmCal scans
# bpscanLists[band]  : scans on source
#
print('---Checking time for ambient and hot load')
timeOFF, timeON, timeAMB, timeHOT, timeTEST, timeREF = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#TEST"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#REFERENCE")
if len(timeAMB) == 0:
    for band_index in list(range(NumBands)):
        atmSPW = atmspwLists[band_index]
        timeXY, Pspec = GetPSpec(msfile, 0, atmSPW[0])
        timeNum, chNum = Pspec.shape[2], Pspec.shape[1]; chRange = list(range(int(0.05*chNum), int(0.95*chNum)))
        chAvgPower = np.mean(Pspec[0][chRange], axis=0)
        onTimeIndex  = indexList(timeON, timeXY)
        onTime, onPower = timeXY[onTimeIndex], chAvgPower[onTimeIndex]
        hot_index, amb_index = np.where(onPower >  np.median(onPower))[0].tolist(), np.where(onPower <  np.median(onPower))[0].tolist()
        hotStart = [hot_index[0]] + np.array(hot_index)[np.where(np.diff(onTime[hot_index]) > 60.0)[0] + 1].tolist()
        ambStart = [amb_index[0]] + np.array(amb_index)[np.where(np.diff(onTime[amb_index]) > 60.0)[0] + 1].tolist()
        if len(hot_index) > len(hotStart):
            hotTimeIndex, ambTimeIndex = np.array(onTimeIndex)[list(set(hot_index) - set(hotStart))], np.array(onTimeIndex)[list(set(amb_index) - set(ambStart))]
        else:
            hotTimeIndex, ambTimeIndex = np.array(onTimeIndex)[hot_index], np.array(onTimeIndex)[amb_index]
        timeAMB = np.append(timeAMB, timeXY[ambTimeIndex])
        timeHOT = np.append(timeHOT, timeXY[hotTimeIndex])
#
azelTime, AntID, AZ, EL = GetAzEl(msfile)
#-------- Get Load Temperature
tempAmb, tempHot  = np.zeros([useAntNum]), np.zeros([useAntNum])
for ant_index in list(range(useAntNum)):
    temp = GetLoadTemp(msfile, useAnt[ant_index], atmspwLists[0][0])
    if len(temp) == 0: continue
    tempAmb[ant_index], tempHot[ant_index] = temp[0], temp[1]
    if tempAmb[ant_index] < 240: tempAmb[ant_index] += 273.15       # Old MS describes the load temperature in Celsius
    if tempHot[ant_index] < 240: tempHot[ant_index] += 273.15       #
# timeOFF : mjd of CALIBRATE_ATMOSPHERE#OFF_SOURCE
# timeON  : mjd of CALIBRATE_ATMOSPHERE#ON_SOURCE (becore Cycle 3, ambient + hot loads
# timeAMB : mjd of CALIBRATE_ATMOSPHERE#AMBIENT (after Cycle 3)
# timeHOT : mjd of CALIBRATE_ATMOSPHERE#HOT (after Cycle 3)
# timeTEST: mjd of CALIBRATE_ATMOSPHERE#TEST (ON_SOURCE since Cycle 10)
#
#
#-------- Trx, TantN, and Tau0
srcDic = GetSourceDic(msfile)
def powerBias(sqld, chav):
    sumX = np.sum(sqld)
    sumY = np.sum(chav)
    sumXX= sqld.dot(sqld)
    sumXY= chav.dot(sqld)
    bias = (sumXX* sumY - sumX* sumXY)/(len(chav)* sumXX - sumX*sumX)
    return 0.0 if np.isnan(bias) else bias
#
for band_index, bandName in enumerate(UniqBands):
    tsysLog = open(prefix + '-' + bandName + '-Tsys.log', 'w')
    #-------- Trx
    atmTimeRef, offSpec, ambSpec, hotSpec, scanList = scanAtmSpec(msfile, useAnt, atmscanLists[band_index], atmspwLists[band_index], timeOFF, timeON, timeAMB, timeHOT)
    scanList.sort()
    atmscanLists[band_index] = scanList
    print('atmCal scans = ', end=''); print(atmscanLists[band_index])
    atmscanNum, spwNum = len(atmscanLists[band_index]), len(atmspwLists[band_index])
    for ant_index, ant in enumerate(antList[useAnt]):
        for spw_index, sqldspw in enumerate(SQLDspwLists[band_index]):
            chavspw = CHAVspwLists[band_index][spw_index]
            timeScan, SQLD = GetPSpecScan(msfile, ant_index, sqldspw, scanList[0])
            powerSQLD = np.array([np.median(SQLD[:,0,indexList(timeOFF, timeScan)], axis=1), np.median(SQLD[:,0,indexList(timeAMB, timeScan)], axis=1), np.median(SQLD[:,0,indexList(timeHOT, timeScan)], axis=1)])
            timeScan, CHAV = GetPSpecScan(msfile, ant_index, chavspw, scanList[0])
            powerCHAV = np.array([np.median(CHAV[:,0,indexList(timeOFF, timeScan)], axis=1), np.median(CHAV[:,0,indexList(timeAMB, timeScan)], axis=1), np.median(CHAV[:,0,indexList(timeHOT, timeScan)], axis=1)])[:,(0,-1)]
            CorrBias = [powerBias(powerSQLD[:,0], powerCHAV[:,0]), powerBias(powerSQLD[:,1], powerCHAV[:,1])]
            SpecAddress = atmscanNum*(ant_index* spwNum + spw_index)
            for index in list(range(SpecAddress, SpecAddress+atmscanNum)):
                hotSpec[index][0] = hotSpec[index][0] - CorrBias[0] 
                hotSpec[index][1] = hotSpec[index][1] - CorrBias[1] 
                ambSpec[index][0] = ambSpec[index][0] - CorrBias[0] 
                ambSpec[index][1] = ambSpec[index][1] - CorrBias[1] 
                offSpec[index][0] = offSpec[index][0] - CorrBias[0] 
                offSpec[index][1] = offSpec[index][1] - CorrBias[1] 
            #
        #
    #
    TrxList, TskyList, scanFlag = TrxTskySpec(useAnt, tempAmb, tempHot, atmspwLists[band_index], atmscanLists[band_index], ambSpec, hotSpec, offSpec)
    #-------- Check sun angle 
    for scan_index, scan in enumerate(scanList):
        sourceID = msmd.sourceidforfield(msmd.fieldsforscan(scan)[0])
        sourceName = srcDic[sourceID]['Name']
        SunAngle = srcDic[sourceID]['SA']
        print('%s : Sun Angle = %.1f' % (sourceName, SunAngle))
        if SunAngle < SunAngleTsysLimit: scanFlag[:,:,:,scan_index] *= 0.01
    for spw_index, spw in enumerate(atmspwLists[band_index]):
        np.save('%s-%s-SPW%d.Trx.npy' % (prefix, bandName, spw), TrxList[spw_index])    # TxList[spw][pol,ch,ant,scan]
    #-------- Az and El position at atmCal and onsource scans
    AtmEL = np.ones([useAntNum, atmscanNum])
    for ant_index, ant in enumerate(useAnt):
        azelTime_index = np.where( AntID == ant )[0].tolist()
        if len(azelTime_index) == 0: azelTime_index = np.where( AntID == useAnt[ant_index + 1] )[0].tolist()
        for scan_index, scan in enumerate(atmscanLists[band_index]):
            AtmEL[ant_index, scan_index] = EL[azelTime_index[np.argmin(abs(azelTime[azelTime_index] - atmTimeRef[scan_index]))]]
    #
    atmsecZ  = 1.0 / np.sin( np.median(AtmEL, axis=0) )
    #-------- Tsky and TantN
    Tau0, Tau0Excess, Tau0Coef, TantN = tau0SpecFit(tempAtm, atmsecZ, useAnt, atmspwLists[band_index], TskyList, scanFlag)
    SPWfreqList, freqList, Tau0med = [], [], []
    for spw_index, spw in enumerate(atmspwLists[band_index]):
        chNum, chWid, freq = GetChNum(msfile, spw)
        freqList = freqList + [freq*1.0e-9]; SPWfreqList = SPWfreqList + [np.median(freq)*1.0e-9]; Tau0med = Tau0med + [np.median(Tau0[spw_index])]
    Tau0med = np.array(Tau0med)
    #-------- Log Tau0 (mean opacity at zenith)
    LogTrx(antList[useAnt], atmspwLists[band_index], SPWfreqList, scanList, atmTimeRef, TrxList, TantN, tsysLog)
    text_sd = 'Tau0 :  '
    for spw_index in list(range(spwNum)): text_sd = text_sd + ' %7.5f | ' % (Tau0med[spw_index])
    tsysLog.write(text_sd + '\n'); print(text_sd)
    tsysLog.close()
    #-------- Save to npy files
    np.save('%s-%s.TrxAnt.npy' % (prefix, bandName), antList[useAnt])     # antList[ant]
    np.save('%s-%s.atmTime.npy' % (prefix, bandName), atmTimeRef)     # antList[ant]
    np.save('%s-%s.TauE.npy' % (prefix, bandName), Tau0Excess)     # antList[ant]
    for spw_index, spw in enumerate(atmspwLists[band_index]):
        np.save('%s-%s-SPW%d.TrxFreq.npy' % (prefix, bandName, spw), freqList[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.Trx.npy' % (prefix, bandName, spw), TrxList[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.TantN.npy' % (prefix, bandName, spw), TantN[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.Tau0.npy' % (prefix, bandName, spw), Tau0[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.Tau0C.npy' % (prefix, bandName, spw), Tau0Coef[spw_index])    # freqList[spw]
    #-------- Violently variable Tau0
    for spw_index, spw in enumerate(CHAVspwLists[band_index]):
        if np.std(Tau0Excess[spw_index]) / Tau0med[spw_index] > 0.15 : ONTAU = True #    variability > 15%
        if len(Tau0Excess[spw_index]) < 2 and Tau0med[spw_index] > 0.05 : ONTAU = True #   single-shot tau > 0.05
        if ONTAU:
            print('SPW%d : sd(Tau0) = %.3f / median(Tau0) = %.3f' % (atmspwLists[band_index][spw_index], np.std(Tau0Excess[spw_index]), Tau0med[spw_index]))
            onSQLD, offSQLD, ambSQLD, hotSQLD, onTime, offTime, ambTime, hotTime = [], [], [], [], [], [], [], []
            for scan_index, scan in enumerate(OnScanLists[band_index]):
                scanOn = []
                for ant_index, ant in enumerate(antList):
                    progress = (scan_index* antNum + ant_index + 1.0) / (antNum* len(OnScanLists[band_index]))
                    sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
                    timeScan, SQLD = GetPSpecScan(msfile, ant_index, spw, scan)
                    scanOn  = scanOn + [SQLD[0,0] + SQLD[1,0]]
                onSQLD = onSQLD + [np.median(np.array(scanOn), axis=0)]
                onTime = onTime + [timeScan]
            TauSQLD = []
            for scan_index, scan in enumerate(atmscanLists[band_index]):
                scanOff, scanAmb, scanHot = [], [], []
                for ant_index, ant in enumerate(antList):
                    timeScan, SQLD = GetPSpecScan(msfile, ant_index, spw, scan)
                    offIndex = indexList(timeOFF, timeScan); scanOff = scanOff + [SQLD[0,0,offIndex] + SQLD[1,0,offIndex]]
                    ambIndex = indexList(timeAMB, timeScan); scanAmb = scanAmb + [SQLD[0,0,ambIndex] + SQLD[1,0,ambIndex]]
                    hotIndex = indexList(timeHOT, timeScan); scanHot = scanHot + [SQLD[0,0,hotIndex] + SQLD[1,0,hotIndex]]
                offSQLD = offSQLD + [np.median(np.array(scanOff), axis=0)]; offTime = offTime + [timeScan[offIndex]]
                ambSQLD = ambSQLD + [np.median(np.array(scanAmb), axis=0)]; ambTime = ambTime + [timeScan[ambIndex]]
                hotSQLD = hotSQLD + [np.median(np.array(scanHot), axis=0)]; hotTime = hotTime + [timeScan[hotIndex]]
            #
            onTimeCont,  onSQLDCont  = concatScans(onTime,  onSQLD)
            offTimeCont, offSQLDCont = concatScans(offTime, offSQLD)
            ambTimeCont, ambSQLDCont = concatScans(ambTime, ambSQLD)
            hotTimeCont, hotSQLDCont = concatScans(hotTime, hotSQLD)
            scaleFact = ATTatm(onTimeCont, onSQLDCont, offTimeCont, offSQLDCont)    # Attenuator between atmCal and onsource
            TskyOff= (offSQLDCont* (np.median(tempHot) - np.median(tempAmb)) + np.median(tempAmb)* np.median(hotSQLDCont) - np.median(tempHot)* np.median(ambSQLDCont)) / (np.median(hotSQLDCont) - np.median(ambSQLDCont))
            TauOff = -np.log( (TskyOff - tempAtm) / (au.Tcmb - tempAtm) )
            az, el = AzElMatch(offTimeCont, azelTime, AntID, ant_index, AZ, EL )
            Tau0Off = TauOff * np.sin(el)
            Tau0Scale = np.median(Tau0Off) / np.median(Tau0[spw_index])
            TskyOn = (onSQLDCont/scaleFact* (np.median(tempHot) - np.median(tempAmb)) + np.median(tempAmb)* np.median(hotSQLDCont) - np.median(tempHot)* np.median(ambSQLDCont)) / (np.median(hotSQLDCont) - np.median(ambSQLDCont))
            TauOn  = -np.log( (TskyOn - tempAtm) / (au.Tcmb - tempAtm) )/Tau0Scale 
            az, el = AzElMatch(onTimeCont, azelTime, AntID, ant_index, AZ, EL )
            TauEOn = TauOn* np.sin(el) - Tau0med[spw_index]
            np.save('%s-%s-SPW%d.TauEon.npy' % (prefix, bandName,atmspwLists[band_index][spw_index]),np.array([onTimeCont,TauEOn]))     # antList[ant]
            if PLOTTAU: plotTauEOn(prefix, bandName, spw, onTimeCont, onSQLDCont, TauEOn+Tau0med[spw_index])
    #---- Plots
    if not 'PLOTFMT' in locals():   PLOTFMT = 'pdf'
    if PLOTTAU:
        plotTauSpec(prefix + '_' + bandName, atmspwLists[band_index], freqList, Tau0) 
        plotTauFit(prefix + '_'  + bandName, antList[useAnt], atmspwLists[band_index], atmsecZ, tempAtm, Tau0, TantN, TskyList, np.min(scanFlag, axis=1)) 
        if len(atmscanLists[band_index]) > 5: plotTau0E(prefix + '_' + bandName, atmTimeRef, atmspwLists[band_index], Tau0, Tau0Excess, np.min(scanFlag, axis=(1,2))) 
    if PLOTTSYS: plotTsys(prefix + '_' + bandName, antList[useAnt], atmspwLists[band_index], freqList, atmTimeRef, TrxList, TskyList)
#
del ONTAU, PLOTTAU, PLOTTSYS
#-------- Plot optical depth
msmd.close()
msmd.done()
