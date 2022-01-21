# TsysCal.py 
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
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
#-------- Check MS file
msfile = wd + prefix + '.ms'
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
BPscan  = GetBPcalScans(msfile)[0]   # BP scan
msmd.open(msfile)
atmScanList = msmd.scansforintent('*ATMOSPHERE*') 
BPspws = msmd.spwsforscan(BPscan)
BBspwList = sort(list( set(BPspws) & set(msmd.almaspws(sqld=True)))).tolist();  BBnum = len(BBspwList)
BBdic = dict(zip(list(range(1,BBnum + 1)), [[]]*BBnum)) # dictionary {BB: [SQLDspw, [atmScans], [onsourceScans]]}
for BBID in list(range(1,BBnum+1)):
    spwList = list(set( msmd.spwsforbaseband(BBID, sqldmode='Only') ) & set(BPspws))
    spwScanList = msmd.scansforspw(spwList[0])
    BBdic[BBID] = [spwList[0], sort(list(set(spwScanList) & set(atmScanList))).tolist(), sort(list(set(spwScanList) - set(atmScanList))).tolist()]
#
msmd.close(); msmd.done()
#--------
def GetPSpecState(msfile, ant, spwID, stateID):
    data_desc_id = SPW2DATA_DESC_ID(msfile, spwID)
    Out='ANTENNA1 == %d && ANTENNA2 == %d && DATA_DESC_ID == %d && STATE_ID == %d' % (ant, ant, data_desc_id, stateID)
    tb.open(msfile)
    MSpointer = tb.query(Out)
    timeStamp = MSpointer.getcol('TIME')
    ACORR  = MSpointer.getcol('DATA')
    tb.close()
    return timeStamp, ACORR.real
#


if 'atmSPWs' not in locals():
    atmSPWs = GetAtmSPWs(msfile); atmSPWs.sort()
atmBandNames = GetBandNames(msfile, atmSPWs); UniqBands = unique(atmBandNames).tolist()
if UniqBands == []: UniqBands = BandList
NumBands = len(UniqBands)
'''
msmd.open(msfile)
atmspwLists, atmscanLists = [], []
for band_index in list(range(NumBands)):
    bandAtmSPWs = np.array(atmSPWs)[indexList(np.array([UniqBands[band_index]]), np.array(atmBandNames))].tolist()
    if atmBandNames == []: bandAtmSPWs = np.array(atmSPWs)
    atmspwLists = atmspwLists + [bandAtmSPWs]
    if 'spwFlag' in locals():
        flagIndex = indexList(np.array(spwFlag), np.array(atmspwLists[band_index]))
        for index in flagIndex: del atmspwLists[band_index][index]
    #
    atmscanList = list(set(msmd.scansforspw(atmspwLists[band_index][0]))& set(msmd.scansforintent("CALIBRATE_ATMOSPHERE*")))
    atmscanList.sort()
    atmscanLists = atmscanLists + [atmscanList]
    print(' %s: atmSPW=' % (UniqBands[band_index]), end=''); print(atmspwLists[band_index])
#
# atmSPWs[band] : SPWs used in atmCal scans
# bpSPWs[band]  : SPWs used in bandpass scan (i.e. SPWs for OBS_TARGET)
# atmscanLists[band] : atmCal scans
# bpscanLists[band]  : scans on source
#
print('---Checking time for ambient and hot load')
timeOFF, timeON, timeAMB, timeHOT = msmd.timesforintent("CALIBRATE_ATMOSPHERE#OFF_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#ON_SOURCE"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#AMBIENT"), msmd.timesforintent("CALIBRATE_ATMOSPHERE#HOT")
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
# timeOFF : mjd of CALIBRATE_ATMOSPHERE#OFF_SOURCE
# timeON  : mjd of CALIBRATE_ATMOSPHERE#ON_SOURCE (becore Cycle 3, ambient + hot loads
# timeAMB : mjd of CALIBRATE_ATMOSPHERE#AMBIENT (after Cycle 3)
# timeHOT : mjd of CALIBRATE_ATMOSPHERE#HOT (after Cycle 3)
#
#-------- Get Load Temperature
tempAmb, tempHot  = np.zeros([useAntNum]), np.zeros([useAntNum])
for ant_index in list(range(useAntNum)):
    tempAmb[ant_index], tempHot[ant_index] = GetLoadTemp(msfile, useAnt[ant_index], atmspwLists[0][0])
    if tempAmb[ant_index] < 240: tempAmb[ant_index] += 273.15       # Old MS describes the load temperature in Celsius
    if tempHot[ant_index] < 240: tempHot[ant_index] += 273.15       #
#
#-------- Trx, TantN, and Tau0
Tau0Max = np.zeros(NumBands)
sourceList, polList = GetSourceList(msfile)
SunAngleSourceList = GetSunAngle(msfile)
for band_index in list(range(NumBands)):
    tsysLog = open(prefix + '-' + UniqBands[band_index] + '-Tsys.log', 'w')
    #-------- Trx
    atmTimeRef, offSpec, ambSpec, hotSpec, scanList = scanAtmSpec(msfile, useAnt, atmscanLists[band_index], atmspwLists[band_index], timeOFF, timeON, timeAMB, timeHOT)
    atmscanLists[band_index] = scanList; scanList.sort()
    print('atmCal scans = ', end=''); print(atmscanLists[band_index])
    atmscanNum, spwNum = len(atmscanLists[band_index]), len(atmspwLists[band_index])
    TrxList, TskyList, scanFlag = TrxTskySpec(useAnt, tempAmb, tempHot, atmspwLists[band_index], atmscanLists[band_index], ambSpec, hotSpec, offSpec)
    #-------- Check sun angle 
    for scan_index in list(range(len(scanList))):
        sourceID = msmd.sourceidforfield(msmd.fieldsforscan(scanList[scan_index])[0])
        SunAngle = SunAngleSourceList[sourceID]
        print('Sun Angle : %.1f' % (SunAngle))
        if SunAngle < SunAngleTsysLimit: scanFlag[:,:,:,scan_index] *= 0.01
    #
    for spw_index in list(range(spwNum)):
        np.save('%s-%s-SPW%d.Trx.npy' % (prefix, UniqBands[band_index], atmspwLists[band_index][spw_index]), TrxList[spw_index])    # TxList[spw][pol,ch,ant,scan]
    #-------- Az and El position at atmCal and onsource scans
    AtmEL = np.ones([useAntNum, atmscanNum])
    for ant_index in list(range(useAntNum)):
        azelTime_index = np.where( AntID == useAnt[ant_index] )[0].tolist()
        if len(azelTime_index) == 0: azelTime_index = np.where( AntID == useAnt[ant_index + 1] )[0].tolist()
        for scan_index in list(range(atmscanNum)): AtmEL[ant_index, scan_index] = EL[azelTime_index[argmin(abs(azelTime[azelTime_index] - atmTimeRef[scan_index]))]]
    #
    atmsecZ  = 1.0 / np.sin( np.median(AtmEL, axis=0) )
    #-------- Tsky and TantN
    Tau0, Tau0Excess, Tau0Coef, TantN = tau0SpecFit(tempAtm, atmsecZ, useAnt, atmspwLists[band_index], TskyList, scanFlag)
    SPWfreqList, freqList = [], []
    Tau0med = np.zeros(spwNum)
    for spw_index in list(range(spwNum)):
        chNum, chWid, freq = GetChNum(msfile, atmspwLists[band_index][spw_index]); freqList = freqList + [freq*1.0e-9]; SPWfreqList = SPWfreqList + [np.median(freq)*1.0e-9]
        Tau0med[spw_index] = np.median(Tau0[spw_index])
    #
    Tau0Max[band_index] = np.max(Tau0med)
    #-------- Log Tau0 (mean opacity at zenith)
    LogTrx(antList[useAnt], atmspwLists[band_index], SPWfreqList, scanList, atmTimeRef, TrxList, TantN, tsysLog)
    text_sd = 'Tau0 :  '
    for spw_index in list(range(spwNum)): text_sd = text_sd + ' %7.5f | ' % (Tau0med[spw_index])
    tsysLog.write(text_sd + '\n'); print(text_sd)
    tsysLog.close()
    #-------- Save to npy files
    np.save('%s-%s.TrxAnt.npy' % (prefix, UniqBands[band_index]), antList[useAnt])     # antList[ant]
    np.save('%s-%s.atmTime.npy' % (prefix, UniqBands[band_index]), atmTimeRef)     # antList[ant]
    np.save('%s-%s.TauE.npy' % (prefix, UniqBands[band_index]), Tau0Excess)     # antList[ant]
    for spw_index in list(range(spwNum)):
        np.save('%s-%s-SPW%d.TrxFreq.npy' % (prefix, UniqBands[band_index], atmspwLists[band_index][spw_index]), freqList[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.Trx.npy' % (prefix, UniqBands[band_index], atmspwLists[band_index][spw_index]), TrxList[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.TantN.npy' % (prefix, UniqBands[band_index], atmspwLists[band_index][spw_index]), TantN[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.Tau0.npy' % (prefix, UniqBands[band_index], atmspwLists[band_index][spw_index]), Tau0[spw_index])    # freqList[spw]
        np.save('%s-%s-SPW%d.Tau0C.npy' % (prefix, UniqBands[band_index], atmspwLists[band_index][spw_index]), Tau0Coef[spw_index])    # freqList[spw]
    #
    #---- Plots
    if not 'PLOTFMT' in locals():   PLOTFMT = 'pdf'
    if PLOTTAU:
        plotTauSpec(prefix + '_' + UniqBands[band_index], atmspwLists[band_index], freqList, Tau0) 
        plotTauFit(prefix + '_' + UniqBands[band_index], antList[useAnt], atmspwLists[band_index], atmsecZ, tempAtm, Tau0, TantN, TskyList, np.min(scanFlag, axis=1)) 
        if len(atmscanLists[band_index]) > 5: plotTau0E(prefix + '_' + UniqBands[band_index], atmTimeRef, atmspwLists[band_index], Tau0, Tau0Excess, np.min(scanFlag, axis=(1,2))) 
    if PLOTTSYS: plotTsys(prefix + '_' + UniqBands[band_index], antList[useAnt], atmspwLists[band_index], freqList, atmTimeRef, TrxList, TskyList)
#
#-------- Plot optical depth
'''
msmd.close()
msmd.done()
