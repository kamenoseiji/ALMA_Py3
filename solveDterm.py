exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Grid.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
from matplotlib.backends.backend_pdf import PdfPages
import pickle
#----------------------------------------- Procedures
def flagOutLier(value, thresh=5.0):
    return np.where(abs(value - np.median(value)) > thresh* np.std(value))[0].tolist()
#
if 'antFlag' not in locals():   antFlag = []
spwNum = len(spwList)
polXindex, polYindex = (arange(4)//2).tolist(), (arange(4)%2).tolist()
#
trkAntSet = set(range(64))
scansFile = []
pattern = r'RB_..'
timeNum = 0
sourceList = []
msfile = wd + prefix + '.ms'
Antenna1, Antenna2 = GetBaselineIndex(msfile, spwList[0], scanList[0])
UseAntList = CrossCorrAntList(Antenna1, Antenna2)
antList = GetAntName(msfile)[UseAntList]
srcDic = GetSourceDic(msfile)
for SSOkey in [key for key in srcDic.keys() if srcDic[key]['RA'] == 0.0]: del srcDic[SSOkey]    # Remove SSO
sourceList = unique([srcDic[ID]['Name'] for ID in srcDic.keys()]).tolist(); numSource = len(sourceList)
sourceScan = []
scanDic   = dict(zip(sourceList, [[]]*numSource)) # Scan list index for each source
timeDic   = dict(zip(sourceList, [[]]*numSource)) # Time index list for each source
StokesDic = dict(zip(sourceList, [[]]*numSource)) # Stokes parameters for each source
scanIndex = 0
print('-- Checking %s ' % (msfile))
#-------- Flagging status of antennas
flagTimeIndex = []
if 'FGprefix' in locals():  # Flag table
    antSPWFlag  = []
    for spw_index, spw in enumerate(spwList):
        TS = np.load('%s-SPW%d.TS.npy' % (FGprefix, spw));
        FG = np.load('%s-SPW%d.FG.npy' % (FGprefix, spw));
        newFlagIndex = np.where( np.median(FG, axis=1) == 0)[0].tolist()
        useAntIndex  = list(set(range(FG.shape[0])) - set(newFlagIndex))
        antSPWFlag = antSPWFlag + [antFlag + antList[newFlagIndex].tolist()]
        flagTimeIndex = flagTimeIndex + np.where(np.min(FG[useAntIndex], axis=0) == 0.0)[0].tolist()
    #
    newAntFlag = antList[newFlagIndex].tolist()
    antFlag = list(set([ant for row in antSPWFlag for ant in row])); antFlag.sort()
    useTimeIndex = list(set(range(len(TS))) - set(flagTimeIndex)); useTimeIndex.sort()
    useTimeStamp = TS[useTimeIndex]
    if len(antFlag) > 0: print('Flagged antennas : %s' % (str(antFlag)))
#
#sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList)
refAntID  = indexList([refantName], antList)
flagAntID = indexList(antFlag, antList)
if len(refAntID) < 1:
    print('Antenna %s didn not participate in this file.' % (refantName))
    sys.exit()
else:
    refAntID = refAntID[0]
#
msmd.open(msfile)
if 'scanList' in locals():
    scanList.sort()
    scanLS = scanList
else:
    scanLS = msmd.scannumbers().tolist()
    scanLS.sort()
#
spwName = msmd.namesforspws(spwList)[0]; BandName = re.findall(pattern, spwName)[0]; bandID = int(BandName[3:5])
BandPA = (BANDPA[bandID] + 90.0)*pi/180.0
#-------- Check source list and Stokes Parameters
for sourceID in srcDic.keys():
    sourceName = srcDic[sourceID]['Name']
    sourceIDscan = msmd.scansforfield(sourceID).tolist()
    scanDic[sourceName] = scanDic[sourceName] + sourceIDscan 
    if len(sourceIDscan) > 0:
        interval, timeStamp = GetTimerecord(msfile, 0, 1, spwList[0], sourceIDscan[0])
        IQU = GetPolQuery(sourceName, timeStamp[0], BANDFQ[bandID], SCR_DIR, R_DIR)
        StokesDic[sourceName] = [IQU[0][sourceName], IQU[1][sourceName], IQU[2][sourceName], 0.0]
        print('---- %s : expected I=%.1f p=%.1f%%' % (sourceName, StokesDic[sourceName][0], 100.0*sqrt(StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2)/StokesDic[sourceName][0]))
    #
#
msmd.done()
antMap = [refAntID] + list(trkAntSet - set([refAntID]))
antMap = indexList(np.array(antMap), np.array(UseAntList))
antNum = len(antMap); blNum = int(antNum * (antNum - 1)/2)
ant0 = ANT0[0:blNum]; ant1 = ANT1[0:blNum]
blMap, blInv= list(range(blNum)), [False]* blNum
for bl_index in list(range(blNum)):
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
#
if not 'bunchNum' in locals(): bunchNum = 1
def bunchVecCH(spec): return bunchVec(spec, bunchNum)
#-------- AZ, EL, PA
azelTime, AntID, AZ, EL = GetAzEl(msfile)
azelTime_index = np.where( AntID == refAntID )[0].tolist()
if len(azelTime_index) == 0: azelTime_index = np.where(AntID == 0)[0].tolist()
timeThresh = np.median( np.diff( azelTime[azelTime_index]))
#-------- Loop for SPW
DxList, DyList, FreqList = [], [], []
for spw_index, spw in enumerate(spwList):
    mjdSec, Az, El, PA, XspecList, timeNum, scanST = [], [], [], [], [], [], []
    #-------- time-independent spectral setups
    chNum, chWid, Freq = GetChNum(msfile, spw); chRange = list(range(int(0.05*chNum/bunchNum), int(0.95*chNum/bunchNum))); FreqList = FreqList + [1.0e-9* bunchVecCH(Freq) ]
    DxSpec, DySpec = np.zeros([antNum, int(ceil(chNum/bunchNum))], dtype=complex), np.zeros([antNum, int(ceil(chNum/bunchNum))], dtype=complex)
    caledVis = np.ones([4,blNum, 0], dtype=complex)
    if 'BPprefix' in locals():  # Bandpass file
        BPantList, BP_ant = np.load(BPprefix + '-REF' + refantName + '.Ant.npy'), np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (BPprefix, refantName, BPscan, spw))
        BP_ant = BP_ant[indexList(antList[antMap], BPantList)]      # BP antenna mapping
    if 'XYprefix' in locals():
        XYspec = np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (XYprefix, refantName, BPscan, spw))
        print('Apply XY phase into Y-pol Bandpass.'); BP_ant[:,1] *= XYspec  # XY phase cal
    #
    BP_bl = BP_ant[ant0][:,polYindex]* BP_ant[ant1][:,polXindex].conjugate()    # Baseline-based bandpass table
    #-------- For visibilities in each scan
    scanST = scanST + [0]
    for scan_index, scan in enumerate(scanList):
        print('-- Loading visibility data %s SPW=%d SCAN=%d...' % (prefix, spw, scan))
        timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)  # Xspec[POL, CH, BL, TIME]
        del Pspec
        if bunchNum > 1: Xspec = np.apply_along_axis(bunchVecCH, 1, Xspec)
        #---- remove flagged records
        flagIndex = indexList(useTimeStamp, timeStamp)
        timeNum = timeNum + [len(flagIndex)]
        if scan_index != 0: scanST = scanST + [scanST[scan_index - 1] + timeNum[scan_index - 1]]
        XspecList = XspecList + [Xspec[:,:,:,flagIndex]]
        del Xspec
        #-------- Expected polarization responses
        scanAz, scanEl = AzElMatch(timeStamp[flagIndex], azelTime, AntID, refAntID, AZ, EL)
        mjdSec, scanPA = mjdSec + timeStamp[flagIndex].tolist(), AzEl2PA(scanAz, scanEl, ALMA_lat) + BandPA
        Az, El, PA = Az + scanAz.tolist(), El + scanEl.tolist(), PA + scanPA.tolist()
    #
    #-------- Combine scans
    VisSpec, chAvgVis = np.zeros([4, int(ceil(chNum/bunchNum)), blNum, np.sum(timeNum)], dtype=complex), np.zeros([4, blNum, np.sum(timeNum)], dtype=complex)
    timeIndex = 0
    if chNum == 1:
        print('  -- Channel-averaged data: no BP and delay cal')
        for scan_index, scan in enumerate(scanList):
            chAvgVis[:, :, timeIndex:timeIndex + timeNum[scan_index]] =  CrossPolBL(XspecList[scan_index][:,:,blMap], blInv)[:,0]
            timeIndex = timeIndex + timeNum[scan_index]
    else:
        print('  -- Apply bandpass cal')
        for scan_index, scan in enumerate(scanList):
            print('   Scan %d : %d records' % (scan, timeNum[scan_index]))
            VisSpec[:,:,:,timeIndex:timeIndex + timeNum[scan_index]] = (CrossPolBL(XspecList[scan_index][:,:,blMap], blInv).transpose(3, 2, 0, 1) / BP_bl).transpose(2,3,1,0)
            timeIndex = timeIndex + timeNum[scan_index]
        #
        chAvgVis = np.mean(VisSpec[:,chRange], axis=1)
    #
    del XspecList
    #-------- Gain solutions
    mjdSec, Az, El, PA, PAnum = np.array(mjdSec), np.array(Az), np.array(El), np.array(PA), len(PA)
    print('---- Antenna-based gain solution using tracking antennas')
    Gain = np.array([ gainComplexVec(chAvgVis[0]), gainComplexVec(chAvgVis[3]) ])   # Parallel-pol gain
    Gamp = np.sqrt(np.mean(abs(Gain)**2, axis=0))
    Gain = Gamp* Gain/abs(Gain)
    #-------- Gain-calibrated visibilities
    print('  -- Apply parallel-hand gain calibration')
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    blWeight = np.mean(caledVis[[0,3]].real, axis=(0,2))**2 / (np.var(caledVis[0].real, axis=1) + np.var(caledVis[3].real, axis=1))
    blWeight = blWeight / np.sum(blWeight)
    Vis = caledVis.transpose(0,2,1).dot(blWeight)
    #-------- Coarse estimation of Q and U using XX and YY
    print('  -- Solution for Q and U')
    QCpUS, UCmQS =  np.zeros(PAnum), np.zeros(PAnum)
    if 'QUmodel' not in locals(): QUmodel = False
    CS, SN = np.cos(2.0* PA), np.sin(2.0* PA)
    for sourceName in scanDic.keys():
        scanLS = scanDic[sourceName]
        if len(scanLS) < 1 : continue
        timeIndex = []
        for scan in scanLS:
            scanIndex = scanList.index(scan)
            timeIndex = timeIndex + list(range(scanST[scanIndex], scanST[scanIndex] + timeNum[scanIndex]))
        timeDic[sourceName] = timeIndex
        if QUmodel:
            QUsol = np.array(StokesDic[sourceName])[[1,2]]/StokesDic[sourceName][0]
        else:
            timeIndex = list(set(timeIndex) - set(flagOutLier(Vis[0].real, 5.0) + flagOutLier(Vis[0].imag, 5.0) + flagOutLier(Vis[3].real, 5.0) + flagOutLier(Vis[3].imag, 5.0)))
            timeDic[sourceName] = timeIndex
            QUsol   = XXYY2QU(PA[timeIndex], Vis[[0,3]][:,timeIndex])             # XX*, YY* to estimate Q, U
            text_sd = '[XX,YY] %s: Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, QUsol[0], QUsol[1], 100.0* np.sqrt(QUsol[0]**2 + QUsol[1]**2), np.arctan2(QUsol[1],QUsol[0])*90.0/pi); print(text_sd)
        #
        QCpUS[timeIndex] = QUsol[0]* CS[timeIndex] + QUsol[1]* SN[timeIndex]
        UCmQS[timeIndex] = QUsol[1]* CS[timeIndex] - QUsol[0]* SN[timeIndex]
    #
    ##-------- XY phase determination
    print('  -- Degenerating pi-ambiguity in XY phase')
    XYphase = XY2Phase(UCmQS, Vis[[1,2]])    # XY*, YX* to estimate X-Y phase
    XYsign = np.sign(np.cos(XYphase))
    text_sd = '  XY Phase = %6.2f [deg]  sign = %3.0f' % (XYphase* 180.0 / pi, XYsign); print(text_sd)
    #-------- Gain adjustment
    print('  -- Polarized gain calibration')
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY* XYsign])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis = caledVis.transpose(0,2,1).dot(blWeight)
    #-------- Fine estimation of Q and U using XY and YX
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        Qsol, Usol = XY2Stokes(PA[timeIndex], Vis[[1,2]][:,timeIndex])
        text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/pi); print(text_sd)
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    Vis = caledVis.transpose(0,2,1).dot(blWeight)
    #-------- XY phase correction
    XYphase, DtotP, DtotM = XY2PhaseVec(mjdSec - np.median(mjdSec), Vis[[1,2]], UCmQS, QCpUS, 1000)
    twiddle = np.exp((1.0j)* XYphase)
    caledVis[1] /= twiddle
    caledVis[2] *= twiddle
    XYvis = Vis[[1,2]]; XYV = 0.5*(XYvis[0] + XYvis[1].conjugate())
    #-------- Display XY cross correlation
    plotMax = 1.2* max(abs(XYV))
    ArrayDx, ArrayDy = 0.5* (DtotP - DtotM), 0.5* (DtotP + DtotM)
    text_Dx, text_Dy = 'Array Dx = %+.4f %+.4fi' % (ArrayDx.real, ArrayDx.imag), 'Array Dy = %+.4f %+.4fi' % (ArrayDy.real, ArrayDy.imag)
    print(text_Dx + ' ' + text_Dy)
    figXY = plt.figure(figsize = (8, 8))
    XYP = figXY.add_subplot( 1, 1, 1 )
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 3 : continue
        plotY = DtotP + DtotM* QCpUS + UCmQS* np.exp((0.0 + 1.0j)*XYphase)
        XYP.plot( UCmQS[timeIndex], plotY[timeIndex].real, '-', label=text_Dx); XYP.plot( UCmQS[timeIndex], plotY[timeIndex].imag, '-', label=text_Dy)
        XYP.plot( UCmQS[timeIndex], XYV[timeIndex].real, '.', label='Re <XY*>'); XYP.plot( UCmQS[timeIndex], XYV[timeIndex].imag, '.', label='Im <XY*>')
    #
    XYP.set_xlabel('U $\cos 2 \psi $ - Q $\sin 2 \psi$'); XYP.set_ylabel('<XY*>'); XYP.set_title('%s-SPW%d-REF%s' % (prefix, spw, refantName))
    XYP.axis([-plotMax, plotMax, -plotMax, plotMax]); XYP.grid()
    XYP.legend(loc = 'best', prop={'size' :12}, numpoints = 1)
    plt.show()
    figXY.savefig('XY_%s-REF%s-SPW%d.pdf' % (prefix, refantName, spw))
    plt.close('all')
    del figXY
    XYvis[0] -= (DtotP + DtotM* QCpUS); XYvis[1] -= (DtotP + DtotM* QCpUS).conjugate()
    #-------- Fine estimation of Q and U using XY and YX
    print('  -- XY phase correction')
    Vis    = np.mean(caledVis, axis=1)
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        Qsol, Usol = XY2Stokes(PA[timeIndex], Vis[[1,2]][:,timeIndex])
        text_sd = '[XY,YX] %s:  Q/I= %6.3f  U/I= %6.3f p=%.2f%% EVPA = %6.2f deg' % (sourceName, Qsol, Usol, 100.0* np.sqrt(Qsol**2 + Usol**2), np.arctan2(Usol,Qsol)*90.0/pi); print(text_sd)
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
    #
    GainX, GainY = polariGain(caledVis[0], caledVis[3], QCpUS)
    GainY *= twiddle
    Gain = np.array([Gain[0]* GainX, Gain[1]* GainY])
    caledVis = chAvgVis / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    # Vis    = np.mean(caledVis, axis=1)
    GainCaledVisSpec = VisSpec.transpose(1,0,2,3) / (Gain[polYindex][:,ant0]* Gain[polXindex][:,ant1].conjugate())
    del VisSpec
    #-------- Antenna-based on-axis D-term (chAvg)
    StokesI = np.ones(PAnum)
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        caledVis[:,:,timeIndex] *= StokesDic[sourceName][0]
        QCpUS[timeIndex] *= StokesDic[sourceName][0]
        UCmQS[timeIndex] *= StokesDic[sourceName][0]
        StokesI[timeIndex] *= StokesDic[sourceName][0]
    #
    Dx, Dy = VisMuiti_solveD(caledVis, QCpUS, UCmQS, np.repeat(ArrayDx, antNum), np.repeat(ArrayDy, antNum), StokesI)
    #-------- D-term-corrected Stokes parameters
    Minv = InvMullerVector(Dx[ant1], Dy[ant1], Dx[ant0], Dy[ant0], np.ones(blNum, dtype=complex))
    print('  -- D-term-corrected visibilities')
    useTimeIndex = []
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        useTimeIndex = useTimeIndex + timeIndex
        srcTimeNum = len(timeIndex)
        if srcTimeNum < 1 : continue
        PS = InvPAVector(PA[timeIndex], np.ones(srcTimeNum))
        StokesVis = PS.reshape(4, 4*srcTimeNum).dot(Minv.reshape(4, 4*blNum).dot(caledVis[:,:,timeIndex].reshape(4*blNum, srcTimeNum)).reshape(4*srcTimeNum)) / (srcTimeNum* blNum)
        Isol, Qsol, Usol = StokesVis[0].real, StokesVis[1].real, StokesVis[2].real
        Ierr, Qerr, Uerr = abs(StokesVis[0].imag), abs(StokesVis[1].imag), abs(StokesVis[2].imag)
        text_sd = '%s: I= %6.3f  Q= %6.3f+-%6.4f  U= %6.3f+-%6.4f EVPA = %6.2f deg' % (sourceName, Isol, Qsol, Qerr, Usol, Uerr, np.arctan2(Usol,Qsol)*90.0/pi); print(text_sd)
        StokesDic[sourceName] = (np.array([Isol, Qsol, Usol, 0.0])).tolist()
        StokesI[timeIndex] = Isol
        QCpUS[timeIndex] = Qsol* CS[timeIndex] + Usol* SN[timeIndex]
        UCmQS[timeIndex] = Usol* CS[timeIndex] - Qsol* SN[timeIndex]
        #for index in timeIndex: GainCaledVisSpec[:,:,:,index] *= Isol
        #GainCaledVisSpec *= StokesI
        GainCaledVisSpec[:,:,:,timeIndex] *= Isol
    #
    useTimeIndex.sort()
    #-------- get D-term spectra
    print('  -- Determining D-term spectra')
    for ch_index in list(range(int(chNum/bunchNum))):
        DxSpec[:,ch_index], DySpec[:,ch_index] = VisMuiti_solveD(GainCaledVisSpec[ch_index][:,:,useTimeIndex], QCpUS[useTimeIndex], UCmQS[useTimeIndex], Dx, Dy, StokesI[useTimeIndex])
        progress = (ch_index + 1.0) / (chNum / bunchNum)
        sys.stderr.write('\r\033[K' + get_progressbar_str(progress)); sys.stderr.flush()
    #
    sys.stderr.write('\n'); sys.stderr.flush()
    if 'Dsmooth' in locals():
        node_index = list(range(3, chNum/bunchNum, Dsmooth))
        bunchedFreq = bunchVecCH(Freq)
        for ant_index in list(range(antNum)):
            DX_real, DX_imag = scipy.interpolate.splrep(bunchedFreq, DxSpec[ant_index].real, k=3, t=bunchedFreq[node_index]), scipy.interpolate.splrep(bunchedFreq, DxSpec[ant_index].imag, k=3, t=bunchedFreq[node_index])
            DY_real, DY_imag = scipy.interpolate.splrep(bunchedFreq, DySpec[ant_index].real, k=3, t=bunchedFreq[node_index]), scipy.interpolate.splrep(bunchedFreq, DySpec[ant_index].imag, k=3, t=bunchedFreq[node_index])
            DxSpec[ant_index] =scipy.interpolate.splev(bunchedFreq, DX_real) + (0.0 + 1.0j)* scipy.interpolate.splev(bunchedFreq, DX_imag)
            DySpec[ant_index] =scipy.interpolate.splev(bunchedFreq, DY_real) + (0.0 + 1.0j)* scipy.interpolate.splev(bunchedFreq, DY_imag)
        #
    #
    #-------- D-term-corrected visibilities (invD dot Vis = PS)
    del chAvgVis, StokesVis
    print('  -- Applying D-term spectral correction')
    M  = InvMullerVector(DxSpec[ant0], DySpec[ant0], DxSpec[ant1], DySpec[ant1], np.ones([blNum,int(chNum/bunchNum)])).transpose(0,3,1,2)
    StokesVis = np.zeros([4, int(chNum/bunchNum), PAnum], dtype=complex )
    for time_index in useTimeIndex: StokesVis[:, :, time_index] = 4.0* np.mean(M* GainCaledVisSpec[:,:,:,time_index], axis=(2,3))
    del GainCaledVisSpec
    chAvgVis = np.mean(StokesVis[:,chRange], axis=1)
    XYC = chAvgVis[[1,2]]
    PS = InvPAVector(PA, np.ones(PAnum))
    for ch_index in list(range(int(chNum/bunchNum))): StokesVis[:,ch_index] = np.sum(PS* StokesVis[:,ch_index], axis=1)
    maxP = 0.0
    for sourceName in sourceList:
        scanLS = scanDic[sourceName]
        colorIndex = lineCmap(sourceList.index(sourceName) / 8.0)
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        QUsol = np.array([StokesDic[sourceName][1], StokesDic[sourceName][2]])
        maxP = max(maxP, sqrt(QUsol.dot(QUsol)))
        EVPA = 0.5* np.arctan2(QUsol[1], QUsol[0])
        ThetaPlot = PA[timeIndex] - EVPA; ThetaPlot = np.arctan(np.tan(ThetaPlot))
        ThetaMin, ThetaMax = min(ThetaPlot), max(ThetaPlot)
        PArange = np.arange(ThetaMin + EVPA, ThetaMax + EVPA, 0.01)
        ThetaRange = np.arange(ThetaMin, ThetaMax, 0.01)
        CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
        UCmQS, QCpUS = QUsol[1]*CSrange - QUsol[0]* SNrange, QUsol[0]*CSrange + QUsol[1]* SNrange
        ThetaRange[ThetaRange >  1.56] = np.inf
        ThetaRange[ThetaRange < -1.56] = -np.inf
        plt.plot(RADDEG* ThetaRange,  QCpUS, '-', color=colorIndex, linestyle='dashed', label=sourceName + ' XX* - I')    # XX* - 1.0
        plt.plot(RADDEG* ThetaRange, -QCpUS, '-', color=colorIndex, linestyle='dashdot',label=sourceName + ' YY* - I')   # YY* - 1.0
        plt.plot(RADDEG* ThetaRange,  UCmQS, '-', color=colorIndex, linestyle='solid')
        plt.plot(RADDEG* ThetaRange,  np.zeros(len(ThetaRange)), '-', color=colorIndex, linestyle='dotted')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[0][timeIndex].real - StokesDic[sourceName][0], 'k.')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[1][timeIndex].real, 'b,', label=sourceName + ' ReXY*')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[1][timeIndex].imag, 'c,', label=sourceName + ' ImXY*')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[2][timeIndex].real, 'b,')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[2][timeIndex].imag, 'c,')
        plt.plot(RADDEG* ThetaPlot, chAvgVis[3][timeIndex].real - StokesDic[sourceName][0], 'k.')
        for scan in scanLS:
            scanIndex = scanList.index(scan)
            scanMJD = mjdSec[scanST[scanIndex]]
            text_sd = 'Scan %d : %s' % (scan, qa.time('%fs' % (scanMJD), form='fits', prec=6)[0][11:21])
            plt.text( RADDEG* np.arctan(np.tan(PA[scanST[scanIndex]] - EVPA)), -1.5* maxP, text_sd, verticalalignment='bottom', fontsize=6, rotation=90)
        #
    #
    plt.xlabel('Linear polarization angle w.r.t. X-Feed [deg]'); plt.ylabel('Cross correlations [Jy]')
    plt.xlim([-90.0, 90.0])
    plt.ylim([-1.5* maxP, 1.5*maxP])
    plt.legend(loc = 'best', prop={'size' :6}, numpoints = 1)
    plt.savefig('%s-SPW%d-%s-QUXY.pdf' % (prefix, spw, refantName))
    #-------- Save Results
    np.save('%s-SPW%d-%s.Ant.npy' % (prefix, spw, refantName), antList[antMap])
    np.save('%s-SPW%d-%s.Azel.npy' % (prefix, spw, refantName), np.array([mjdSec, Az, El, PA]))
    np.save('%s-SPW%d-%s.TS.npy' % (prefix, spw, refantName), mjdSec )
    np.save('%s-SPW%d-%s.GA.npy' % (prefix, spw, refantName), Gain )
    np.save('%s-SPW%d-%s.XYPH.npy' % (prefix, spw, refantName), XYphase )
    np.save('%s-SPW%d-%s.XYV.npy' % (prefix, spw, refantName), XYvis )
    np.save('%s-SPW%d-%s.XYC.npy' % (prefix, spw, refantName), XYC )
    for ant_index, ant in enumerate(antList[antMap]):
        DtermFile = np.array([FreqList[spw_index], DxSpec[ant_index].real, DxSpec[ant_index].imag, DySpec[ant_index].real, DySpec[ant_index].imag])
        np.save('%s-SPW%d-%s.DSpec.npy' % (prefix, spw, ant), DtermFile)
    #
    plt.close('all')
    #-------- Plot Stokes spectra
    polLabel, Pcolor = ['I', 'Q', 'U', 'V'], ['black', 'blue', 'red', 'green']
    for sourceName in sourceList:
        timeIndex = timeDic[sourceName]
        if len(timeIndex) < 1 : continue
        figSP = plt.figure(figsize = (11, 8))
        figSP.suptitle(prefix + ' ' + sourceName)
        figSP.text(0.45, 0.05, 'Frequency [GHz]')
        figSP.text(0.03, 0.45, 'Stokes visibility amplitude [Jy]', rotation=90)
        #
        StokesI_SP = figSP.add_subplot( 2, 1, 1 )
        StokesP_SP = figSP.add_subplot( 2, 1, 2 )
        StokesSpec, StokesErr =  np.mean(StokesVis[:,:,timeIndex], axis=2).real, abs(np.mean(StokesVis[:,:,timeIndex], axis=2).imag)
        np.save('%s-REF%s-%s-SPW%d.StokesSpec.npy' % (prefix, refantName, sourceName, spw), StokesSpec)
        np.save('%s-REF%s-%s-SPW%d.StokesErr.npy' % (prefix, refantName, sourceName, spw), StokesErr)
        np.save('%s-REF%s-%s-SPW%d.Freq.npy' % (prefix, refantName, sourceName, spw), Freq)
        StokesDic[sourceName] = np.mean(StokesSpec, axis=1).tolist()
        #
        IMax = np.max(StokesSpec[0])
        StokesI_SP.step(Freq[chRange], StokesSpec[0][chRange], where='mid', label=polLabel[0], color=Pcolor[0])
        StokesP_SP.step(Freq[chRange], StokesSpec[1][chRange], where='mid', label=polLabel[1], color=Pcolor[1])
        StokesP_SP.step(Freq[chRange], StokesSpec[2][chRange], where='mid', label=polLabel[2], color=Pcolor[2])
        StokesP_SP.step(Freq[chRange], StokesSpec[3][chRange], where='mid', label=polLabel[3], color=Pcolor[3])
        StokesI_SP.tick_params(axis='both', labelsize=6)
        StokesP_SP.tick_params(axis='both', labelsize=6)
        StokesI_SP.axis([np.min(Freq[chRange]), max(Freq[chRange]), 0.0, 1.25*IMax])
        StokesP_SP.axis([np.min(Freq[chRange]), max(Freq[chRange]), -0.15*IMax, 0.15*IMax])
        StokesI_SP.text(min(Freq[chRange]), IMax*1.35, sourceName)
        StokesI_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        StokesP_SP.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        figSP.savefig('SP_%s-REF%s-%s-SPW%d.pdf' % (prefix, refantName, sourceName, spw))
        plt.close('all')
    #
    DxList, DyList = DxList + [DxSpec], DyList + [DySpec]
    #---- Save Stokes parameters of the calibraors
    fileDic = open('Stokes.%s-SPW%d.dic' % (prefix, spw), mode='wb')
    pickle.dump(StokesDic, fileDic)
    fileDic.close()
    StokesTextFile = open('Stokes.%s-SPW%d.txt' % (prefix, spw), mode='w')
    text_sd = 'Source       I      Q      U      V      p%    EVPA'
    print(text_sd); StokesTextFile.write(text_sd + '\n')
    for sourceName in sourceList:
        if StokesDic[sourceName] == []: continue
        text_sd = '%s %6.3f %6.3f %6.3f %6.3f %5.2f %5.2f' % (sourceName, StokesDic[sourceName][0], StokesDic[sourceName][1], StokesDic[sourceName][2], StokesDic[sourceName][3], 100.0* sqrt(StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2)/StokesDic[sourceName][0], 90.0* np.arctan2(StokesDic[sourceName][2], StokesDic[sourceName][1]) / np.pi)
        print(text_sd); StokesTextFile.write(text_sd + '\n')
    #
    StokesTextFile.close()
#
