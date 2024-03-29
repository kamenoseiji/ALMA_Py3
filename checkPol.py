import sys
import analysisUtils as au
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Grid.py').read())
fileNum = len(prefixList)
QAresult = ['Fail', 'Pass']
QAresult = ['Fail', 'Pending QA0+', 'Pass']
det_thresh, XY_thresh = 20, 0.05   # Determinant > 20, XY cross correlation > 50 mJy
if 'PHASECAL' not in locals(): PHASECAL = False
#-------- Check SPWs for polarization
bpSPWList, bandNames, BandPA = [], [], []
for prefix in prefixList:
    msfile = wd + prefix + '.ms'
    bpSPWs = GetBPcalSPWs(msfile); bpSPWList = bpSPWList + [bpSPWs]
    print('---Checking spectral windows for %s : ' % (prefix), end=''); print(bpSPWs)
    msmd.open(msfile)
    bpspwNames = msmd.namesforspws(bpSPWs)
    bandNamePattern, bandNamesInEB = r'RB_..', []
    for spwName in bpspwNames : bandNamesInEB = bandNamesInEB + re.findall(bandNamePattern, spwName)
    bandNames = bandNames + [bandNamesInEB]
    msmd.close(); msmd.done()
#
UniqBands = unique(bandNames).tolist(); NumBands = len(UniqBands)
for band_index in list(range(NumBands)):
    BandPA = BandPA + [(BANDPA[int(UniqBands[band_index][3:5])] + 90.0)*pi/180.0]
    #print('%s : BandPA = %.2f' % (UniqBands[band_index], BandPA[band_index]))
#
#-------- Loop for MS files
for band_index in list(range(NumBands)):
    bpscanLists = []
    bandID = int(UniqBands[band_index][3:5])
    for file_index in list(range(fileNum)):
        prefix = prefixList[file_index]; msfile = wd + prefix + '.ms'
        bpSPWs = bpSPWList[file_index]
        msmd.open(msfile)
        #-------- Check SPWs for polarization
        bpspwNames = bandNames[file_index]
        bpspwList  = np.array(bpSPWs)[indexList( np.array([UniqBands[band_index]]), np.array(bpspwNames))].tolist()
        bpscanLists = bpscanLists + [msmd.scansforspw(bpspwList[0]).tolist()]
        msmd.close(); msmd.done()
    #
    interval, refTime = GetTimerecord(msfile, 0, 0, bpSPWs[0], bpscanLists[file_index][0])
    #-------- Check polcal source list
    PolScanList, polSourceList, polScanSource = [], [], []
    for file_index in list(range(fileNum)):
        prefix = prefixList[file_index]; msfile = wd + prefix + '.ms'
        msmd.open(msfile)
        #print('---Checking source list for %s' % (prefix))
        srcDic = GetSourceDic(msfile)
        sourceList = list(dict.fromkeys([ srcDic[ID]['Name'] for ID in srcDic.keys() ]))
        #sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList)
        PolScans = msmd.scansforintent("*CALIBRATE_BANDPASS*")
        PolScans = np.append(PolScans, msmd.scansforintent("*CALIBRATE_POLARIZATION*"))
        if PHASECAL: PolScans = np.append(PolScans, msmd.scansforintent("*CALIBRATE_PHASE*"))
        PolScanList = PolScanList + [PolScans[indexList(np.array(bpscanLists[file_index]), PolScans)].tolist()]
        polSourceID, polSourceName = [], []
        for PolScan in PolScanList[file_index]:
            sourceID = msmd.sourceidforfield(msmd.fieldsforscan(PolScan))
            polSourceID = polSourceID + [sourceID]
            polSourceName = polSourceName + [sourceList[sourceID]]
        #
        polSourceIDList = unique(polSourceID).tolist()
        for sourceID in polSourceIDList: polSourceList = polSourceList + [sourceList[sourceID]]
        polScanSource = polScanSource + [polSourceName]
        msmd.close(); msmd.done()
    #
    polSourceList = unique(polSourceList).tolist(); numPolSource = len(polSourceList)
    if numPolSource == 0: continue
    StokesDic = dict(zip(polSourceList, [[]]*numPolSource))
    ScansDic  = dict(zip(polSourceList, [[]]*numPolSource))
    for source in polSourceList: ScansDic[source] = [[]]*fileNum
    #---- Relate files and scans for each source
    for file_index in list(range(fileNum)):
        for scan_index in list(range(len(PolScanList[file_index]))):
            print('Scan %3d : %s ' % (PolScanList[file_index][scan_index], polScanSource[file_index][scan_index]))
            ScansDic[polScanSource[file_index][scan_index]][file_index] = ScansDic[polScanSource[file_index][scan_index]][file_index] + [PolScanList[file_index][scan_index]]
        #
    #-------- Check QU catalog
    os.system('rm -rf CalQU.data')
    text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, qa.time('%fs' % (refTime[0]), form='ymd')[0], BANDFQ[bandID])
    for source in polSourceList: text_sd = text_sd + ' ' + source
    os.system(text_sd)
    fp = open('CalQU.data')
    lines = fp.readlines()
    fp.close()
    for eachLine in lines:
        sourceName = eachLine.split()[0]
        StokesDic[sourceName] = [float(eachLine.split()[1]), float(eachLine.split()[2]), float(eachLine.split()[3]), 0.0]
    #
    for sourceName in polSourceList:
        if len(StokesDic[sourceName]) < 4:
            StokesDic.pop(sourceName)
            polSourceList.remove(sourceName)
    #
    numPolSource = len(StokesDic)
    #-------- for each polcal source
    combPA, combQCpUS, combUCmQS = [], [], []
    print('Source     |  I[Jy]  Q[Jy]  U[Jy] p[%] EVPA[d]|   min   max    det  QA')
    print('-----------+----------------------------------+-------------------------')
    figPL = plt.figure(figsize = (11, 8))
    lineCmap = plt.get_cmap('Set1')
    figPL.suptitle('Session %s : Band %d : Expected Cross Polarization' % (prefixList[0].split('_')[4], bandID) )
    PolPL = figPL.add_subplot( 1, 1, 1 )
    maxP = 0.01
    for src_index in list(range(numPolSource)):
        QA = 0    # 0:fail, 1:Pending QA0+, 2:pass
        sourceName = polSourceList[src_index]
        PolAZ, PolEL, PolPA, refTime, textPA, textSD = [], [], [], [], [], []
        colorIndex = plt.rcParams['axes.prop_cycle'].by_key()['color'][src_index % 10]
        for file_index in list(range(fileNum)):
            prefix = prefixList[file_index]; msfile = wd + prefix + '.ms'
            #-------- AZEL in the MS
            azelTime, AntID, AZ, EL = GetAzEl(msfile)
            for trialID in list(range(64)):
                azelTime_index = np.where( AntID == trialID )[0].tolist() 
                if len(azelTime_index) > 1: break
            #
            for scanID in ScansDic[sourceName][file_index]:
                interval, timeStamp = GetTimerecord(msfile, trialID, trialID, bpSPWList[file_index][0], scanID)
                AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, trialID, AZ, EL)
                PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]
                PolAZ, PolEL, PolPA, refTime = PolAZ + AzScan.tolist(), PolEL + ElScan.tolist(),  PolPA + PA.tolist(), refTime + timeStamp.tolist()
                text_sd = 'Scan %d : %s %s' % (scanID, qa.time('%fs' % (timeStamp[0]), form='fits', prec=6)[0][11:21], sourceName)
                textSD = textSD + [text_sd]
                textPA = textPA + [PA[0]]
            #
        #
        CS, SN = np.cos(2.0* np.array(PolPA)), np.sin(2.0* np.array(PolPA))
        QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
        UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
        combPA    = combPA + PolPA
        combQCpUS = combQCpUS + QCpUS.tolist()
        combUCmQS = combUCmQS + UCmQS.tolist()
        polDeg, EVPA = np.sqrt( StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2 ) / StokesDic[sourceName][0], 0.5* np.arctan2(StokesDic[sourceName][2],StokesDic[sourceName][1])
        det_D = (UCmQS.dot(UCmQS)* len(UCmQS) - np.sum(UCmQS)**2) / UCmQS.dot(UCmQS) # Determinant for D-term
        maxUCmQS, minUCmQS, maxXY = np.max(UCmQS), np.min(UCmQS), np.max(abs(UCmQS)) 
        if det_D > det_thresh: QA += 1
        if maxXY > XY_thresh:  QA += 1
        print('%s | %6.3f %6.3f %6.3f %4.1f %6.1f | %5.2f %5.2f %6.1f  %s' % (sourceName, StokesDic[sourceName][0], StokesDic[sourceName][1], StokesDic[sourceName][2], 100.0*polDeg, EVPA*180.0/np.pi, minUCmQS, maxUCmQS, det_D, QAresult[QA]))
        plotPA = np.array(PolPA) - EVPA; plotPA = np.arctan(np.tan(plotPA))
        ThetaPlot = np.array(PolPA) - EVPA; ThetaPlot = np.arctan(np.tan(ThetaPlot))
        ThetaMin, ThetaMax = min(ThetaPlot), max(ThetaPlot)
        PArange = np.arange(ThetaMin + EVPA, ThetaMax + EVPA, 0.01)
        ThetaRange = np.arange(ThetaMin, ThetaMax, 0.01)
        CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
        UCMQS, QCPUS = StokesDic[sourceName][2]*CSrange - StokesDic[sourceName][1]* SNrange, StokesDic[sourceName][1]*CSrange + StokesDic[sourceName][2]* SNrange
        maxP = max(maxP, np.max(abs(QCPUS)), np.max(abs(UCMQS)))
        ThetaRange[ThetaRange >  1.56] = np.inf
        ThetaRange[ThetaRange < -1.56] = -np.inf
        PolPL.plot(RADDEG* ThetaRange,  QCPUS, '-', color=colorIndex, linestyle='dashed', linewidth=0.5, label=sourceName + ' XX* - I')     # XX* - 1.0
        PolPL.plot(RADDEG* ThetaRange,  UCMQS, '-', color=colorIndex, linestyle='solid', label=sourceName + ' Re(XY*)')     # Real part of XY*
        PolPL.plot(RADDEG* plotPA,  QCpUS, '.', color=colorIndex , markersize=0.5, label=sourceName + '(XX* - YY*)/2')     # Real part of XY*
        PolPL.plot(RADDEG* plotPA,  UCmQS, 'o', color=colorIndex, label=sourceName + ' Re(XY*)')     # Real part of XY*
        for index in list(range(len(textPA))):
            plt.text(RADDEG* np.arctan(np.tan(textPA[index] - EVPA)), -1.09* maxP, textSD[index], verticalalignment='bottom', fontsize=6, rotation=90)
        #
    #
    plt.ylim([-1.1* maxP, 1.1*maxP])
    if numPolSource < 10: PolPL.legend(loc = 'best', prop={'size' :8}, numpoints = 1)
    QA = 0
    combUCmQS = np.array(combUCmQS)
    det_D = (combUCmQS.dot(combUCmQS)* len(combUCmQS) - np.sum(combUCmQS)**2) / combUCmQS.dot(combUCmQS) # Determinant for D-term
    if det_D > det_thresh: QA += 1
    if np.max(abs(combUCmQS)) > XY_thresh : QA += 1
    print('----------------------------------------------+-------------------------')
    print('Combined solution                             | %5.2f %5.2f %6.1f  %s' % (np.min(combUCmQS), np.max(combUCmQS), det_D, QAresult[QA]))
    PolPL.set_xlabel('Linear polarization angle w.r.t. X-Feed [deg]'); PolPL.set_ylabel('Cross correlations [Jy]'); PolPL.set_title(prefixList); PolPL.grid()
    figPL.savefig('%s.Band%d.QUXY.png' % (prefixList[0].split('_')[4], bandID))
    plt.close('all')
#
