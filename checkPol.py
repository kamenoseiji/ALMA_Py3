import sys
import analysisUtils as au
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Grid.py').read())
fileNum = len(prefixList)
QAresult = ['Fail', 'Pass']
det_thresh, XY_thresh = 400, 0.05   # Determinant > 400, XY cross correlation > 50 mJy
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
    QA = 1    # pass QA0
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
        sourceList, posList = GetSourceList(msfile); sourceList = sourceRename(sourceList)
        PolScans = msmd.scansforintent("*CALIBRATE_POLARIZATION*")
        PolScanList = PolScanList + [PolScans[indexList(PolScans, np.array(bpscanLists[file_index]))].tolist()]
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
    if numPolSource == 0: conrinue
    StokesDic = dict(zip(polSourceList, [[]]*numPolSource))
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
    #-------- for each polcal source
    print('Source     : I[Jy]  Q[Jy]  U[Jy] p[%] EVPA[d]:   min   max   det   QA')
    figPL = plt.figure(figsize = (11, 8))
    lineCmap = plt.get_cmap('Set1')
    figPL.suptitle('Session %s : Band %d : Expected Cross Polarization' % (prefixList[0].split('_')[4], bandID) )
    PolPL = figPL.add_subplot( 1, 1, 1 )
    for src_index in list(range(numPolSource)):
        PolAZ, PolEL, PolPA, refTime = [], [], [], []
        colorIndex = lineCmap(int(src_index / 8.0))
        for file_index in list(range(fileNum)):
            prefix = prefixList[file_index]; msfile = wd + prefix + '.ms'
            #-------- AZEL in the MS
            azelTime, AntID, AZ, EL = GetAzEl(msfile)
            for trialID in list(range(64)):
                azelTime_index = np.where( AntID == trialID )[0].tolist() 
                if len(azelTime_index) > 1: break
            #
            for scan_index in list(range(len(PolScanList[file_index]))):
                if polScanSource[file_index][scan_index] != polSourceList[src_index] : continue
                interval, timeStamp = GetTimerecord(msfile, trialID, trialID, bpSPWList[file_index][0], PolScanList[file_index][scan_index])
                AzScan, ElScan = AzElMatch(timeStamp, azelTime, AntID, trialID, AZ, EL)
                PA = AzEl2PA(AzScan, ElScan) + BandPA[band_index]
                PolAZ, PolEL, PolPA, refTime = PolAZ + AzScan.tolist(), PolEL + ElScan.tolist(),  PolPA + PA.tolist(), refTime + timeStamp.tolist()
            #
        #
        CS, SN = np.cos(2.0* np.array(PolPA)), np.sin(2.0* np.array(PolPA))
        QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
        UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
        polDeg, EVPA = np.sqrt( StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2 ) / StokesDic[sourceName][0], 0.5* np.arctan2(StokesDic[sourceName][2],StokesDic[sourceName][1])
        det_D = np.sum(UCmQS**2)*len(UCmQS) - (np.sum(UCmQS))**2      # Determinant for D-term
        maxUCmQS, minUCmQS, maxXY = np.max(UCmQS), np.min(UCmQS), np.max(abs(UCmQS)) 
        if det_D < det_thresh: QA = int(QA*0)
        if maxXY < XY_thresh:  QA = int(QA*0)
        print('%s : %.3f %.3f %.3f %4.1f %6.1f : %5.2f %5.2f %.1f  %s' % (sourceName, StokesDic[sourceName][0], StokesDic[sourceName][1], StokesDic[sourceName][2], 100.0*polDeg, EVPA*180.0/np.pi, minUCmQS, maxUCmQS, det_D, QAresult[QA]))
        plotPA = np.array(PolPA) - EVPA; plotPA = np.arctan(np.tan(plotPA))
        ThetaPlot = np.array(PolPA) - EVPA; ThetaPlot = np.arctan(np.tan(ThetaPlot))
        ThetaMin, ThetaMax = min(ThetaPlot), max(ThetaPlot)
        PArange = np.arange(ThetaMin + EVPA, ThetaMax + EVPA, 0.01)
        ThetaRange = np.arange(ThetaMin, ThetaMax, 0.01)
        CSrange, SNrange = np.cos(2.0*PArange), np.sin(2.0*PArange)
        UCMQS, QCPUS = StokesDic[sourceName][2]*CSrange - StokesDic[sourceName][1]* SNrange, StokesDic[sourceName][1]*CSrange + StokesDic[sourceName][2]* SNrange
        ThetaRange[ThetaRange >  1.56] = np.inf
        ThetaRange[ThetaRange < -1.56] = -np.inf
        PolPL.plot(RADDEG* ThetaRange,  QCPUS, '-', color=lineCmap(src_index/8), linestyle='dashed', label=polSourceList[src_index] + ' XX* - I')     # XX* - 1.0
        PolPL.plot(RADDEG* ThetaRange,  UCMQS, '-', color=lineCmap((src_index + 2)/8), linestyle='solid', label=polSourceList[src_index] + ' Re(XY*)')     # Real part of XY*
        PolPL.plot(RADDEG* plotPA,  QCpUS, 'o', color=lineCmap((src_index + 4)/8) , label=polSourceList[src_index] + '(XX* - YY*)/2')     # Real part of XY*
        PolPL.plot(RADDEG* plotPA,  UCmQS, 'o', color=lineCmap((src_index + 6)/8), label=polSourceList[src_index] + ' Re(XY*)')     # Real part of XY*
    #
    PolPL.legend(loc = 'best', prop={'size' :8}, numpoints = 1)
    PolPL.set_xlabel('Linear polarization angle w.r.t. X-Feed [deg]'); PolPL.set_ylabel('Cross correlations [Jy]'); PolPL.set_title(prefixList); PolPL.grid()
    figPL.savefig('%s.Band%d.QUXY.png' % (prefixList[0].split('_')[4], bandID))
    plt.close('all')
#
