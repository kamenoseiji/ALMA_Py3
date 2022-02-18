import sys
import analysisUtils as au
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Grid.py').read())
exec(open(SCR_DIR + 'ASDM_XML.py').read())
fileNum = len(prefixList)
QAresult = ['Fail', 'Pass']
det_thresh, XY_thresh = 100, 0.05   # Determinant > 400, XY cross correlation > 50 mJy
if 'PHASECAL' not in locals(): PHASECAL = False
#-------- Check SPWs for polarization
bandNames, BandPA = [], []
scanList, SrcList, STList, ETList, sourceListXML, PosListXML = [], [], [], [], [], []
for prefix in prefixList:
    bandNames = bandNames + BandList( prefix )
    scanL, srcL, STL, ETL = PolScan( prefix )
    scanList, SrcList, STList, ETList = scanList + scanL, SrcList + srcL, STList + STL, ETList + ETL
    sourceList, PosList = SourceList( prefix )
    sourceListXML = sourceListXML + sourceList
    PosListXML    = PosListXML + PosList
#
#-------- Check Source Coordinate
uniqSrcList = list(set(sourceListXML)); srcNum = len(uniqSrcList)
SrcDic  = dict(zip(uniqSrcList, [[]]*srcNum))
for source in uniqSrcList:
    index = np.where(np.array(sourceListXML) == source)[0][0]
    SrcDic[source] = PosListXML[index]
#    
del sourceListXML, PosListXML
UniqBands = list(set(bandNames)); NumBands = len(UniqBands)
for band_index in list(range(NumBands)):
    bandID = int(UniqBands[band_index].split('_')[-1])
    BandPA = BandPA + [(BANDPA[bandID] + 90.0)*pi/180.0]
    print('%s : BandPA = %.2f' % (UniqBands[band_index], BandPA[band_index]))
    #
    polSourceList = unique(SrcList).tolist(); numPolSource = len(polSourceList)
    if numPolSource == 0: continue
    StokesDic = dict(zip(polSourceList, [[]]*numPolSource))
    ScansDic  = dict(zip(polSourceList, [[]]*numPolSource))
    for source in polSourceList: ScansDic[source] = []
    #---- Relate files and scans for each source
    for scan_index in list(range(len(scanList))):
        print('Scan %3d : %s %s - %s' % (scanList[scan_index], SrcList[scan_index], qa.time('%fs' % (STList[scan_index]*1.0e-9), form='ymd')[0], qa.time('%fs' % (ETList[scan_index]*1.0e-9), form='ymd')[0]))
        ScansDic[SrcList[scan_index]] = ScansDic[SrcList[scan_index]] + [scanList[scan_index]]
    #-------- Check QU catalog
    os.system('rm -rf CalQU.data')
    text_sd = R_DIR + 'Rscript %spolQuery.R -D%s -F%f' % (SCR_DIR, qa.time('%fs' % (np.median(STList)*1.0e-9), form='ymd')[0], BANDFQ[bandID])
    print(text_sd)
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
        QA = 1    # pass QA0
        sourceName = polSourceList[src_index]
        refTime, textPA, textSD = [], [], []
        colorIndex = plt.rcParams['axes.prop_cycle'].by_key()['color'][src_index % 10]
        scan_index = np.where(np.array(SrcList) == sourceName)[0].tolist()
        for index in scan_index: refTime = refTime + np.arange(1.0e-9* STList[index],  1.0e-9* ETList[index], 10).tolist()
        HA = gst2ha( mjd2gmst( np.array(refTime)/86400.0, 0.0), ALMA_long, SrcDic[sourceName][0])
        AZ, EL, PA = ha2azel( HA, SrcDic[sourceName][1] )
        #
        PA = PA + BandPA[band_index]
        CS, SN = np.cos(2.0* np.array(PA)), np.sin(2.0* np.array(PA))
        QCpUS = StokesDic[sourceName][1]*CS + StokesDic[sourceName][2]*SN   # Qcos + Usin
        UCmQS = StokesDic[sourceName][2]*CS - StokesDic[sourceName][1]*SN   # Ucos - Qsin
        combPA    = combPA + PA.tolist()
        combQCpUS = combQCpUS + QCpUS.tolist()
        combUCmQS = combUCmQS + UCmQS.tolist()
        polDeg, EVPA = np.sqrt( StokesDic[sourceName][1]**2 + StokesDic[sourceName][2]**2 ) / StokesDic[sourceName][0], 0.5* np.arctan2(StokesDic[sourceName][2],StokesDic[sourceName][1])
        det_D = np.sum(UCmQS**2)*len(UCmQS) - (np.sum(UCmQS))**2      # Determinant for D-term
        maxUCmQS, minUCmQS, maxXY = np.max(UCmQS), np.min(UCmQS), np.max(abs(UCmQS)) 
        if det_D < det_thresh: QA = int(QA*0)
        if maxXY < XY_thresh:  QA = int(QA*0)
        print('%s | %6.3f %6.3f %6.3f %4.1f %6.1f | %5.2f %5.2f %6.1f  %s' % (sourceName, StokesDic[sourceName][0], StokesDic[sourceName][1], StokesDic[sourceName][2], 100.0*polDeg, EVPA*180.0/np.pi, minUCmQS, maxUCmQS, det_D, QAresult[QA]))
        plotPA = np.array(PA) - EVPA; plotPA = np.arctan(np.tan(plotPA))
        ThetaPlot = np.array(PA) - EVPA; ThetaPlot = np.arctan(np.tan(ThetaPlot))
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
    combUCmQS = np.array(combUCmQS)
    det_D = np.sum(combUCmQS**2)*len(combUCmQS) - (np.sum(combUCmQS))**2      # Determinant for D-term
    if (det_D > det_thresh) & (np.max(abs(combUCmQS)) > XY_thresh) : QA = 1
    print('----------------------------------------------+-------------------------')
    print('Combined solution                             | %5.2f %5.2f %6.1f  %s' % (np.min(combUCmQS), np.max(combUCmQS), det_D, QAresult[QA]))
    PolPL.set_xlabel('Linear polarization angle w.r.t. X-Feed [deg]'); PolPL.set_ylabel('Cross correlations [Jy]'); PolPL.set_title(prefixList); PolPL.grid()
    figPL.savefig('%s.Band%d.QUXY.png' % (prefixList[0].split('_')[4], bandID))
    plt.close('all')
#
