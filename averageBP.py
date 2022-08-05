#---- Script for Band-3 Astroholograpy Data
import sys
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
#
XYSPW, BPSPW = [], []
if 'chTrim' not in locals(): chTrim = 0.06
chRange = list(range(int(chTrim*chNum), int((1.0 - chTrim)*chNum)))
#-------- Procedures
for spw in spwList:
    print('SPW %2d:---------------------------------------------------------------' % (spw))
    BPList, XYList = [], []
    for scan in scanList:
        BPList = BPList + [np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refant, scan, spw))]
        XYList = XYList + [np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, refant, scan, spw))]
    #
    BPant  = np.array(BPList)
    XYspec = np.array(XYList)
    scanNum, antNum, parapolNum  = BPant.shape[0], BPant.shape[1], BPant.shape[2]
    BPweight = np.zeros([scanNum, antNum, parapolNum], dtype=complex)
    #-------- Channel binning
    if 'bunchNum' not in locals(): bunchNum = 1
    if bunchNum > 1:
        def bunchVecCH(spec): return bunchVec(spec, bunchNum)
        BPant  = np.apply_along_axis(bunchVecCH, 3, BPant)
        XYspec = np.apply_along_axis(bunchVecCH, 1, XYspec)
    #---- Reference scan
    if 'BPscan' in locals():
        bpScanIndex = scanList.index(BPscan)
    else:
        bpScanIndex = np.argmax(abs(np.mean(XYspec, axis=1)))
    BPmean = np.mean(BPant, axis=0)
    XYmean = XYspec[bpScanIndex]; chNum = len(XYmean)
    for iter in list(range(10)):
        #-------- BP table 
        for ant_index in list(range(antNum)):
            for pol_index in list(range(parapolNum)):
                BPpower = np.sum(BPant[:, ant_index, pol_index]* BPant[:, ant_index, pol_index].conjugate(), axis=1).real
                BPcorr = BPant[:, ant_index, pol_index].dot(BPmean[ant_index, pol_index].conjugate()) / sqrt(BPpower* (BPmean[ant_index, pol_index].dot(BPmean[ant_index, pol_index].conjugate()).real))
                BPvar  = -np.log(abs(BPcorr))
                BPweight[:, ant_index, pol_index]  = BPcorr.conjugate() / (BPvar + np.percentile(BPvar, 100/scanNum))
                BPweight[:, ant_index, pol_index] = BPweight[:, ant_index, pol_index] / np.sum(abs(BPweight[:, ant_index, pol_index]))
            #
        #
        BPmean = np.sum(BPant.transpose(3,0,1,2)* BPweight, axis=1).transpose(1,2,0)
        BPscale = np.mean(abs(BPmean[:,:,chRange]), axis=2)
        BPmean = (BPmean.transpose(2,0,1) / BPscale).transpose(1,2,0)
        #-------- XY phase 
        XYcorr = XYspec.dot(XYmean.conjugate()) / chNum
        #text_amp, text_phs = '',''
        #for index in list(range(len(XYcorr))):
        #    text_amp = text_amp + '%5.2f ' % (abs(XYcorr[index]))
        #    text_phs = text_phs + '%5.1f ' % (180.0*np.angle(XYcorr[index])/np.pi)
        #
        #print(text_amp )
        #print(text_phs )
        XYsign = np.sign(XYcorr.real)
        if 'XYwgt' in locals():
            XYweight = XYsign* np.array(XYwgt)
        else:
            XYvar  = -np.log(abs(XYcorr))
            XYweight =  XYsign / (XYvar + np.percentile(XYvar, 100/scanNum))
        #
        XYmean   = (XYspec.T).dot(XYweight); XYmean = XYmean / abs(XYmean)
    #
    text_BPwgt, text_XYwgt, text_scan = 'BP wgt:', 'XY wgt:', 'Scan  :'
    for scan_index, scan in enumerate(scanList):
        text_scan   = text_scan   + '    %3d ' % (scan)
        text_BPwgt  = text_BPwgt + '%7.3f ' % (np.median(abs(BPweight), axis=(1,2))[scan_index])
        text_XYwgt  = text_XYwgt + '%7.1f ' % (XYweight[scan_index])
    #
    print(text_scan)
    print(text_BPwgt)
    print(text_XYwgt)
    XYSPW = XYSPW + [XYmean]
    BPSPW = BPSPW + [BPmean]
    np.save('%s-REF%s-SC0-SPW%d-BPant.npy'  % (prefix, refant, spw), BPmean) 
    np.save('%s-REF%s-SC0-SPW%d-XYspec.npy' % (prefix, refant, spw), XYmean) 
#
#-------- Plots
if 'BPPLOT' not in locals(): BPPLOT = False
if BPPLOT:
    pp = PdfPages('XYP_%s_REF%s_Scan0.pdf' % (prefix, refant))
    plotXYP(pp, prefix, spwList, XYSPW, bunchNum) 
    pp = PdfPages('BP_%s_REF%s_Scan0.pdf'  % (prefix, refant))
    if 'plotMin' not in locals(): plotMin = 0.0
    if 'plotMax' not in locals(): plotMax = 1.2
    antList = np.load('%s-REF%s.Ant.npy' % (prefix, refant))
    FreqList = []
    if bunchNum > 1:
        for spw in spwList: FreqList = FreqList + [bunchVecCH(np.load('%s-SPW%d-Freq.npy' % (prefix, spw)))]
    else:
        for spw in spwList: FreqList = FreqList + [np.load('%s-SPW%d-Freq.npy' % (prefix, spw))]
    plotSP(pp, prefix, antList, spwList, FreqList, BPSPW, plotMin, plotMax) 
#
del chRange
