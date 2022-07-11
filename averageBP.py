#---- Script for Band-3 Astroholograpy Data
import sys
exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
#
XYSPW, BPSPW = [], []
#-------- Procedures
for spw in spwList:
    BPList, XYList = [], []
    for scan in scanList:
        BPList = BPList + [np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refant, scan, spw))]
        XYList = XYList + [np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, refant, scan, spw))]
    #
    BPant  = np.array(BPList)
    XYspec = np.array(XYList)
    #---- Reference scan
    if 'BPscan' in locals():
        bpScanIndex = scanList.index(BPscan)
    else:
        bpScanIndex = np.argmax(abs(np.mean(XYspec, axis=1)))
    #
    Bpweight = 1.0 / np.var(BPant, axis=3)  # BPweight[scan, ant, pol]
    BPmean = (np.sum(BPant.transpose(3,0,1,2)* Bpweight, axis=1) / np.sum(Bpweight, axis=0)).transpose(1,2,0)
    XYmean = XYspec[bpScanIndex]; chNum = len(XYmean)
    for iter in list(range(10)):
        XYcorr = XYspec.dot(XYmean.conjugate()) / chNum
        #text_amp, text_phs = '',''
        #for index in list(range(len(XYcorr))):
        #    text_amp = text_amp + '%5.2f ' % (abs(XYcorr[index]))
        #    text_phs = text_phs + '%5.1f ' % (180.0*np.angle(XYcorr[index])/np.pi)
        #
        #print(text_amp )
        #print(text_phs )
        XYsign = np.sign(XYcorr.real)
        XYvar  = -np.log(abs(XYcorr))
        XYweight =  XYsign / (XYvar + np.percentile(XYvar, 100/len(scanList)))
        XYmean   = (XYspec.T).dot(XYweight); XYmean = XYmean / abs(XYmean)
    #
    XYSPW = XYSPW + [XYmean]
    BPSPW = BPSPW + [BPmean]
    np.save('%s-REF%s-SC0-SPW%d-BPant.npy'  % (prefix, refant, spw), BPmean) 
    np.save('%s-REF%s-SC0-SPW%d-XYspec.npy' % (prefix, refant, spw), XYmean) 
#
#-------- Plots
if 'BPPLOT' not in locals(): BPPLOT = False
if BPPLOT:
    pp = PdfPages('XYP_%s_REF%s_Scan0.pdf' % (prefix, refant))
    plotXYP(pp, prefix, spwList, XYSPW) 
    pp = PdfPages('BP_%s_REF%s_Scan0.pdf'  % (prefix, refant))
    if 'plotMin' not in locals(): plotMin = 0.0
    if 'plotMax' not in locals(): plotMax = 1.2
    antList = np.load('%s-REF%s.Ant.npy' % (prefix, refant))
    FreqList = []
    for spw in spwList: FreqList = FreqList + [np.load('%s-SPW%d-Freq.npy' % (prefix, spw))]
    plotSP(pp, prefix, antList, spwList, FreqList, BPSPW, plotMin, plotMax) 
#
