import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from interferometry import delay_search
from Plotters import plotSP, plotXYP
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-b', dest='bunchNum', metavar='bunchNum',
    help='Channel binning', default='1')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan ID  e.g. 3,5,7', default='')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA42', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
parser.add_option('-P', dest='PLOTPDF', metavar='PLOTPDF',
    help='Plot PDF', action="store_true")
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
refant  = options.refant
bunchNum=  int(options.bunchNum)
spwList = [int(spw) for spw in options.spwList.split(',')]
scanList= [int(scan) for scan in options.scanList.split(',')]
PLOTPDF = options.PLOTPDF
#-------- Procedures
XYSPW, BPSPW, XYdelayList = [], [], []
FreqList = [np.load('%s-SPW%d-Freq.npy' % (prefix, spw)) for spw in spwList]
for spw_index, spw in enumerate(spwList):
    print('SPW %2d:---------------------------------------------------------------' % (spw))
    BPant   = np.array([np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refant, scan, spw)) for scan in scanList])
    XYspec  = np.array([np.load('%s-REF%s-SC%d-SPW%d-XYspec.npy' % (prefix, refant, scan, spw)) for scan in scanList])
    scanNum, antNum, parapolNum  = BPant.shape[0], BPant.shape[1], BPant.shape[2]
    #-------- Channel binning
    if 'bunchNum' not in locals(): bunchNum = 1
    if bunchNum > 1:
        def bunchVecCH(spec): return bunchVec(spec, bunchNum)
        BPant  = np.apply_along_axis(bunchVecCH, 3, BPant)
        XYspec = np.apply_along_axis(bunchVecCH, 1, XYspec)
        FreqList[spw_index] = bunchVecCH(FreqList[spw_index])
    BW = max(FreqList[spw_index]) - min(FreqList[spw_index])
    #-------- Weight and phase using the reference SPW
    if spw_index == 0:
        BPmean = np.mean(BPant, axis=0)
        chNum = BPmean.shape[2]
        chRange = list(range(int(0.1*chNum), int(0.95*chNum)))
        #-------- Iteration to optimize BPweight and XYweight
        for iter in list(range(10)):
            #-------- BP table 
            BPcorr = BPant * BPmean.conjugate()
            BPweight = np.mean(BPcorr[:,:,:,chRange], axis=3).conjugate()
            BPmean =  np.mean(BPant.transpose(3,0,1,2)* BPweight, axis=1).transpose(1,2,0)
        #-------- XY phase
        bpScanIndex = 0
        if 'BPscan' in locals(): bpScanIndex = scanList.index(BPscan)   # reference scan for XYspec 
        XYmean = XYspec[bpScanIndex]
        XYcorr = XYspec.dot(XYmean.conjugate()) / chNum
        XYsign = np.sign(XYcorr.real)
        if 'XYwgt' in locals():
            XYweight = XYsign* np.array(XYwgt)
        else:
            XYvar  = -np.log(abs(XYcorr) - 1.0e-3)
            XYweight =  XYsign / (XYvar + np.percentile(XYvar, 100/scanNum))
        #
    else:
        BPmean =  np.mean(BPant.transpose(3,0,1,2)* BPweight, axis=1).transpose(1,2,0)
        chNum = BPmean.shape[2]
        chRange = list(range(int(0.1*chNum), int(0.95*chNum)))
    #
    BPmean = (BPmean.transpose(2,0,1) /  np.mean(abs(BPmean[:,:,chRange]), axis=2)).transpose(1,2,0)
    XYmean   = (XYspec.T).dot(XYweight); XYmean = XYmean / abs(XYmean)
    XYdelay, snr = delay_search(XYmean); XYdelayList = XYdelayList + [0.5e9*XYdelay/BW] 
    text_BPwgt, text_XYwgt, text_scan = 'BP phs:', 'XY wgt:', 'Scan  :'
    for scan_index, scan in enumerate(scanList):
        text_scan   = text_scan   + '    %3d ' % (scan)
        text_BPwgt  = text_BPwgt + '%7.3f ' % (np.median(180* np.angle(BPweight)/np.pi, axis=(1,2))[scan_index])
        text_XYwgt  = text_XYwgt + '%7.1f ' % (XYweight[scan_index])
    #
    print(text_scan)
    print(text_BPwgt)
    print(text_XYwgt)
    #
    XYSPW = XYSPW + [XYmean]
    BPSPW = BPSPW + [BPmean]
    np.save('%s-REF%s-SC0-SPW%d-BPant.npy'  % (prefix, refant, spw), BPmean) 
    np.save('%s-REF%s-SC0-SPW%d-XYspec.npy' % (prefix, refant, spw), XYmean) 
#
#-------- Plots
if 'PLOTPDF' not in locals(): PLOTPDF = False
if PLOTPDF:
    pp = PdfPages('XYP_%s_REF%s_Scan0.pdf' % (prefix, refant))
    plotXYP(pp, prefix, spwList, XYSPW, XYdelayList, bunchNum) 
    pp = PdfPages('BP_%s_REF%s_Scan0.pdf'  % (prefix, refant))
    if 'plotMin' not in locals(): plotMin = 0.0
    if 'plotMax' not in locals(): plotMax = 1.2
    antList = np.load('%s-REF%s.Ant.npy' % (prefix, refant))
    plotSP(pp, prefix, antList, spwList, FreqList, BPSPW, plotMin, plotMax) 
#
del chRange
