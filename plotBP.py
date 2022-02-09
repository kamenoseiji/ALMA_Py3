import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
#exec(open(SCR_DIR + 'interferometry.py').read())
exec(open(SCR_DIR + 'Plotters.py').read())
#
#-------- Bandpass Table
if 'scanList' not in locals(): scanList = []
antList = np.load('%s-REF%s.Ant.npy' % (prefix, refantName))
BPList, FreqList = [], []
spwNum, BPscanNum, scanNum  = len(spwList), len(BPscanList), len(scanList)
if 'bunchNum' not in locals(): bunchNum = 1
SPplot = False
if scanNum > 0: SPplot = True   # show bandpass-corrected spectral shape
for spw_index in list(range(spwNum)):
    Freq   = np.load('%s-SPW%d-Freq.npy' % (prefix, spwList[spw_index]))
    #-- First BP scan
    BP_ant = np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refantName, BPscanList[0], spwList[spw_index]))
    #-- 2nd BP scan and later
    for scan_index in list(range(1, BPscanNum)):
        tempBP = np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refantName, BPscanList[scan_index], spwList[spw_index]))
        BP_ant = BP_ant + tempBP
    BP_ant = BP_ant / BPscanNum
    if SPplot:
        #-- First Target scan
        SP_ant = np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refantName, scanList[0], spwList[spw_index]))
        for scan_index in list(range(1, scanNum)):
            tempBP = np.load('%s-REF%s-SC%d-SPW%d-BPant.npy' % (prefix, refantName, scanList[scan_index], spwList[spw_index]))
            SP_ant = SP_ant + tempBP
        SP_ant = SP_ant / scanNum
        BP_ant = SP_ant / BP_ant
    #
    BPList, FreqList = BPList + [BP_ant], FreqList + [Freq]
#
ppolNum = BP_ant.shape[1]
PolList = ['X', 'Y']
#
#-------- Plots
if SPplot:
    pp = PdfPages('SP_%s_REF%s_Scan%d.pdf' % (prefix, refantName, BPscan))
else:
    pp = PdfPages('BP_%s_REF%s_Scan%d.pdf' % (prefix, refantName, BPscan))
#if 'spurRFLists' in locals():
#        plotSP(pp, prefix, antList[antMap], spwList, BPscan, BPList, bunchNum, 1.2, spurRFLists) 
#    else:
if 'plotMin' not in locals(): plotMin = 0.0
if 'plotMax' not in locals(): plotMax = 1.2
plotSP(pp, prefix, antList, spwList, FreqList, BPList, plotMin, plotMax) 
#
