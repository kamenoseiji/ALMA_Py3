####
# Script for Bispectra (closure phase)
####
import numpy as np
import datetime
#import matplotlib as plt
#exec(open(SCR_DIR + 'interferometry.py').read())
from interferometry import ANT0, ANT1, Ant2Bl, Ant2BlD, indexList, bestRefant, GetAntName, GetUVW, GetVisAllBL, ParaPolBL, specBunch, bunchVec
from Plotters import plotBispec
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-a', dest='antFlag', metavar='antFlag',
    help='Antennas to flag e.g. DA41,DV08', default='')
parser.add_option('-b', dest='BPprefix', metavar='BPprefix',
    help='Bandpass uid and scan to apply e.g. uid___A002_X10dadb6_X18e63', default='')
parser.add_option('-c', dest='scanList', metavar='scanList',
    help='Scan ID  e.g. 3,5,7', default='')
parser.add_option('-s', dest='spw', metavar='spw',
    help='SPW to plot e.g. 17', default='')
parser.add_option('-t', dest='timeBunch', metavar='timeBunch',
    help='Time average', default='1')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA42', default='')
parser.add_option('-S', dest='startTime', metavar='startTime',
    help='Start time e.g. 2020-03-03T14:11:25', default='')
#
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
antFlag = [ant for ant in options.antFlag.split(',')]
scanList= [int(scan) for scan in options.scanList.split(',')]
spw     = int(options.spw)
timeBunch = int(options.timeBunch)
if options.startTime != '': startMJD = qa.convert(options.startTime, 's')['value']
#if options.refant != '': refant = options.refant
refant = options.refant
#-------- Definitions
msfile = prefix + '.ms'
antList = GetAntName(msfile)
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
#blMap = list(range(blNum))
#if 'startTime' in locals(): startMJD = qa.convert(startTime, 's')['value']
#if 'BPscan' not in locals(): BPscan = 3 # default bandpass scan
#-------- Array Configuration
print('---Checking array configulation in scan %d' % (scanList[0]))
flagAnt = indexList(antFlag, antList)
UseAnt = list(set(range(antNum)) - set(flagAnt)); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
if refant not in antList[UseAnt]: refant = ''
if refant == '':
    timeStamp, UVW = GetUVW(msfile, spw, scanList[0])
    uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    refantID = bestRefant(uvDist)
else:
    refantID = np.where(antList[UseAnt] == refant)[0][0]
print( '  Use %s as the refant.' % (antList[UseAnt[refantID]]))
#-------- Baseline Mapping
print('---Baseline Mapping')
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
ant0 = ANT0[0:UseBlNum]; ant1 = ANT1[0:UseBlNum]
blMap, blInv = Ant2BlD(np.array(antMap)[ant0], np.array(antMap)[ant1])
print( '  %d baselines are inverted.' % (len(np.where( blInv )[0])))
'''
print('---Checking array configuration')
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
text_sd = '  Usable antennas (%d) : ' % (len(UseAnt))
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print(text_sd)
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
blMap = Ant2Bl(np.array(antMap)[ant0], np.array(antMap)[ant1])
#blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
#for bl_index in list(range(UseBlNum)): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spw, scanList[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
if 'refant' in locals():    refantID = indexList(np.array([refant]), antList[UseAnt])[0]
else: refantID = bestRefant(uvDist)
print('  Use ' + antList[UseAnt[refantID]] + ' as the refant.')
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print('  %d baselines are inverted.' % (len(np.where( blInv )[0])))
'''
figInch  = max(16,UseAntNum)
fontSize = min(32, figInch)
#-------- Bandpass Table
#if 'BPprefix' in locals():
if options.BPprefix != '':
    BPfileName = '%s-REF%s-SC%d-SPW%d-BPant.npy' % (options.BPprefix, antList[UseAnt[refantID]], BPscan, spw)
    print('---Loading bandpass table : ' + BPfileName)
    BP_ant = np.load(BPfileName)
#-------- Loop for Scan
DT = []
for scan_index, scan in enumerate(scanList):
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]; print
    polIndex = [0,3] if polNum == 4 else [0, 1] if polNum == 2 else [0]
    if chNum == 1:  chRange = [0]
    if 'chRange' not in locals(): chRange = list(range(int(0.05*chNum), int(0.95*chNum)))
    tempSpec = ParaPolBL(Xspec[polIndex][:,:,blMap], blInv).transpose(3,2,0,1)  # Parallel Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    if 'BP_ant' in locals():
        print('Applying bandpass calibration...')
        BPCaledXspec = (tempSpec / (BP_ant[ant0]* BP_ant[ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    else:
        BPCaledXspec = tempSpec.transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    if timeBunch > 1:
        chAvgVis = np.array([specBunch(chAvgVis[0], 1, timeBunch), specBunch(chAvgVis[1], 1, timeBunch)])   # chAvgVis[pol, bl, time]
        timeNum = chAvgVis.shape[2]
        timeStamp = bunchVec(timeStamp, timeBunch)
    #
    for mjdSec in timeStamp.tolist(): #DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
        DT = DT + [datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f')]
    if scan_index == 0:
        scanVis = chAvgVis
        ST_text = qa.time('%fs' % (timeStamp[0]), form='fits', prec=6)[0]
    else:
        scanVis = np.r_['2', scanVis, chAvgVis]
    #
#
polNum = len(polIndex)
polColor, polName = ['b','g'], ['X', 'Y']
ET_text = qa.time('%fs' % (timeStamp[-1]), form='fits', prec=6)[0]
#-------- Prepare Plots
labelList = ['UTC on %s' % (DT[0].strftime('%Y-%m-%d')), 'Closure Phase [rad]', 'Amplitude', '%s SPW=%d %s - %s' % (prefix, spw, ST_text, ET_text)]
plotFile  = 'BS_%s_SPW%d' % (prefix, spw)
pMax = np.percentile(abs(scanVis), 98) if 'plotMax' not in locals() else plotMax
plotBispec(antList[antMap], scanVis, DT, plotFile, labelList, pMax)
os.system('pdftoppm -png %s %s' % (plotFile + '.pdf', plotFile))
