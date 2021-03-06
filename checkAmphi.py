####
# Script for all-baseline time - visibility amplitue/phase plot
####
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
exec(open(SCR_DIR + 'interferometry.py').read())
#-------- Definitions
if 'antFlag' not in locals(): antFlag = []
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
blMap = list(range(blNum))
if 'startTime' in locals(): startMJD = qa.convert(startTime, 's')['value']
if 'BPscan' not in locals(): BPscan = 3 # default bandpass scan
#
print('---Checking array configuration')
flagAnt = np.ones([antNum]); flagAnt[indexList(antFlag, antList)] = 0.0
UseAnt = np.where(flagAnt > 0.0)[0].tolist(); UseAntNum = len(UseAnt); UseBlNum  = int(UseAntNum* (UseAntNum - 1) / 2)
text_sd = '  Usable antennas (%d) : ' % (len(UseAnt))
for ants in antList[UseAnt].tolist(): text_sd = text_sd + ants + ' '
print(text_sd)
blMap, blInv= list(range(UseBlNum)), [False]* UseBlNum
ant0, ant1 = ANT0[0:UseBlNum], ANT1[0:UseBlNum]
for bl_index in list(range(UseBlNum)): blMap[bl_index] = Ant2Bl(UseAnt[ant0[bl_index]], UseAnt[ant1[bl_index]])
timeStamp, UVW = GetUVW(msfile, spw, scanList[0])
uvw = np.mean(UVW[:,blMap], axis=2); uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
if 'refant' in locals():    refantID = indexList(np.array([refant]), antList[UseAnt])[0]
else: refantID = bestRefant(uvDist)
print('  Use ' + antList[UseAnt[refantID]] + ' as the refant.')
antMap = [UseAnt[refantID]] + list(set(UseAnt) - set([UseAnt[refantID]]))
for bl_index in list(range(UseBlNum)): blMap[bl_index], blInv[bl_index]  = Ant2BlD(antMap[ant0[bl_index]], antMap[ant1[bl_index]])
print('  %d baselines are inverted.' % (len(np.where( blInv )[0])))
figInch  = max(16,UseAntNum)
fontSize = min(32, figInch)
#-------- Bandpass Table
if 'BPprefix' in locals():
    BPfileName = '%s-REF%s-SC%d-SPW%d-BPant.npy' % (BPprefix, antList[UseAnt[refantID]], BPscan, spw)
    print('---Loading bandpass table : ' + BPfileName)
    BP_ant = np.load(BPfileName)
#
#-------- Loop for Scan
scanNum = len(scanList)
DT = []
for scan_index in range(scanNum):
    scan = scanList[scan_index]
    #-------- Baseline-based cross power spectra
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spw, scan)
    timeNum, polNum, chNum = Xspec.shape[3], Xspec.shape[0], Xspec.shape[1]; print
    polIndex = [0,3] if polNum == 4 else [0, 1] if polNum == 2 else [0]
    if chNum == 1:  chRange = [0]
    if 'chRange' not in locals(): chRange = list(range(int(0.05*chNum), int(0.95*chNum)))
    tempSpec = ParaPolBL(Xspec[polIndex][:,:,blMap], blInv).transpose(3,2,0,1)  # Parallel Polarization Baseline Mapping : tempSpec[time, blMap, pol, ch]
    Pspec = Pspec[polIndex][:,:,antMap]
    if 'BP_ant' in locals():
        print('Applying bandpass calibration...')
        BPCaledXspec = (tempSpec / (BP_ant[ant0]* BP_ant[ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    else:
        BPCaledXspec = tempSpec.transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    chAvgAC  = np.mean(Pspec[:,chRange], axis=1)
    if 'timeBunch' in locals():
        chAvgVis = np.array([specBunch(chAvgVis[0], 1, timeBunch), specBunch(chAvgVis[1], 1, timeBunch)])   # chAvgVis[pol, bl, time]
        chAvgAC  = np.array([specBunch(chAvgAC[0], 1, timeBunch),  specBunch(chAvgAC[1], 1, timeBunch)])    # chAvgAC[pol, ant, time]
        timeNum = chAvgVis.shape[2]
        timeStamp = bunchVec(timeStamp, timeBunch)
    #
    for mjdSec in timeStamp.tolist(): #DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
        DT = DT + [datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f')]
    if scan_index == 0:
        scanVis = chAvgVis
        scanAC  = chAvgAC
        ST_text = qa.time('%fs' % (timeStamp[0]), form='fits', prec=6)[0]
    else:
        scanVis = np.r_['2', scanVis, chAvgVis]
        scanAC  = np.r_['2', scanAC,  chAvgAC]
    #
#
polNum = len(polIndex)
polColor, polName = ['b','g'], ['X', 'Y']
ET_text = qa.time('%fs' % (timeStamp[-1]), form='fits', prec=6)[0]
#-------- Prepare Plots
figSPW = plt.figure(figsize=(figInch, figInch))
figSPW.text(0.475, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')), fontsize=fontSize)
figSPW.text(0.05, 0.5, 'Phase [rad]', rotation=90, fontsize=fontSize)
figSPW.text(0.95, 0.5, 'Amplitude', rotation=-90, fontsize=fontSize)
text_timerange = '%s - %s' % (ST_text, ET_text)
figSPW.suptitle('%s SPW=%d %s' % (prefix, spw, text_timerange), fontsize=fontSize)
pMax = np.percentile(abs(scanVis), 98) if 'plotMax' not in locals() else plotMax
aMax = np.percentile(abs(scanAC), 98)
#
for bl_index in list(range(UseBlNum)):
    ants = Bl2Ant(bl_index)
    #-------- Plot cross correlations
    BLamp = figSPW.add_subplot(UseAntNum, UseAntNum, ants[1]*UseAntNum + ants[0] + 1)
    BLphs = figSPW.add_subplot(UseAntNum, UseAntNum, ants[0]*UseAntNum + ants[1] + 1)
    for pol_index in list(range(polNum)):
        plotVis = scanVis[pol_index, bl_index]
        BLamp.step(DT, abs(plotVis), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
        BLphs.plot(DT, np.angle(plotVis), '.', color=polColor[pol_index], label = 'Pol=' + polName[pol_index])
    #
    BLamp.axis([np.min(DT), np.max(DT), 0.0, 1.25*pMax])
    BLamp.xaxis.set_major_locator(plt.NullLocator())
    BLphs.axis([np.min(DT), np.max(DT), -math.pi, math.pi])
    BLphs.tick_params(axis='x', labelsize=int(fontSize*0.25), labelrotation=-90)
    #
    if ants[1] == 0:    # Antenna label in the top and leftside
        BLamp.set_title( antList[antMap[ants[0]]] )
        BLphs.set_ylabel( antList[antMap[ants[0]]] )
    else:
        BLphs.yaxis.set_major_locator(plt.NullLocator())
    #
    if ants[0] == UseAntNum - 1:    # Antenna at rightside
        BLamp.yaxis.tick_right()
    else:
        BLamp.yaxis.set_major_locator(plt.NullLocator())
    if ants[0] < UseAntNum - 1:    # except bottom panel : skip drawing X-axis
        BLphs.xaxis.set_major_locator(plt.NullLocator())
    #
#-------- Plot autocorrelations
for ant_index in list(range(UseAntNum)):
    BLamp = figSPW.add_subplot(UseAntNum, UseAntNum, ant_index*UseAntNum + ant_index + 1)
    BLamp.patch.set_facecolor('lavenderblush')
    BLamp.tick_params(axis='x', labelsize=int(fontSize*0.25), labelrotation=-90)
    for pol_index in list(range(polNum)):
        BLamp.step(DT, abs(scanAC[pol_index, ant_index]), color=polColor[pol_index], where='mid', label = polName[pol_index])
    #
    BLamp.axis([np.min(DT), np.max(DT), 0.0, 1.25*aMax])
    if ant_index < UseAntNum - 1:
        BLamp.xaxis.set_major_locator(plt.NullLocator())
    else:
        BLamp.yaxis.tick_right()
    if ant_index > 0:
        BLamp.yaxis.set_major_locator(plt.NullLocator())
    else:
        BLamp.set_title( antList[antMap[0]])
        BLamp.set_ylabel( antList[antMap[0]])
        BLamp.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
    #
#
plt.show()
pngFile = 'VT_%s_SPW%d' % (prefix, spw)
pdfFile = pngFile + '.pdf'
figSPW.savefig(pdfFile, format='pdf', dpi=144)
plt.close('all')
os.system('pdftoppm -png %s %s' % (pdfFile, pngFile))
