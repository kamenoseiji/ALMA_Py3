####
# Script for Bispectra (closure phase)
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
#-------- Array Configuration
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
    if 'BP_ant' in locals():
        print('Applying bandpass calibration...')
        BPCaledXspec = (tempSpec / (BP_ant[ant0]* BP_ant[ant1].conjugate())).transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    else:
        BPCaledXspec = tempSpec.transpose(2,3,1,0) # Bandpass Cal ; BPCaledXspec[pol, ch, bl, time]
    #
    #-------- Antenna-based Gain correction
    chAvgVis = np.mean(BPCaledXspec[:, chRange], axis=1)
    if 'timeBunch' in locals():
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
Pclamp = recursiveP(UseAntNum)
CLnum = Pclamp.shape[0]
CLList = P2CL(Pclamp)
#-------- Prepare Plots
figSPW = plt.figure(figsize=(figInch, figInch))
figSPW.text(0.475, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')), fontsize=fontSize)
figSPW.text(0.05, 0.5, 'Closure Amplitude', rotation=90, fontsize=fontSize)
text_timerange = '%s - %s' % (ST_text, ET_text)
figSPW.suptitle('%s SPW=%d %s' % (prefix, spw, text_timerange), fontsize=fontSize)
pMax = 2.0 if 'plotMax' not in locals() else plotMax
for cl_index in list(range(CLnum)):
    CLampAntList = CLList[cl_index]
    text_sd = '[%s-%s][%s-%s] / [%s-%s][%s-%s]' % (antList[antMap[CLampAntList[0]]], antList[antMap[CLampAntList[1]]], antList[antMap[CLampAntList[2]]], antList[antMap[CLampAntList[3]]],
         antList[antMap[CLampAntList[0]]], antList[antMap[CLampAntList[2]]], antList[antMap[CLampAntList[1]]], antList[antMap[CLampAntList[3]]])
    print(text_sd)
    if UseAntNum % 2 == 0:
        NumPlotY, NumPlotX = int(UseAntNum/2), int(UseAntNum -3)
    else:
        NumPlotY, NumPlotX = int(UseAntNum), int((UseAntNum -3)/2)
    #
    BLamp = figSPW.add_subplot(NumPlotY, NumPlotX, cl_index + 1)
    for pol_index in list(range(polNum)):
        Numer_bl = np.where(Pclamp[cl_index] == 1)[0].tolist()
        Denom_bl = np.where(Pclamp[cl_index] ==-1)[0].tolist()
        plotVis = abs(scanVis[pol_index][Numer_bl[0]])* abs(scanVis[pol_index][Numer_bl[1]]) / (abs(scanVis[pol_index][Denom_bl[0]])* abs(scanVis[pol_index][Denom_bl[1]]))
        BLamp.step(DT, abs(plotVis), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
    #
    BLamp.axis([np.min(DT), np.max(DT), 0.0, pMax])
    BLamp.tick_params(axis='x', labelsize=int(fontSize*0.5), labelrotation=-90)
    BLamp.set_title(text_sd, fontsize=0.5*fontSize)
    if cl_index < CLnum - NumPlotX: BLamp.xaxis.set_major_locator(plt.NullLocator())
#
plt.show()
pngFile = 'CA_%s_SPW%d' % (prefix, spw)
pdfFile = pngFile + '.pdf'
figSPW.savefig(pdfFile, format='pdf', dpi=144)
plt.close('all')
os.system('pdftoppm -png %s %s' % (pdfFile, pngFile))
