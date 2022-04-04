####
# Script for Bispectra (closure phase)
####
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
exec(open(SCR_DIR + 'interferometry.py').read())
#-------- Definitions
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
blMap = list(range(blNum))
spwNum  = len(spwList)
figInch = max(16,antNum)
fontSize = min(32, figInch)
if 'startTime' in locals(): startMJD = qa.convert(startTime, 's')['value']
if 'BPscan' not in locals(): BPscan = 3 # default bandpass scan
#
#-------- Baseline mapping
#-------- Loop for SPW
for spw_index in range(spwNum):
    interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], BPscan)
    integDuration = np.median(interval)
    timeNum = len(timeStamp)
    if 'integTime' in locals(): timeNum = int(np.ceil(integTime / integDuration))
    timeNum = min(timeNum, len(timeStamp))
    DT = []
    for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(qa.time('%fs' % (mjdSec), form='fits', prec=9)[0], '%Y-%m-%dT%H:%M:%S.%f'))
    #-------- Prepare Plots
    figSPW = plt.figure(figsize=(figInch, figInch))
    figSPW.text(0.475, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')), fontsize=fontSize)
    figSPW.text(0.05, 0.5, 'Closure Phase [rad]', rotation=90, fontsize=fontSize)
    figSPW.text(0.95, 0.5, 'Amplitude', rotation=-90, fontsize=fontSize)
    #
    #-------- Plot VisAmp/Phs
    print(' Loading SPW = %d' % (spwList[spw_index]))
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPscan)
    #---- integration timerange
    st_index = 0
    if 'startMJD' in locals(): st_index = np.argmin( abs(timeStamp - startMJD) )
    if st_index + timeNum > len(timeStamp): st_index = len(timeStamp) - timeNum
    timeRange = list(range(st_index, st_index + timeNum))
    text_timerange = qa.time('%fs' % (timeStamp[st_index]), form='fits', prec=6)[0] + ' - ' + qa.time('%fs' % (timeStamp[timeRange[-1]]), form='fits', prec=6)[0]
    print('Integration in %s (%.1f sec)' % (text_timerange, timeNum* integDuration))
    figSPW.suptitle('%s SPW=%d Scan=%d %s (%.1f sec)' % (prefix, spwList[spw_index], BPscan, text_timerange, (timeNum* integDuration)), fontsize=fontSize)
    #---- polarization format
    polNum = Xspec.shape[0]
    if polNum == 4: pol = [0,3]; polName = ['X', 'Y']         # parallel pol in full-pol correlations
    if polNum == 2: pol = [0,1]; polName = ['X', 'Y']         # XX and YY
    if polNum == 1: pol = [0]; polName = ['X']           # Only XX
    polColor = ['b','g']
    polNum = len(pol)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap][:,:,:,timeRange]
    Pspec = Pspec[pol]; Pspec = Pspec[:,:,:,timeRange]
    if chNum > 1:
        chRange = list(range(int(chNum* 0.05), int(chNum* 0.95)))
        tempVis = np.mean(Xspec[:,chRange], axis=1)    # tempVis[pol, bl, time]
        tempAC  = np.mean(Pspec[:,chRange], axis=1)    # tempVis[pol, ant, time]
    else:
        tempVis = Xspec[:,0]    # tempVis[pol, bl, time]
        tempAC  = Pspec[:,0]    # tempVis[pol, ant, time]
    #
    pMax = np.percentile(abs(tempVis), 98) if 'plotMax' not in locals() else plotMax
    polColor = ['b', 'g']
    for bl_index in list(range(blNum)):
        ants = Bl2Ant(bl_index)
        #-------- Plot visibility amplitude
        BLamp = figSPW.add_subplot(antNum-1, antNum-1, ants[1]*(antNum -1) + ants[0])
        for pol_index in list(range(polNum)):
            plotVis = tempVis[pol_index, bl_index]
            BLamp.step(DT, abs(plotVis), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
        #
        BLamp.axis([np.min(DT), np.max(DT), 0.0, 1.25*pMax])
        BLamp.xaxis.set_major_locator(plt.NullLocator())
        if bl_index == 0: BLamp.set_ylabel(antList[0])
        if ants[1] == 0:    # Antenna label in the top and leftside
            BLamp.set_title( antList[ants[0]] )
        if ants[0] == antNum - 1:    # Antenna at rightside
            BLamp.yaxis.tick_right()
        else:
            BLamp.yaxis.set_major_locator(plt.NullLocator())
        #-------- Plot closure phase
        if ants[1] > 0:         # plot Closure phase
            BLphs = figSPW.add_subplot(antNum-1, antNum-1, (ants[0] - 1)*(antNum-1) + ants[1])
            tri0, tri1, tri2 = Ant2Bl(ants[0], 0), Ant2Bl(ants[1], 0), Ant2Bl(ants[0], ants[1])
            for pol_index in list(range(polNum)):
                plotVis = tempVis[pol_index, tri0].conjugate()* tempVis[pol_index, tri1]* tempVis[pol_index, tri2]
                BLphs.plot(DT, np.angle(plotVis), '.', color=polColor[pol_index], label = 'Pol=' + polName[pol_index])
            #
            print('%d : %d - %d - %d (ant %s, %s, %s)' % (bl_index, tri0, tri1, tri2, antList[0], antList[ants[1]], antList[ants[0]]))
            BLphs.set_title('%s-%s-%s' % (antList[0], antList[ants[1]], antList[ants[0]] ))
            BLphs.axis([np.min(DT), np.max(DT), -math.pi, math.pi])
            BLphs.tick_params(axis='x', labelsize=int(fontSize*0.125), labelrotation=-90)
            if ants[1] == 1: BLphs.set_ylabel(antList[ants[0]] )
            if ants[1] > 1: BLphs.yaxis.set_major_locator(plt.NullLocator())
            if ants[0] < antNum - 1:    # except bottom panel : skip drawing X-axis
                BLphs.xaxis.set_major_locator(plt.NullLocator())
            #
        #
    #
    plt.show()
    pngFile = 'BS_%s_Scan%d_SPW%d' % (prefix, BPscan, spwList[spw_index])
    pdfFile = pngFile + '.pdf'
    figSPW.savefig(pdfFile, format='pdf', dpi=144)
    plt.close('all')
    os.system('pdftoppm -png %s %s' % (pdfFile, pngFile))
#
