#---- Script for Band-3 Astroholograpy Data
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
exec(open(SCR_DIR + 'interferometry.py').read())
#-------- Definitions
msfile = wd + prefix + '.ms'
antList = GetAntName(msfile)
if 'refant' not in locals(): refant = range(len(antList))
antList = antList[refant]
antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
spwNum  = len(spwList)
if 'chBunch' not in locals(): chBunch = 1
if 'startTime' in locals(): startMJD = qa.convert(startTime, 's')['value']
#
#-------- Procedures
timeRange = np.zeros([2])	# Start and End time period 
blMap = list(range(blNum))
blInv = [False]* blNum		# True -> inverted baseline
for bl_index in list(range(blNum)):
    ants = Bl2Ant(bl_index)
    blMap[bl_index], blInv[bl_index]  = Ant2BlD(refant[ants[0]], refant[ants[1]])
#
#-------- Procedures
interval, timeStamp = GetTimerecord(msfile, 0, 0, spwList[0], BPscan)
integDuration = np.median(interval)
timeNum = len(timeStamp)
if 'integTime' in locals(): timeNum = int(np.ceil(integTime / integDuration))
timeNum = min(timeNum, len(timeStamp))
#
#-------- Prepare Plots
for spw_index in range(spwNum):
    figSPW = plt.figure(spw_index, figsize = (64, 64), dpi=72, frameon=False, tight_layout=True, constrained_layout=True)
    figSPW.text(0.45, 0.05, 'Frequency [GHz]')
    figSPW.text(0.05, 0.5, 'Phase [rad]', rotation=90)
    figSPW.text(0.95, 0.5, 'Amplitude', rotation=-90)
    #
    #-------- Plot BP
    #for spw_index in range(spwNum):
    print(' Loading SPW = %d' % (spwList[spw_index]))
    #figSPW = plt.figure(spw_index)
    chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* Freq  # GHz
    chRange = list(range(int(round(chNum/chBunch * 0.05)), int(round(chNum/chBunch * 0.95))))
    timeStamp, Pspec, Xspec = GetVisAllBL(msfile, spwList[spw_index], BPscan)
    #---- integration timerange
    st_index = 0
    if 'startMJD' in locals(): st_index = np.argmin( abs(timeStamp - startMJD) )
    if st_index + timeNum > len(timeStamp): st_index = len(timeStamp) - timeNum
    timeRange = list(range(st_index, st_index + timeNum))
    text_timerange = qa.time('%fs' % (timeStamp[st_index]), form='fits', prec=6)[0] + ' - ' + qa.time('%fs' % (timeStamp[timeRange[-1]]), form='fits', prec=6)[0]
    print('Integration in %s (%.1f sec)' % (text_timerange, timeNum* integDuration))
    figSPW.suptitle('%s SPW=%d Scan=%d Integration in %s (%.1f sec)' % (prefix, spwList[spw_index], BPscan, text_timerange, (timeNum* integDuration)))
    #---- polarization format
    polNum = Xspec.shape[0]
    if polNum == 4: pol = [0,3]; polName = ['XX', 'YY']         # parallel pol in full-pol correlations
    if polNum == 2: pol = [0,1]; polName = ['XX', 'YY']         # XX and YY
    if polNum == 1: pol = [0]; polName = ['XX']           # Only XX
    polColor = ['b','g']
    polNum = len(pol)
    Xspec = Xspec[pol]; Xspec = Xspec[:,:,blMap][:,:,:,timeRange]
    Pspec = Pspec[pol]; Pspec = Pspec[:,:,:,timeRange]
    tempVis = np.mean(Xspec, axis=3)    # tempVis[pol, ch, bl]
    tempAC  = np.mean(Pspec, axis=3)    # tempVis[pol, ch, ant]
    pMax = np.percentile(abs(tempVis), 98) if 'plotMax' not in locals() else plotMax
    polColor = ['b', 'g']
    for bl_index in list(range(blNum)):
        ants = Bl2Ant(bl_index)
        print('Preparing baseline %s - %s' % (antList[ants[1]], antList[ants[0]]))
        BLampPL = figSPW.add_subplot( antNum, antNum, antNum* ants[1] + ants[0] + 1 )
        BLphsPL = figSPW.add_subplot( antNum, antNum, antNum* ants[0] + ants[1] + 1 )
        for pol_index in list(range(polNum)):
            plotVis = tempVis[pol_index, :, bl_index]
            BLampPL.step( Freq, abs(plotVis), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
            BLphsPL.plot( Freq, np.angle(plotVis), '.', color=polColor[pol_index], label = 'Pol=' + polName[pol_index])
        #
        BLampPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25* pMax])
        BLampPL.xaxis.set_major_locator(plt.NullLocator())
        BLphsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
        if ants[1] == 0:    # Antenna label in the top and leftside
            BLampPL.set_title( antList[ants[0]] )
            BLphsPL.set_ylabel( antList[ants[0]] )
        #
        if ants[0] < antNum - 1: BLphsPL.xaxis.set_major_locator(plt.NullLocator())
        if ants[1] > 0: BLphsPL.yaxis.set_major_locator(plt.NullLocator())
        if ants[0] < antNum - 1: BLampPL.yaxis.set_major_locator(plt.NullLocator())
        #
        BLphsPL.tick_params(axis='x', which='major', labelsize=5)   # Frequency axis label
        BLampPL.tick_params(axis='y', which='major', labelsize=5, labelright=True, labelleft=False, right=True, left=False) # Amplitude axis label
    #
    for ant_index in list(range(antNum)):
        ACampPL = figSPW.add_subplot( antNum, antNum, (antNum + 1)* ant_index + 1 )
        ACampPL.patch.set_facecolor('pink')
        for pol_index in list(range(polNum)):
            ACampPL.step( Freq, abs(tempAC[pol_index, :, ant_index]), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
            ACampPL.xaxis.set_major_locator(plt.NullLocator())
            ACampPL.yaxis.set_major_locator(plt.NullLocator())
        #
        if ant_index == 0:
            ACampPL.set_title( antList[0] )
            ACampPL.set_ylabel( antList[0] )
            ACampPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
        #
    #
    figSPW.savefig('PS_%s_Scan%d_SPW%d.png' % (prefix, BPscan, spwList[spw_index]))
#
