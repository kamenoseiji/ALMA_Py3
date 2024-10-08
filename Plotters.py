import math
import analysisUtils as au
import numpy as np
import scipy
import datetime
from interferometry import GetChNum, bunchVec, delay_search, Bl2Ant, Ant2Bl
from Grid import tauSMTH
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
polColor = ['b', 'g', 'm', 'y']
polName = ['X', 'Y', 'XY', 'YX']
#
#-------- Set Color Map
lineCmap = plt.get_cmap('Set1')
#-------- Plot optical depth
def plotTauSpec(prefix, spwList, freqList, Tau0spec):
    figTauSP = plt.figure(0, figsize = (11,8))
    figTauSP.suptitle(prefix + ' Zenith Opacity')
    figTauSP.text(0.45, 0.05, 'Frequency [GHz]')
    figTauSP.text(0.03, 0.45, 'Optical Depth', rotation=90)
    spwNum = len(spwList)
    plotMax = 0.01
    for spw_index in range(spwNum): plotMax = max(plotMax, np.max(Tau0spec[spw_index]))
    for spw_index in range(spwNum):
        chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
        TauPL = figTauSP.add_subplot(1, spwNum, spw_index + 1 )
        TauPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, 1.05* plotMax])
        TauPL.tick_params(axis='both', labelsize=6)
        TauPL.step(freqList[spw_index][chRange], Tau0spec[spw_index][chRange], color='k', where='mid')
        text_sd = 'SPW = %d' % (spwList[spw_index])
        TauPL.text(np.min(freqList[spw_index]), 1.01* plotMax, text_sd, fontsize='8')
    #
    figTauSP.savefig('TAUS_' + prefix + '.pdf')
    plt.close('all')
    return
#
#-------- Plot Tau fit
def plotTauFit(prefix, antList, spwList, secZ, tempAtm, Tau0, TantN, TskyList, scanFlag):
    antNum = len(antList)
    spwNum = len(spwList)
    airmass = np.arange( 1.0, 1.25*np.max(secZ), 0.01)
    figTauFit = plt.figure(0, figsize = (11,8))
    figTauFit.suptitle(prefix + ' Optical Depth')
    figTauFit.text(0.45, 0.05, 'Airmass')
    figTauFit.text(0.03, 0.45, 'Sky Temperature [K]', rotation=90)
    for spw_index in range(spwNum):
        chAvgTsky = np.median(TskyList[spw_index].transpose(2, 1, 0) - TantN[spw_index], axis=2).T  # chAvgTsky[ant, scan]
        chAvgTau0 = np.median(Tau0[spw_index])
        plotMax = 1.2 * np.max(chAvgTsky)
        TskyPL = figTauFit.add_subplot(1, spwNum, spw_index + 1 )
        TskyPL.axis([1.0, 2.5, 0.0, plotMax])
        TskyPL.plot( airmass, au.Tcmb* np.exp(-chAvgTau0* airmass) + tempAtm* (1.0 - np.exp(-chAvgTau0* airmass)), '-', color='k', alpha=0.5)
        for ant_index in range(antNum):
            rgb = lineCmap(float(ant_index) / antNum )
            plotTsky = chAvgTsky[ant_index]
            TskyPL.scatter( secZ, plotTsky, s=10.0* scanFlag[spw_index, ant_index], color=rgb, label = antList[ant_index])
        #
        text_sd = 'Tau(zenith)=%6.4f' % (chAvgTau0)
        TskyPL.text(1.01, 0.95* plotMax, text_sd, fontsize='9')
        TskyPL.set_title('SPW %d' % (spwList[spw_index]))
    #
    TskyPL.legend(loc = 'lower right', prop={'size' :7}, numpoints = 1)
    figTauFit.savefig('TAUF_' + prefix + '.pdf')
    plt.close('all')
    return
#
#-------- Plot plotTauEOn
def plotTauEOn(prefix, bandName, spw, onTime, onData, TauEOn):
    DT = []
    for mjdSec in onTime.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    figTau = plt.figure(0, figsize = (11,8))
    figTau.suptitle('%s %s SPW %d CHAV data and zenith optical depth' % (prefix, bandName, spw))
    figTau.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')));
    CHAVPL = figTau.add_subplot(2, 1, 1)
    CHAVPL.plot(DT, onData, 'k.')
    CHAVPL.set_ylabel('Channel-Averaded Power')
    TAU0PL = figTau.add_subplot(2, 1, 2)
    TAU0PL.plot(DT, TauEOn, 'b.')
    TAU0PL.set_ylabel('Estimated Zenith Optical Depth')
    figTau.savefig('TauEOn_%s-%s-SPW%d.pdf' % (prefix, bandName, spw))
    plt.close('all')
    return
#-------- Plot Tau0-Excess
def plotTau0E(prefix, atmTime, spwList, Tau0, Tau0Excess, scanFlag):
    spwNum = len(spwList)
    figTauE = plt.figure(0, figsize = (11,8))
    figTauE.suptitle(prefix + ' Zenith Optical Depth')
    DT, DTSpl = [], []
    for mjdSec in atmTime.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    mjdSpl = np.arange(atmTime[0], atmTime[-1], 1)
    for mjdSec in mjdSpl.tolist(): DTSpl.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    figTauE.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')));
    figTauE.text(0.03, 0.45, 'Zenith Optical Depth', rotation=90)
    for spw_index in range(spwNum):
        Tau0E = np.median(Tau0[spw_index]) + Tau0Excess[spw_index]
        SP = tauSMTH( atmTime-atmTime[0], Tau0E )
        Tau0ESpl = scipy.interpolate.splev(mjdSpl-atmTime[0], SP)
        TauEPL = figTauE.add_subplot(1, spwNum, spw_index + 1 )
        TauEPL.plot( DTSpl, Tau0ESpl, '-', color='k')
        TauEPL.scatter( DT, Tau0E, s=10.0* scanFlag[spw_index], color='c')
        TauEPL.tick_params(axis='x', labelsize=6)
        TauEPL.set_title('SPW %d' % (spwList[spw_index]))
    figTauE.savefig('TAUE_' + prefix + '.pdf')
    plt.close('all')
    return
#-------- Plot Tsys and Trx Spectrum
def plotTsys(prefix, antList, spwList, freqList, atmTime, TrxList, TskyList):
    pp = PdfPages('TSYS_' + prefix + '.pdf')
    #-------- Plots for Tsys spectra
    antNum, spwNum, scanNum, polNum  = len(antList), len(spwList), len(atmTime), TrxList[0].shape[0]
    #-------- Prepare Plots
    figAnt = plt.figure(figsize = (8, 11))
    figAnt.suptitle(prefix)
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.42, 'Tsys (X:blue, Y:green) and Trec (X:magenta, Y:yellow) [K]', rotation=90)
    #-------- Plot BP
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in TsysPL: figAnt.delaxes(PL)
        #
        TsysPL, TsysMax = [], []
        for spw_index in range(spwNum):
            TsysMax = TsysMax + [1.7* np.percentile(TrxList[spw_index], 75, axis=(0,1,3))[ant_index] + 1.5* np.median(TskyList[spw_index])]
            #TsysMax = TsysMax + [1.7* np.percentile(TrxList[spw_index], 75, axis=(0,1,3))[ant_index] + np.max(TskyList[spw_index])]
        plotMax = max(TsysMax)
        for spw_index in range(spwNum):
            chNum = len(freqList[spw_index]); chRange = range(int(0.05*chNum), int(0.95*chNum))
            for scan_index in range(scanNum):
                currentPL = figAnt.add_subplot(scanNum, spwNum, spwNum* scan_index + spw_index + 1 )
                TsysPL = TsysPL + [currentPL]
                timeLabel = au.call_qa_time('%fs' % (atmTime[scan_index]), form='fits')
                for pol_index in range(polNum):
                    plotTrx  = TrxList[spw_index][pol_index, chRange, ant_index, scan_index]
                    plotTsys = TskyList[spw_index][chRange, ant_index, scan_index] + plotTrx
                    currentPL.step( freqList[spw_index][chRange], plotTsys, where='mid', color=polColor[pol_index], label = 'Tsys Pol '+ polName[pol_index])
                    currentPL.step( freqList[spw_index][chRange], plotTrx,  color=polColor[pol_index+2], where='mid', label = 'Trec Pol ' + polName[pol_index])
                #
                currentPL.axis([np.min(freqList[spw_index]), np.max(freqList[spw_index]), 0.0, plotMax])
                currentPL.tick_params(axis='both', labelsize=6)
                if scan_index == 0: currentPL.set_title('SPW %d' % (spwList[spw_index]))
                if scan_index == 0 and spw_index == 0: currentPL.text(1.2* np.min(freqList[0]) - 0.2* np.max(freqList[0]), (0.9 + 0.075*scanNum)*plotMax, antList[ant_index], fontsize='16')
                if scan_index < scanNum - 1: currentPL.set_xticklabels([])
                if spw_index == 0: currentPL.text(np.min(freqList[spw_index]), 0.8* plotMax, timeLabel, fontsize='8')
                else: currentPL.set_yticklabels([])
                #
            #
        #
        figAnt.savefig(pp, format='pdf')
        #
        #
    #
    for PL in TsysPL: figAnt.delaxes(PL)
    plt.close('all')
    pp.close()
    del(TsysPL)
    del(figAnt)
    return
#
#-------- Plot baseline-AllanVariance map
def plotBLAV(prefix, antList, spw, AV_bl):
    pp = PdfPages('AV_%s-SPW%d.pdf' % (prefix, spw))
    figAV = plt.figure(figsize = (8, 8))
    figAV.suptitle('Allan Variance %s SPW%d' % (prefix, spw))
    axAV = figAV.add_subplot(1, 1, 1)
    antNum, blNum = len(antList), AV_bl.shape[1]
    AVmap = np.zeros([2, antNum, antNum])
    for bl_index in list(range(blNum)):
        ants = Bl2Ant(bl_index)
        AVmap[0, ants[1], ants[0]] = AV_bl[0, bl_index]  # pol-X
        AVmap[1, ants[0], ants[1]] = AV_bl[1, bl_index]  # pol-Y
    #
    AVmap = np.sum(AVmap, axis=0)
    for ant_index, ant in enumerate(antList): AVmap[ant_index, ant_index] = 0.0
    imAV = axAV.imshow(AVmap, cmap='coolwarm', vmin=0.0, vmax=0.01)
    figAV.colorbar(imAV, ax=axAV)
    axAV.plot([0,antNum-1], [0,antNum-1], color='w')
    axAV.set_xticks(list(range(antNum))); axAV.set_yticks( list(range(antNum)))
    axAV.set_xticklabels(antList.tolist()); axAV.set_yticklabels(antList.tolist())
    figAV.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot autocorrelation power spectra
def plotAC(prefix, antList, spwList, freqList, AC):
    pp = PdfPages('AC_' + prefix + '.pdf')
    antNum, spwNum, polNum = len(antList), len(spwList), AC[0].shape[2]
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Power Spectra')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.5, 'Median amplitude and variation [dB]', rotation=90)
    #-------- Plot AC
    for ant_index in range(antNum):
        if ant_index > 0:
            for PL in ACList: figAnt.delaxes(PL)
            for PL in SDList: figAnt.delaxes(PL)
        #
        ACList, SDList = [], []
        ACMAX = 0.0
        for spw_index in range(spwNum): ACMAX = max(ACMAX, np.max(AC[spw_index][ant_index]))
        for spw_index in range(spwNum):
            Freq = freqList[spw_index]
            for pol_index in range(polNum):
                ACPL = figAnt.add_subplot(4, spwNum, spwNum* pol_index + spw_index + 1)
                SDPL = figAnt.add_subplot(4, spwNum, spwNum* (2+pol_index) + spw_index + 1)
                ACList = ACList + [ACPL]; SDList = SDList + [SDPL]
                avgAC = np.mean(AC[spw_index][ant_index, :, pol_index], axis=0)
                plotAC = 10.0* np.log10(avgAC )
                chNum = len(plotAC)
                chRange = range( int(chNum*0.06), int(chNum* 0.96))
                maxAC, minAC, maxFreq = np.max(plotAC), np.min(plotAC), Freq[np.argmax(plotAC)]
                text_sd = 'Peak = %.1f dB at %.2f GHz' % (maxAC, maxFreq)
                #plotMax, plotMin = 10.0* np.log10(ACMAX), 10.0* np.log10(ACMAX) - 20.0
                #plotMax, plotMin = 10.0* np.log10(ACMAX), 10.0* np.log10(minAC) - 1.0
                plotMax, plotMin = 10.0* np.log10(ACMAX), 118.0
                if spw_index == 0 and pol_index == 0: ACPL.text(1.2* np.min(Freq) - 0.2* np.max(Freq), 1.3*plotMax - 0.3* plotMin, antList[ant_index], fontsize=10)
                ACPL.axis([np.min(Freq), np.max(Freq), plotMin, 1.1* plotMax - 0.1* plotMin])
                ACPL.tick_params(axis='both', labelsize=6)
                ACPL.set_xticklabels([])
                ACPL.text( np.min(Freq), 1.2*plotMax - 0.2*plotMin, 'AC SPW=%d Pol-%s' % (spwList[spw_index], polName[pol_index]), fontsize=7)
                ACPL.text( np.min(Freq), 1.12*plotMax - 0.12*plotMin, text_sd, fontsize=7)
                ACPL.step(Freq, plotAC, where='mid')
                #
                plotSD = 10.0* np.log10(np.std(AC[spw_index][ant_index, :, pol_index]/avgAC, axis=0))
                maxSD, minSD, maxFreq = np.max(plotSD[chRange]), np.min(plotSD[chRange]), Freq[np.argmax(plotSD[chRange])]
                text_sd = 'Peak = %.1f dB at %.2f GHz' % (maxSD, maxFreq)
                bgcolor = 'green'
                if maxSD > -30.0: bgcolor = 'orange'
                if maxSD > -20.0: bgcolor = 'red'
                #plotMax, plotMin = max(-20.0, maxSD), min(-40.0, minSD)
                plotMax, plotMin = max(-20.0, maxSD), min(-20.0, minSD)
                SDPL.axis([np.min(Freq), np.max(Freq), 1.05*plotMin - 0.05*plotMax, 1.05*plotMax - 0.05*plotMin])
                SDPL.axhspan(ymin=-30.0, ymax=plotMax, color=bgcolor, alpha=0.1) 
                SDPL.tick_params(axis='both', labelsize=6)
                SDPL.get_xaxis().get_major_formatter().set_useOffset(False)
                if pol_index == 0: SDPL.set_xticklabels([])
                SDPL.text( np.min(Freq), 1.16* plotMax - 0.16*plotMin, 'SD SPW=%d Pol-%s' % (spwList[spw_index], polName[pol_index]), fontsize=7)
                SDPL.text( np.min(Freq), 1.08* plotMax - 0.08*plotMin, text_sd, fontsize=7)
                SDPL.step(Freq, plotSD, where='mid')
            #
        #
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all')
    pp.close()
    del(ACList)
    del(SDList)
    del(ACPL)
    del(SDPL)
    return
#
#-------- Plot Cross power spectrum
def plotSP(pp, prefix, antList, spwList, freqList, BPList, plotMin=0.0, plotMax=1.5, delayMessage=False):
    antNum, spwNum = len(antList), len(spwList)
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Bandpass')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #-------- Plot BP
    text_delay = 'Ant  '
    for spw_index, spw in enumerate(spwList):
        text_delay = text_delay + ' |SPW%02d'%spw
        text_delay = text_delay + '  ' + polName[0]
        if BPList[spw_index].shape[1] == 2: text_delay = text_delay + '         ' + polName[1]
    if delayMessage: print(text_delay)
    text_delay = '------+'
    for spw_index, spw in enumerate(spwList): text_delay = text_delay + '--- delay [ns] ----+'
    if delayMessage: print(text_delay)
    for ant_index, antName in enumerate(antList):
        if ant_index > 0:
            for PL in AmpList: figAnt.delaxes(PL)
            for PL in PhsList: figAnt.delaxes(PL)
        #
        AmpList, PhsList = [], []
        text_delay = antName + '  |'
        for spw_index, spw in enumerate(spwList):
            Freq = freqList[spw_index]
            plotFreq = 1.0e-9 * Freq
            chNum = len(Freq); chRange = list(range(int(0.10*chNum), int(0.95*chNum)))
            BW = Freq[chRange[-1]] - Freq[chRange[0]]
            AmpPL = figAnt.add_subplot(2, spwNum, spw_index + 1 )
            PhsPL = figAnt.add_subplot(2, spwNum, spwNum + spw_index + 1 )
            AmpList = AmpList + [AmpPL]
            PhsList = PhsList + [PhsPL]
            ppolNum = BPList[spw_index].shape[1]
            for pol_index in list(range(ppolNum)):
                plotBandpass = BPList[spw_index][ant_index,pol_index]
                delayAnt, delaySNR = delay_search(plotBandpass[chRange])
                # text_sd = '%s %s : %+f SNR=%.1f' % (antName, polName[pol_index], 0.5e9* delayAnt/ BW, delaySNR)
                AmpPL.fill_between(plotFreq[chRange], plotMin, plotMax, color='yellow', alpha=0.1)
                AmpPL.step(plotFreq, abs(plotBandpass), color=polColor[pol_index], where='mid', label = 'Pol-%s' % (polName[pol_index]))
                text_delay = text_delay + '%+f ' % (0.5e9* delayAnt/ BW)
                PhsPL.fill_between(plotFreq[chRange], -np.pi, np.pi, color='yellow', alpha=0.1)
                PhsPL.plot(plotFreq, np.angle(plotBandpass), '.', color=polColor[pol_index], label = 'Pol-%s %+f [ns]' % (polName[pol_index], 0.5e9* delayAnt/ BW))
            #
            if spw_index == 0: AmpPL.set_title(antList[ant_index])
            AmpPL.axis([np.min(plotFreq), np.max(plotFreq), plotMin, plotMax])
            AmpPL.tick_params(axis='both', labelsize=6)
            AmpPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            AmpPL.text( np.min(plotFreq), 0.9* plotMax + 0.1* plotMin, 'SPW=%d Amp' % (spwList[spw_index]))
            PhsPL.axis([np.min(plotFreq), np.max(plotFreq), -np.pi, np.pi])
            PhsPL.tick_params(axis='both', labelsize=6)
            PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            PhsPL.text( np.min(plotFreq), 2.5, 'SPW=%s Phs' % (str(spw)))
        #
        if delayMessage: print(text_delay)
        figAnt.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot Bandpass
def plotBP(pp, prefix, antList, spwList, BPscan, BPList, bunchNum=1, plotMax=1.2, plotMarker=[[]]):
    msfile = prefix + '.ms'
    antNum, spwNum = len(antList), len(spwList)
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' Scan %d' % (BPscan))
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'Bandpass Amplitude and Phase', rotation=90)
    #-------- Plot BP
    for ant_index in list(range(antNum)):
        if ant_index > 0:
            for PL in AmpList: figAnt.delaxes(PL)
            for PL in PhsList: figAnt.delaxes(PL)
        #
        AmpList, PhsList = [], []
        for spw_index in list(range(spwNum)):
            chNum, chWid, Freq = GetChNum(msfile, spwList[spw_index]); Freq = 1.0e-9* bunchVec(Freq, bunchNum)  # GHz
            AmpPL = figAnt.add_subplot(2, spwNum, spw_index + 1 )
            PhsPL = figAnt.add_subplot(2, spwNum, spwNum + spw_index + 1 )
            AmpList = AmpList + [AmpPL]
            PhsList = PhsList + [PhsPL]
            for pol_index in list(range(BPList[spw_index].shape[1])):
                plotBandpass = BPList[spw_index][ant_index,pol_index]
                '''
                delayAnt, delaySNR = delay_search(plotBandpass)
                text_sd = '%s %s %+f %.1f' % (antList[ant_index], polName[pol_index], 0.5*delayAnt / (chNum*chWid), delaySNR)
                print(text_sd)
                '''
                AmpPL.step(Freq, abs(plotBandpass), color=polColor[pol_index], where='mid', label = 'Pol-' + polName[pol_index])
                PhsPL.plot( Freq, np.angle(plotBandpass), '.', color=polColor[pol_index], label = 'Pol-' + polName[pol_index])
            #
            if len(plotMarker[0]) > 0: 
                for spurIndex in list(range(len(plotMarker[spw_index]))): AmpPL.vlines(x=1.0e-9 * plotMarker[spw_index][spurIndex], ymin=0.0, ymax=1.25*plotMax, color='gray') 
            if spw_index == 0: AmpPL.set_title(antList[ant_index])
            AmpPL.axis([np.min(Freq), np.max(Freq), 0.0, 1.25*plotMax])
            AmpPL.tick_params(axis='both', labelsize=6)
            AmpPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            AmpPL.text( np.min(Freq), 1.1* plotMax, 'SPW=%d Amp' % (spwList[spw_index]))
            PhsPL.axis([np.min(Freq), np.max(Freq), -math.pi, math.pi])
            PhsPL.tick_params(axis='both', labelsize=6)
            PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
            PhsPL.text( np.min(Freq), 2.5, 'SPW=%d Phs' % (spwList[spw_index]))
        #
        #plt.show()
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all')
    pp.close()
    del(AmpList)
    del(PhsList)
    del(AmpPL)
    del(PhsPL)
    del(figAnt)
    return
#
#-------- Plot Gain
def plotGain(prefix, spw):
    #-------- Load tables
    antFile = prefix + '.Ant.npy'; antList = np.load(antFile); antNum = len(antList)
    timeFile = '%s-SPW%d.TS.npy' % (prefix, spw)
    GainFile = '%s-SPW%d.GA.npy' % (prefix, spw)    # Gain[ant, pol, time]
    FlagFile = '%s-SPW%d.FG.npy' % (prefix, spw)    # Flag[ant, time] 
    pp = PdfPages('GA_%s-SPW%d.pdf' %  (prefix, spw))
    DT, timeStamp, Gain, FG = [], np.load(timeFile), np.load(GainFile), np.load(FlagFile)
    polNum = Gain.shape[1]
    for mjdSec in timeStamp.tolist(): DT.append(datetime.datetime.strptime(au.call_qa_time('%fs' % (mjdSec), form='fits', prec=9), '%Y-%m-%dT%H:%M:%S.%f'))
    #-------- Prepare Plots
    figAmp, figPhs = plt.figure(figsize = (8, 11)), plt.figure(figsize = (8, 11))
    figAmp.suptitle(GainFile + ' Gain Amplitude'); figPhs.suptitle(GainFile + ' Gain Phase')
    figAmp.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d'))); figPhs.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')))
    figAmp.text(0.03, 0.45, 'Gain Amplitude = sqrt(correlated flux / SEFD)', rotation=90); figPhs.text(0.03, 0.45, 'Gain Phase [deg]', rotation=90)
    plotMin, plotMax = 0.0, 1.1* np.percentile(np.max(abs(Gain), axis=1)* FG, 90)
    #-------- Plot Gain
    gainXmed, gainYmed = [], []
    for ant_index in list(range(antNum)):
        ampMed = [0.0, 0.0]
        flag_index = np.where(FG[ant_index] > 0.01)[0].tolist()
        AmpPL = figAmp.add_subplot( int(np.ceil(antNum/2.0)), 2, ant_index + 1 )
        PhsPL = figPhs.add_subplot( int(np.ceil(antNum/2.0)), 2, ant_index + 1 )
        for pol_index in list(range(polNum)): AmpPL.plot( np.array(DT)[flag_index], abs(Gain[ant_index, pol_index, flag_index]), '.', markersize=3, color=polColor[pol_index])
        for pol_index in list(range(polNum)): PhsPL.plot( np.array(DT)[flag_index], np.angle(Gain[ant_index, pol_index, flag_index])*180.0/math.pi, '.', markersize=3, color=polColor[pol_index])
        if len(flag_index) > 0: ampMed = np.median(abs(Gain[ant_index][:, flag_index]), axis=1)
        gainXmed, gainYmed = gainXmed + [ampMed[0]], gainYmed + [ampMed[1]]
        AmpPL.yaxis.set_major_formatter(ptick.ScalarFormatter(useMathText=True))
        AmpPL.yaxis.offsetText.set_fontsize(3)
        PhsPL.yaxis.offsetText.set_fontsize(3)
        AmpPL.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
        AmpPL.tick_params(labelsize=4)
        PhsPL.tick_params(labelsize=4)
        AmpPL.axis([np.min(DT), np.max(DT), plotMin, plotMax])
        PhsPL.axis([np.min(DT), np.max(DT), -180.0, 180.0])
        if polNum == 2: text_sd = '%s : Gain(median) = (%.2f%% %.2f%%)' % (antList[ant_index], 100.0* ampMed[0], 100.0* ampMed[1])
        else: text_sd = '%s : Gain(median) = (%.2f%%)' % (antList[ant_index], 100.0* ampMed[0])
        print(text_sd)
        AmpPL.text( 0.05, 1.02, text_sd, transform=AmpPL.transAxes, fontsize=5)
        PhsPL.text( 0.05, 1.02, antList[ant_index], transform=PhsPL.transAxes, fontsize=5)
    #
    text_sd = 'all: Gain(median) = (%.2f%% %.2f%%)' % (100.0* np.median(gainXmed), 100.0* np.median(gainYmed)); print(text_sd)
    figAmp.savefig(pp, format='pdf'); figPhs.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot visibility amplitude and closure phase
def plotBispec(antList, scanVis, DT, plotFile, labelList, pMax):
    antNum = len(antList)
    polNum, blNum  = scanVis.shape[0], scanVis.shape[1]
    figInch  = max(16,antNum)
    fontSize = min(32, figInch)
    figSPW = plt.figure(figsize=(figInch, figInch))
    figSPW.text(0.475, 0.05, labelList[0], fontsize=fontSize)
    figSPW.text(0.05, 0.5, labelList[1], rotation=90, fontsize=fontSize)
    figSPW.text(0.95, 0.5, labelList[2], rotation=-90, fontsize=fontSize)
    figSPW.suptitle(labelList[3], fontsize=fontSize)
    #
    for bl_index in list(range(blNum)):
        ants = Bl2Ant(bl_index)
        #-------- Plot visibility amplitude
        BLamp = figSPW.add_subplot(antNum-1, antNum-1, ants[1]*(antNum -1) + ants[0])
        for pol_index in list(range(polNum)):
            plotVis = scanVis[pol_index, bl_index]
            BLamp.step(DT, abs(plotVis), color=polColor[pol_index], where='mid', label = 'Pol=' + polName[pol_index])
        #
        BLamp.axis([np.min(DT), np.max(DT), 0.0, 1.25*pMax])
        BLamp.xaxis.set_major_locator(plt.NullLocator())
        if bl_index == 0:
            BLamp.legend(loc = 'best', prop={'size' :7}, numpoints = 1)
        if ants[0] - ants[1] == 1:
            BLamp.set_ylabel(antList[ants[1]])
        if ants[1] == 0:    # Antenna label in the top and leftside
            BLamp.set_title( antList[ants[0]] )
        if ants[0] == antNum - 1:    # Antenna at rightside
            BLamp.yaxis.tick_right()
        else:
            BLamp.yaxis.set_major_locator(plt.NullLocator())
        #-------- Plot closure phase
        if bl_index < 2: continue
        if ants[1] > 0:         # plot Closure phase
            BLphs = figSPW.add_subplot(antNum-1, antNum-1, (ants[0] - 1)*(antNum-1) + ants[1])
            BLphs.patch.set_facecolor('lightyellow')
            tri0, tri1, tri2 = Ant2Bl(ants[0], 0), Ant2Bl(ants[1], 0), Ant2Bl(ants[0], ants[1])
            for pol_index in list(range(polNum)):
                plotVis = scanVis[pol_index, tri0].conjugate()* scanVis[pol_index, tri1]* scanVis[pol_index, tri2]
                BLphs.plot(DT, np.angle(plotVis), '.', color=polColor[pol_index], label = 'Pol=' + polName[pol_index])
            #
            print('%d : %d - %d - %d (ant %s, %s, %s)' % (bl_index, tri0, tri1, tri2, antList[0], antList[ants[1]], antList[ants[0]]))
            BLphs.set_title('%s-%s-%s' % (antList[0], antList[ants[1]], antList[ants[0]] ), fontsize=0.5*fontSize)
            BLphs.axis([np.min(DT), np.max(DT), -math.pi, math.pi])
            BLphs.tick_params(axis='x', labelsize=int(fontSize*0.25), labelrotation=-90)
        if ants[1] == 1: BLphs.set_ylabel(antList[ants[0]] )
        if ants[1] > 1: BLphs.yaxis.set_major_locator(plt.NullLocator())
        if ants[0] < antNum - 1:    # except bottom panel : skip drawing X-axis
            BLphs.xaxis.set_major_locator(plt.NullLocator())
        #
    #
    figSPW.savefig(plotFile + '.pdf', format='pdf', dpi=144)
    plt.close('all')
    return
#
#-------- Plot XY-phase spectra
def plotXYP(pp, prefix, spwList, XYspec, bunchNum=1):
    spwNum = len(spwList)
    figXYP = plt.figure(figsize = (11, 8))
    figXYP.suptitle(prefix + ' XY Phase')
    figXYP.text(0.45, 0.05, 'Frequency [GHz]')
    figXYP.text(0.03, 0.45, 'XY Phase [deg]', rotation=90)
    for spw_index in list(range(spwNum)):
        spw = spwList[spw_index]
        chNum, chWid, Freq = GetChNum(prefix + '.ms', spw); Freq = 1.0e-9* bunchVec(Freq, bunchNum)  # GHz
        PhsPL = figXYP.add_subplot(1, spwNum, spw_index + 1)
        XYP  = XYspec[spw_index]
        PhsPL.plot( Freq, np.angle(XYP)*180.0/math.pi, '.', label = 'SPW %d' % (spw))
        PhsPL.axis([np.min(Freq), np.max(Freq), -180.0, 180.0])
        PhsPL.tick_params(axis='both', labelsize=6)
        PhsPL.legend(loc = 'lower left', prop={'size' :7}, numpoints = 1)
    #
    #plt.show()
    figXYP.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
#-------- Plot D-term spectra
def plotDSpec(pp, prefix, antList, spwList, FreqList, DxList, DyList):
    plotMax = 0.12
    antNum, spwNum = len(antList), len(spwList)
    figAnt = plt.figure(figsize = (11, 8))
    figAnt.suptitle(prefix + ' D-term spectra')
    figAnt.text(0.45, 0.05, 'Frequency [GHz]')
    figAnt.text(0.03, 0.45, 'D-term Spectra (Real and Imaginary)', rotation=90)
    for ant_index in list(range(antNum)):
        if ant_index > 0:
            for PL in DxPList: figAnt.delaxes(PL)
            for PL in DyPList: figAnt.delaxes(PL)
        #
        DxPList, DyPList = [], []
        for spw_index in list(range(spwNum)):
            DxPL = figAnt.add_subplot( 2, spwNum, spw_index + 1 )
            DyPL = figAnt.add_subplot( 2, spwNum, spw_index + spwNum + 1 )
            DxPList = DxPList + [DxPL]
            DyPList = DyPList + [DyPL]
            #
            plotDx, plotDy = DxList[ant_index + antNum* spw_index], DyList[ant_index + antNum* spw_index]
            DxPL.step( FreqList[spw_index], plotDx.real, where='mid', label = 'reDx')
            DxPL.step( FreqList[spw_index], plotDx.imag, where='mid', label = 'imDx')
            DxPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -plotMax, plotMax])
            DyPL.step( FreqList[spw_index], plotDy.real, where='mid', label = 'reDy')
            DyPL.step( FreqList[spw_index], plotDy.imag, where='mid', label = 'imDy')
            DyPL.axis([np.min(FreqList[spw_index]), np.max(FreqList[spw_index]), -plotMax, plotMax])
            #
            if spw_index == 0: DxPL.set_title(antList[ant_index])
            DxPL.tick_params(axis='both', labelsize=6)
            DyPL.tick_params(axis='both', labelsize=6)
        #
        DxPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        DyPL.legend(loc = 'upper right', prop={'size' :7}, numpoints = 1)
        #plt.show()
        figAnt.savefig(pp, format='pdf')
    #
    plt.close('all'); pp.close()
    del(DxPList); del(DyPList); del(DxPL); del(DyPL); del(figAnt)
    return
#
#-------- Plot Flagging map
def plotFlag(pp, prefix, antList, DT, spwList, FGList):
    antNum = len(antList)
    figFG = plt.figure(figsize = (8, 11))
    figFG.suptitle('%s Flag Map' % (prefix));
    figFG.text(0.45, 0.05, 'UTC on %s' % (DT[0].strftime('%Y-%m-%d')));
    for spw_index, spw in enumerate(spwList):
        FGPL = figFG.add_subplot(len(spwList), 1, spw_index + 1 )
        FG = FGList[spw_index]
        for ant_index, ant in enumerate(antList):
            FGPL.scatter(DT, np.linspace(ant_index, ant_index, len(DT)), vmin=0, vmax=1, marker='s', s=1, c= 1.0 - FG[ant_index], cmap=cm.seismic)
        FGPL.set_title('SPW%d' % (spw))
        FGPL.tick_params(labelsize=3)
        FGPL.set_yticks(list(range(antNum)))
        FGPL.set_yticklabels(antList.tolist())
    #
    figFG.savefig(pp, format='pdf')
    plt.close('all')
    pp.close()
    return
#
