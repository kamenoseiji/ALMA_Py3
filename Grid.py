import os
import math
import numpy as np
import scipy
#SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
SSOCatalog = ['Uranus', 'Neptune', 'Callisto', 'Ganymede', 'Titan', 'Io', 'Europa', 'Ceres', 'Pallas', 'Vesta', 'Juno', 'Mars', 'Mercury', 'Venus']
SSOscore   = [
[ 5.0,     4.0,       1.0,        1.0,        0.1,     0.2,  0.3,      0.2,     0.1,     0.1,    10.0,   10.0,     10.0],   # Band 1
[ 6.0,     5.0,       1.0,        1.0,        0.1,     0.6,  0.5,      0.3,     0.2,     0.2,    10.0,   10.0,     10.0],   # Band 2
[ 7.0,     6.0,       1.0,        1.0,        0.2,     0.7,  0.7,      0.5,     0.3,     0.3,    10.0,   10.0,     10.0],   # Band 3
[ 8.0,     7.0,       2.0,        2.0,        0.3,     0.9,  0.9,      0.6,     0.4,     0.4,     8.0,   10.0,      8.0],   # Band 4
[ 9.0,     8.0,       3.0,        3.0,        0.5,     1.0,  1.0,      0.8,     0.5,     0.5,     4.0,    4.0,      4.0],   # Band 5
[10.0,     9.0,       4.0,        4.0,        5.0,     3.0,  3.0,      1.0,     0.6,     0.6,     2.0,    2.0,      2.0],   # Band 6
[10.0,     9.0,       5.0,        5.0,        7.0,     4.0,  4.0,      3.0,     0.7,     0.7,     1.0,    1.0,      1.0],   # Band 7
[ 5.0,     6.0,       7.0,        7.0,        7.0,     8.0,  4.0,      3.0,     0.7,     0.7,     1.0,    1.0,      1.0],   # Band 8
[ 2.0,     3.0,       7.0,        7.0,        7.0,     9.0,  4.0,     10.0,     7.0,     4.0,     1.0,    1.0,      1.0]]   # Band 9
ELshadow = math.pi* 35.0 / 180.0
#ELshadow = math.pi* 65.0 / 180.0
SunAngleThresh = 5.0
#-- Acceptable opacity       B3   B4    B5   B6   B7   B8,  B9, B10
TauLimit = [0.0,  0.1, 0.1, 0.2, 0.2, 0.25, 0.5, 0.5, 0.5, 0.5, 0.6]
#-- Source name alias
sourceDic = {
'J0006-063':'J0006-0623',
'0006-063':'J0006-0623',
'J0237+288':'J0237+2848',
'NGC_1052': 'J0241-0815',
'3c84':'J0319+4130',
'0334-401':'J0334-4008',
'J0334-401':'J0334-4008',
'J0423-013':'J0423-0120',
'J0510+180':'J0510+1800',
'J0519-454':'J0519-4546',
'J0522-364':'J0522-3627',
'J0538-440':'J0538-4405',
'J0750+125':'J0750+1231',
'J0854+201':'J0854+2006',
'J1037-295':'J1037-2934',
'J1058+015':'J1058+0133',
'J1107-448':'J1107-4449',
'J1146+399':'J1146+3958',
'3c273':'J1229+0203',
'3C273':'J1229+0203',
'3c279':'J1256-0547',
'3C279':'J1256-0547',
'J1337-129':'J1337-1257',
'1337-129':'J1337-1257',
'3C286':'J1331+3030',
'1427-421':'J1427-4206',
'J1427-421':'J1427-4206',
'J1517-243':'J1517-2422',
'1517-243':'J1517-2422',
'J1550+054':'J1550+0527',
'1550+054':'J1550+0527',
'J1613-586':'J1617-5848',
'1613-586':'J1617-5848',
'J1625-254':'J1625-2527',
'3c345':'J1642+3948',
'J1733-130':'J1733-1304',
'1733-130':'J1733-1304',
'nrao530':'J1733-1304',
'J1751+096':'J1751+0939',
'1751+096':'J1751+0939',
'J1924-292':'J1924-2914',
'1924-292':'J1924-2914',
'J2025+337':'J2025+3343',
'J2056-472':'J2056-4714',
'2056-472':'J2056-4714',
'2157-694':'J2157-6941',
'J2148+069':'J2148+0657',
'J2232+117':'J2232+1143',
'3c454.3':'J2253+1608',
'J2258-279':'J2258-2758',
'BL_Lac':'J2202+4216',
'VLA_J1650+0824':'J1650+0824',
'VLA_J1652+0618':'J1652+0618',
'VLA_J1656+1826':'J1656+1826',
'VLA_J1658+0741':'J1658+0741',
'VLA_J1658+0515':'J1658+0515',
'VLA_J1700+0522':'J1700+0522',
'VLA_J1707+0148':'J1707+0148',
'VLA_J1707+1331':'J1707+1331',
'VLA_J1707+1846':'J1707+1846',
'VLA_J1708+0035':'J1708+0035',
'VLA_J1716+2152':'J1716+2152',
'VLA_J1719+0658':'J1719+0658',
'VLA_J1719+1745':'J1719+1745',
'VLA_J1719+0817':'J1719+0817',
'VLA_J1722+1013':'J1722+1013',
'VLA_J1726+0639':'J1726+0639',
'VLA_J1728+1215':'J1728+1215',
'VLA_J1728+0427':'J1728+0427',
'VLA_J1730+0024':'J1730+0024',
'VLA_J1734+0926':'J1734+0926',
'VLA_J1737+0621':'J1737+0621',
'VLA_J1740+2211':'J1740+2211',
'VLA_J1730+0024':'J1730+0024',
'VLA_J1734+0926':'J1734+0926',
'VLA_J1737+0621':'J1737+0621',
'VLA_J1740+2211':'J1740+2211'}

def sourceRename(sourceList):
    renameList = []
    for srcname in sourceList:
        srcname = srcname.replace('*','')
        renameList = renameList + [sourceDic.get(srcname, srcname)]
    return renameList
#
def AeNominal(msfile, antList):
    msmd.open(msfile)
    antDia = [msmd.antennadiameter(antName)['value'] for antName in antList]
    msmd.close()
    msmd.done()
    return 0.7* 0.25* np.pi* antDia**2      # Nominal Collecting Area
#
def nearestValue():
    TauEonList[spw_index][0]
#-------- Apply Tsys calibration to visibilities
def applyTsysCal(prefix, BandName, BandbpSPW, scanDic, SSODic, XspecList):
    from interferometry import ANT0, ANT1, Ant2Bl, kb, Tcmb, GetAntName, GetAntD, GetTemp, indexList, smoothValue
    #---- Check antenna list
    antList = GetAntName(prefix + '.ms')
    antNum = len(antList); blNum = int(antNum* (antNum - 1)/2)
    TrxAntList = np.load('%s-%s.TrxAnt.npy' % (prefix, BandName))
    antDia = GetAntD(TrxAntList)
    useAnt = indexList(TrxAntList,antList); useAntNum = len(useAnt); useBlNum = int(useAntNum* (useAntNum - 1)/2)
    ant0, ant1 = ANT0[0:useBlNum], ANT1[0:useBlNum]
    useBlMap = [Ant2Bl(useAnt[ant0[bl_index]], useAnt[ant1[bl_index]])  for bl_index in list(range(useBlNum))]
    nominalAe = 0.72
    #---- Load Tsys data
    tempAtm = GetTemp(prefix + '.ms')
    TauE = np.load('%s-%s.TauE.npy' % (prefix, BandName))   #  TauE[spw,scan]: time-variable excexs of zenith optical depth
    atmTime = np.load('%s-%s.atmTime.npy' % (prefix, BandName))#  atmTime[scan] : mjdSed at TauE measurements
    atmReltime = atmTime - atmTime[0]
    polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
    Tau0List, Tau0CList, TrxList, TaNList, TrxFreq, TauEonList = [], [], [], [], [], []
    for spw_index, spw in enumerate(BandbpSPW['spw']):
        Tau0List = Tau0List + [np.load('%s-%s-SPW%d.Tau0.npy' % (prefix, BandName, spw))]   # Tau0List[spw] [ch]
        Tau0CList= Tau0CList+ [np.load('%s-%s-SPW%d.Tau0C.npy'% (prefix, BandName, spw))]   # Tau0CList[spw] [intercept,slope]
        TrxList  = TrxList  + [np.load('%s-%s-SPW%d.Trx.npy'  % (prefix, BandName, spw))]   # TrxList[spw] [pol, ch, ant, scan]
        TaNList  = TaNList  + [np.load('%s-%s-SPW%d.TantN.npy'% (prefix, BandName, spw))]   # TaNList[spw] [ant, ch]
        TrxFreq  = TrxFreq  + [np.load('%s-%s-SPW%d.TrxFreq.npy'% (prefix, BandName, spw))] # TrxFreq[spw] [ch]
        TauEonPath = '%s-%s-SPW%d.TauEon.npy' % (prefix, BandName, spw)
        if os.path.isfile(TauEonPath): 
            TauEonList = TauEonList + [np.load(TauEonPath)]
        else:
            TauEonList = TauEonList + [np.array([])]
    for scan_index, scan in enumerate(scanDic.keys()):
        scanTau = []
        TsysScanDic = dict(zip(TrxAntList, [[]]* len(TrxAntList)))
        source = scanDic[scan]['source']
        for spw_index, spw in enumerate(BandbpSPW['spw']):
            chNum = BandbpSPW['chNum'][spw_index]
            TrxAnt = (np.median(TrxList[spw_index], axis=3) + TaNList[spw_index].T).transpose(1, 2, 0)  # [ch, ant, pol]
            StokesI = SSODic[source][1][spw_index] if source in SSOCatalog else scanDic[scan]['I'] 
            Tant = StokesI* nominalAe* np.pi* antDia**2 / (8.0* kb)                     # Antenna temperature of SSO
            Tau0SP = np.outer(Tau0List[spw_index], np.ones(len(scanDic[scan]['mjdSec'])))
            secZ = 1.0 / np.sin(scanDic[scan]['EL'])                           # Airmass
            if TauEonList[spw_index].shape[0] == 2:
                #Tau0SP = Tau0SP + TauEonList[spw_index][1][indexList(scanDic[scan]['mjdSec'], TauEonList[spw_index][0])]
                Tau0SP = Tau0SP + smoothValue(TauEonList[spw_index][0], TauEonList[spw_index][1], scanDic[scan]['mjdSec']).tolist()
            else: 
                SP = tauSMTH(atmReltime, TauE[spw_index] )
                Tau0SP = Tau0SP + scipy.interpolate.splev(scanDic[scan]['mjdSec'] - atmTime[0], SP)
            zenithTau = Tau0SP + Tau0CList[spw_index][0] + Tau0CList[spw_index][1]*secZ   # Smoothed zenith optical depth
            scanTau = scanTau + [zenithTau * secZ]  # Optical depth at the elevation
            exp_Tau = np.exp(-zenithTau * secZ )    # Atmospheric attenuation
            atmCorrect = np.mean(1.0 / exp_Tau, axis=1)              # Correction for atmospheric attenuation
            TsysScan = (Tcmb + Tant + (TrxAnt.transpose(2,1,0)* atmCorrect + tempAtm* (atmCorrect - 1.0)).transpose(0,2,1)).transpose(2,0,1)
            #-------- Tsys correction
            Xspec = XspecList[spw_index][scan_index][:,:,useBlMap].transpose(3,2,0,1)* np.sqrt(TsysScan[ant0][:,polXindex]* TsysScan[ant1][:,polYindex])
            XspecList[spw_index][scan_index][:,:,useBlMap] = Xspec.transpose(2,3,1,0)
            for ant_index, ant in enumerate(TrxAntList): TsysScanDic[ant] = TsysScanDic[ant] + [TsysScan[ant_index]]
        scanDic[scan]['Tau']  = scanTau
        scanDic[scan]['Tsys'] = TsysScanDic
    #
    return scanDic, XspecList
#
def aprioriSEFD(Ae, EL, TrxSpec, Tau0Spec):
    secZ = 1.0 / np.sin(EL)
    zenithTau = Tau0spec + scipy.interpolate.splev(np.median(timeStamp), exTauSP) + Tau0Coef[spw_index][0] + Tau0Coef[spw_index][1]*secZ
    exp_Tau = np.exp(-zenithTau * secZ )
    TsysEQScan = np.mean(TrxList[spw_index].transpose(2,0,1)[:,:,chRange] + Tcmb*exp_Tau[chRange] + tempAtm* (1.0 - exp_Tau[chRange]), axis=2)[Trx2antMap] # [antMap, pol]
    return 2.0* kb* TsysEQScan.T / Ae
#
#-------- Disk Visibility
def diskVis(diskRadius, u):
    # diskRadius : radius of planet disk [rad]
    # u          : spatial frequency (= baseline / wavelength)
    argument = 2.0* np.pi* u* diskRadius
    return 2.0* scipy.special.jn(1, argument) / argument
#
#-------- Disk Visibility with primary beam correction, u must be smaller than 0.3/diskRadius
def diskVisBeam(diskShape, u, v, primaryBeam):
    # diskShape  : planet disk diameter [MajorAxis, MinorAxis, PA] (rad)
    # u,v        : spatial frequency (= baseline / wavelength)
    # primaryBeam: FWHM of primary beam [rad]
    from interferometry import beamF
    cs, sn = np.cos(diskShape[2]), np.sin(diskShape[2])
    diskRadius = 0.5* np.sqrt(diskShape[0]* diskShape[1])
    DSmaj = 1.0 / np.sqrt( (0.30585 / diskShape[0])**2 + 2.0* np.log(2.0)/(np.pi* primaryBeam)**2 )    # Primary-beam correction
    DSmin = 1.0 / np.sqrt( (0.30585 / diskShape[1])**2 + 2.0* np.log(2.0)/(np.pi* primaryBeam)**2 )    # Primary-beam correction
    uvDisp = (DSmin*(u* cs - v* sn))**2 + (DSmaj*(u* sn + v* cs))**2
    return beamF(diskRadius/primaryBeam)* np.exp(-0.5* uvDisp)
#
#-------- Apertue effciency measurements using Solar System Objects
def SSOAe(antList, spwDic, uvw, scanDic, SSODic, XSList):
    # antList   : List of antenna name
    # spwDic    : SPW dictionary ['spw', 'freq', 'chNum', 'chRange', 'BW']
    # uvw       : baseline vector [m]
    # scanDic   : scan dictionary ['msfile', 'source', 'mjdSec', 'EL', 'PA', 'I', 'QCpUS', 'Tau', 'Tsys', 'Gain']
    # SSODic    : SSO dictionary
    # XSList    : Cross Correlation XspecList[spw][pol, ch, bl, time]
    from interferometry import GetAntD, Bl2Ant, ANT0, ANT1, kb, gainComplex
    SSOname = scanDic['source']
    text_sd = ' Flux Calibrator : %10s EL=%.1f' % (SSOname, 180.0*np.median(scanDic['EL'])/np.pi)
    uvDist = np.sqrt(uvw[0]**2 + uvw[1]**2)
    UVlimit = 0.5 / SSODic[SSOname][2][0]          # Maximum usable uv distance [lambda]
    text_sd = text_sd + ' uv limit = %5.0f klambda' % (UVlimit*1.0e-3); print(text_sd)
    antNum, blNum, antDia = len(antList), len(uvDist), GetAntD(antList)
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    #-------- Check usable antenna/baselines
    uvFlag = np.ones(blNum)
    for spw_index, spw in enumerate(spwDic['spw']):
        chRange = spwDic['chRange'][spw_index]
        centerFreq = np.mean(spwDic['freq'][spw_index][chRange])
        uvFlag[np.where( uvDist > UVlimit* 299792458 / centerFreq)[0].tolist()] *= 0.0
    #
    SAant = [0] + [Bl2Ant(bl_index)[0] for bl_index in np.where(uvFlag > 0.1)[0].tolist() if Bl2Ant(bl_index)[1] == 0]
    SAbl  = [bl_index for bl_index in list(range(blNum)) if Bl2Ant(bl_index)[0] in SAant and Bl2Ant(bl_index)[1] in SAant]
    if len(SAant) < 4: return #  Too few antennas
    #-------- Baseline map
    print('Subarray : ',); print(antList[SAant])
    for ant_index in list(range(1, antNum)):
        text_sd = antList[ant_index] + ' : '; print(text_sd, end='')
        blList = [bl_index for bl_index in list(range(blNum)) if Bl2Ant(bl_index)[0] == ant_index]
        for bl_index in blList:
            if uvFlag[bl_index] < 1.0: text_sd = '\033[91m%4.0f\033[0m' % (uvDist[bl_index])
            else: text_sd = '%4.0f' % (uvDist[bl_index])
            print(text_sd, end=' ')
        print('')
    print('       ', end='')
    for ant_index, ant in enumerate(antList[0:antNum-1]): print(ant, end=' ')
    print('')
    AeSPW, WgSPW, SSOmodel = [], [], []
    #-------- Aperture Efficiency
    for spw_index, spw in enumerate(spwDic['spw']):
        uvWave = uvw[0:2,:] * centerFreq / 299792458    # UV distance in wavelength
        primaryBeam = 1.13* 299792458 / (np.pi * antDia* centerFreq)
        SSOmodelVis = SSODic[SSOname][1][spw_index]*  diskVisBeam(SSODic[SSOname][2], uvWave[0], uvWave[1], primaryBeam[ant0]* primaryBeam[ant0]* np.sqrt(2.0 / (primaryBeam[ant0]**2 + primaryBeam[ant1]**2))) + 1.0e-9
        VisChav = np.mean(XSList[spw_index][:,chRange][:,:,SAbl], axis=(1,3)) / SSOmodelVis[SAbl]
        Gain = np.apply_along_axis(gainComplex, 1, VisChav)  # Gain[pol, ant]
        Aeff = 8.0* kb* abs(Gain)**2 / (np.pi* antDia[SAant]**2)
        Ae, Wg = np.zeros([antNum, 2]), np.zeros([antNum, 2])
        for ant_index, SA in enumerate(SAant):
            Ae[SA] = Aeff[:, ant_index]
            Wg[SA] = np.sign(Aeff[:, ant_index])* np.median(abs(SSOmodelVis))
        if np.median(Ae) > 0.99 : Wg *= 1.0e-6
        if np.median(Ae) < 0.19 : Wg *= 1.0e-2
        if np.median(scanDic['EL']) < ELshadow : Wg *= 1.0e-1
        AeSPW = AeSPW + [Ae]
        WgSPW = WgSPW + [Wg]
        SSOmodel = SSOmodel + [SSOmodelVis]
    FscaleDic = {
        'Ae'    : AeSPW,
        'Wg'    : WgSPW,
        'model' : SSOmodel}
    return FscaleDic
#-------- Average Ae among multiple SSOs
def averageAe(FscaleDic, spwList):
    AeList, WgList = [], []
    for SSO_index, SSOname in enumerate(FscaleDic.keys()):
        if FscaleDic[SSOname] is None: continue
        if np.median(np.array(FscaleDic[SSOname]['Wg'])) < 1.0e-9: continue
        AeSPW = FscaleDic[SSOname]['Ae']; AeList = AeList + [np.array(AeSPW)]
        WgSPW = FscaleDic[SSOname]['Wg']; WgList = WgList + [np.array(WgSPW)]
    return  (np.sum(np.array(WgList) * np.array(AeList), axis=0)/(np.sum(np.array(WgList), axis=0)+1.0e-9)).transpose(1,2,0)   # Aeff[ant, pol, spw]
#-------- Transfer and equalize aperture efficiencies
def AeTransfer(VisChav, Aeff, antDia):
    from interferometry import Bl2Ant, gainComplexVec
    blNum = VisChav.shape[1]
    SAant = np.where(np.min(Aeff, axis=1) > 0.25)[0].tolist()
    SAbl  = [bl_index for bl_index in list(range(blNum)) if Bl2Ant(bl_index)[0] in SAant and Bl2Ant(bl_index)[1] in SAant]
    GainSA  = gainComplexVec(VisChav[:,SAbl].T)
    GainAll = gainComplexVec(VisChav.T)
    scaleFlux = np.median(abs(GainSA**2).T / (antDia[SAant]**2 * Aeff[SAant].T), axis=1)
    return ((abs(GainAll)**2 / scaleFlux).T / (antDia**2)).T
#
#-------- Gain scaling
def GainScale(Aeff, antDia, polVis):
    # Aeff :    Aperure Efficiency Table
    # antDia :  Antenna diameter
    # polVis[pol,bl,time] :    cross power spectra
    from interferometry import kb, ANT0, ANT1
    polXindex, polYindex = (np.arange(4)//2).tolist(), (np.arange(4)%2).tolist()
    blNum = polVis.shape[1]
    ant0, ant1 = ANT0[0:blNum], ANT1[0:blNum]
    fluxScale = np.sqrt( 2.0* kb / (0.25* np.pi* antDia**2* Aeff.T))
    return (polVis.transpose(2,0,1)* fluxScale[polYindex][:,ant0]* fluxScale[polXindex][:,ant1]).transpose(1,2,0)
#-------- Stokes parameters
def Vis2Stokes(VisChav, Dcat, PA):
    # VisChav   : channel-averaged visibiliities [pol, bl, time]
    # Dcat      : D-term [ant, pol]
    # PA        : Parallactic Angle with respect to X-feed [time]
    from interferometry import ANT0, ANT1, InvMullerMatrix, InvPAVector
    PAnum = len(PA)
    PS = InvPAVector(PA, np.ones(PAnum))
    blNum = VisChav.shape[1]
    Stokes = np.zeros([4,blNum], dtype=complex)
    for bl_index in list(range(blNum)):
        Minv = InvMullerMatrix(Dcat[ANT1[bl_index], 0], Dcat[ANT1[bl_index], 1], Dcat[ANT0[bl_index], 0], Dcat[ANT0[bl_index], 1])
        Stokes[:,bl_index] = PS.reshape(4, 4*PAnum).dot(Minv.dot( VisChav[:,bl_index]).reshape(4*PAnum)) / PAnum
    #
    return Stokes
#-------- Linear regression for visibility-baseline relation
def lmStokes(StokesVis, uvDist):
    # StokesVis  : Stokes visibilities [Stokes, bl]
    # uvDist     : Projected baseline length [bl]
    StokesFlux, StokesSlope, StokesErr = np.zeros([4]), np.zeros([4]), np.ones([4])
    percent80, sdvis = np.percentile(StokesVis[0].real, 80),  np.std(StokesVis[0].real)
    visFlag = np.where(abs(StokesVis[0].real - percent80) < 2.0* sdvis )[0]      # 3-sigma critesion
    weight = np.zeros(len(uvDist)); weight[visFlag] = 1.0/(np.var(StokesVis[0][visFlag].real)* uvDist[visFlag])
    P, W = np.c_[np.ones(len(weight)), uvDist], np.diag(weight)
    PtWP_inv = scipy.linalg.inv(P.T.dot(W.dot(P)))
    solution, solerr = PtWP_inv.dot(P.T.dot(weight* StokesVis[0].real)),  np.sqrt(np.diag(PtWP_inv)) # solution[0]:intercept, solution[1]:slope
    if abs(solution[1]) < 2.0* solerr[1]: solution[0], solution[1] = np.median(StokesVis[0][visFlag].real), 0.0
    StokesFlux[0], StokesSlope[0], StokesErr[0] = solution[0], solution[1], solerr[0]
    for pol_index in [1,2,3]:
        StokesFlux[pol_index] = StokesSlope[0] * np.median(StokesVis[pol_index].real)/StokesFlux[0]
        solution[0] = (weight.dot(StokesVis[pol_index].real) - StokesSlope[pol_index]* weight.dot(uvDist))/(np.sum(weight))
        StokesFlux[pol_index] = solution[0]
        resid = StokesVis[pol_index].real - StokesSlope[pol_index]* uvDist - solution[0]
        StokesErr[pol_index] = np.sqrt(weight.dot(resid**2)/np.sum(weight))
    return StokesFlux, StokesSlope, StokesErr
#-------- Smooth time-variable Tau
def tauSMTH( timeSample, TauE ):
    if len(timeSample) > 5:
        SplineWeight = np.ones(len(timeSample) + 4)
        flagIndex = (np.where(abs(TauE - np.median(TauE))/np.std(TauE) > 3.0)[0] + 2).tolist()
        SplineWeight[flagIndex] = 0.01
        tempTime = np.append([timeSample[0]-500.0, timeSample[0]-300.0], np.append(timeSample, [timeSample[-1]+300.0, timeSample[-1]+500.0]))
        tempTauE = np.append([TauE[0], TauE[0]], np.append(TauE, [TauE[-1], TauE[-1]]))
        smthTau = scipy.interpolate.splrep(tempTime, tempTauE, k=3, w=SplineWeight, t=tempTime[list(range(1, len(tempTime), 2))] - 60.0 )
    else:
        tempTime = np.arange(np.min(timeSample) - 3600.0,  np.max(timeSample) + 3600.0, 300.0)
        tempTauE = np.repeat(np.median(TauE), len(tempTime))
        smthTau = scipy.interpolate.splrep(tempTime, tempTauE, k=3)
    #
    return smthTau
#
