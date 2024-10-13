import numpy as np
import os
SCR_DIR = '/users/skameno/ALMA_Py3/'
R_DIR = '/usr/bin/'
wd = './'
QUMODEL = True
INTList = 'POL,BANDPASS'
flagAnt = ''
sessionFile = open('SessionList', 'r')
sessionList = sessionFile.readlines()
sessionFile.close()
for sessionEntry in sessionList:
    #-------- Step1 : load ASDM, make MS, split, and concatinate
    if sessionEntry[0] == '#': continue
    entry = sessionEntry.split()
    Session = entry[0]
    UIDList = entry[1:]
    freqRes = 15 # [MHz]
    text_sd = 'casa -c %sExportPol.py -I %s -S %s -f %d' % (SCR_DIR, INTList, Session, freqRes)
    if flagAnt != '': text_sd = text_sd + ' -a %s' % (flagAnt)
    text_sd = text_sd + ' -u '
    for UID in UIDList:
        text_sd = text_sd + UID.replace("/", "_").replace(":","_").replace(" ","") + ','
    print(text_sd[:-1])
    #os.system(text_sd[:-1])
    #-------- Step2 : Check a priori properties of polarization calibrators
    prefix = Session
    text_sd = 'casa -c %scheckPolCalScans.py -u %s -Q' % (SCR_DIR, Session)
    print(text_sd)
    os.system(text_sd)
    fp = open('%s-PolQuery.log' % (prefix))
    polLines = fp.readlines()
    fp.close()
    scanList, spwList = [], []
    for eachLine in polLines:
        if eachLine.split()[1] == 'SPW':
            for spw_index, spw in enumerate( eachLine.split()[2:] ): spwList = spwList + [int(spw)]
        if eachLine.split()[0].isdecimal(): scanList = scanList + [int(eachLine.split()[0])]
    #
    print(scanList)
    #-------- Step3 : Check parallel-hand antenna-based gain and generate flag table
    for spw in spwList:
        text_sd = 'casa -c %scheckGain.py -u %s -s %d -T 0.2 -P -c ' % (SCR_DIR, prefix, spw)
        for scan in scanList: text_sd = text_sd + '%d,' % (scan)
        print(text_sd[:-1])
        os.system(text_sd[:-1])
    text_sd = 'casa -c %splotFG.py -u %s -s ' % (SCR_DIR, prefix)
    for spw in spwList:
        text_sd = text_sd + '%s,' % (spw)
    print(text_sd[:-1])
    os.system(text_sd[:-1])
    #-------- Step4 : Bandpass table
    refantName = np.load('%s.Ant.npy' % (prefix))[0]
    for scan in scanList:
        text_sd = 'casa -c %scheckBP.py -u %s -c %d -P -s ' % (SCR_DIR, prefix, scan)
        for spw in spwList: text_sd = text_sd + '%d,' % (spw)
        fp = open('%s-PolFlag.log' % (prefix))
        fgLines = fp.readlines()
        flagAntList = []
        if len(fgLines) > 0:
            for eachLine in fgLines: flagAntList = flagAntList + [eachLine.split()[1]]
            flagAntList = list(set(flagAntList))
            text_sd = text_sd[:-1] + ' -a '
            for flagAnt in flagAntList: text_sd = text_sd + '%s,' % (flagAnt)
        #
        print(text_sd[:-1])
        os.system(text_sd[:-1])
    #
    #-------- Step5 : Average bandpass table
    text_sd = 'casa -c %saverageBP.py -u %s -P -s ' % (SCR_DIR, prefix)
    for spw in spwList: text_sd = text_sd + '%d,' % (spw)
    text_sd = text_sd[:-1] + ' -c '
    for scan in scanList: text_sd = text_sd + '%d,' % (scan)
    text_sd = text_sd[:-1] + ' -R ' + refantName
    print(text_sd)
    os.system(text_sd)
    #-------- Step6 : Solve for D-term
    text_sd = 'casa -c %ssolveDterm.py -u %s -Q -R %s -s ' % (SCR_DIR, prefix, refantName)
    for spw in spwList: text_sd = text_sd + '%d,' % (spw)
    text_sd = text_sd[:-1] + ' -c '
    for scan in scanList: text_sd = text_sd + '%d,' % (scan)
    if len(flagAntList) > 0:
        text_sd = text_sd[:-1] + ' -a '
        for flagAnt in flagAntList: text_sd = text_sd + '%s,' % (flagAnt)
    print(text_sd[:-1])
    os.system(text_sd[:-1])
    #-------- Step7 : Plots D-term and XY phase 
    text_sd = 'casa -c %splotDterm.py -u %s -R %s -s ' % (SCR_DIR, prefix, refantName)
    for spw in spwList: text_sd = text_sd + '%d,' % (spw)
    print(text_sd[:-1])
    os.system(text_sd[:-1])
    text_sd = 'casa -c %splotXYC.py -u %s -R %s -s ' % (SCR_DIR, prefix, refantName)
    for spw in spwList: text_sd = text_sd + '%d,' % (spw)
    print(text_sd[:-1])
    os.system(text_sd[:-1])
    #-------- Step8 : Store reports into a directory
    os.system('rm -rf ' + Session)
    os.system('mkdir ' + Session)
    os.system('mkdir ' + Session + '/Dterm')
    os.system('mkdir ' + Session + '/Scripts')
    os.system('mkdir ' + Session + '/PDF')
    os.system('mkdir ' + Session + '/NPY')
    os.system('mkdir ' + Session + '/LOG')
    os.system('mv *DSpec.npy ' + Session + '/Dterm/')
    os.system('mv *.pdf ' + Session + '/PDF/')
    os.system('cp *.py ' + Session + '/Scripts/')
    os.system('mv *.npy ' + Session + '/NPY/')
    os.system('mv *.log ' + Session + '/LOG/')
    os.system('mv *.listobs ' + Session + '/LOG/')
    os.system('mv *.dic ' + Session + '/LOG/')
    os.system('mv *.txt ' + Session + '/LOG/')
#
