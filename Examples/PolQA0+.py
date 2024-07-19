SCR_DIR = '/users/skameno/ALMA_Py3/'
R_DIR = '/usr/bin/'
wd = './'
QUMODEL = True
INTList = ['BANDPASS', 'POL']
sessionFile = open('SessionList', 'r')
sessionList = sessionFile.readlines()
sessionFile.close()
for sessionEntry in sessionList:
    if sessionEntry[0] == '#': continue
    entry = sessionEntry.split()
    Session = entry[0]
    UIDList = entry[1:]
    chBunch = 1
    prefixList = []
    for UID in UIDList:
        prefixList = prefixList + [UID.replace("/", "_").replace(":","_").replace(" ","")]
    exec(open(SCR_DIR + 'ExportPol.py').read())
    if 'spwNames' in locals(): del(spwNames)
    prefix = Session
    exec(open(SCR_DIR + 'checkPolCalScans.py').read())
    scanList = onsourceScans
    spwList = bpspwLists[-1]
    UseSPWList = spwList
    BPPLOT = True
    prefix = Session
    SNR_THRESH = 2
    antFlag = []
    antSPWFlag = []
    for spw in spwList:
        exec(open(SCR_DIR + 'checkGain.py').read())
        antSPWFlag = antSPWFlag + [antFlag + newAntFlag]
    antFlag = antFlag + list(set([ant for row in antSPWFlag for ant in row]))
    refantName = antList[UseAnt[refantID]]
    refant = refantName
    for BPscan in scanList:
        exec(open(SCR_DIR + 'checkBP.py').read())
    exec(open(SCR_DIR + 'averageBP.py').read())
    FGprefix = Session
    BPprefix = Session
    XYprefix = Session
    BPscan = 0
    exec(open(SCR_DIR + 'solveDterm.py').read())
    del antList
    exec(open(SCR_DIR + 'plotDterm.py').read())
    exec(open(SCR_DIR + 'plotXYC.py').read())
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

