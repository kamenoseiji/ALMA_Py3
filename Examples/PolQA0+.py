SCR_DIR = '/users/skameno/ALMA_Py3/'
R_DIR = '/usr/bin/'
wd = './'
QUMODEL = True
sessionFile = open('SessionList', 'r')
sessionList = sessionFile.readlines()
sessionFile.close()
for sessionEntry in sessionList:
    if sessionEntry[0] == '#': continue
    entry = sessionEntry.split()
    Session = entry[0]
    UIDList = entry[1:]
    prefixList = []
    for UID in UIDList:
        prefixList = prefixList + [UID.replace("/", "_").replace(":","_").replace(" ","")]
    #
    exec(open(SCR_DIR + 'ExportPol.py').read())
    if 'spwNames' in locals(): del(spwNames)
    prefix = Session
    exec(open(SCR_DIR + 'checkPolCalScans.py').read())
    scanList = onsourceScans
    spwList = bpspwLists[-1]
    antFlag = ['DV04']
    for spw in spwList:
        exec(open(SCR_DIR + 'checkGain.py').read())
    BPscan = BPScan
    exec(open(SCR_DIR + 'checkBP.py').read())
    refantName = antList[UseAnt[refantID]]
    BPprefix = Session
    XYprefix = Session
    exec(open(SCR_DIR + 'solveDterm.py').read())
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
    os.system('mv *.py ' + Session + '/Scripts/')
    os.system('mv *.npy ' + Session + '/NPY/')
    os.system('mv *.log ' + Session + '/LOG/')
    os.system('mv *.listobs ' + Session + '/LOG/')
    os.system('mv *.dic ' + Session + '/LOG/')
    os.system('scp -r ' + Session + '/ skameno@ssh.alma.cl:/home/skameno/public_html/POL/')
#
