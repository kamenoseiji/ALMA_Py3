SCR_DIR = '/users/skameno/ALMA_Py3/'
R_DIR = '/usr/bin/'
wd = './'
sessionFile = open('SessionList', 'r')
sessionList = sessionFile.readlines()
sessionFile.close()
for sessionEntry in sessionList:
    if sessionEntry[0] == '#': continue
    entry = sessionEntry.split()
    sessionName = entry[0]
    EBList = entry[1:]
    prefixList = []
    for UID in EBList:
        prefixList = prefixList + [UID.replace("/", "_").replace(":","_").replace(" ","")]
    #
    for prefix in prefixList:
        if not os.path.isdir(prefix): os.system('asdmExport -m ' + prefix )
    exec(open(SCR_DIR + 'ASDMPol.py').read())
#
