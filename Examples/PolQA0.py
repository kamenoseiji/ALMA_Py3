SCR_DIR = '/users/skameno/ALMA_Py3/'
R_DIR = '/usr/bin/'
wd = './'
prefixList = ['uid___A002_Xc1e2be_X39e9', 'uid___A002_Xc1e2be_X3d13', 'uid___A002_Xc1e2be_X3f7a']
for prefix in prefixList:
     if not os.path.isdir(prefix): os.system('asdmExport -m ' + prefix )
exec(open(SCR_DIR + 'ASDMPol.py').read())
