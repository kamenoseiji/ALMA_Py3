import os
import glob
import numpy as np
from datetime import datetime
def generetaCheckEach( prefix ):
    return('sed -e "s/uid___A002_X[0-9a-f]*_X[0-9a-f]*/%s/g" template.py > checkEach.py' % (prefix))
#
def replaceList( UID, status ):
    UID_text = UID.replace(':','\:').replace('/','\/')
    return('sed -e "s/new %s/%s %s/g" UID/UIDList > UID/UIDnew' % (UID_text, status, UID_text))
#
UIDfile = open('./UID/UIDList', 'r')
UIDentry = UIDfile.readlines()
UIDfile.close()
UIDentry = list(set(UIDentry))
newEntry = [entry for entry in UIDentry if entry.split()[0] == 'new']
FSUIDs = list(set([entry.split('/')[3] for entry in newEntry]))     # Each FSR
for FS in FSUIDs:
    FSentry = [entry for entry in UIDentry if FS in entry]
    ArrayIDs = list(set([entry.split()[2] for entry in FSentry]))   # Each Array
    for Array in ArrayIDs:
        ARentry = [entry for entry in UIDentry if Array in entry]
        DT     = [datetime.strptime(entry.split()[4], "%Y-%m-%dT%H:%M:%S").timestamp() for entry in ARentry]
        EBList = [entry.split()[1].replace("/", "_").replace(":","_").replace(" ","") for entry in ARentry]
        sort_index = np.argsort(DT).tolist()
        for prefix in EBList:
            if not os.path.isdir(prefix + '.ms'):
                if not os.path.isdir(prefix):
                    if not os.path.isdir('UID/' + prefix):
                        os.system('asdmExport %s' % (prefix))
                        os.system('mv %s UID/' % (prefix))
                    os.system('ln -s UID/' + prefix + ' .')
                importasdm(prefix)
        if len(ARentry) > 1:
            text_sd = 'casa -c ~/ALMA_Py3/splitMerge.py -u '
            for index in sort_index: text_sd = text_sd + EBList[index] + ','
            print(text_sd[:-1])
            os.system(text_sd[:-1])
            for index in sort_index: os.system('rm -rf %s.ms*' % (EBList[index]))
            prefixElement = EBList[0].split('_X'); newPrefix = prefixElement[0] + '_X' + prefixElement[1]
            MSList = glob.glob(newPrefix + '.RB_*.ms')
            for msfile in MSList:
                text_sd = 'casa -c ~/ALMA_Py3/checkGridSurvey.py -u %s' % (msfile.split('.ms')[0])
                os.system(text_sd)
                print(text_sd)
                text_sd = 'casa -c ~/ALMA_Py3/splitFluxLog.py -u '
                for index in sort_index: text_sd = text_sd + EBList[index] + ','
                os.system(text_sd[:-1])
                print(text_sd[:-1])
                for index in sort_index:
                    print('scp ' + EBList[index] + '*Flux.log skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/')
                    os.system('scp ' + EBList[index] + '*Flux.log skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/')
                    text_sd = replaceList(EBList[index], 'done')
                    os.system(text_sd)
                    os.system('mv ./UID/UIDnew ./UID/UIDList')
                os.system('rm -rf %s' % (msfile))
                os.system('mv *.npy NPY/')
                os.system('mv *.pdf PDF/')
                os.system('mv *.log LOG/')
                os.system('rm -rf *.cl')
        else:
            text_sd = 'casa -c ~/ALMA_Py3/checkGridSurvey.py -u %s' % (EBList[0])
            os.system(text_sd)
            print(text_sd)
            print('scp ' + EBList[0] + '*Flux.log skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/')
            text_sd = replaceList(EBList[0], 'done')
            os.system(text_sd)
            os.system('mv ./UID/UIDnew ./UID/UIDList')
            os.system('mv *.npy NPY/')
            os.system('mv *.pdf PDF/')
            os.system('mv *.log LOG/')
            os.system('rm -rf *.cl')
#
