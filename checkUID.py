import os
import glob
import numpy as np
from datetime import datetime
from interferometry import GetAntName, GetBandNames
def replaceList( UID, status ):
    UID_text = UID.replace("___", r"\:\/\/").replace("_",r"\/")
    return('sed -e "s/new %s/%s %s/g" UID/UIDList > UID/UIDnew' % (UID_text, status, UID_text))
#
UIDfile = open('./UID/UIDList', 'r')
UIDentry = UIDfile.readlines()
UIDfile.close()
UIDentry = list(set(UIDentry))
newEntry = [entry for entry in UIDentry if (entry.split()[0] == 'new') and ('B9' not in entry.split()[-1])]
FSUIDs = list(set([entry.split('/')[3] for entry in newEntry]))     # Each FSR
for FS in FSUIDs:
    FSentry = [entry for entry in newEntry if FS in entry]
    ArrayIDs = list(set([entry.split()[2] for entry in FSentry]))   # Each Array
    for Array in ArrayIDs:
        ARentry = [entry for entry in FSentry if Array in entry]
        DT     = [datetime.strptime(entry.split()[4], "%Y-%m-%dT%H:%M:%S").timestamp() for entry in ARentry]
        EBList = [entry.split()[1].replace("/", "_").replace(":","_").replace(" ","") for entry in ARentry]
        #-------- Generage MS
        for prefix in EBList:                                       # Each EB
            if not os.path.isdir(prefix + '.ms'):
                if not os.path.isdir(prefix):
                    if not os.path.isdir('UID/' + prefix):
                        os.system('asdmExport %s' % (prefix))
                        os.system('mv %s UID/' % (prefix))
                    os.system('ln -s UID/' + prefix + ' .')
                importasdm(prefix)
        #-------- Classification by band
        bandEB = [np.unique(GetBandNames(prefix + '.ms')).tolist() for prefix in EBList]
        bandList = np.unique(np.array(bandEB)).tolist()
        #-------- Check Number of Antennas in the array
        antList = GetAntName(prefix + '.ms')
        if len(antList) < 4: 
            for index in sort_index:
                text_sd = replaceList(EBList[index], 'fail')
                os.system(text_sd)
                os.system('mv ./UID/UIDnew ./UID/UIDList')
            continue
        #-------- Each band
        for bandName in bandList:
            EBindexList = [index for index,band in enumerate(bandEB) if bandName in band]
            if len(EBindexList) > 1:
                DTband, EBband = np.array(DT)[EBindexList], np.array(EBList)[EBindexList]
                sort_index = np.argsort(DTband).tolist()
                #-------- Concatinate multiple EBs with the same array
                text_sd = 'casa -c ~/ALMA_Py3/splitMerge.py -u '
                for index in sort_index: text_sd = text_sd + EBband[index] + ','
                print(text_sd[:-1])
                os.system(text_sd[:-1])
                for index in sort_index: os.system('rm -rf %s.ms*' % (EBband[index]))
                #-------- Run AMAPOLA reduction for concatenated MS
                prefixElement = EBband[0].split('_X'); newPrefix = prefixElement[0] + '_X' + prefixElement[1]
                prefix = '%s.%s' % (newPrefix, bandName)
                text_sd = 'casa -c ~/ALMA_Py3/checkGridSurvey.py -u %s' % (prefix)
                print(text_sd)
                os.system(text_sd)
                #-------- Separate Flux.log
                text_sd = 'casa -c ~/ALMA_Py3/splitFluxLog.py -u '
                for index in sort_index: text_sd = text_sd + EBband[index] + ','
                print(text_sd[:-1])
                os.system(text_sd[:-1])
                for index in sort_index:
                    print('scp ' + EBband[index] + '-' + bandName + '-Flux.log skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/')
                    os.system('scp ' + EBband[index] + '-' + bandName + '-Flux.log skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/')
                    text_sd = replaceList(EBband[index], 'done')
                    os.system(text_sd)
                    os.system('mv ./UID/UIDnew ./UID/UIDList')
                    os.system('rm -rf %s.%s.ms' % (EBband[index], bandName))
                os.system('rm -rf %s.ms' % (prefix))
                os.system('rm -rf casa*.log')
                os.system('mv *.npy NPY/')
                os.system('mv *.pdf PDF/')
                os.system('mv *.log LOG/')
                os.system('rm -rf *.cl')
            else:
                prefix = EBList[EBindexList[0]]
                text_sd = 'casa -c ~/ALMA_Py3/checkGridSurvey.py -u %s' % (prefix)
                print(text_sd)
                os.system(text_sd)
                print('scp ' + prefix + '-' + bandName + '-Flux.log skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/')
                os.system('scp ' + prefix + '-' + bandName + '-Flux.log skameno@ssh.alma.cl:/home/skameno/public_html/Grid/Stokes/')
                text_sd = replaceList(prefix, 'done')
                os.system(text_sd)
                os.system('mv ./UID/UIDnew ./UID/UIDList')
                os.system('rm -rf casa*.log')
                os.system('mv *.npy NPY/')
                os.system('mv *.pdf PDF/')
                os.system('mv *.log LOG/')
                os.system('rm -rf *.cl')
#
