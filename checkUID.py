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
        '''
        for prefix in EBList:
            if not os.path.isdir(prefix + '.ms'):
                if not os.path.isdir(prefix):
                    if not os.path.isdir('UID/' + prefix):
                        os.system('asdmExport %s' % (prefix))
                        os.system('mv %s UID/' % (prefix))
                os.system('ln -s UID/' + prefix + ' .')
                importasdm(prefix)
        '''
        if len(ARentry) > 1:
            text_sd = 'casa -c ~/ALMA_Py3/splitMerge.py -u '
            for index in sort_index: text_sd = text_sd + EBList[index] + ','
            #os.system(text_sd[:-1])
            print(text_sd[:-1])
            prefixElement = EBList[0].split('_X'); newPrefix = prefixElement[0] + '_X' + prefixElement[1]
            MSList = glob.glob(newPrefix + '*.ms')
            for msfile in MSList:
                text_sd = 'casa -c ~/ALMA_Py3/checkGridSurvey.py -u %s' % (msfile.split()[0])
                #os.system(text_sd)
                print(text_sd)
                text_sd = 'casa -c ~/ALMA_Py3/splitFluxLog.py -u '
                for index in sort_index: text_sd = text_sd + EBList[index] + ','
                #os.system(text_sd[:-1])
                print(text_sd[:-1])
        else:
            text_sd = 'casa -c ~/ALMA_Py3/checkGridSurvey.py -u %s' % (EBList[0])
            #os.system(text_sd)
            print(text_sd)
#
                

'''



newList, doneList, failList = [],[],[]
for entry in UIDentry:
    if entry.split()[0][0] == '#':   continue   # Comment out byt #
    if entry.split()[0] == 'done': doneList = doneList + [entry.split()[1]]
    if entry.split()[0] == 'fail': failList = failList + [entry.split()[1]]
    if entry.split()[0] == 'new':  newList  = newList  + [entry.split()[1]]
#
for newUID in newList:
    prefix = newUID.replace("/", "_").replace(":","_").replace(" ","")



if len(newList) == 0: os.system('sleep 3600')
newUID = newList[0]
newUID = newUID.rstrip('\n')
prefix = newUID.replace("/", "_").replace(":","_").replace(" ","")
errList = glob.glob('%s.*.trial' % (prefix))
errList.sort()
errCount = len(errList)
if errCount > 0: errCount = int(errList[-1].split('.')[1]) + 1
print('%s error count = %d' % (prefix, errCount))
os.system('touch %s.%d.trial' % (prefix, errCount))
if errCount > 3:    # Failed
    text_sd = replaceList(newUID, 'fail')
    os.system(text_sd)
    os.system('mv ./UID/UIDnew ./UID/UIDList')
#
#if errCount > 1:   # Clean ASDM and MS
#    os.system('rm -rf %s*' % (prefix))
#    os.system('rm -rf UID/%s*' % (prefix))
#
if errCount <= 3:  # Prepare ASDM and MS
    if not os.path.isdir(prefix + '.ms'):
        if not os.path.isdir(prefix):
            if not os.path.isdir('UID/' + prefix):
                os.system('asdmExport %s' % (prefix))
                os.system('mv %s UID/' % (prefix))
            #
            os.system('ln -s UID/' + prefix + ' .')
        #
        importasdm(prefix)
    #
if not os.path.isdir(prefix + '.ms'):
    print('No available ASDM %s' % (prefix))
    text_sd = replaceList(newUID, 'nada')
    os.system(text_sd)
    os.system('mv ./UID/UIDnew ./UID/UIDList')
else:
    text_sd = generetaCheckEach(prefix)
    os.system(text_sd)
    os.system('touch %s.go' % (prefix))
#
'''
