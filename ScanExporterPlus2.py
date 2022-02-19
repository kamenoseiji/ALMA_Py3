#!/usr/bin/env python

from __future__ import print_function
import os, sys, email, struct
from optparse import OptionParser
from xml.etree import ElementTree

def getSpectralWindowDic(asdm):
    """
    Return {spwid: spwtype} dictionary for Spectral Window table of given ASDM.
    """
    tree = ElementTree.parse(os.path.join(asdm, 'SpectralWindow.xml'))
    elem = tree.getroot()
    spws = {}
    for row in elem.findall('row'):
        for spwEntry in row.findall('spectralWindowId'): spwID = spwEntry.text.strip()
        for spwEntry in row.findall('name'): name = spwEntry.text.strip()
        spwtype = 'UNKNOWN'
        for t in ['SQLD', 'WVR', 'FULL_RES', 'CH_AVG']:
            if t in name: spwtype = t
        #
        spws[spwID] = spwtype
    #
    return spws

def getDataDescriptionDic(asdm):
    """
    Return {ddid: spwid} dictionary for Data Description table of given ASDM.
    """
    tree = ElementTree.parse(os.path.join(asdm, 'DataDescription.xml'))
    elem = tree.getroot()
    dds = {}
    for row in elem.findall('row'):
        for entry in row.findall('dataDescriptionId'): ddid = entry.text.strip()
        for entry in row.findall('spectralWindowId'): spwID = entry.text.strip()
        dds[ddid] = spwID
    return dds

def getConfigDescriptionDic(asdm):
    """
    Return {cdid: spwid} dictionary for Config Description table of given ASDM.
    """
    tree = ElementTree.parse(os.path.join(asdm, 'ConfigDescription.xml'))
    elem = tree.getroot()
    cds = {}
    for row in elem.findall('row'):
        for entry in row.findall('configDescriptionId'): cdid = entry.text.strip()
        for entry in row.findall('dataDescriptionId'): ddid = entry.text.strip()
        cds[cdid] = ddid.split()[2:]
    return cds

def getScanDic(asdm):
    """
    Return {scan: intent} dictionary for Sscan table of given ASDM.
    """
    tree = ElementTree.parse(os.path.join(asdm, 'Scan.xml'))
    elem = tree.getroot()
    scans = {}
    for row in elem.findall('row'):
        for entry in row.findall('scanNumber'): sn = int(entry.text.strip())
        for entry in row.findall('scanIntent'): si = entry.text.strip()
        if sn not in list(scans.keys()): scans[sn] = ''
        scans[sn] = si
    return scans

def findBDF(asdm, scans, confdescs):
    """
    Find and return the BDF name for given ASDM, scan, intent, & configdescid.
    """
    tree = ElementTree.parse(os.path.join(asdm, 'Main.xml'))
    elem = tree.getroot()
    bdfs = []
    for row in elem.findall('row'):
        for entry in row.findall('scanNumber'): sn = int(entry.text.strip())
        for entry in row.findall('configDescriptionId'): cdid = entry.text.strip()
        for entry in row.findall('dataUID'):
            for dataUIDentry in entry.findall('EntityRef'): entityId = dataUIDentry.get('entityId')
        #
        if sn in scans and cdid in confdescs: bdfs.append(entityId)
    return bdfs

def downloadBinaries(asdm, uids):
    for uid in uids:
        print('Exporting: %s ' % uid)
        url = 'http://ngas01.sco.alma.cl:7777/RETRIEVE?file_id=%s' % (uid.replace('uid://', ''))
        outfile = '%s/ASDMBinary/%s' % (asdm, uid.replace(':', '_').replace('/', '_'))
        if not os.path.isfile(outfile):
            command = 'wget %s -O %s' % (url, outfile)
            print(command)
            os.system(command)


parser = OptionParser()
parser.add_option('-u', dest='asdm', metavar='ASDM',
    help='EB UID.  Default: none', default='')
parser.add_option('-t', dest='datatype', metavar='DataType',
    help='Data type (WVR, SQLD, FULL_RES, CH_AVG), or their comma-concatenation.  Default: "" (=all)', default='')
parser.add_option('-s', dest='scan', metavar='Scan',
    help='Scan number, or their comma-concatenation.  Default: "" (=all)', default='')
parser.add_option('-i', dest='intent', metavar='Intent',
    help='Scan intent (or its substring), or their comma-concatenation.  Default: "" (=all)', default='')
(options, args) = parser.parse_args()

if options.asdm=='':
    print('-f (ASDM) must be specified.')
    sys.exit()

asdm = options.asdm.replace(':', '_').replace('/', '_')

if not os.path.exists(asdm):
    os.system('asdmExport -m %s' % asdm)
if not os.path.exists(os.path.join(asdm, 'ASDM.xml')):
    print('"%s" does not seem like an ASDM.' % asdm)
    sys.exit()
spwdic = getSpectralWindowDic(asdm)
datadescdic = getDataDescriptionDic(asdm)
confdescdic = getConfigDescriptionDic(asdm)
scandic = getScanDic(asdm)

if len(options.scan)==0:
    scans_tmp = []
else:
    scans_tmp = [int(i) for i in options.scan.strip().split(',')]

if len(options.intent)==0:
    intents = []
else:
    intents = options.intent.strip().split(',')

if len(options.datatype)==0:
    datatypes = []
else:
    datatypes = options.datatype.strip().split(',')

confdescs = []
for cdi in list(confdescdic.keys()):
    for ddi in confdescdic[cdi]:
        if spwdic[datadescdic[ddi]] in datatypes or len(datatypes)==0:
            confdescs.append(cdi)

scans = []
for scannum in list(scandic.keys()):
    if len(intents)==0:
        intentmatches = True
    else:
        intentmatches = False
        for ii in intents:
            intentmatches = intentmatches or ii in scandic[scannum]
    if (scannum in scans_tmp or len(scans_tmp)==0) and intentmatches:
        scans.append(scannum)

bdfs = findBDF(asdm, scans, confdescs)
downloadBinaries(asdm, bdfs)
