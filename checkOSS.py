#!/usr/local/bin/python
#-------- Usage 
# Step.1 copy acslog from ape1-gns.osf.alma.cl:/alma/acslogs/logs/APE1/SYSTEM/[YYYY-MM-DD]/log-APE1-[ONLINE SYSTEM VERSION]-YYYY-MM-DD-hh-mm-ss_YYYY-MM-DD-hh-mm-ss.xml.gz
# Step.2 extract using gzip -d *.gz
# Step.3 Run this script 
import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-t', dest='LogFile', metavar='LogFile',
    help='Log File Name   e.g. log-APE1-ONLINE-CYCLE10-B-34-2024-04-17-05-00-00_2024-08-10T09:33:48.004_2024-08-10T09:38:08.647.xml', default='')
parser.add_option('-a', dest='Array', metavar='Array',
    help='Array name      e.g. Array8-ACA', default='')
#
(options, args) = parser.parse_args()
#-------- Read log file
LogFileName = options.LogFile
ArrayName   = options.Array
keyword = '%s/CalibratorCatalog' % (ArrayName)
os.system("grep -e '%s' %s > temp.xml" % (keyword, LogFileName))
logfile = open('temp.xml')
logLines = logfile.readlines()
logfile.close()
for logLine in logLines:
    entryLine = logLine.rstrip(']]></Info>\n').rsplit(keyword + ']')[-1]
    print(entryLine)
#
