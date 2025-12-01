import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.backends.backend_pdf import PdfPages
from Plotters import plotStokesSpectra
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-u', dest='prefix', metavar='prefix',
    help='EB UID   e.g. uid___A002_X10dadb6_X18e6', default='')
parser.add_option('-S', dest='sourceList', metavar='sourceList',
    help='Source List', default='')
parser.add_option('-s', dest='spwList', metavar='spwList',
    help='SPW List e.g. 0,1,2,3', default='')
parser.add_option('-R', dest='refant', metavar='refant',
    help='Reference antenna e.g. DA42', default='')
(options, args) = parser.parse_args()
prefix  = options.prefix.replace("/", "_").replace(":","_").replace(" ","")
sourceList = [source for source in  options.sourceList.split(',')]
spwList = [int(spw) for spw in  options.spwList.split(',')]
refantName = options.refant
print('  Source      Frequency    Stokes I [Jy]  Stokes Q [Jy]    Stokes U [Jy]    Stokes V [Jy]    p [%]             EVPA [deg]')
for sourceName in sourceList:
    for spw in spwList:
        plotStokesSpectra(prefix, refantName, sourceName, spw)
    #
#
