# ALMA_Py3
Scripts used for ALMA array diagnostics, mainly in Python 3, executable on CASA version 6 or later.

- How to download the codes:
  At preferred directory, type 

  git clone https://github.com/kamenoseiji/ALMA_Py3.git

  Then, scripts will be extracted in ALMA_Py3 directory.

- Requirements
  + CASA version 6 or later
  + analysisUtils
  + numpy
  + matplotlib

- Setup PYTHONPATH
  To import modules, add following line in your ~/.bashrc (tweak the real path to ALMA_Py3) and activate by source ~/.bashrc
  ```
  export PYTHONPATH=${PYTHONPATH}:/users/skameno/ALMA_Py3/
  ```

- How to use
  + Run with casa -c
  + Run with option -h to show help.
  + For example, running

casa -c ~/ALMA_Py3/checkBP.py -h

  will show:

Usage: checkBP.py [options]

Options:
  -h, --help   show this help message and exit
  -u prefix    EB UID   e.g. uid___A002_X10dadb6_X18e6
  -a antFlag   Antennas to flag e.g. DA41,DV08
  -b bunchNum  Channel binning
  -c scanID    Scan ID  e.g. 2
  -f           Apply flagging
  -m plotMin   Plot range minimum
  -M plotMax   Plot range minimum
  -s spwList   SPW List e.g. 0,1,2,3
  -P           Plot PDF
  -R refant    reference antenna

  + See Wiki page https://github.com/kamenoseiji/ALMA_Py3/wiki for delay.



