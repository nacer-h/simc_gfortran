#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-07-26 13:23:57 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys,os,math

KIN = sys.argv[1]
SPEC = sys.argv[2]

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''
# Import package for cuts
from ltsep import Root

lt=Root(os.path.realpath(__file__))

# Add this to all files for more dynamic pathing
UTILPATH=lt.UTILPATH
SIMCPATH=lt.SIMCPATH

################################################################################################################################################

InputSIMC = "Heep_%s_%s" % (SPEC,KIN)

# Open inp_f file to grab prescale values and tracking efficiency
inp_f = SIMCPATH+"/input/%s.inp" % InputSIMC

with open(inp_f, 'r') as f:
    f_data = f.read()

with open(inp_f, 'r') as f:
    # Search for keywords, then save as value in dictionary
    for line in f:
        data = line.split('=')
        if 'Ebeam' in data[0]:
            if not 'dEbeam' in data[0]:
                ebeam = data[1].split(";")[0]
        if 'spec%e%P' in data[0]:
            eP = data[1].split(";")[0]
        if 'spec%e%theta' in data[0]:
            eTh = data[1].split(";")[0]
        if 'spec%p%P' in data[0]:
            pP = data[1].split(";")[0]
        if 'spec%p%theta' in data[0]:
            pTh = data[1].split(";")[0]

inpDict = {
    "ebeam" : float(ebeam)/1000, # Convert to GeV
    "eP" : float(eP)/1000,
    "eTh" : eTh,
    "pP" : float(pP)/1000,
    "pTh" : pTh,
}

BashPathEntry=("%s,%s,%s,%s,%s" % (inpDict["ebeam"],inpDict["eTh"],inpDict["eP"],inpDict["pTh"],inpDict["pP"]))
print(BashPathEntry)
