#!/usr/arch/bin/python
"""
does the ligand dbgen collection in case you messed it up the first time.
temp script that shows how to call the db compilation step.
"""

import os
import string
import sys
import disk_server

dbsplitsize = 3000 #temp
topdir = os.getcwd()

for dudtarget in os.listdir(topdir + "/dbgen/ligands"):
  gendir = topdir + "/dbgen/ligands/" + dudtarget
  dbname = topdir + "/dud_" + dudtarget + "_lig"
  disk_server.collect_results(dbname, gendir, dbsplitsize)
