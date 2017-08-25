#!/usr/bin/env python

"""
Run a filter over omega fits and remove the "bad" sites, see definitions
below

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (05.11.2015)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import codecs

f = codecs.open("data/processed/omega_fluxnet_PM.csv", "r",
                encoding='utf-8', errors='ignore')
df1 = pd.read_csv(f)

f = codecs.open("data/processed/omega_fluxnet_screened_PM.csv", "r",
                encoding='utf-8', errors='ignore')
df2 = pd.read_csv(f)

print(len(np.unique(df1.site)), len(np.unique(df2.site)))
print(len(df1), len(df2))

print("But actual total excludes mixed forests & wetlands")

pfts = ['ENF','EBF','DBF','TRF','SAV','SHB','GRA','C3C','C4C']

print(len(pfts))
print(len(np.unique(df2.PFT)))

print(np.unique(df2.PFT), len(df2[df2.PFT == "WET"]), len(df2[df2.PFT == "MF"]), len(df2[df2.PFT == " TBD"]))
df2 = df2[df2.PFT.isin(pfts)]
print(np.unique(df2.PFT), len(np.unique(df2.PFT)))

print(len(df2), len(np.unique(df2.site)))

print("Exclude sites where omega was bunk")
df2 = df2[df2.omega >= 0]
print(len(df2), len(np.unique(df2.site)))
