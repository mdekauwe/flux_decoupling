#!/usr/bin/env python

"""
Summarise the lit review

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (15.08.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import codecs
import numpy as np

def main():

    fdir = "data/lit_review"
    fname = os.path.join(fdir, 'decoupling_database.xlsx')
    df = pd.read_excel(fname)

    # Actual total excludes mixed forests & wetlands
    pfts = ['ENF','EBF','DBF','TRF','SAV','SHB','GRA','C3C','C4C']
    df = df[df.PFT.isin(pfts)]

    print(len(df))
    print(df['PFT'].value_counts())
    

if __name__ == "__main__":

    main()
