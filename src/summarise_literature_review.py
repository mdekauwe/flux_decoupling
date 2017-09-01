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

    print(len(np.unique(df['reference'])))

    print("\nInformation for table\n")

    for pft in pfts:
        p = df[df.PFT == pft]
        print(pft, "&", round(p.omega.mean(),2), "&", round(p.omega.std(),2),
              "&", round(p.omega.min(),2), "&", round(p.omega.max(),2), "&",  len(p), "\\")


if __name__ == "__main__":

    main()
