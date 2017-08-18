#!/usr/bin/env python

"""
For R2, of the sites with SWC data, figure out how many match the sites we used
for our analysis.

We are also going to figure out how many sites had more than X% of actual data.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (18.08.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import brewer2mpl
import codecs
import matplotlib

def main():

    fdir = "data/processed"
    f = codecs.open(os.path.join(fdir, "omega_fluxnet_screened_PM.csv"), "r",
                    encoding='utf-8', errors='ignore')
    df1 = pd.read_csv(f)

    f = codecs.open(os.path.join(fdir, "sites_with_soil_water.csv"), "r",
                    encoding='utf-8', errors='ignore')
    df2 = pd.read_csv(f, names=["site", "year", "nSWC1", "nSWC2", "percent_SW1",
                                "percent_SW2"])

    df2 = df2[df2.percent_SW1 > 20.0]
    sites = np.unique(df2.site)

    x = []
    for site in np.unique(df1.site):
        if site in sites:
            print(site)
            x.append(site)
    print(len(x))

    x = []
    df2 = df2[df2.percent_SW2 > 20.0]
    sites = np.unique(df2.site)

    for site in np.unique(df1.site):
        if site in sites:
            print(site)
            x.append(site)
    print(len(x))


if __name__ == "__main__":

    main()
