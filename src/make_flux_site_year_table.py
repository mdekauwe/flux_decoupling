#!/usr/bin/env python

"""
Compile a list of sites & site years used in the final analysis

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
from tabulate import tabulate

def main():

    fdir = "data/processed"
    f = codecs.open(os.path.join(fdir, "omega_fluxnet_screened_PM.csv"), "r",
                    encoding='utf-8', errors='ignore')
    df1 = pd.read_csv(f)

    # Actual total excludes mixed forests & wetlands
    pfts = ['ENF','EBF','DBF','TRF','SAV','SHB','GRA','C3C','C4C']
    df1 = df1[df1.PFT.isin(pfts)]

    # Exclude sites where omega was bunk"
    df1 = df1[df1.omega >= 0]

    # Need to get long names as well!
    fdir = "data/raw_data/LaThuile_fluxnet_data/ancillary_files/csv/raw"
    f = codecs.open(os.path.join(fdir, "SummaryAnc_txt.csv"), "r",
                    encoding='utf-8', errors='ignore')
    df2 = pd.read_csv(f)


    for site in np.unique(df1.site):
        site_yrs =  df1[df1.site == site].year.values
        site_yrs = " ".join(str(x) for x in site_yrs)
        #site_yrs = site_yrs.replace(" ", " & ")

        long_name =  df2[df2.sitename == site].longname.values
        long_name = " ".join(str(x) for x in long_name)
        long_name = long_name.replace("&", "and")
        long_name = long_name.replace("_", " ")
        long_name = long_name.replace("#", "")

        #re.sub(r"&(?!#\d{4};)", "and", long_name)
        #print(site, "&", long_name, "&", site_yrs, "\\\\")
        print(site, "&", site_yrs, "\\\\")

        #site_yrs = site_yrs.replace(" ", ",")
        #print("%s,%s" % (site, site_yrs))


if __name__ == "__main__":

    main()
