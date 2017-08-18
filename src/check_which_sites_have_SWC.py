#!/usr/bin/env python
"""
For R2, check which sites have SWC?

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (15.08.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import calendar
import datetime as dt
from rmse import rmse

def main():

    fdir = "data/raw_data/LaThuile_fluxnet_data/raw_data"
    flist = glob.glob(os.path.join(fdir, "*.csv"))

    for i, fname in enumerate(flist):
        site = os.path.basename(fname).split(".")[0]
        yr = os.path.basename(fname).split(".")[1]

        df = pd.read_csv(fname)

        df1 = df.copy()
        df2 = df.copy()

        df1 = df1[ (df1['SWC1_fqc'] == 1) ]
        df2 = df2[ (df2['SWC2_fqc'] == 1) ]

        df1 = df1[ (df1['SWC1_f'] > 0.0) ]
        df2 = df2[ (df2['SWC2_f'] > 0.0) ]

        if len(df1.SWC1_f) > 50 or len(df2.SWC2_f) > 50:
            print "%s,%s,%d,%d,%.1f,%.1f" % \
                    (site, yr, len(df1.SWC1_f), len(df2.SWC2_f),
                     len(df1.SWC1_f) / float(len(df)) * 100.0,
                     len(df2.SWC2_f) / float(len(df)) * 100.0)




def date_converter(self, *args):
    year = int(float(args[0]))
    doy = int(float(args[1]))
    # in leap years the rogue date from the following year will be 367
    # as we are correctin this below it really doesn't matter what we set
    # it to but lets make it 366 for now so that the code doesn't barf
    if doy == 367:
        doy = 366

    hour = int(args[2].split(".")[0])
    minutes = int((float(args[2]) - hour) * 60.)
    date = "%s %s %s:%s" % (year, doy, hour, minutes)

    return pd.datetime.strptime(date, '%Y %j %H:%M')

def fix_rogue_date(self, df, drop=False):
    files_year = np.median(df.index.year)

    if drop:
        # drop 30 min slot from following year
        df = df[df.index.year == files_year]
    else:
        # files contain a rouge date from the following year, fix it.
        dates = pd.Series([pd.to_datetime(date) \
                           for date in df.index]).values
        fixed_date = "%s %s %s:%s" % (int(files_year + 1), 1, 0, 0)
        dates[-1] = pd.datetime.strptime(fixed_date, '%Y %j %H:%M')
        df = df.reindex(dates)

    return df




if __name__ == "__main__":

    main()
