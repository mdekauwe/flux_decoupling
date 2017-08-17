#!/usr/bin/env python
"""
Fit omega from flux data.

NB. this code uses the python MPI package to speed up fits, so expect all your
cores to go to work when you run this...

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (06.03.2017)"
__email__ = "mdekauwe@gmail.com"

import matplotlib
matplotlib.use('agg') # stop windows popping up

import os
import sys
import glob
import pandas as pd
import numpy as np
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt
import calendar
import datetime as dt
from scipy.stats import pearsonr
from rmse import rmse
import scipy.stats as stats
import multiprocessing as mp
import re
import codecs

from penman_monteith import PenmanMonteith
from estimate_pressure import estimate_pressure

class FitOmega(object):
    """
    Fit Jarvis & McNaughton decoupling factor by inverting PM against
    fluxnet data using the lmfit package
    """

    def __init__(self, fdir, adir, co2dir, ofdir, site_fname, lai_fname,
                 global_co2_fname, ofname):

        self.flist = glob.glob(os.path.join(fdir, "*.csv"))
        self.site_fname = os.path.join(adir, site_fname)
        self.lai_fname = os.path.join(adir, lai_fname)
        self.global_co2_fname = os.path.join(co2dir, global_co2_fname)
        self.ofname = os.path.join(ofdir, ofname)

        self.out_cols = ["site","name","country","year","latitude",\
                         "longitude","PFT","n","ebr","omega","wind",\
                         "ga","gs","lai","canopy_ht","tower_ht","map",
                         "sprecip"]

        # W/m2 = 1000 (kg/m3) * 2.45 (MJ/kg) * 10^6 (J/kg) * 1 mm/day * \
        #        (1/86400) (day/s) * (1/1000) (mm/m)
        # 2.45 * 1E6 W/m2 = kg/m2/s or mm/s
        self.WM2_TO_KG_M2_S = 1.0 / ( 2.45 * 1E6 )
        self.KG_TO_G = 1000.0
        self.MOL_TO_MMOL = 1000.0
        self.G_TO_MOL_H20 = 1.0 / 18.0
        self.HPA_TO_KPA = 0.1
        self.KPA_TO_PA = 1000.0

    def main(self):

        #f = codecs.open(self.site_fname, "r",
        #                encoding='utf-8', errors='ignore')
        #df_site = pd.read_csv(f)
        #
        #f = codecs.open(self.lai_fname, "r",
        #                encoding='utf-8', errors='ignore')
        #df_lai = pd.read_csv(f)
        df_site = pd.read_csv(self.site_fname)
        df_lai = pd.read_csv(self.lai_fname)

        results = self.mpi_wrapper(df_site, df_lai)

        # merge all the processor DF fits into one big dataframe
        df = pd.concat(results, ignore_index=True)

        if os.path.exists(self.ofname):
            os.remove(self.ofname)
        df.to_csv(self.ofname, index=False)

    def mpi_wrapper(self, df_site, df_lai):

        # setup multiprocessor stuff
        num_processors = mp.cpu_count()
        chunk_size = int(np.ceil(len(self.flist) / float(num_processors)))
        pool = mp.Pool(processes=num_processors)
        queue = mp.Queue() # define an output queue

        # break up the files list equally between prcoessors, of course it won't
        # quite fit eqaully so account for this
        processes = []
        for i in range(num_processors):

            start = chunk_size * i
            end = chunk_size * (i + 1)
            if end > len(self.flist):
                end = len(self.flist)

            # setup a list of processes that we want to run
            p = mp.Process(target=self.worker,
                           args=(queue, self.flist[start:end], df_site,
                                 df_lai, i))
            processes.append(p)

        # Run processes
        for p in processes:
            p.start()

        # OS pipes are not infinitely long, so the queuing process can get
        # blocked when using the put function - a "deadlock". The following
        # code gets around this issue

        # Get process results from the output queue
        results = []
        while True:
            running = any(p.is_alive() for p in processes)
            while not queue.empty():
                results.append(queue.get())
            if not running:
                break

        # Exit the completed processes
        for p in processes:
            p.join()

        return results

    def worker(self, output, flist, df_site, df_lai, processor_number):

        df_out = pd.DataFrame(columns=self.out_cols)

        for i, fname in enumerate(flist):

            d = self.get_site_info(df_site, df_lai, fname)
            print("%s %s" % (d['site'], d['yr']))

            #f = codecs.open(fname, "r",
            #                encoding='utf-8', errors='ignore')
            #df = pd.read_csv(f, index_col='date',
            #                 parse_dates={'date': ["Year","DoY","Time"]},
            #                 date_parser=self.date_converter)

            df = pd.read_csv(fname, index_col='date',
                             parse_dates={'date': ["Year","DoY","Time"]},
                             date_parser=self.date_converter)

            # files contain a rouge date from the following year, fix it.
            df = self.fix_rogue_date(df, drop=True)

            (months, missing_nee) = self.get_three_most_productive_months(df)
            if missing_nee:
                df_out = self.write_bad_row(df_out, d)
                continue

            # kPa
            df.loc[:, 'VPD_f'] *= self.HPA_TO_KPA

            # mol m-2 s-1
            conv = self.WM2_TO_KG_M2_S * self.KG_TO_G * self.G_TO_MOL_H20
            df["ET"] = df['LE_f'] * conv

            # Calculate gs using pressure estimated from elevation if
            # possible, if not we will use a standard pressure assumption
            df,d = self.filter_dataframe(df, d, months)

            if 'gs_est' not in df:
                df_out = self.write_bad_row(df_out, d)
            else:
                # Filter extreme omega are ridiculous
                extreme = df['omega'].mean() + (3.0 * df['omega'].std())
                df = df[df['omega'] < extreme]

                if (len(df['gs_est']) > 20 and
                    d['pft'] != "TBD" and
                    d['lat'] != "TBD" and
                    d['pft'] != "WET"):

                    # some columns have strings, force to float otherwise we
                    # will produce a bug when we attempt to fit this bad data
                    # that will be screend anyway
                    df['VPD_f'] = df['VPD_f'].astype('float64')
                    df['CO2'] = df['CO2'].astype('float64')
                    df['GPP_f'] = df['GPP_f'].astype('float64')

                    self.make_plot(F, d, df)

                    # Add a tropics class
                    if self.is_tropics(d):
                        d['pft'] = "TropRF"
                        df_out = self.write_row(df_out, d)
                    else:
                        df_out = self.write_row(df_out, d)

                else:
                    df_out = self.write_bad_row(df_out, d)

        output.put(df_out)

    def is_tropics(self, d):
        tropics = False
        if ((d['clim_grp'] == "Tropical" and d['pft'] == "DBF") or
            (d['clim_grp'] == "Tropical" and d['pft'] == "EBF") or
            (d['clim_grp'] == "Tropical" and d['pft'] == "ENF") or
            (d['lat'] >= -23.43 and d['lat'] <= 23.43 and d['pft'] == "DBF") or
            (d['lat'] >= -23.43 and d['lat'] <= 23.43 and d['pft'] == "EBF") or
            (d['lat'] >= -23.43 and d['lat'] <= 23.43 and d['pft'] == "ENF")):
            tropics = True

        return tropics

    def get_site_info(self, df_site, df_lai, fname):

        d = {}
        s = os.path.basename(fname).split(".")[0]
        d['site'] = s
        d['yr'] = os.path.basename(fname).split(".")[1]
        d['lat'] = df_site.loc[df_site.Site_ID == s,'Latitude'].values[0]
        d['lon'] = df_site.loc[df_site.Site_ID == s,'Longitude'].values[0]
        d['pft'] = df_site.loc[df_site.Site_ID == s,'IGBP_class'].values[0]
        d['clim_cl'] = df_site.loc[df_site.Site_ID == s,
                                   'Climate_class'].values[0]
        d['clim_grp'] = df_site.loc[df_site.Site_ID == s,
                                    'Climate_group'].values[0]

        # remove commas from country tag as it messes out csv output
        name = df_site.loc[df_site.Site_ID == s,'Name'].values[0]
        d['name'] = name.replace("," ,"")
        d['country'] = df_site.loc[df_site.Site_ID == s,'Country'].values[0]
        d['elev'] = df_site.loc[df_site.Site_ID == s,'Elevation'].values[0]
        d['lai'] = df_lai.loc[df_lai.sitename == s,'LAI'].values[0]
        d['lai_max'] = df_lai.loc[df_lai.sitename == s,'LAI_MAX'].values[0]
        d['lai_min'] = df_lai.loc[df_lai.sitename == s,'LAI_MIN'].values[0]

        fdir = "data/raw_data/anna_meta"
        f = codecs.open(os.path.join(fdir, "site_metadata.csv"), "r",
                        encoding='utf-8', errors='ignore')
        df_m = pd.read_csv(f)

        d['tower_ht'] = -999.9
        try:
            ht = df_m.loc[df_m.SiteCode == s, 'TowerHeight'].values[0]
            if ~np.isnan(ht):
                d['tower_ht'] = ht
        except IndexError:
            pass

        d['canopy_ht'] = -999.9
        try:
            ht = df_m.loc[df_m.SiteCode == s, 'CanopyHeight'].values[0]
            if ~np.isnan(ht):
                d['canopy_ht'] = ht
        except IndexError:
            pass

        return (d)

    def write_row(self, df_out, d):

        row = pd.Series([d['site'], d['name'], d['country'], d['yr'],
                         d['lat'], d['lon'], d['pft'], d['num_pts'],
                         d['EBR'], d['omega'], d['wind'], d['ga'],
                         d['gs'], d['lai'], d['tower_ht'], d['canopy_ht'],
                         d['ppt'], d['sprecip']],
                         index=self.out_cols)


        result = df_out.append(row, ignore_index=True)
        return result

    def write_bad_row(self, df_out, d):
        row = pd.Series([d['site'], d['name'], d['country'], d['yr'],
                         d['lat'], d['lon'], d['pft'], -999.9,
                         -999.9, -999.9, -999.9, -999.9, -999.9,
                         -999.9, -999.9, -999.9, -999.9, -999.9],
                         index=self.out_cols)

        result = df_out.append(row, ignore_index=True)
        return result

    def make_plot(self, F, d, df):

        fig = plt.figure(figsize=(16,6))
        fig.subplots_adjust(wspace=0.05)
        fig.subplots_adjust(hspace=0.3)
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = "sans-serif"
        plt.rcParams['font.sans-serif'] = "Helvetica"
        plt.rcParams['axes.labelsize'] = 12
        plt.rcParams['font.size'] = 12
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12

        almost_black = '#262626'
        # change the tick colors also to the almost black
        plt.rcParams['ytick.color'] = almost_black
        plt.rcParams['xtick.color'] = almost_black

        # change the text colors also to the almost black
        plt.rcParams['text.color'] = almost_black

        # Change the default axis colors from black to a slightly lighter black,
        # and a little thinner (0.5 instead of 1)
        plt.rcParams['axes.edgecolor'] = almost_black
        plt.rcParams['axes.labelcolor'] = almost_black

        ax1 = fig.add_subplot(111)


        ax1.plot(df["omega"], "ko", alpha=0.3)
        ax1.set_ylim(0, 1)

        ax1.set_ylabel("$\Omega$ (-)")

        fig.suptitle("%s:%s: $\Omega$ = %.2f; N = %d" %
                     (d['site'], d['yr'], d['omega'], d['num_pts']),
                     x=0.25)

        ax1.locator_params(nbins=6)

        fig.savefig("plots/omega_fits/%s_%s_omega.pdf" % (d['site'], d['yr']),
                    bbox_inches='tight', pad_inches=0.1)
        #plt.show()
        plt.close(fig)

    def get_three_most_productive_months(self, df):
        # filter three most productive months
        df_m = df.resample("M").mean()
        missing_nee = False

        df_m["NEP"] = np.where(df_m.NEE_f > -1000.,
                               df_m.NEE_f * -1., df_m.NEE_f)
        try:
            df_m = df_m.sort_values("NEP", ascending=False)#[:3]
            months = df_m.index.month
        except KeyError:
            missing_nee = True
            months = None

        return (months, missing_nee)

    def filter_dataframe(self, df, d, months):
        # Only trusted half hourly data are used, defined as original data or
        # those empirically modeled with a high degree of confidence
        # (quality control flag fqcOK = 1, see www.fluxdata.org for details
        # and definition, and Reichstein et al. [2005, Appendix A]). From
        # Williams WATER RESOURCES RESEARCH, VOL. 48, W06523,
        # doi:10.1029/2011WR011586, 2012

        # For the fully coupled this is just filler, it gets set in the PM
        # version
        d['omega'] = -999.9
        no_G = False

        df_p = df.copy()
        total_length = len(df_p)
        df_p = df_p[df_p['Precip_fqcOK'] == 1]
        df_y = df_p.groupby(df_p.index.year).sum()

        if len(df_p) == 0:
            pptx = -999.9
        else:
            pptx = df_y.Precip_f.values[0]
            if pptx < 0.0:
                pptx = -999.9
            elif len(df_p) < total_length * 0.9:
                #print("%.2f " % (len(df_p)/float(total_length)*100.0))
                pptx = -999.9

        d['ppt'] = pptx

        # filter three most productive months
        df = df[(df.index.month == months[0]) |
                (df.index.month == months[1]) |
                (df.index.month == months[2])]

        # save months
        d['most_prod_mths'] = str(months).strip('[]')

        # NB, I'm not slicing the df here, setting this to a copy
        df_bal = df[(df['LE_fqcOK'] == 1) &
                    (df['H_fqcOK'] == 1) &
                    (df['Rn_fqcOK'] == 1) &
                    (df['G_fqcOK'] == 1)] # soil heat flux

        top = np.sum(df_bal.LE_f + df_bal.H_f)
        bot = np.sum(df_bal.Rn_f - df_bal.G_f)
        #if bot > 0.0:
        #    d['EBR'] = top / bot
        #else:
        #    d['EBR'] = -999.9

        # figure out 3 mth precip
        df_sp = df.copy()
        df_sp = df_sp[df_sp['Precip_fqcOK'] == 1]
        df_sp = df_sp.groupby(df_sp.index.year).sum()

        if len(df_sp) == 0:
            pptxx = -999.9
        else:
            pptxx = df_sp.Precip_f.values[0]
            if pptx < 0.0:
                pptxx = -999.9
            elif len(df_sp) < total_length * 0.9:
                pptxx = -999.9

        d['sprecip'] = pptxx

        # filter daylight hours, good LE data, GPP, CO2
        #
        # If we have no ground heat flux, just use Rn
        if len(df[(df['G_fqcOK'] == 1)]) == 0:
            df = df[(df.index.hour >= 8) &
                    (df.index.hour <= 16) &
                    (df['LE_fqcOK'] == 1) &
                    (df['H_fqcOK'] == 1) &
                    (df['VPD_fqcOK'] == 1) &
                    (df['Precip_fqcOK'] == 1) &
                    (df['WS_fqcOK'] == 1) &
                    (df['ET'] > 0.01 / 1000.) & # check in mmol, but units are mol
                    (df['VPD_f'] > 0.05) &
                    (df['Rn_fqcOK'] == 1) &
                    (df['Ta_fqcOK'] == 1)]
            no_G = True
            d['EBR'] = -999.9
        else:
            df = df[(df.index.hour >= 8) &
                    (df.index.hour <= 16) &
                    (df['LE_fqcOK'] == 1) &
                    (df['H_fqcOK'] == 1) &
                    (df['VPD_fqcOK'] == 1) &
                    (df['Precip_fqcOK'] == 1) &
                    (df['WS_fqcOK'] == 1) &
                    (df['VPD_f'] > 0.05) &
                    (df['G_fqcOK'] == 1) &
                    (df['Rn_fqcOK'] == 1) &
                    (df['Ta_fqcOK'] == 1)]


            # Correct based on method 4 from Wohlfahrt et al. Agricultural and
            # Forest Meteorology 149
            #if top > 0.0:
            #    corection_factor = bot/top
            #    df.LE_f *= corection_factor
            #    df.H_f *= corection_factor
            #
            #    if bot > 0.0:
            #        d['EBR'] = top / bot
            #else:
            #    d['EBR'] = -999.9

            if bot > 0.0:
                d['EBR'] = top / bot
            else:
                d['EBR'] = -999.9
        #df_rain = df[df.Precip_f > 0.0]
        idx = df[df.Precip_f > 0.0].index.tolist()

        bad_dates = []
        for rain_idx in idx:
            bad_dates.append(rain_idx)
            for i in range(48):
                new_idx = rain_idx + dt.timedelta(minutes=30)
                bad_dates.append(new_idx)
                rain_idx = new_idx

        # There will be duplicate dates most likely so remove these.
        bad_dates = np.unique(bad_dates)

        # remove rain days...
        df = df[~df.index.isin(bad_dates)]

        # Estimate gs from inverting the penman-monteith
        df,d = self.penman_montieth_wrapper(d, df, no_G)

        # screen for bad data, or data I've set to bad
        df = df[(df['gs_est'] > 0.0) & (np.isnan(df['gs_est']) == False)]

        # screen data with low wind speed
        #df = df[df['WS_f'] > 1.0]

        return (df,d)

    def penman_montieth_wrapper(self, d, df, no_G):

        # screen by low u*, i.e. conditions which are often indicative of
        # poorly developed turbulence, after Sanchez et al. 2010, HESS, 14,
        # 1487-1497. Some authors use 0.3 m s-1 (Oliphant et al. 2004) or
        # 0.35 m s-1 (Barr et al. 2006) as a threshold for u*
        df = df[df.ustar>= 0.25]

        # convert all units...
        vpd = df['VPD_f'] * self.KPA_TO_PA # Pa
        wind = df['WS_f'] # m s-1
        rnet = df['Rn_f'] # W m-2
        tair = df['Ta_f'] # deg C
        ustar = df['ustar'] # frictional velocity, m s-1

        if len(df.G_f) > 0:
            G = df['G_f'] # W m-2
        else:
            G = None
        """
        if no_G:
            if len(rnet) > 0 and len(df['LE_f']) > 0 and len(df['H_f']) > 0:
                #Soil heat flux (G) assumed as the difference between Rn and
                #the sum of latent and sensible heat fluxes
                G = df['Rn_f'] - df['LE_f'] - df['H_f']
            else:
                G = None
        """
        if no_G:
            G = None



        #us = ustar[ustar> -5000.]
        #print(np.nanmean(us), np.max(us), np.min(us))



        # This is a pain, but elevation could be a float or a string
        # and we need to try and catch both possibilities
        if isinstance(d['elev'], float) and d['elev'] > -500.0:
            press = estimate_pressure(df['Ta_f'], float(d['elev']))

        # some elevations have dates or a dash e.g. 450-570
        elif (isinstance(d['elev'], float) == False and
              d['elev'] != "TBD" and
              any(c.isalpha() for c in d['elev']) == False and
              re.search(r"-", d['elev']) == False): # some elevations
            press = estimate_pressure(df['Ta_f'], float(d['elev']))

        else:
            press = 101135.0

        PM = PenmanMonteith(use_ustar=True)
        new_press = estimate_pressure(tair, float(d['elev']))
        lambdax = PM.calc_latent_heat_of_vapourisation(tair)

        # W m-2 -> mol m-2 s-1
        #trans = df['LE_f'] / lambdax
        conv = self.WM2_TO_KG_M2_S * self.KG_TO_G * self.G_TO_MOL_H20
        trans = df['LE_f'] * conv
        df['gs_est']  = PM.invert_penman(vpd, wind, rnet, tair,
                                         press, trans,
                                         ustar=ustar, G=G)

        (omega,ga,gs) = PM.calc_decoupling_coefficent(wind, tair, press,
                                                      df["gs_est"], ustar=ustar)
        omega = omega[(omega >= 0.0) & (omega <= 1.0)]
        if len(omega) > 1:
            d['omega'] = np.nanmean(omega)
            #df['omega'] = omega
            df.loc[:, "omega"] = omega
        else:
            d['omega'] = -999.9
            #df['omega'] = np.ones(len(df)) * -999.9
            df.loc[:, "omega"] = np.ones(len(df)) * -999.9

        d['num_pts'] = len(df["gs_est"])
        d['wind'] = np.mean(wind[~np.isnan(wind)])

        if len(ga) > 1:
            d['ga'] = np.nanmean(ga)
        else:
            d['ga'] = -999.9

        if len(gs) > 1:
            d['gs'] = np.nanmean(gs)
        else:
            d['gs'] = -999.9

        return df, d

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

    def remove_from_list(self, flist, bad_files):
        return [i for i in flist if i not in bad_files]



if __name__ == "__main__":

    F = FitOmega(fdir="data/raw_data/LaThuile_fluxnet_data/raw_data",
                 adir="data/raw_data/LaThuile_fluxnet_data/ancillary_files/csv/raw/",
                 ofdir="data/processed/",
                 co2dir="data/raw_data/global_CO2_data/",
                 site_fname="CommonAnc_LATEST.csv",
                 lai_fname = "LAI_values.csv",
                 global_co2_fname="Global_CO2_mean_NOAA.csv",
                 ofname="omega_fluxnet_PM.csv")
    F.main()
