#!/usr/bin/env python

"""
Plot decoupling factor as a function of 3 most productive months of rainfall

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (28.04.2017)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import brewer2mpl
import seaborn as sns
import codecs
import matplotlib.patches as mpatches

def main():

    fdir = "data/processed"
    f = codecs.open(os.path.join(fdir, "omega_fluxnet_screened_PM.csv"), "r",
                    encoding='utf-8', errors='ignore')
    df = pd.read_csv(f)
    df = df[(df.omega >= 0) & (df.lai > 0.0)]

    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    #colour_list = ["#E69F00","#56B4E9", "#009E73", "#CC79A7"]
    colour_list = brewer2mpl.get_map('Set2', 'qualitative', 6).mpl_colors

    golden_mean = 0.6180339887498949
    width = 9
    height = width * golden_mean
    #width = 6*2*(1/golden_mean)

    fig = plt.figure(figsize=(width, height))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
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
    forest_PFTs = ['ENF', 'EBF', 'DBF', 'TRF']
    nonforest_PFTs = ['SHB', 'GRA','C3C', 'C4C']
    grasses_PFTs = ['GRA']


    forest = df[df.PFT.isin(forest_PFTs)]
    nonforest = df[df.PFT.isin(nonforest_PFTs)]
    grasses = df[df.PFT.isin(grasses_PFTs)]


    x = []
    y = []
    grasses = grasses[grasses.sprecip > 0.0]
    for site in np.unique(grasses.site):
        s = grasses[grasses.site == site]

        ax1.errorbar(np.mean(s.sprecip), np.mean(s.omega), ls=" ", marker="s",
                 color=colour_list[0], alpha=0.9, yerr=np.std(s.omega))
        x.append(np.mean(s.sprecip))
        y.append(np.mean(s.omega))
    x = np.asarray(x)
    y = np.asarray(y)

    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(x,y)
    (m, b) = np.polyfit(x,y, 1)
    print(p_value, r_value, r_value**2)
    if p_value <= 0.05:
        #ax1.plot(x, m*x+b, ls="-", color=colour_list[0], label="r = %.2f" % (r_value))
        ax1.plot(x, m*x+b, ls="-", color=colour_list[0])


    #forest = forest[forest.sprecip > 0.0]
    #for site in np.unique(forest.site):
    #    s = forest[forest.site == site]
    #    ax1.errorbar(np.mean(s.sprecip), np.mean(s.omega), ls=" ", marker="o",
    #             color="black", alpha=0.9, yerr=np.std(s.omega),
    #             xerr=np.std(s.sprecip))


    x = []
    y = []
    enf_PFTs = ['ENF']
    enf = df[df.PFT.isin(enf_PFTs)]
    enf = enf[enf.sprecip > 0.0]
    for site in np.unique(enf.site):
        s = enf[enf.site == site]
        ax1.errorbar(np.mean(s.sprecip), np.mean(s.omega), ls=" ", marker="o",
                 color=colour_list[1], alpha=0.9, yerr=np.std(s.omega))
        x.append(np.mean(s.sprecip))
        y.append(np.mean(s.omega))
    x = np.asarray(x)
    y = np.asarray(y)

    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(x,y)
    (m, b) = np.polyfit(x,y, 1)
    print(p_value, r_value, r_value**2)
    if p_value <= 0.05:
        #ax1.plot(x, m*x+b, ls="-", color=colour_list[1], label="r = %.2f" % (r_value))
        ax1.plot(x, m*x+b, ls="-", color=colour_list[1])


    x = []
    y = []
    ebf_PFTs = ['EBF']
    ebf = df[df.PFT.isin(ebf_PFTs)]
    ebf = ebf[ebf.sprecip > 0.0]
    for site in np.unique(ebf.site):
        s = ebf[ebf.site == site]
        ax1.errorbar(np.mean(s.sprecip), np.mean(s.omega), ls=" ", marker="p",
                 color=colour_list[2], alpha=0.9, yerr=np.std(s.omega))
        x.append(np.mean(s.sprecip))
        y.append(np.mean(s.omega))
    x = np.asarray(x)
    y = np.asarray(y)
    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(x,y)
    (m, b) = np.polyfit(x,y, 1)
    print(p_value, r_value, r_value**2)
    if p_value <= 0.05:
        #ax1.plot(x, m*x+b, ls="-", color=colour_list[2], label="r = %.2f" % (r_value))
        ax1.plot(x, m*x+b, ls="-", color=colour_list[2])


    x = []
    y = []
    dbf_PFTs = ['DBF']
    dbf = df[df.PFT.isin(dbf_PFTs)]
    dbf = dbf[dbf.sprecip > 0.0]
    for site in np.unique(dbf.site):
        s = dbf[dbf.site == site]
        ax1.errorbar(np.mean(s.sprecip), np.mean(s.omega), ls=" ", marker="^",
                 color=colour_list[3], alpha=0.9, yerr=np.std(s.omega))
        x.append(np.mean(s.sprecip))
        y.append(np.mean(s.omega))
    x = np.asarray(x)
    y = np.asarray(y)

    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(x,y)
    (m, b) = np.polyfit(x,y, 1)
    print(p_value, r_value, r_value**2)
    if p_value <= 0.05:
        #ax1.plot(x, m*x+b, ls="-", color=colour_list[3], label="r = %.2f" % (r_value))
        ax1.plot(x, m*x+b, ls="-", color=colour_list[3])



    x = []
    y = []
    trf_PFTs = ['TRF']
    trf = df[df.PFT.isin(trf_PFTs)]
    trf = trf[trf.sprecip > 0.0]
    for site in np.unique(trf.site):
        s = trf[trf.site == site]
        ax1.errorbar(np.mean(s.sprecip), np.mean(s.omega), ls=" ", marker="v",
                 color=colour_list[4], alpha=0.9, yerr=np.std(s.omega), zorder=0)
        x.append(np.mean(s.sprecip))
        y.append(np.mean(s.omega))
    x = np.asarray(x)
    y = np.asarray(y)

    #from scipy.stats import linregress
    #(slope, intercept,
    # r_value, p_value, std_err) = linregress(x,y)
    #(m, b) = np.polyfit(x,y, 1)
    #print(p_value, r_value, r_value**2)
    #if p_value <= 0.05:
    #    ax1.plot(x, m*x+b, ls="-", color=colour_list[3])

    # Add legned, the above is in a loop so we can't do it then
    ax1.errorbar(np.zeros(1)*np.nan, np.zeros(1)*np.nan, ls=" ", marker="s",
             color=colour_list[0], label='GRA ; r = 0.46', alpha=0.9, yerr=np.std(s.omega))


    ax1.errorbar(np.zeros(1)*np.nan, np.zeros(1)*np.nan, ls=" ", marker="o",
             color=colour_list[1], label='ENF ; r = 0.40', alpha=0.9, ms=8, yerr=np.std(s.omega))
    ax1.errorbar(np.zeros(1)*np.nan, np.zeros(1)*np.nan, ls=" ", marker="^",
             color=colour_list[3], label='DBF ; r = 0.64', alpha=0.9, ms=8, yerr=np.std(s.omega))
    ax1.errorbar(np.zeros(1)*np.nan, np.zeros(1)*np.nan, ls=" ", marker="p",
             color=colour_list[2], label='EBF', alpha=0.9, ms=8, yerr=np.std(s.omega))
    ax1.errorbar(np.zeros(1)*np.nan, np.zeros(1)*np.nan, ls=" ", marker="v",
             color=colour_list[4], label='TRF', alpha=0.9, ms=8, yerr=np.std(s.omega))

    """
    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(forest.lai,forest.omega)
    (m, b) = np.polyfit(forest.lai,forest.omega, 1)
    print(p_value, r_value, r_value**2)
    if p_value <= 0.05:
        ax1.plot(forest.lai, m*forest.lai+b, ls="-",
                 color="black")

    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(nonforest.lai,nonforest.omega)
    (m, b) = np.polyfit(nonforest.lai,nonforest.omega, 1)
    print(p_value, r_value, r_value**2)
    if p_value <= 0.05:
        ax1.plot(forest.lai, m*forest.lai+b, ls="-",
                 color="black")
    """
    #ax1.legend(numpoints=1, ncol=1, loc="best", frameon=False)


    ax1.locator_params(nbins=7, axis="x")

    ax1.set_xlabel("Precipitation (mm 3-month$^{-1}$)")
    ax1.set_ylabel("$\Omega$ (-)")
    #plt.setp(ax1.get_yticklabels(), visible=False)

    ax1.set_ylim(0, 0.8)
    ax1.set_xlim(0, 600)
    ax1.locator_params(nbins=3, axis="x")
    ax1.locator_params(nbins=6, axis="y")
    ax1.legend(numpoints=1, ncol=1, loc="upper left", frameon=False)

    odir = "/Users/%s/Dropbox/Decoupling_paper/figures/figs" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir, "omega_vs_sprecip.pdf"),
                    bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
