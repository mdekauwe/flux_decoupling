#!/usr/bin/env python

"""
Plot decoupling factor as a function of wind speed.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (17.04.2017)"
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

def main():

    fdir = "data/processed"
    f = codecs.open(os.path.join(fdir, "omega_fluxnet_screened_PM.csv"), "r",
                    encoding='utf-8', errors='ignore')
    df = pd.read_csv(f)
    df = df[df.omega >= 0]

    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    golden_mean = 0.6180339887498949
    width = 9
    height = width * golden_mean
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
    nonforest_PFTs = ['SAV', 'SHB', 'C3G','C3C']

    forest = df[df.PFT.isin(forest_PFTs)]
    nonforest = df[df.PFT.isin(nonforest_PFTs)]

    # screen by the odd very high wind speed
    forest = forest[forest.wind<6]
    nonforest = nonforest[nonforest.wind<6]

    ax1.plot(forest.wind, forest.omega, ls=" ", marker="o", color="black",
             label='Forest', alpha=0.7)

    #ax1.plot(nonforest.wind, nonforest.omega, ls=" ", marker="o",
    #         color="royalblue", label='Non-Forest', alpha=0.8)

    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(forest.wind, forest.omega)
    (m, b) = np.polyfit(forest.wind, forest.omega, 1)
    print(p_value, r_value, r_value**2)
    if p_value <= 0.05:
        ax1.plot(forest.wind, m*forest.wind+b, ls="-",
                 color="black")

    #lab = "y=%.2fx + %.2f\nr = %.2f; p < 0.001" % (m, b, r_value)
    lab = "r = %.2f" % (r_value)
    props = dict(boxstyle='round', facecolor='white', alpha=1.0,
                 ec="white")
    ax1.text(0.02, 0.95, lab, transform=ax1.transAxes,
             fontsize=12, verticalalignment='top',
             bbox=props)
    """
    from scipy.stats import linregress
    (slope, intercept,
     r_value, p_value, std_err) = linregress(nonforest.wind, nonforest.omega)
    (m, b) = np.polyfit(nonforest.wind, nonforest.omega, 1)
    print(p_value, r_value)
    if p_value <= 0.05:
        ax1.plot(nonforest.wind, m*nonforest.wind+b, ls="-",
                 color="black")
    """
    #ax1.legend(numpoints=1, ncol=1, loc="best", frameon=False)


    ax1.locator_params(nbins=7, axis="x")

    ax1.set_xlabel("Wind (m s$^{-1}$)")
    ax1.set_ylabel("$\Omega$ (-)")
    ax1.set_xlim(0, 6)

    odir = "/Users/%s/Dropbox/Decoupling_paper/figures/figs" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir, "omega_vs_wind.pdf"),
                    bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
