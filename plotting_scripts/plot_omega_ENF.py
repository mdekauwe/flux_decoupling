#!/usr/bin/env python

"""
Explore if data variability explains ENF decoupling factor range?

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
    df = df[df.PFT == 'ENF']
    sites = np.unique(df.site)

    good = []
    bad = []
    for s in sites:
        site = df[df.site == s]

        mu =  np.mean(site.omega)
        # coefficient of variation
        cov = np.std(site.omega) / mu * 100.

        if cov < 20:
            good.append(s)
        else:
            bad.append(s)

    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

    fig = plt.figure(figsize=(12,6))
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

    colour_list = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

    #colour_list = sns.palplot(sns.color_palette("colorblind", 10))
    #colour_list = sns.color_palette("Set2", 10)
    # CB palette  with grey:
    # from http://jfly.iam.u-tokyo.ac.jp/color/image/pallete.jpg
    #colour_list = ["#56B4E9", "#009E73", "#0072B2", "#F0E442",\
    #               "#E69F00", "#D55E00", "#CC79A7", "#999999"]
    colour_list = sns.color_palette("Accent", 10)

    ax1 = fig.add_subplot(111)

    good = []
    bad = []
    poor_samples = []
    for s in sites:
        site = df[df.site == s]
        nsamp = len(site.omega)

        mu =  np.mean(site.omega)
        # coefficient of variation
        cov = np.std(site.omega) / mu * 100.

        if nsamp <= 2:
            poor_samples.append(s)
        elif cov < 10:
            good.append(s)
        else:
            bad.append(s)
    xlabs = []
    xpos = []
    vals = np.zeros(0)
    xindex = np.arange(len(good))
    for i,s in enumerate(good):
        site = df[df.site == s]
        x = xindex[i]
        ax1.errorbar(x, np.mean(site.omega),
                     yerr=np.std(site.omega), ls=" ", marker="o",
                     color="black",
                     markeredgecolor="lightgrey",
                     alpha=0.8, capsize=False)

        #ax1.annotate(s, (x-0.3, 0.07),
        #             rotation=90, fontsize=8)

        vals = np.append(vals, site.omega.values)

        xlabs.append(s)
        xpos.append(x)

    ax1.axhline(y=np.mean(vals), xmin=0.04, xmax=xindex[-1]/60+0.035,
                ls="-", color="lightgrey")

    ax1.axvline(x=x+2.5, ls="--", color="grey")

    vals = np.zeros(0)
    xindex = (xindex[-1]+5) + np.arange(len(bad))
    for i,s in enumerate(bad):
        site = df[df.site == s]
        x = xindex[i]
        ax1.errorbar(x, np.mean(site.omega),
                     yerr=np.std(site.omega), ls=" ", marker="o",
                     color="black",
                     markeredgecolor="lightgrey",
                     alpha=0.8, capsize=False)
        #ax1.annotate(s, (x-0.3, 0.07),
        #             rotation=90, fontsize=8)
        vals = np.append(vals, site.omega.values)

        xlabs.append(s)
        xpos.append(x)

    ax1.axvline(x=x+2.5, ls="--", color="grey")

    ax1.annotate("(a)", (-1, 0.56), fontsize=12, fontweight="bold")
    ax1.annotate("(b)", (12, 0.56), fontsize=12, fontweight="bold")
    ax1.annotate("(c)", (45, 0.56), fontsize=12, fontweight="bold")
    ax1.axhline(y=np.mean(vals), xmin=xindex[0]/60+0.015, xmax=xindex[-1]/60-0.04,
                ls="-", color="lightgrey")


    vals = np.zeros(0)
    xindex = (xindex[-1]+5) + np.arange(len(poor_samples))
    for i,s in enumerate(poor_samples):
        site = df[df.site == s]
        x = xindex[i]
        ax1.errorbar(x, np.mean(site.omega),
                     yerr=np.std(site.omega), ls=" ", marker="o",
                     color="black",
                     markeredgecolor="lightgrey",
                     alpha=0.8, capsize=False)
        #ax1.annotate(s, (x-0.3, 0.1),
        #             rotation=90, fontsize=8, annotation_clip=True)

        xlabs.append(s)
        xpos.append(x)



        #trans = ax1.get_xaxis_transform() # x in data untis, y in axes fraction
        #ax1.annotate(s, xy=(10, -0.1), xycoords=trans)
        #plt.text(-0.3, 0.1, s, fontsize=8, transform=plt.gcf().transFigure)

        vals = np.append(vals, site.omega.values)

    ax1.set_xticklabels(xlabs, rotation=90, fontsize=8)
    ax1.set_xticks(xpos)

    ax1.axhline(y=np.mean(vals), xmin=xindex[0]/60-0.06, xmax=xindex[-1]/60-0.09,
                ls="-", color="lightgrey")
    ax1.set_ylim(0, 0.6)

    #plt.setp(ax1.get_xticklabels(), visible=False)

    #ax.set_xlim(-60, 90)

    #ax3.set_xlabel("Latitude (degrees)")
    ax1.set_ylabel("$\Omega$ (-)")
    #for ax in [ax1, ax2]:
    #    plt.setp(ax.get_xticklabels(), visible=False)

    odir = "/Users/%s/Dropbox/Decoupling_paper/figures/figs/" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir, "omega_ENF.pdf"),
                bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
