#!/usr/bin/env python

"""
Plot decoupling factor by PFT

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
import seaborn as sns
import brewer2mpl
import codecs
import matplotlib

def main():

    fdir = "data/processed"
    f = codecs.open(os.path.join(fdir, "omega_fluxnet_screened_PM.csv"), "r",
                    encoding='utf-8', errors='ignore')
    df = pd.read_csv(f)
    #print(len(df))
    #print(len(np.unique(df.site)))
    df = df[df.omega >= 0]
    #print(len(df))
    #print(len(np.unique(df.site)))

    #pfts = ['ENF','EBF','DBF','TRF','SAV','SHB','C3G','C3C', 'C4C']
    pfts = ['ENF','EBF','DBF','TRF','SAV','SHB','GRA','C3C', 'C4C']

    # Get literature review data.
    fdir = "data/lit_review"
    fname = os.path.join(fdir, 'decoupling_database.xlsx')
    df_lit = pd.read_excel(fname)
    df_lit = df_lit[df_lit.PFT.isin(pfts)]

    df['type']='FLUXNET'
    df_lit['type']='Literature'
    dfs = pd.concat([df, df_lit])


    sns.set_style("ticks")
    sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
    set2 = sns.color_palette("Set3", 9)

    #colours = ["#97adee"]*8
    colours = ["#E4E6D7"]*8

    fig = plt.figure(figsize=(9,4.5))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"

    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    m1 = df.groupby(['PFT'])['omega'].median()
    print( m1 )

    print("\n")

    m2 = df_lit.groupby(['PFT'])['omega'].median()
    print( m2 )

    ax = fig.add_subplot(111)

    # define outlier properties
    flierprops = dict(marker='o', markersize=3, markerfacecolor="black")
    ax = sns.boxplot(x="PFT", y="omega", hue="type", data=dfs, palette="PRGn",
                     order=pfts, flierprops=flierprops)

    for i,p in enumerate(pfts):

        if i == 0:
            offset = 0


        if p != "C4G":
            plt.text(-0.4+offset, 1.0-0.1,
                     "n=%d/%d\n" % (len(df[df.PFT  == p].omega), len(df_lit[df_lit.PFT  == p].omega)),
                     horizontalalignment='left', size='small', color="black")

        offset += 1.0


    ax.set_ylabel("$\Omega$ (-)")
    ax.set_xlabel("Plant functional type")

    ax.set_ylim(0, 0.9)

    ax.set_xticklabels(pfts)

    #ax.legend_.remove()
    #plt.legend(loc="upper right")

    l = ax.legend(loc="upper left")
    l.set_title('')

    odir = "/Users/%s/Dropbox/Decoupling_paper/figures/figs" % \
            (os.getlogin())
    plt.savefig(os.path.join(odir,
                #"Fluxnet_decoupling_boxplot_EBR_correction_excluding_sites_with_G.pdf"),
                "Fluxnet_decoupling_boxplot.pdf"),
                bbox_inches='tight', pad_inches=0.1)

    #plt.show()

if __name__ == "__main__":

    main()
