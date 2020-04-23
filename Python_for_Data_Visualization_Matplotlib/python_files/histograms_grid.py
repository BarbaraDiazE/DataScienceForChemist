"""
Plot many histograms as figure's subplots.
Generaly, histograms for suplementary material.
"""
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
from scipy.stats import uniform
import scipy.special


def grid_histogram(n_bins):
    # define sigure and subplots
    fig, axes = plt.subplots(
        nrows=2,
        ncols=3,
        # sharex=True, sharey=True
    )
    #### histogram_one_descriptor###
    data = pd.read_csv(
        "/home/babs/Documents/DIFACQUIM/admetox_research/Libraries/Master_results.csv"
    )
    l = ["FDA", "BIOFACQUIM", "TCM", "AfroDB", "NuBBEDB"]
    colors = ["blueviolet", "green", "mediumvioletred", "darkorange", "dodgerblue"]
    ####
    # descriptor1 Silicos-IT LogSw
    ###
    X1 = list()
    X2 = list()
    X3 = list()
    X4 = list()
    X5 = list()
    # X6 = list()
    for i in l:
        d = data[data["Library"] == i]
        X1.append(d["Silicos-IT LogSw"].tolist())
        X2.append(d["Consensus Log P"].tolist())
        X3.append(d["Intestinal absorption"].tolist())
        X4.append(d["Fraction unbound"].tolist())
        X5.append(d["Total Clearance"].tolist())
        # X6.append(d["Silicos-IT LogSw"].tolist())
    label = ["FDA", "BIOFACQUIM", "TCM", "AfroDB", "NuBBEDB"]
    axes[0, 0].grid(color="silver", linestyle="-", linewidth=1.5)
    axes[0, 0].hist(X1, n_bins, histtype="bar", density=True, label=label, color=colors)
    axes[0, 0].legend()
    axes[0, 0].set_xlabel("Silicos-IT LogSw")
    axes[0, 0].set_ylabel("Fraction")
    axes[0, 0].set_ylim(0, 0.15)
    ###
    # descriptor2
    ###
    axes[0, 1].grid(color="silver", linestyle="-", linewidth=1.5)
    axes[0, 1].hist(X2, n_bins, histtype="bar", density=True, label=label, color=colors)
    axes[0, 1].legend()
    axes[0, 1].set_xlabel("Consensus Log P")
    axes[0, 1].set_ylabel("Fraction")
    ###
    # descriptor3
    ###
    axes[0, 2].grid(color="silver", linestyle="-", linewidth=1.5)
    axes[0, 2].hist(X3, n_bins, histtype="bar", density=True, label=label, color=colors)
    axes[0, 2].legend()
    axes[0, 2].set_xlabel("Intestinal absorption")
    axes[0, 2].set_ylabel("Fraction")
    ###
    # descriptor4
    ###
    axes[1, 0].grid(color="silver", linestyle="-", linewidth=1.5)
    axes[1, 0].hist(X4, n_bins, histtype="bar", density=True, label=label, color=colors)
    axes[1, 0].legend()
    axes[1, 0].set_xlabel("Fraction unbound")
    axes[1, 0].set_ylabel("Fraction")
    ###
    # descriptor5
    ###
    axes[1, 1].grid(color="silver", linestyle="-", linewidth=1.5)
    axes[1, 1].hist(X5, n_bins, histtype="bar", density=True, label=label, color=colors)
    axes[1, 1].legend()
    axes[1, 1].set_xlabel("Fraction unbound")
    axes[1, 1].set_xlabel("Total Clearance")

    # ###
    # descriptor6
    ###
    axes[1, 2].axis("off")

    ##titles
    axes[0, 0].set_title("A)", loc="left")
    axes[0, 1].set_title("B)", loc="left")
    axes[0, 2].set_title("C)", loc="left")
    axes[1, 0].set_title("D)", loc="left")
    axes[1, 1].set_title("E)", loc="left")
    # axes[1, 2].set_title("F)", loc="left")
    fig.tight_layout()
    plt.show()


# execute function
grid_histogram(10)
