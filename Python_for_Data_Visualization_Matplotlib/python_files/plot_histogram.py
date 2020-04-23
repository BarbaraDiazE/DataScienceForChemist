"""
plot histogram using matplotlib.
"""
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
from scipy.stats import uniform
import scipy.special


def grid_histogram(n_bins, descriptor):
    #### histogram_one_descriptor###
    data = pd.read_csv(
        "/home/babs/Documents/DIFACQUIM/admetox_research/Libraries/Master_results.csv"
    )
    plt.figure()
    # define libraries
    l = ["FDA", "BIOFACQUIM", "TCM", "AfroDB", "NuBBEDB"]
    # define legends and colors
    color_palette = {
        "FDA": "blueviolet",
        "BIOFACQUIM": "green",
        "TCM": "mediumvioletred",
        "AfroDB": "darkorange",
        "NuBBEDB": "dodgerblue",
    }

    # descriptor
    x_multi = list()
    for i in l:
        d = data[data["Library"] == i]
        x_multi.append(d[descriptor].tolist())
    plt.grid(color="silver", linestyle="-", linewidth=0.5)
    plt.hist(
        x_multi,
        n_bins,
        histtype="bar",
        density=True,
        label=color_palette.keys(),
        color=color_palette.values(),
    )
    plt.legend()
    plt.xlabel(descriptor)
    plt.ylabel("Fraction")
    # plt.ylim(0, 0.31)

    ##plot title
    plt.title("A)", loc="left")
    # show
    plt.show()


# execute function
grid_histogram(10, "Consensus Log P")
