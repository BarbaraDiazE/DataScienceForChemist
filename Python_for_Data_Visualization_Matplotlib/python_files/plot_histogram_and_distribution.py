"""
plot histogram and distribution on a grid using matplotlib.
"""
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
from scipy.stats import uniform
import scipy.special


def grid_histogram(n_bins, descriptor):
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
    #### histogram_one_descriptor###
    data = pd.read_csv(
        "/home/babs/Documents/DIFACQUIM/admetox_research/Libraries/Master_results.csv"
    )
    # define libraries
    l = ["FDA", "BIOFACQUIM", "TCM", "AfroDB", "NuBBEDB"]
    # define colors
    colors = ["blueviolet", "green", "mediumvioletred", "darkorange", "dodgerblue"]
    # descriptor1
    x_multi = list()
    for i in l:
        d = data[data["Library"] == i]
        x_multi.append(d[descriptor].tolist())
    label = ["FDA", "BIOFACQUIM", "TCM", "Afrodb", "NuBBEdb"]
    # subplot histogram
    axes[0].grid(color="silver", linestyle="-", linewidth=0.5)
    axes[0].hist(
        x_multi, n_bins, histtype="bar", density=True, label=label, color=colors
    )
    axes[0].legend()
    axes[0].set_xlabel(descriptor)
    axes[0].set_ylabel("Fraction")
    axes[0].set_ylim(0, 0.31)
    ### subplot normal distribution ###
    label = ["FDA", "BIOFACQUIM", "TCM", "AfroDB", "NuBBEDB"]
    colors_dict = {
        "FDA": "blueviolet",
        "BIOFACQUIM": "green",
        "TCM": "mediumvioletred",
        "AfroDB": "darkorange",
        "NuBBEDB": "dodgerblue",
    }
    axes[1].grid(color="silver", linestyle="-", linewidth=0.5)
    for i in label:
        d = data[data["Library"] == i]
        x = np.sort(d["Consensus Log P"].to_numpy())
        mu = np.mean(x)
        sigma = np.std(x)
        hist, edges = np.histogram(x, density=True)
        y = np.array(
            1
            / (sigma * np.sqrt(2 * np.pi))
            * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
        )
        axes[1].plot(x, y, label=i, color=colors_dict[i])
    axes[1].legend()
    axes[1].set_xlabel(descriptor)
    axes[1].set_ylabel("Probability density")
    ##titles
    axes[0].set_title("A)", loc="left")
    axes[1].set_title("B)", loc="left")

    fig.tight_layout()
    plt.show()


# execute function
grid_histogram(10, "Consensus Log P")
