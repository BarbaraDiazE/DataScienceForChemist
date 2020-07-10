import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

import statistics as st
import itertools as it

"""Plot ROC by differents Kernels"""

root = {
    "root": "/home/barbara/Documents/DIFACQUIM/PPI_classifier/phase-1/Databases/morgan2/"
}


class PlotSimilarity:
    def __init__(self, files):
        self.files = files
        print(self.files)

    def data(self):
        libraries = list()
        plot_data = dict()
        for i in files:
            db = pd.read_csv(i, index_col="Unnamed: 0")
            # print(db.head())
            # print(db.columns)
            sim = np.array(db.sim)
            y = np.array(db.y)
            _ = db.library.loc[0]
            print(db.library.loc[0])
            libraries.append(_)
            print(libraries)
            plot_data[_] = {"sim": sim, "y": y}
        self.libraries = libraries
        print(plot_data)
        return plot_data

    def plot_sim(self, colors):
        plot_data = self.data()
        fig = plt.figure()
        lw = 2
        libraries = self.libraries
        for i in range(len(libraries)):
            print(libraries[i])
            plt.plot(
                plot_data[libraries[i]]["sim"],
                plot_data[libraries[i]]["y"],
                color=colors[i],
                lw=lw,
                linestyle="-",
                label=libraries[i],
            )
        plt.xlim([0.0, 1.01])
        plt.ylim([0.0, 1.01])
        plt.xlabel("Similarity")
        plt.ylabel("CDF")
        plt.title("Diversity Analysis")
        plt.legend(loc="lower right", ncol=1, shadow=False, fancybox=False)
        plt.show()
        # plt.savefig("div_analysis.png")
        fig.savefig("plot.png")


###Define variables ###
# files is a list list with individual database files
# colors is a list with the nessesary number of colors for each database
files = ["MACCS_Tanimoto_BIOFACQUIM2V_.csv", "MACCS_Tanimoto_NUBBE2V_.csv"]
colors = ["mediumvioletred", "forestgreen"]

# Execute plot
a = PlotSimilarity(files)
a.plot_sim(colors)

