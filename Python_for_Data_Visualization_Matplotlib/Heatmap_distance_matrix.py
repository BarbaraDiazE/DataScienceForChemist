"""Compute similarity matrix from MaccsKeys"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.DataManip.Metric.rdMetricMatrixCalc import (
    GetTanimotoSimMat,
    GetTanimotoDistMat,
)


def compute_maccskeys_fp(Data):
    smiles = list(Data["canonical_smiles"])  # set smile colum
    smi = [Chem.MolFromSmiles(x) for x in smiles]
    fps = [MACCSkeys.GenMACCSKeys(x) for x in smi]
    tanimoto_sim_mat_lower_triangle = GetTanimotoSimMat(fps)
    n_mol = len(fps)
    similarity_matrix = np.ones([n_mol, n_mol])
    i_lower = np.tril_indices(n=n_mol, m=n_mol, k=-1)
    i_upper = np.triu_indices(n=n_mol, m=n_mol, k=1)
    similarity_matrix[i_lower] = tanimoto_sim_mat_lower_triangle
    similarity_matrix[i_upper] = similarity_matrix.T[i_upper]
    distance_matrix = np.subtract(1, similarity_matrix)
    distance_matrix = np.round(distance_matrix, 2)
    print(distance_matrix)
    ids = list(Data["chembl_id"])  # set id columns
    return distance_matrix, ids


def plot(distance_matrix, ids):
    plt.subplots(figsize=(10, 10))
    ax = sns.heatmap(
        data=distance_matrix,
        yticklabels=ids,
        xticklabels=ids,
        linewidths=0.2,
        # vmin= int,
        # vmax=int,
        # annot=True,
    )
    ax.set_title("Tanimoto distance matrix", fontsize=30)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    # plt.xlabel("x-axis",fontsize = 24)
    # plt.ylabel("y-axis", fontsize = 24)
    plt.show()
    plt.savefig("heatmap.png")
    return plt


Data = pd.read_csv("sample1_ChEMBL.csv", sep=",", index_col="Unnamed: 0")
Data = Data.sample(10, random_state=1992)
print(Data.columns)
distance_matrix, ids = compute_maccskeys_fp(Data)
plot = plot(distance_matrix, ids)
