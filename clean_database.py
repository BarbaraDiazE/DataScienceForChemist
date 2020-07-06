import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import SaltRemover
from molvs.standardize import Standardizer
from molvs.tautomer import TautomerCanonicalizer

"""
script to genetare  a new csv as product of cleaned a database
"""

"""
funcions
"""
# delete organometalics
def check_atoms(mol):
    """
    filter compounds that contains just desired atoms
    """
    desired_elements = {"H", "B", "C", "N", "O", "F", "P", "S", "Cl"}
    elements = set([i.GetSymbol() for i in mol.GetAtoms()])
    atom_difference = len(elements - desired_elements)
    return atom_difference


def remove_salts(mol):
    remover = SaltRemover.SaltRemover()
    res = remover.StripMol(mol)
    # return Chem.MolToSmiles(res)
    return res


def standarization(mol):
    s = Standardizer()
    s_mol = s.standardize(mol)
    return s_mol


def get_tautomer(mol):
    TC = TautomerCanonicalizer()
    return TC(mol)


"""
execution
"""


def clean_database(input_file, smiles_column, output_file):
    data = pd.read_csv(input_file, sep=",", index_col="Unnamed: 0")
    data = data.reset_index()
    data = data.drop("index", axis=1)
    print(data.columns)
    smi = data[smiles_column].tolist()
    #
    cleaned_smiles = list()
    index = list()
    for i in range(len(smi)):
        print(smi[i])
        try:
            mol = Chem.MolFromSmiles(smi[i])
            if mol == None:
                print("smile not valid")
            else:

                print("mol exist")
                # check atoms
                r = check_atoms(mol)
                if r == 0:
                    # remove salts
                    ns_mol = remove_salts(mol)
                    # std
                    s_mol = standarization(ns_mol)
                    # tautomer
                    t_mol = get_tautomer(s_mol)
                    # Draw.MolToFile(t_mol, str(i) + "_tmol.png")
                    index.append(i)
                    cleaned_smiles.append(Chem.MolToSmiles(t_mol))
                else:
                    print("a not allowed atom is present")

        except:
            return "function not valid"
    # print("index", index)
    cleaned_data = data.loc[index, :]
    cleaned_data = cleaned_data.reset_index()
    cleaned_data = cleaned_data.drop("index", axis=1)
    cleaned_data["cleaned_smiles"] = cleaned_smiles
    cleaned_data.to_csv(output_file, sep=",", index=True)
    print(
        f'{"initial compounds "}{data.shape[0]}{", "} {"final compounds "}{cleaned_data.shape[0]}'
    )
    return cleaned_data


cleaned_data = clean_database(
    "sample1_ChEMBL.csv", "canonical_smiles", "cleaned_ChEMBL.csv"
)
print(cleaned_data)
print(cleaned_data.shape)
