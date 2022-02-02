about = """
======================================================================================================
Description: csv2sdf.py converts a csv file containing SMILES column with smiles strings to 2D sdf file.
Usage: python csv2sdf.py in.csv out.sdf
======================================================================================================
"""
print(about)

print("Info: Make sure you have installed RDKit and it is working.")
print("Info: Best way to get it installed on Windows machine is through Anaconda")
print("Info: Importing required packages")
print("============================================================================")

print("Importing required packages.")
import os
import pandas as pd
import sys
import argparse

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit import RDConfig

print("Required packages imported.")


print("RDKit Version: {}".format(rdBase.rdkitVersion))
print("RDKit Boost Version: {}".format(rdBase.boostVersion))

def csv2sdf(in_file,out_file):
    df = pd.read_csv(in_file)
    PandasTools.AddMoleculeColumnToFrame(df,'SMILES','Molecule')
    PandasTools.WriteSDF(df, out_file, molColName='Molecule', properties=list(df.columns))
 
if (len(sys.argv) == 3):
    csv2sdf(sys.argv[1],sys.argv[2])
else:
    '''Usage instructions'''
    print("=====================")
    print("Something went wrong!")
    print("=====================")
    print("Usage instructions: python csv2sdf.py in.csv out.sdf")
    print("Make sure column containing smiles is named 'SMILES'")