
about = """
======================================================================================================
Description: sdf2csv.py converts a sdf file to 2D csv file.
Usage: python sdf2csv.py in.sdf out.csv
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

def sdf2csv(in_file,out_file):
    frame = PandasTools.LoadSDF(in_file,smilesName='SMILES')
    frame = frame.drop(['ID','ROMol'], axis=1)
    frame.to_csv(out_file,index=False)
 
if (len(sys.argv) == 3):
    sdf2csv(sys.argv[1],sys.argv[2])
else:
    '''Usage instructions'''
    print("=====================")
    print("Something went wrong!")
    print("=====================")
    print("Usage instructions: python sdf2csv.py in.sdf out.csv")
    