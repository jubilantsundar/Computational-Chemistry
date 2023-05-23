#!/usr/bin/env python

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit import RDConfig

import click
import os
import sys
import numpy as np

from rdkit import rdBase
print("RDKit_Version {}".format(rdBase.rdkitVersion))
#print(rdBase.boostVersion)

print("Author: Sundar Jubilant")
print("Email: jubilantsundar@gmail.com")
print("##############################################")
print(">>>Following is required to run this script. <<<")
print(">>>This script takes only csv as query input. <<<")
print(">>>This script takes either SDF file of a text file")
print("where the smiles and ID are separated either by") 
print("space or comma as DB input. <<<")
print("##############################################")

@click.command()
@click.option('--querycsv', '-qin', "qcsv", help='csv file with query SMILES and names', required=True)
@click.option('--querysmilesheader', '-qsh', "qsmilesheader", help='column name containing SMILES', default="SMILES", required=True)
@click.option('--querynamesheader', '-qnh', "qnamesheader", help='column name containing names', default="ID", required=True)
@click.option('--databasefile', '-din', "dbfile", help='csv or sdf or cxsmiles file with database SMILES and names', required=True)
@click.option('--databasesmilesheader', '-dsh', "dbsmilesheader", help='database column name containing SMILES', default="SMILES", required=True)
@click.option('--databasenamesheader', '-dnh', "dbnamesheader", help='database column name containing names', default="ID", required=True)
@click.option('--separator', '-sep', "separator", help='symbol separating the smiles from ID such as " " or ","', required=True)
@click.option('--output', '-o', "outfilename", help='output file path', default='out.csv')
##@click.option('--rseed', '-s', default=12345, type=int, help='Use to obtain same results')
##@click.option('--randomcoords', '-c', default=False, type=bool, help='Use for random coordinate embedding')

def substrSearch(qcsv,qsmilesheader,qnamesheader,dbfile,dbsmilesheader,dbnamesheader,separator,outfilename):

    ## Reading input in csv format
    if qcsv.endswith('.csv'):
        query_file = pd.read_table(qcsv,sep=",")
        query_file = query_file[[qsmilesheader,qnamesheader]]
        print("{} queries were successfully loaded!".format(query_file.shape[0]))
    else:
        print("Input file is not in CSV format")
        exit()

    ## Reading DB file in either csv or sdf format
    if dbfile.endswith('.sdf'):
        db = PandasTools.LoadSDF(dbfile)
        db[dbsmilesheader] = db.ROMol.apply(Chem.MolToSmiles)
        db = db[[dbsmilesheader,dbnamesheader]]
        print("{} DB entries were successfully loaded!".format(db.shape[0]))
    else:
        db = pd.read_table(dbfile,sep=separator)
        db = db[[dbsmilesheader,dbnamesheader]]
        print("{} DB entries were successfully loaded!".format(db.shape[0]))


    # the list for the dataframe
    query, target, query_s, target_s, simi = [], [], [], [], []

    for i in range(len(query_file)):
        query_mol = Chem.MolFromSmiles(query_file[qsmilesheader][i])
        query_id = query_file[qnamesheader][i]
    
        for j in range(len(db)):
            try:
                db_smi = db[dbsmilesheader][j]
                db_id = db[dbnamesheader][j]
                if (db_smi is not None):
                    db_mol = Chem.MolFromSmiles(db[dbsmilesheader][j])
                    db_id = db[dbnamesheader][j]
            
                    if (db_mol.HasSubstructMatch(query_mol)):
                        print("{} {} a match!".format(query_id,db_id))
                        query.append(query_id)
                        target.append(db_id)
                        query_s.append(query_file[qsmilesheader][i])
                        target_s.append(db[dbsmilesheader][j])
                else:
                    print("please check SMILES {}".format(db_id))
            except Exception:
                pass

    # build the dataframe and sort it
    d = {'query_ID':query, 'query_smiles':query_s, 'db_hit_ID':target, 'db_hit_smiles':target_s}
    df_final = pd.DataFrame(data=d)

    print("No. of cpds with Substructure {}".format(df_final.shape[0]))
    # save as csv
    df_final.to_csv(outfilename,index=False)
    print("Final output csv is saved!")

if __name__=='__main__':
    substrSearch()
