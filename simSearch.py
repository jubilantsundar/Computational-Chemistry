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

@click.command()
@click.option('--querycsv', '-qin', "qcsv", help='csv file with query SMILES and names', required=True)
@click.option('--querysmilesheader', '-qsh', "qsmilesheader", help='SMILES column (SMILES)', default="SMILES")
@click.option('--querynamesheader', '-qnh', "qnamesheader", help='Name column (ID)', default="ID")
@click.option('--databasefile', '-din', "dbfile", help='database in csv or sdf format', required=True)
@click.option('--databasesmilesheader', '-dsh', "dbsmilesheader", help='SMILES column (SMILES)', default="SMILES")
@click.option('--databasenamesheader', '-dnh', "dbnamesheader", help='Name column (ID)', default="ID")
@click.option('--tanimotocutoff', '-c', "cutoff", help='Similarity cutoff (0.7)', default=0.7)
@click.option('--output', '-o', "outfilename", help='output file path (out.csv)', default='out.csv')
##@click.option('--rseed', '-s', default=12345, type=int, help='Use to obtain same results')
##@click.option('--randomcoords', '-c', default=False, type=bool, help='Use for random coordinate embedding')

def simSearch(qcsv,qsmilesheader,qnamesheader,dbfile,dbsmilesheader,dbnamesheader,cutoff,outfilename):

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
        db = pd.read_csv(dbfile)
        db = db[[dbsmilesheader,dbnamesheader]]
        print("{} DB entries were successfully loaded!".format(db.shape[0]))

    # the list for the dataframe
    query, target, query_s, target_s, simi = [], [], [], [], []

    for i in range(len(query_file)):
        query_smi = Chem.MolFromSmiles(query_file[qsmilesheader][i])
        query_id = query_file[qnamesheader][i]
        try:
            fp1 = Chem.RDKFingerprint(query_smi)
        except:
            print("Error occurred 1.")
    
        for j in range(len(db)):
            try:
                db_smi = Chem.MolFromSmiles(db[dbsmilesheader][j])
                if (db_smi is not None):
                    db_id = db[dbnamesheader][j]
                    fp2 = Chem.RDKFingerprint(db_smi)
                    tani = DataStructs.TanimotoSimilarity(fp1,fp2)
                    if (tani >= cutoff):
                        print("{} {} {}".format(query_id,db_id,tani))
                        query.append(query_id)
                        target.append(db_id)
                        query_s.append(query_file[qsmilesheader][i])
                        target_s.append(db[dbsmilesheader][j])
                        simi.append(tani)
                else:
                    print("please check SMILES {}".format(db_id))
            except:
                print("Exception Error!")


    # build the dataframe and sort it
    d = {'query_ID':query, 'query_smiles':query_s, 'db_hit_ID':target, 'db_hit_smiles':target_s, 'similarity':simi}
    df_final = pd.DataFrame(data=d)
    df_final = df_final.sort_values('similarity',ascending=False)
    print("No. of cpds with Tanimoto greater than {}: {}".format(cutoff,df_final.shape[0]))

    # save as csv
    df_final.to_csv(outfilename,index=False)
    print("Final output csv is saved as {}".format(outfilename))

if __name__=='__main__':
    simSearch()
