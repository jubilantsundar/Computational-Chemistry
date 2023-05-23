# Install RDKit with 
#"conda install -c rdkit rdkit"

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit import RDConfig
import os
import sys
import numpy as np

from rdkit import rdBase
#print("RDKit_Version {}".format(rdBase.rdkitVersion))
#print(rdBase.boostVersion)

if len(sys.argv) < 8:
	print("Usage: python3 similarity.py input.csv header_for_SMILES header_for_ID database_file.csv DB_header_for_SMILES DB_header_for_ID output.csv tanimoto_cutoff")
    #print("Example: python3 similarity.py input.csv SMILES Name database_file.csv SMILES ID output.csv 0.8")
	exit()

in_filename = sys.argv[1]
smi_colname = sys.argv[2]
mol_colname = sys.argv[3]

db_filename = sys.argv[4]
db_smi_colname = sys.argv[5]
db_mol_colname = sys.argv[6]

out_filename = sys.argv[7]
sim_cutoff = sys.argv[8]
sim_cutoff = float(sim_cutoff)

query_file = pd.read_table(in_filename,sep=",")
query_file = query_file[[smi_colname,mol_colname]]
print("{} queries were successfully loaded!".format(query_file.shape[0]))

db = pd.read_table(db_filename,sep=",")
db = db[[db_smi_colname,db_mol_colname]]
db = db.replace(np.nan, '', regex=True)
print("{} DB entries were successfully loaded!".format(db.shape[0]))

# the list for the dataframe
query, target, query_s, target_s, simi = [], [], [], [], []

for i in range(len(query_file)):
    #print(query_file['SMILES'][i])
    query_smi = Chem.MolFromSmiles(query_file['SMILES'][i])
    query_id = query_file['Name'][i]
    try:
        fp1 = Chem.RDKFingerprint(query_smi)
    except:
        print("Error occurred 1.")
    
    for j in range(len(db)):
        db_smi_test = db['SMILES'][j]
        if (db_smi_test is not None):
            try:
                #print(db['SMILES'][j])
                db_smi = Chem.MolFromSmiles(db['SMILES'][j])
                db_id = db['Name'][j]
                fp2 = Chem.RDKFingerprint(db_smi)
                tani = DataStructs.TanimotoSimilarity(fp1,fp2)
                #tracking the process
                #if (j%100000 == 0):
                print("{} {} {}".format(query_id,db_id,tani))
                if (tani >= sim_cutoff):
                    print("{} {} {}".format(query_id,db_id,tani))
                    query.append(query_id)
                    target.append(db_id)
                    query_s.append(query_file['SMILES'][i])
                    target_s.append(db['SMILES'][j])
                    simi.append(tani)
            except:
                print("Error")
        else:
            print("please check SMILES {}".format(db_id))


# build the dataframe and sort it
d = {'query_ID':query, 'query_smiles':query_s, 'db_hit_ID':target, 'db_hit_smiles':target_s, 'similarity':simi}
df_final = pd.DataFrame(data=d)
df_final = df_final.sort_values('similarity',ascending=False)
print("No. of cpds with Tanimoto greater than {}: {}".format(sim_cutoff,df_final.shape[0]))

# save as csv
df_final.to_csv(out_filename,index=False)
