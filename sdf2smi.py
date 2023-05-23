#!/usr/bin/env python

## Importing requred Python packages
import click
from rdkit import Chem
from rdkit.Chem import PandasTools
import os

from rdkit import rdBase
print("RDKit_Version {}".format(rdBase.rdkitVersion))
#print(rdBase.boostVersion)

print("##############################################")
print("Author: Sundar Jubilant")
print("Email: jubilantsundar@gmail.com")
print("##############################################\n")

@click.command()
@click.option('--inputsdf', '-i', "infile", help='Input sdf file', required=True)
@click.option('--outsmi', '-o', "outfile", help='Output smi file. Default = out.smi', default="out.smi")
##@click.option('--idprefix', '-p', "id_prefix", help="Prefix for new IDs for the cpds. Default = mol_", default="mol_")
##@click.option('--rseed', '-s', default=12345, type=int, help='Use to obtain same results')
##@click.option('--randomcoords', '-c', default=False, type=bool, help='Use for random coordinate embedding')

def sdf2smi(infile,outfile):

    user_input = input('Question: Would you like to generate new IDs for the compounds or keep existing IDs? Type new/keep: ')

    if (user_input.lower() == 'new'):
        prefix = input('Please provide a prefix such as "mol" for new ID series: ')
        print("New IDs will be {}_1 {}_2 and so on.".format(prefix,prefix))

        df = PandasTools.LoadSDF(infile)
        df['SMILES'] = df.ROMol.apply(Chem.MolToSmiles)
        
        #Generating new IDs
        n_cpds = df.shape[0]
        print("Number of compounds: {}".format(n_cpds))
        id_list = [prefix + '{}'.format(i+1) for i in range(n_cpds)] 

        df['ID'] = id_list
        df = df[["SMILES","ID"]]
        df.to_csv(outfile,index=False)
        print('{} with newly generated IDs written as output.'.format(outfile))

    elif (user_input.lower() == 'keep'):
        print('Keeping the existing IDs.')
        df = PandasTools.LoadSDF(infile)
        df['SMILES'] = df.ROMol.apply(Chem.MolToSmiles)
        df = df[["SMILES","ID"]]
        df.to_csv(outfile,index=False)
        print('{} written as output.'.format(outfile))

    elif (user_input == None):
        print('You did not answer. Exiting NOW!!!')
        exit()

if __name__=='__main__':
    sdf2smi()
