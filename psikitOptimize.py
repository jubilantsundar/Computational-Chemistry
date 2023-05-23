#!/usr/bin/env python

## QM optimization using Psi4
## Using Psikit as a python wrapper
## Inspired from https://iwatobipen.wordpress.com/2019/01/24/draw-homo-lumo-with-psikit-rdkit-psi4-pymol/

import click
from rdkit import Chem
from psikit import Psikit
import pandas as pd
from biopandas.pdb import PandasPdb
import os
import shutil

import rdkit
print(rdkit.__version__)

@click.command()
@click.option('--input', '-i', "insmiles", help='input SMILES', required=True)
@click.option('--inputid', '-n', "inname", help='input ARK ID', required=True)
##@click.option('--output', '-o', help='output file path', default='gamess_out.sdf')
##@click.option('--rseed', '-s', default=12345, type=int, help='Use to obtain same results')
##@click.option('--randomcoords', '-c', default=False, type=bool, help='Use for random coordinate embedding')

def psikitOpt(insmiles,inname):
    pk = Psikit()
    pk.read_from_smiles(insmiles)
    # It is better to optimize although it is single point calculation due to time
    pk.optimize(basis_sets="scf/sto-3g")
    
    print("############################################")
    print("Optimizer: Optimization complete!")
    print("Compoud: {} Energy: {}".format(inname,pk.energy()))
    print("############################################")
        
    # get MO view! to visualize in PyMol
    pk.getMOview()
    # Makes a frontier.py files
    # run frontier.py to visualize orbitals in Pymol
    # You may want to make sure the paths in the file are correct
    pk.save_frontier() 

    ## Reading in all the cube files generated in a tmp folder.
    ## Please confirm if the folder is correct.
    if os.path.exists(inname):
        print("{} directory exists. Please check for its content and make sure to delete!".format(inname))
    else:
        os.mkdir(inname)
        tmpdir = pk.tempdir
        print("Moving cube files from temp directory: {}".format(tmpdir))

        ## Moving cube files from tmp directory
        cube_files = [f for f in os.listdir(tmpdir)]
        for cube_file in cube_files:
            old_path = tmpdir + "/" + cube_file
            new_path = inname + "/" + cube_file
            shutil.move(old_path,inname)

    ## Get ready for RESP charges
    print("##############################################################################")
    print("Getting RESP charges. Changing directory to write PDB files with RESP charges.")
    print("##############################################################################")
    os.chdir(inname)

    mol = pk.mol
    atoms = mol.GetAtoms()
    RESP = pk.calc_resp_charges()
    data = {'RESP': [float(atom.GetProp('RESP')) for atom in atoms]}
    df = pd.DataFrame(data)

    Chem.MolToPDBFile(pk.mol, '{}.pdb'.format(inname))
    pdbobj = PandasPdb().read_pdb('{}.pdb'.format(inname))
    
    # Following step is editing state, you can see the approach like a pandas method ;)
    pdbobj.df['HETATM']['b_factor'] = df['RESP']
    pdbobj.to_pdb('{}.pdb'.format(inname))
    print("{}.pdb was written with RESP charges in B-factor column.".format(inname))
    print("##############################################################################")

if __name__=='__main__':
    psikitOpt()
