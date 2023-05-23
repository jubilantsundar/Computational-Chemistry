
import click
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

from rdkit.Chem import rdDistGeom
from rdkit.Chem import TorsionFingerprints
import numpy as np

import rdkit
print(rdkit.__version__)


@click.command()
@click.option('--input', '-i', help='inputfile MOL', required=True)
@click.option('--output', '-o', help='output file path', default='gen_confs.sdf')
@click.option('--rseed', '-s', default=12345, type=int, help='Use to obtain same results')
@click.option('--randomcoords', '-c', default=False, type=bool, help='Use for random coordinate embedding')
@click.option('--forcetol', '-f', default=0.0135, type=float, help='Default changed from 0.001 to 0.0135 based on Greg Landrum work')
@click.option('--prunermsthresh', '-t', default=0.1, type=float, help='Retain only the conformations out of ‘numConfs’')
@click.option('--numconf', default=50, type=int)
@click.option('--add_ref', '-r', default=False, type=bool)

def confgen(input, output, rseed, randomcoords, forcetol, prunermsthresh, numconf, add_ref):
    mol = Chem.AddHs(Chem.MolFromMolFile(input), addCoords=True)
    refmol = Chem.AddHs(Chem.Mol(mol))
    param = rdDistGeom.ETKDGv3()
    param.randomSeed = rseed
    param.useRandomCoords = randomcoords
    param.optimizeForceTol = forcetol
    param.pruneRmsThresh = prunermsthresh
    param.verbose = False
    param.numThreads = 8

    cids = rdDistGeom.EmbedMultipleConfs(mol, numconf, param)
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant='MMFF94s')
    w = Chem.SDWriter(output)
    if add_ref:
        refmol.SetProp('CID', '-1')
        refmol.SetProp('Energy', '')
        w.write(refmol)
    res = []

    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))
    sorted_res = sorted(res, key=lambda x:x[1])
    rdMolAlign.AlignMolConformers(mol)
    for cid, e in sorted_res:
        mol.SetProp('CID', str(cid))
        mol.SetProp('Energy', str(e))
        w.write(mol, confId=cid)
    w.close()

if __name__=='__main__':
    confgen()
