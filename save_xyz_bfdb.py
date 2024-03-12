from generate_csv import create_newxyzs
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

just_unique = True

def get_csv_bfdb(path_bfdb_xyzs):
    os.chdir(path_bfdb_xyzs)
    items = os.listdir()
    all_xyz = []

    for dimer in items[:3]:
        if dimer[-9:] == "dimer.xyz": # only dimers
            raw_mol = Chem.MolFromXYZFile(dimer)
            mol = Chem.Mol(raw_mol)
            rdDetermineBonds.DetermineConnectivity(mol)
            m = mol
            # get charge and multiplicity from comment line
            f = open(dimer, "r") 
            cl = f.readlines()[1]
            charge, mult = cl.split(" ")[:2]

            newxyzs = create_newxyzs(m, str(dimer)[:-4], 'NA', just_unique, charge, mult)  
            all_xyz.extend(newxyzs)  

    all_xyz = pd.DataFrame(all_xyz, columns=['conformer', 'molecule_name', 'xyz',
                                        'xyz_reduced', 'atom_index', 'attached_atom',
                                        'neighboring_atoms', 'fragment', 'netcharge', 'multiplicity', 'smiles_orig', 'smiles_reduced']) 
    os.chdir("../../")

    return all_xyz


dfs = []
bfdb_datasets = ["BBI_xyzfiles/BBI25", "BBI_xyzfiles/BBI_less_BBI25",
                 "SSI_xyzfiles/SSI_less_SSI500", "SSI_xyzfiles/SSI500_less_SSI100", "SSI_xyzfiles/SSI100"]

for pth in bfdb_datasets:
    dfs.append(get_csv_bfdb(pth))

all_xyz = pd.concat(dfs)

all_xyz.to_csv('all_xyz_bfdb.csv', index = None)

