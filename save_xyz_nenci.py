from generate_csv import create_newxyzs
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

just_unique = True

def save_csv_nenci(path_nenci_xyzs):
    os.chdir(path_nenci_xyzs)
    items = os.listdir()
    all_xyz = []

    for dimer in items[:3]:
        raw_mol = Chem.MolFromXYZFile(dimer)
        mol = Chem.Mol(raw_mol)
        rdDetermineBonds.DetermineConnectivity(mol)
        m = mol
        newxyzs = create_newxyzs(m, str(dimer)[4:-9], str(dimer)[-8:-4], just_unique)  
        all_xyz.extend(newxyzs)  

    all_xyz = pd.DataFrame(all_xyz, columns=['conformer', 'molecule_name', 'xyz',
                                        'xyz_reduced', 'atom_index', 'attached_atom',
                                        'neighboring_atoms', 'fragment']) 
    os.chdir("../")
    all_xyz.to_csv('all_xyz_nenci.csv', index=None)

    return 0

save_csv_nenci("xyznenci")

