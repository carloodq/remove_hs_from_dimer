from create_newxyzs import create_newxyzs
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

just_unique = True

def save_csv_nenci(path_nenci_xyzs):
    os.chdir(path_nenci_xyzs)
    items = os.listdir()
    all_xyz = []

    for dimer in items[:30]:
        raw_mol = Chem.MolFromXYZFile(dimer)
        mol = Chem.Mol(raw_mol)
        rdDetermineBonds.DetermineConnectivity(mol)
        m = mol
        # get charge and multiplicity from comment line
        f = open(dimer, "r") 
        cl = f.readlines()[1]
        charge, mult = cl.replace("  ", " ").split(" ")[:2]
        
        newxyzs = create_newxyzs(m, str(dimer)[4:-9], str(dimer)[-8:-4], just_unique, charge, mult)  
        all_xyz.extend(newxyzs)  

    all_xyz = pd.DataFrame(all_xyz, columns=['conformer', 'molecule_name', 'xyz',
                                        'xyz_reduced', 'atom_index', 'attached_atom',
                                        'neighboring_atoms', 'fragment', 'netcharge',
                                        'multiplicity', 'smiles_orig', 'smiles_reduced']) 
    
    # add molecule index
    all_xyz['item_index' ] = all_xyz.groupby("molecule_name").ngroup()

    os.chdir("../")
    all_xyz.to_csv('all_xyz_nenci.csv', index=None)

    return 0


save_csv_nenci("xyznenci")

