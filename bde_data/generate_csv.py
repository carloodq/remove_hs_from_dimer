import pandas as pd
from extract_xyz import extract_xyz
import os
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
df = pd.read_csv("Donchev et al DES370K.csv")
# extract_xyz(df.head(100))

def create_newxyzs(m, item):
    ''' 
    it takes all the Hs that could be removed, it records the 1st connection and the neighbor's connections
    takes as input a Mol RDkit object of a molecule, returns a list of broken xyzs with respective parameters
    '''
    newxyzs = []
    for atom in m.GetAtoms():
        if atom.GetSymbol() == 'H':
            h = atom.GetIdx()
            ns = atom.GetNeighbors()
            if len(ns) > 1: # if not outside
                continue
            else:
                for n in ns:
                    if n.GetSymbol() in 'CNOPS':
                        organic = n.GetIdx()
                        neighbors =  [n.GetSymbol() for n in n.GetNeighbors()] 
                        # removing atom to create mol 
                        mw = Chem.RWMol(m)
                        mw.RemoveAtom(h)
                        xyz = Chem.MolToXYZBlock(mw)
                        newxyzs.append([item , xyz , h , n.GetSymbol() , neighbors])
    return newxyzs

cwd = os.getcwd()
items = os.listdir()
all_xyz = []
# iterates through all generated folders and process the first xyz
for item in items:
    if os.path.isdir(item) and item != '__pycache__':
        all_xyzs = os.listdir(os.path.join(cwd, item, 'qm_opt_dimer'))
        # path to the xyz file
        xyzf = os.path.join(cwd, item, 'qm_opt_dimer', all_xyzs[0])
        xyz = os.listdir(xyzf)
        xyzp = os.path.join(cwd, item, 'qm_opt_dimer', all_xyzs[0], xyz[0])
        # open the xyz and add the bonds
        raw_mol = Chem.MolFromXYZFile(xyzp)
        mol = Chem.Mol(raw_mol)
        rdDetermineBonds.DetermineConnectivity(mol)
        m = mol
        newxyzs = create_newxyzs(m, item)  
        all_xyz.extend(newxyzs)   

all_xyz = pd.DataFrame(all_xyz, columns=['molecule_name', 'xyz', 'atom_index', 'attached_atom', 'neighboring_atoms'])
# save it 
all_xyz.to_csv('all_xyz.csv')
