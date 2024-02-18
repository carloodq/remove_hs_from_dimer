import pandas as pd
from extract_xyz import extract_xyz
import os
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
# df = pd.read_csv("Donchev et al DES370K.csv")
# extract_xyz(df.head(100))

just_unique = True

def create_newxyzs(m, item, conformer, just_unique):
    ''' 
    it takes all the Hs that could be removed, it records the 1st connection and the neighbor's connections
    takes as input a Mol RDkit object of a molecule, returns a list of broken xyzs with respective parameters
    '''
    newxyzs = []
    # checking if dimer isonly made f CNOPS and H
    if len(m.GetAtoms()) == len([a.GetSymbol() for a in m.GetAtoms() if a.GetSymbol() in 'CNOPSH' ]):
        for atom in m.GetAtoms():
            if atom.GetSymbol() == 'H':
                h = atom.GetIdx()
                fragment = str(atom.GetOwningMol())
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
                            xyz_reduced = Chem.MolToXYZBlock(mw)
                            xyz = Chem.MolToXYZBlock(m)
                            newxyzs.append([conformer, item , xyz , xyz_reduced,  h , 
                                            n.GetSymbol() , neighbors, fragment])

    if just_unique:                      
        newxyzs = pd.DataFrame(newxyzs, columns = ['conformer', 'item' , 'xyz' , 'xyz_reduced',  'h' , 
                        '1st' , 'neighbors', 'fragment'])
        
        newxyzs['neighbors'] = newxyzs['neighbors'].apply(lambda x: str(x))

        newxyzs = newxyzs.groupby(['1st', 'neighbors']).apply(lambda x: x.sample(1)).reset_index(drop=True)
        newxyzs = [list(x) for x in newxyzs.values]

    return newxyzs


os.chdir("xyzs")
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
        # add all the conformers for each dimer 
        for conformer in xyz:
            xyzp = os.path.join(cwd, item, 'qm_opt_dimer', all_xyzs[0], conformer)
            # open the xyz and add the bonds
            raw_mol = Chem.MolFromXYZFile(xyzp)
            mol = Chem.Mol(raw_mol)
            rdDetermineBonds.DetermineConnectivity(mol)
            m = mol
            newxyzs = create_newxyzs(m, item, str(conformer), just_unique)  
            all_xyz.extend(newxyzs)   

all_xyz = pd.DataFrame(all_xyz, columns=['conformer', 'molecule_name', 'xyz',
                                          'xyz_reduced', 'atom_index', 'attached_atom',
                                            'neighboring_atoms', 'fragment'])
# save it 
os.chdir("../")
all_xyz.to_csv('all_xyz.csv')
