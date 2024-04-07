import pandas as pd
from extract_xyz import extract_xyz
import os
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
# df = pd.read_csv("Donchev et al DES370K.csv")
# extract_xyz(df.head(100))


def create_newxyzs(m, item, conformer, just_unique, net_charge = None, multiplicity = None):
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
                            smiles_reduced = Chem.MolToSmiles(mw)
                            smiles_orig = Chem.MolToSmiles(m)
                            newxyzs.append([conformer, item , xyz , xyz_reduced,  h , 
                                            n.GetSymbol() , neighbors, fragment, net_charge, multiplicity, smiles_orig, smiles_reduced ])

    if just_unique:                      
        newxyzs = pd.DataFrame(newxyzs, columns = ['conformer', 'item' , 'xyz' , 'xyz_reduced',  'h' , 
                        '1st' , 'neighbors', 'fragment', 'netcharge', 'multiplicity', 'smiles_orig', 'smiles_reduced' ])
        
        newxyzs['neighbors'] = newxyzs['neighbors'].apply(lambda x: str(x))

        newxyzs = newxyzs.groupby(['1st', 'neighbors']).apply(lambda x: x.sample(1)).reset_index(drop=True)
        newxyzs = [list(x) for x in newxyzs.values]


    

    return newxyzs

