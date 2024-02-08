# Removing Hs connected to CNOPS

## Part 1
Set up the directory with the xyz file by running this command

```python
df = pd.read_csv("Donchev et al DES370K.csv")
extract_xyz(df)
```

## Part 2
You can run `generate_csv.py` to create the new xyzs, which correspond to the original xyz for each dimer, minus one hydrogen.

For the hydrogen to be removed it must be connected to one of C,N,O,P or S.

The resulting file, `all_xyz.csv`, contains all new generated xyzs. 
Each row includes:
- **molecule_name**: smiles of the dimer
- **xyz**: xyz of the original dimer
- **xyz_reduced**: xyz of the new dimer (without one H)
- **atom_index**: index of the removed atom in the original XYZ
- **attached_atom**: 1st connection to the removed H
- **neighboring_atoms**: neighbouring atoms to the removed H
- **fragment**: monomer the abstracted H belonged to 

Core subset 15k:
- QM up to 4 confromers
- MD up to 10 conformers

Here are some examples of one original and one modified dimer. You can notice that the second dimer is missing one H.

![](whole.png)
![](1.png)

## Part 3
Once the .csv has been generated, you make create folders and .xyz files.
Make sure to be at the same path as the `make_folders.py` script and `all_xyz.xyz`, then run the following command 
```
python make_folders.py
```

The folder structure generated will be
```
├── <dimer name>
│   ├── <conformer_k_index>
│   │   ├── original.xyz
│   │   ├── <index of abstracted h>.xyz
│   │   ├── <index of abstracted h>.xyz
│   │   ├── ..
│   │   ├── ..
```





