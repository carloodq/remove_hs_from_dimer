# Removing Hs connected to CNOPS

These scripts take as input some .xyz files (of dimers) and generate new .xyz files which are many possible abstractions of Hydrogens in the original structure, starting with all possible abstractions and then only keeping a diverse subset, based on criteria such as 1st and 2nd connectivity to the abstracted hydrogen. 

## Part 1
Set up the directory with the xyz file by running this command

```python
df = pd.read_csv("Donchev et al DES15K.csv")
extract_xyz(df)
```

## Part 2.a
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
- QM up to 4 conformers
- MD up to 10 conformers

After generating a new row for each H in the dimer, we drop duplicates. The logic adopted is to randomly keep just one row when there are more than 1 H abstractions with identical 1st connectivity and 2nd connectivity. We can instead decide to keep all rows by setting ```just_unique = False``` in  `generate_csv.py`.

Here are some examples of one original and one modified dimer. You can notice that the second dimer is missing one H.

![](whole.png)
![](1.png)


## Part 2.b
In addition to generating abstracted molecules from **Donchev et al DES15K**, we decided to extend the generation to another dimer dataset, **NENCI-2021**. It includes a total of 7741 dimer configurations, generated from 141 unique dimers at 7 intermediate distances and 9 angles.

The initial .xyz files are found in the `xyznenci` folder. Running the script `save_xyz_nenci.py` generates the file `all_xyz_nenci.csv`, which contains the same columns as `all_xyz.csv`.

## Part 2.c
We will add SSI and BBI from the BFDb.
You can generate the structures with the script `save_xyz_bfdb.py`, which will be saved in the file `all_xyz_bfdb.csv` (in the same columns format as other .csvs).

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