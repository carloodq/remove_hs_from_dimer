import pandas as pd
import os

def write_xyz(i, df, dataset):
    try:
        smiles  , k_index , abstracted_atom_index , netcharge,  multiplicity, bond_type = df['smiles_reduced']  , df['conformer'] , df['atom_index'],df['netcharge'], 2, df['attached_atom'] # multiplicity is 2 for fragment
        if dataset == "des_15k":
            k_index = round(int(k_index.replace(".xyz", "")) * 0.1, 3)     
        my_name =  dataset + "_frag" + str(i) + "_" + bond_type + "-H_" + str(k_index) + "_" + str(abstracted_atom_index) + ".xyz"
        content = df['xyz_reduced']
        comment = str(netcharge) + " " + str(multiplicity) + " " + smiles 

        first_line = content.split('\n')[0]
        after_comment = "\n".join(content.split('\n')[2:])

        f = open(os.path.join('first_batch', my_name), "w")
        f.write(first_line + "\n" + comment + "\n" +  after_comment)
        f.close()
    except:
        print("Writing failed for " + str(df) )

    return 0


def write_xyz_original(i, df, dataset):
    
    smiles , netcharge,  multiplicity = df['smiles_orig']  , df['netcharge'],  1 # multiplicity is 1 for original 
    my_name = dataset + "_orig"  + str(i) +    ".xyz"
    content = df['xyz']
    comment = str(netcharge) + " " + str(multiplicity) + " " + smiles 

    first_line = content.split('\n')[0]
    after_comment = "\n".join(content.split('\n')[2:])

    f = open(os.path.join('first_batch', my_name), "w")
    f.write(first_line + "\n" + comment + "\n" +  after_comment)
    f.close()


    return 0



# make des15k fragment and original
all_xyz = r"C:\Users\carlo\Dropbox\big projects dropbox\ml chem\bde\clean\remove_hs_from_dimer\all_xyz.csv"
df = pd.read_csv(all_xyz)
for i in range(5):
    write_xyz(i, df.iloc[i,:], 'des_15k')

df_originals = df.groupby(['xyz']).first().reset_index()


for i in range(1):
    write_xyz_original(i, df_originals.iloc[i,:], 'des15k')



# make nenci fragment and original
all_xyz = r"C:\Users\carlo\Dropbox\big projects dropbox\ml chem\bde\clean\remove_hs_from_dimer\all_xyz_nenci.csv"
df = pd.read_csv(all_xyz)
for i in range(3):
    write_xyz(i, df.iloc[i,:], 'nenci')

df_originals = df.groupby(['xyz']).first().reset_index()


for i in range(1):
    write_xyz_original(i, df_originals.iloc[i,:], 'nenci')



# make bfdb fragment and original
all_xyz = r"C:\Users\carlo\Dropbox\big projects dropbox\ml chem\bde\clean\remove_hs_from_dimer\all_xyz_bfdb.csv"
df = pd.read_csv(all_xyz)
for i in range(3):
    write_xyz(i, df.iloc[i,:], 'bfdb')

df_originals = df.groupby(['xyz']).first().reset_index()


for i in range(1):
    write_xyz_original(i, df_originals.iloc[i,:], 'bfdb')






