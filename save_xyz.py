import pandas as pd
import os

def write_xyz(df, dataset):
    try:
        item_index, smiles  , k_index , abstracted_atom_index , netcharge,  multiplicity, bond_type = df['item_index'], df['smiles_reduced']  , df['conformer'] , df['atom_index'],df['netcharge'], 2, df['attached_atom'] # multiplicity is 2 for fragment
        if dataset == "des15k":
            k_index = round(int(k_index.replace(".xyz", "")) * 0.1, 3)     
        my_name =  bond_type + "-H_" + dataset + "_frag" + str(item_index) + "_" + str(k_index) + "_" + str(abstracted_atom_index) + ".xyz"
        content = df['xyz_reduced']
        comment = str(netcharge) + " " + str(multiplicity) + " " + smiles 

        first_line = content.split('\n')[0]
        after_comment = "\n".join(content.split('\n')[2:])

        
        xyztxt = first_line + "\n" + comment + "\n" +  after_comment
        f = open(os.path.join('first_batch', my_name), "w")
        f.write(xyztxt)
        f.close()
    except:
        print("Writing failed for " + str(df) )
        return 0, 0

    return my_name, xyztxt


def write_xyz_original(df, dataset):
    
    item_index, smiles , netcharge,  multiplicity = df['item_index'], df['smiles_orig']  , df['netcharge'],  1 # multiplicity is 1 for original 
    my_name = dataset + "_orig"  + str(item_index) +    ".xyz"
    content = df['xyz']
    comment = str(netcharge) + " " + str(multiplicity) + " " + smiles 

    first_line = content.split('\n')[0]
    after_comment = "\n".join(content.split('\n')[2:])

    xyztxt = first_line + "\n" + comment + "\n" +  after_comment
    f = open(os.path.join('first_batch', my_name), "w")
    f.write(xyztxt)
    f.close()


    return my_name, xyztxt



all_dfs = []

# make des15k fragment and original
all_xyz = r"C:\Users\carlo\Dropbox\big projects dropbox\ml chem\bde\clean\remove_hs_from_dimer\all_xyz_des15k.csv"
df = pd.read_csv(all_xyz)


def add_col(df):
    frag_my_name, frag_xyztxt = write_xyz(df, 'des15k')
    df['frag_filename'] = frag_my_name
    df['frag_xyztxt'] = frag_xyztxt
    orig_my_name, orig_xyztxt = write_xyz_original(df, 'des15k')
    df['orig_filename'] = orig_my_name
    df['orig_xyztxt'] = orig_xyztxt
    return df
    
df = df.apply(add_col,  axis='columns')
all_dfs.append(df)

print(len(all_dfs))



# make nenci fragment and original
all_xyz = r"C:\Users\carlo\Dropbox\big projects dropbox\ml chem\bde\clean\remove_hs_from_dimer\all_xyz_nenci.csv"
df = pd.read_csv(all_xyz)

def add_col(df):
    frag_my_name, frag_xyztxt = write_xyz(df, 'nenci')
    df['frag_filename'] = frag_my_name
    df['frag_xyztxt'] = frag_xyztxt
    orig_my_name, orig_xyztxt = write_xyz_original(df, 'nenci')
    df['orig_filename'] = orig_my_name
    df['orig_xyztxt'] = orig_xyztxt
    return df
    
df = df.apply(add_col,  axis='columns')
all_dfs.append(df)
print(len(all_dfs))



# make bfdb fragment and original
all_xyz = r"C:\Users\carlo\Dropbox\big projects dropbox\ml chem\bde\clean\remove_hs_from_dimer\all_xyz_bfdb.csv"
df = pd.read_csv(all_xyz)

def add_col(df):
    frag_my_name, frag_xyztxt = write_xyz(df, 'bfdb')
    df['frag_filename'] = frag_my_name
    df['frag_xyztxt'] = frag_xyztxt
    orig_my_name, orig_xyztxt = write_xyz_original(df, 'bfdb')
    df['orig_filename'] = orig_my_name
    df['orig_xyztxt'] = orig_xyztxt
    return df
    
df = df.apply(add_col,  axis='columns')
all_dfs.append(df)

print(len(all_dfs))


all_xyz = pd.concat(all_dfs)
all_xyz.to_csv('all_datasets.csv', index = None)


