import pandas as pd
import os

def writexyz(content, path_and_name):
    f = open(path_and_name, "w")
    f.write(content)
    f.close()

def folders_from_csv(path = "all_xyz.csv"):
    df = pd.read_csv(path)
    latest_molecule = ""
    latest_conformer = ""
    for x in df.values:
        if x[2] != latest_molecule:
            os.mkdir(x[2])
        # if change in molecule or confermer I create new conformer folder (eg if confrmer is same but molecule if different)
        if "c_" + str(x[1]).split('.')[0] != latest_conformer or x[2] != latest_molecule:
            fdr = "c_" + str(x[1]).split('.')[0]
            os.mkdir(f"{x[2]}/{fdr}")
            writexyz(x[3], f"{x[2]}/{fdr}/original.xyz")
        latest_conformer = "c_" + str(x[1]).split('.')[0]
        latest_molecule = x[2]
        writexyz(x[4], f"{latest_molecule}/{latest_conformer}/{x[5]}.xyz")
    return 0


folders_from_csv()



