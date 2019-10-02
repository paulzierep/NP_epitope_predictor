from rdkit import Chem
import pandas as pd
import os
from rdkit.Chem import AllChem

from io import StringIO
import sys

def simple_sanitizer(smiles_df):
    for index, row in smiles_df.iterrows():

        mol = Chem.MolFromSmiles(row["smiles"])
        fp = AllChem.GetMorganFingerprint(
                                mol,
                                3,
                                useFeatures = False, 
                                useChirality = True,
                                ).GetNonzeroElements()

#df = pd.read_csv("ML_data/chebi_san_2.csv")
# simple_sanitizer(df)
# exit()


def Mol_Sanitizer(IN_CSV_PATH, OUT_CSV_PATH, ERROR_PATH):

    """
    Only check if the smiles can be converted to mol and back, 
    additional functionalities might get added (removal of salts,
    tautomeric conversion, neutralization ... )
    """

    in_csv_df = pd.read_csv(IN_CSV_PATH, index_col = "index")
    out_csv_df = pd.DataFrame()

    error_df = pd.DataFrame()

    #taken form http://rdkit.blogspot.com/2016/03/capturing-error-information.html
    Chem.WrapLogs()
    sio = sys.stderr = StringIO()

    for counter, (index, row) in enumerate(in_csv_df.iterrows()):

        #show process 
        if counter % 1000 == 0:
            print(counter)

        try:
            # row["smiles"] = Chem.MolToSmiles(Chem.MolFromSmiles(row["smiles"]))
            
            #do not actally convert the smiles, just check if it is possible,
            #the reason is, that the Chem.ForwardSDMolSupplier(SDF_PATH) covnerst them slightly different
            # which leads to inconsistancy with (paul notebook examples)
            # in general it would be OK to do so:
            # check differnt smils and fps for: CHEBI:51736 (C1C2C3*24567C1~C4~C5~C6~C7~3, C1C2*34567C1~C3~C4~C5~C6~C7~2)
            Chem.MolToSmiles(Chem.MolFromSmiles(row["smiles"])) 
            
            out_csv_df = out_csv_df.append(row)
            sio = sys.stderr = StringIO() # reset the error logger


        except Exception as exe:
            row["error"] = exe
            row["rdkit_error"] = sio.getvalue()
            error_df = error_df.append(row)
            sio = sys.stderr = StringIO() # reset the error logger

    #store the dataframes
    out_csv_df.index.name = "index"
    out_csv_df.to_csv(OUT_CSV_PATH)
    error_df.to_csv(ERROR_PATH)

    #redirect stderr to normal interpreter stderr after this function
    sys.stderr = sys.__stderr__

###################
#Test the Mol_Sanitizer
###################

# DATA_PATH = "ML_data"
# IN_CSV_PATH = os.path.join(DATA_PATH, "chebi_san.csv")
# OUT_CSV_PATH = os.path.join(DATA_PATH, "chebi_san2.csv")
# ERROR_PATH = os.path.join(DATA_PATH, "chebi_parse_errors2.csv")

# Mol_Sanitizer(IN_CSV_PATH, OUT_CSV_PATH, ERROR_PATH)

# exit()



def SDF_Parser(SDF_PATH, CSV_PATH, ERROR_PATH, prop_map = {"ChEBI ID":"index", 'ChEBI Name':"name"}, isomericSmiles=True):
    """
    Parser which creates a csv from a sdf file 
    and a error csv for file which cannot be parsed. 
    """


    mols = Chem.ForwardSDMolSupplier(SDF_PATH) #parse with rdkit

    mapped_props = [] #store smiles and props from the SDF
    parse_errors = [] #catch error for unparseable smiles

    #taken form http://rdkit.blogspot.com/2016/03/capturing-error-information.html
    Chem.WrapLogs()
    sio = sys.stderr = StringIO()

    for index, mol in enumerate(mols):
        
        #prop_mol = mol
        #mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))


        #show process 
        if index % 1000 == 0:
            print(index)

        if mol == None:

            pe = {}
            pe["index"] = index
            pe["mol_error"] = sio.getvalue()
            #pe["fp_error"] = None
            sio = sys.stderr = StringIO() # reset the error logger
            parse_errors.append(pe)

        else:

            smi = Chem.MolToSmiles(mol,isomericSmiles=isomericSmiles)

            #store all specified props, if they are not given an error is stored and the mol not saved
            mp = {}
            for mol_prop, csv_prop in prop_map.items():
                mp[csv_prop] = mol.GetProp(mol_prop)

            mp["smiles"] = smi

            mapped_props.append(mp)
            sio = sys.stderr = StringIO() # reset the error logger

    #store the dataframes
    csv_df = pd.DataFrame(mapped_props)
    csv_df.to_csv(CSV_PATH, index = False)

    error_df = pd.DataFrame(parse_errors)
    error_df.to_csv(ERROR_PATH, index = False)

    #redirect stderr to normal interpreter stderr after this function
    sys.stderr = sys.__stderr__

###################
#Test the SDF parser
###################

# DATA_PATH = "ML_data"
# SDF_PATH = os.path.join(DATA_PATH, "ChEBI_lite_3star.sdf")
# CSV_PATH = os.path.join(DATA_PATH, "chebi_san.csv")
# ERROR_PATH = os.path.join(DATA_PATH, "chebi_parse_errors.csv")

# SDF_parser(SDF_PATH, CSV_PATH, ERROR_PATH)

# exit()

def map_wrong_chabi_entries(SDF_PATH, ERROR_PATH):
    """
    This is not needed for the np_epitope pipeline,
    but it maps the parsing errors to the actual Chebi ID,
    which is nice for ChEBI if they want to check that. 
    """


    # error_df = pd.read_csv(ERROR_PATH, index_col = "index")
    # print(error_df)
    # exit()

    with open(SDF_PATH, "r") as i_file:
        lines = i_file.readlines()

    id_mapper = {}
    id_counter = 0

    for counter, line in enumerate(lines):

        if counter % 10000 == 0:
            print(counter)

        if "CHEBI:" in line:
            id_mapper[id_counter] = line.strip("\n")
            id_counter += 1

    error_df = pd.read_csv(ERROR_PATH, index_col = "index")
    error_df["new_id"] = error_df.index
    error_df["new_id"].replace(id_mapper, inplace = True)
    error_df.set_index("new_id", drop = True, inplace = True)
    error_df.to_csv(ERROR_PATH)

###################
#Test map_wrong_chabi_entries
###################

# DATA_PATH = "ML_data"
# SDF_PATH = os.path.join(DATA_PATH, "ChEBI_lite_3star.sdf")
# ERROR_PATH = os.path.join(DATA_PATH, "chebi_parse_errors.csv")

# map_wrong_chabi_entries(SDF_PATH, ERROR_PATH)

# exit()