from mol_sanitizer import SDF_Parser, Mol_Sanitizer, filter_duplicates
import os
import pandas as pd

class NP_Epitope_Data_Conversion:
    """
    creates the data for the np_epitope_prediction pipeline
    final output should be a CSV with format of

    ,index,smiles,b_cell,t_cell
    0,CHEBI:90,Oc1cc(O)c2c(c1)O[C@H](c1ccc(O)c(O)c1)[C@H](O)C2,0.0,0.0
    1,CHEBI:165,CC1(C)C(=O)[C@@]2(C)CC[C@@H]1C2,0.0,0.0
    2,CHEBI:598,*C(=O)OC(CO)CO[1*],0.0,0.0
    """


    def __init__(self, OUTPUT_FOLDER):
        """
        Set up all the path for the files to create
        """

        self.OUTPUT_FOLDER = OUTPUT_FOLDER

        # since the SDF parser cannot catch all sanitizing errors a second
        # step is needed, hence two sanitizing and error files
        self.SDF_CSV_PATH_1 = os.path.join(OUTPUT_FOLDER, "chebi_san_1.csv")
        self.SDF_CSV_PATH_2 = os.path.join(OUTPUT_FOLDER, "chebi_san_2.csv")

        self.SDF_ERROR_PATH_1 = os.path.join(OUTPUT_FOLDER, "chebi_parse_errors_1.csv")
        self.SDF_ERROR_PATH_2 = os.path.join(OUTPUT_FOLDER, "chebi_parse_errors_2.csv")
        #self.SDF_UNCOMPLETE  =  os.path.join(OUTPUT_FOLDER, "chebi_uncomplete.csv")
        self.SDF_DUBS_PATH = os.path.join(OUTPUT_FOLDER, "chebi_doublicates.csv")
        self.OUTPUT_PATH  =  os.path.join(OUTPUT_FOLDER, "chebi_san_assigned.csv")

    def create_ML_data(self, SDF_PATH, B_CELL_PATH, T_CELL_PATH, skip_sdf = False, **kwargs):
        """
        Run the data creation steps as a pipeline
        SDF should look like the one downloaded from ChEBI 
        website: ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_lite_3star.sdf.gz
        kwargs can be passed to the SDF_Parser in order to 
        change stuff like isomericSmiles parsing (default 
        kwargs should be good !)

        B_CELL and T_CELL should be a CSV created by query of the IEDB, (basically 
        only one row is needed which has the ChEBI ID, default is: "Non-peptidic epitope Accession",
        can be changed with the id_column kwarg)
        """

        #ChEBI SDF to CSV (only index and smiles are needed), (store the parsing errors)
        #can be skipped when the CSV data was created outside of this pipe
        if not skip_sdf:
            self.convert_chebi_sdf(SDF_PATH, **kwargs)

        #load the csv, drop duplicates (store them)
        self.load_chebi_csv(**kwargs)

        #assign each essay type
        self.assign_epitopes(B_CELL_PATH, "b_cell", **kwargs)
        self.assign_epitopes(T_CELL_PATH, "t_cell", **kwargs)

        #store the final df
        #can be passed like that to the np_epitope_prediction pipeline
        self.store_ML_df()

    def convert_chebi_sdf(self, SDF_PATH, **kwargs):
        SDF_Parser(SDF_PATH, self.SDF_CSV_PATH_1, self.SDF_ERROR_PATH_1, **kwargs)
        Mol_Sanitizer(self.SDF_CSV_PATH_1, self.SDF_CSV_PATH_2, self.SDF_ERROR_PATH_2, **kwargs)

    def load_chebi_csv(self, **kwargs):
        csv_df = pd.read_csv(self.SDF_CSV_PATH_2, index_col = "index")

        print("Input rows: {0}".format(csv_df.shape[0]))

        no_doubs, doub_stats = filter_duplicates(csv_df)

        #store doublicates
        # csv_df["is_doublicated"] = csv_df.duplicated(subset = "smiles")
        doub_stats.to_csv(self.SDF_DUBS_PATH, index = True)

        #################################
        #Optional: Remove smiles with a *
        #################################

        print("Input rows without duplicates: {0}".format(no_doubs.shape[0]))

        self.csv_df_trim = no_doubs.loc[:,["smiles"]]

    def assign_epitopes(self, CSV_PATH, epitope_type, id_column = "Non-peptidic epitope Accession"):

        csv_df = pd.read_csv(CSV_PATH)

        #epitope ID when queried from IEDB are of type CHEBI_90, convert to CHEBI:90
        accs = [acc.replace("_",":") for acc in csv_df[id_column]]

        #assing the epitopes binary in a new column
        self.csv_df_trim.loc[self.csv_df_trim.index.isin(accs), epitope_type] = 1
        self.csv_df_trim[epitope_type].fillna(0, inplace=True)

    def store_ML_df(self):
        #self.csv_df_trim.set_index("index",  drop=True, inplace = True)
        self.csv_df_trim.to_csv(self.OUTPUT_PATH)


###################
#Test the np_epitope_data_conversion
###################

# DATA_PATH = "ML_data"
# SDF_PATH = os.path.join(DATA_PATH, "ChEBI_lite_3star.sdf")
# B_CELL = os.path.join(DATA_PATH, "epitope_table_b_cell_pos.csv")
# T_CELL = os.path.join(DATA_PATH, "epitope_table_t_cell_pos.csv")

# converter = NP_Epitope_Data_Conversion(DATA_PATH)
# converter.create_ML_data(SDF_PATH, B_CELL, T_CELL, skip_sdf = True)

# exit()
# converter.load_chebi_csv()    
# converter.assign_epitopes(B_CELL, "b_cell")
# converter.assign_epitopes(B_CELL, "t_cell")
# print(converter.csv_df_trim)