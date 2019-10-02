import os
import shutil
import pickle

import pandas as pd

from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier

from fp_generation import Smiles_DF_To_Folded_FP_DF, Smiles_DF_To_Unfolded_FP_DF, Single_Smiles_To_Unfolded_FP_DF
from fp_drawing import DrawBitFromSmiles

from feature_stats import Compute_Feature_Enrichment_DF

from similarity_utils import Euclidiean_Distance_Query, Add_ChEBI_IP, RDKit_Mol_Similarity

from output_utils import Results_To_Html, Results_To_Json

#important for the df output (truncuation kills the SVG image)
#pd.set_option('display.max_colwidth', -1)

class Molecule_Group_Classifier:

    def __init__(self, n_clusters):

        self.n_clusters = n_clusters

        self.clf = KMeans(n_clusters=n_clusters, random_state=0, n_jobs = -1)

        #store specific informations for each class
        self.member_counts = {}
        self.warnings = {}
        self.names = {}
        self.ontos = {}

    def fit_predict(self, X):

        clusters = self.clf.fit_predict(X)

        self.member_counts = pd.Series(clusters).value_counts().to_dict()

        #assign empty dicts with the ids as keys
        self.warnings = {k:None for (k, v) in self.member_counts.items()}
        self.names = {k:None for (k, v) in self.member_counts.items()}
        self.ontos = {k:None for (k, v) in self.member_counts.items()}

        return(clusters)

    def assing_onto_mapping(self, ONTO_FOLDER):
        """
        Searches for ontology output from ChEBI in the specified folder and 
        adds the df to the clf information
        """

        for file in os.listdir(ONTO_FOLDER):
            for cluster_id in range(len(os.listdir(ONTO_FOLDER))):
                if str(cluster_id) in file:

                    f_path = os.path.join(ONTO_FOLDER, file) 

                    df = pd.read_csv(f_path, sep = "\t")
                    df = df.loc[(df["Corr-PValue"] <= 0.05),:]
                    df = df.sort_values(by = ["PValue", "Fold"], ascending = [True, False])
                    df = df.iloc[:10,:]

                    self.ontos[cluster_id] = df


    # def assing_warnings(self, warn_dict):
    #     self.warnings = warn_dict

    # def assing_names(self, name_dict):
    #     self.names = name_dict

    def predict(self, X):

        prediction = self.clf.predict(X)[0]

        cluster_infos = {}
        cluster_infos["Cluster"] = int(prediction)
        cluster_infos["Members"] = self.member_counts[prediction]
        cluster_infos["Warning"] = self.warnings[prediction]
        cluster_infos["Name"] = self.names[prediction]
        cluster_infos["Ontology"] = self.ontos[prediction]

        return(cluster_infos)


###################
#Test Epitope_Clf
###################

# DATA_PATH = "ML_data"
# FP_PATH = os.path.join(DATA_PATH, "fp_folded.csv")
# ONTO_PATH = os.path.join(DATA_PATH, "onto")
# fp_df = pd.read_csv(FP_PATH, index_col = "index")
# fp_df = fp_df.iloc[:100,:]

# fp_df.to_csv("test_data/fp_folded_test.csv")
# exit()
# mg_clf = Molecule_Group_Classifier(8)
# mg_clf.fit_predict(fp_df)
# mg_clf.assing_onto_mapping(ONTO_PATH)
# print(mg_clf.predict(fp_df.iloc[[1],:]))

# exit()

# mg_clf.

# print(mg_clf.cluster_infos)
# print(mg_clf.onto_dict)



class Epitope_Predictor:
    """
    Extends the normal sci-kit learn clf by adding extra attributes.

    Important: the fps_idx is used to store the fingerprint,
    they are needed to generate a matching DataFrame !
    """

    def __init__(self, clf, X, y, y_smiles):
        """
        Set the initial attributes
        """

        self.clf = clf
        self.cannot_be_fitted = False

        #although storing the fitting data takes a lot of space for
        #the pickle file, the similarity lookup and the FP
        #seach need the actual data to work with

        self.X = X
        self.y = y

        self.y_smiles = y_smiles

        #this info is returned, when the clf is used for prediction
        #fitting info : can the data be fitted (enough samples)
        #custom info : is the prediction meaningfull or based on strange data collection
        #class support
        self.clf_information = {"fitting_info":None, 
                                "p_support":None, 
                                "n_support":None,}

        #here the enrichment for this class is stored
        self.enrichment_df = None

        #store the fps indices (the column of X)
        #this allows to match the indices for new samples
        #as there is a cutoff for FPs which are found in less then 10 samples 
        #a new molecule can potentially posses a different fingerprint.
        self._fps_idx = self.X.columns.to_list()

    def fit(self):

        #store the sample support
        self.clf_information["p_support"] = self.y.loc[(self.y == 1)].shape[0]
        self.clf_information["n_support"] = self.y.loc[(self.y == 0)].shape[0]

        #currently we assume, that the negative samples are always enough
        #store info text for clf with to little positive samples
        p_samples = self.y.loc[(self.y == 1)].shape[0] 

        if p_samples <= 10:
            info_text = """There are not enough known epitopes ({0}) for this molecule group""".format(p_samples) 

            self.clf_information["fitting_info"] = info_text
            self.cannot_be_fitted = True

        else:
            self.clf.fit(self.X, self.y)


    def fit_FP_importance(self):
        """
        Computes the DF which stores the importance of the FPs for this cluster,
        this allows to compute a table of FPs for a predicted smiles
        """

        if not self.cannot_be_fitted:
            self.enrichment_df = Compute_Feature_Enrichment_DF(self.X,self.y)

    def get_FP_importance(self, smiles):
        """
        Adds the FPs of the query smiles to the FP df, allows to directly observe the 
        important FPs for this smiles
        """

        if self.cannot_be_fitted:
            return(pd.DataFrame())

        #get all fingerprints for the smiles (if they are part of _fps_idx) 
        X = Single_Smiles_To_Unfolded_FP_DF(smiles, self._fps_idx)
        X_sel = X.loc[0,:]
        X_sel = X_sel[(X_sel != 0)]
        #X_sel.rename(columns={0:"observed"}, inplace = True)

        #the function returns the FPs as int, but need to be str for concat with enrichment df 
        #either conversion would work
        X_sel.index = X_sel.index.astype("str")

        #combine the gernerated fingerprints with the FP info table
        FP_imp = pd.concat([X_sel, self.enrichment_df], join = "inner", axis = 1)
        FP_imp.rename(columns={0:"Observed"}, inplace = True)
        FP_imp.index.name = "Fingerprint ID"

        #add the figures
        FP_imp["Figure"] = FP_imp.index.map(lambda fp_id: DrawBitFromSmiles(smiles, fp_id))

        #the sorting gets lost due to concat (?)
        FP_imp = FP_imp.sort_values(by = ["Fold enrich."], ascending = False)

        FP_imp = FP_imp.T

        return(FP_imp)

    def get_FP_importance_with_sim(self, 
        smiles, 
        only_epitopes = False, 
        as_html = False, 
        compute_k_best = None, 
        show_k_best = 5,
        sort_order = "E",
        **kwargs,
        ):
        """
        Like get_FP_importance, but adds also similar
        compounds to the df, based on euclidean distance
        """

        if self.cannot_be_fitted:
            return(pd.DataFrame())

        FP_imp = self.get_FP_importance(smiles)

        query_fps = FP_imp.loc[["Observed"],:]

        if only_epitopes:
            target_fps = self.X.loc[(self.y == 1),query_fps.columns]
        else:
            target_fps = self.X.loc[:,query_fps.columns]

        similar_mols = Euclidiean_Distance_Query(query_fps, target_fps, k_best = compute_k_best)

        #add the real tanimoto sim for the top similar mols
        #get the corresponding smiles
        similar_smiles = self.y_smiles.loc[similar_mols.index]
        #use the basic rdkit function (no need to reinvent the wheel !)
        similar_mols["Tanimoto"] = similar_smiles.apply(lambda target: RDKit_Mol_Similarity(smiles, target))
        # #print(self.y_smiles.loc[similar_mols.index].apply(lambda target: RDKit_Mol_Similarity(smiles, target)))

        #lookup if the mol is a epitope
        similar_mols["Epitope"] = self.y.loc[similar_mols.index]

        similar_mols["IDs"] = similar_mols.index

        similar_mols.index = ["Target {0}".format(x) for x in range(similar_mols.shape[0])]

        #add link to CheBI
        similar_mols["IDs"] = similar_mols["IDs"].apply(Add_ChEBI_IP)

        if sort_order == "E":
            similar_mols = similar_mols.sort_values(by = ["E Dist.", "Tanimoto"], ascending = [True,False])
        if sort_order == "T":
            similar_mols = similar_mols.sort_values(by = ["Tanimoto", "E Dist."], ascending = [False,True])

        similar_mols = similar_mols.iloc[:show_k_best,:]

        imp_with_sim = pd.concat([FP_imp, similar_mols], sort=False)

        imp_with_sim.fillna("", inplace = True)

        if as_html:
            imp_with_sim = imp_with_sim.to_html(escape=False)

        return(imp_with_sim)

    def predict(self, smiles):
        """
        Predicts the probability of the smiles to be an epitope
        """

        if self.cannot_be_fitted:
            return(None)
        else:
            #create the FPs (X) for the smiles
            X = Single_Smiles_To_Unfolded_FP_DF(smiles, self._fps_idx)

            #predict the proba to be a epitope
            proba = self.clf.predict_proba(X)
            return(float(proba[0][1]))


###################
#Test Epitope_Clf
###################

# DATA_PATH = "ML_data"

# cluster_id, cell_type = 1, "b_cell"

# FP_PATH = os.path.join(DATA_PATH, "unfolded_fps")
# EXAMPLE_FP_PATH = os.path.join(FP_PATH, "{0}.csv".format(cluster_id))
# X = pd.read_csv(EXAMPLE_FP_PATH, index_col = "index")

# SMILES_PATH  = os.path.join(DATA_PATH, "target.csv")
# smiles_df = pd.read_csv(SMILES_PATH, index_col = "index")
# y = smiles_df.loc[X.index,:]

# clf = RandomForestClassifier(n_estimators = 1000, n_jobs = -1)
# e_clf = Epitope_Clf(clf, X, y[cell_type], y["smiles"], cluster_id, cell_type)

# print("fit 1")
# e_clf.fit()

# print("fit 2")
# e_clf.fit_FP_importance()

# ONTO_FOLDER = os.path.join(DATA_PATH, "onto")

# e_clf.assing_onto_mapping(ONTO_FOLDER)

# print(e_clf.clf_information)

# exit()

# smiles = 'CCCCCCCC/C=C\CCCCCCCC(=O)OC[C@@H]1COP(=O)(O)O1'

# print("get imp")
# results = e_clf.get_FP_importance_with_sim(smiles, only_epitopes = False)

# # print("sim mols")
# # results = e_clf.add_similar_mols2FP(FP_imp)

# pd.set_option('display.max_colwidth', -1)

# with open("test.html", "w") as out_file:
#     results.to_html(out_file, escape=False)

# proba = e_clf.predict(smiles)
# print(proba)
# exit()


class NP_Epitope_Prediction:

    def __init__(self, data_storage):
        """
        Set up all the paths needed for storage of the data,
        they are hard coded, only the main folder needs to be set
        """

        self.data_storage = data_storage

        #smiles, name, b_cell, t_cell assignment
        self.target_storage = os.path.join(data_storage, "target.csv")

        #folded FPs for clustering
        self.fp_folded_storage = os.path.join(data_storage, "fp_folded.csv")

        #folder used to store the FPs for classification
        self.fp_unfolded_storage = os.path.join(data_storage, "unfolded_fps")

        #folder for the pickled classifiers, only load and read to predict the cluster
        self.cluster_clf_storage = os.path.join(data_storage, "cluster_clf.pickle")

        #folder for the pickled classifiers, only load and read to predict the class
        self.classifier_clf_storage = os.path.join(data_storage, "classifier_clf")

        #dump the ontology IDs here
        self.ontology_folder_out = os.path.join(data_storage, "onto_out")

        #assign the onotogy from here
        self.ontology_folder_in = os.path.join(data_storage, "onto_in")

        #assign the onotogy from there
        self.custom_cluster_descri = os.path.join(data_storage, "cluster_descri.csv")

    def load_target(self, smiles_df):
        """
        Copies the smiles df to the storage folder, 
        good to keep it all together, in this df the clusters will be assigned
        """
        smiles_df.to_csv(self.target_storage)

    def compile_pre_cluster_FPs(self):
        """
        Generates the foldede FPs needed for clustering
        """

        smiles_df = pd.read_csv(self.target_storage, index_col = "index")
        Smiles_DF_To_Folded_FP_DF(smiles_df, self.fp_folded_storage, count_vector = False)


    def sanitize_pre_cluster_data(self):
        """
        Sortes the smiles and fps, Sklearn does not work with name indices but position indices
        sorting makes sure, that clustering and classification works correctly, here is also a checkup
        if all the indices match.
    
        Also keep in mind, that KMeans clustering is different for different sorted data, the sorting, therefore 
        allows to get consistent output for different runs.
        """

        target = pd.read_csv(self.target_storage, index_col = "index")
        fps = pd.read_csv(self.fp_folded_storage, index_col = "index")

        target.sort_index(inplace = True)
        fps.sort_index(inplace = True)

        print("Target shape: \t", target.shape)
        print("FPs shape: \t", fps.shape)
        print("All indices match: \t", (target.index == fps.index).all())

        target.to_csv(self.target_storage)
        fps.to_csv(self.fp_folded_storage)

    def compile_cluster_FPs(self):
        """
        Generates unfolded FPs for each cluster
        """

        smiles_df = pd.read_csv(self.target_storage, index_col = "index")
        
        #create the fp folder (do not overwrite a pre existing one)
        os.makedirs(self.fp_unfolded_storage, exist_ok= True)

        #creates one FP file for each cluster
        groups = smiles_df.groupby("clusters")

        for cluster, group in groups:

            path = os.path.join(self.fp_unfolded_storage, "{0}.csv".format(cluster))
            Smiles_DF_To_Unfolded_FP_DF(group, path)

    def get_cluster_stats(self):
        """
        Utililty function, prints a table which gives you the main infos,
        about mol and epitope types in each class 
        """

        target = pd.read_csv(self.target_storage)
        groups = target.groupby(by = "clusters")

        stats_df = pd.DataFrame()

        for cluster, group in groups:

            stats_df.loc[cluster, "Molecules"] = group.shape[0]
            stats_df.loc[cluster, "B cell"] = group.loc[(group["b_cell"] == 1),:].shape[0]
            stats_df.loc[cluster, "T cell"] = group.loc[(group["t_cell"] == 1),:].shape[0]

        stats_df = stats_df.astype(int)
        print(stats_df)

    def get_IDs_for_onto_mapping(self):
        """
        Utililty function to extract the ChEBI IDs needed to create
        the ontology enrichment.
        The ontology enrichment needs to be computed at:
        https://www.ebi.ac.uk/chebi/tools/binche/ by hand.
        """

        os.makedirs(self.ontology_folder_out, exist_ok= True)

        target = pd.read_csv(self.target_storage, index_col = "index")
        groups = target.groupby(by = "clusters")

        for cluster, group in groups:

            f_path = os.path.join(self.ontology_folder_out,"{0}.txt".format(cluster))

            with open(f_path, "w") as out_file:
                for acc in group.index.to_list():
                    out_file.write("{0}\n".format(acc))

    def fit_clustering(self, n_clusters = 8):
        """
        Fits the clustering algorithm based on the FPs, assings the
        FPs in the smiles_df (target_storage)
        """

        self.clusters = n_clusters

        target = pd.read_csv(self.target_storage, index_col = "index")
        fps = pd.read_csv(self.fp_folded_storage, index_col = "index")

        cluster_clf = Molecule_Group_Classifier(n_clusters=n_clusters)

        # cluster_clf =  KMeans(n_clusters=n_clusters, random_state=0, n_jobs = -1)
        target["clusters"] = cluster_clf.fit_predict(fps)
        target.to_csv(self.target_storage)

        #cluster_clf.assing_onto_mapping(self.ontology_folder_in)

        with open(self.cluster_clf_storage,'wb') as clf_f:
            pickle.dump(cluster_clf,clf_f)

    def update_onto_mapping(self):

        with open(self.cluster_clf_storage,'rb') as clf_f:
            cluster_clf = pickle.load(clf_f)

            cluster_clf.assing_onto_mapping(self.ontology_folder_in)

        with open(self.cluster_clf_storage,'wb') as clf_f:
            pickle.dump(cluster_clf,clf_f)

    def update_cluster_descri(self):

        descri_df = pd.read_csv(self.custom_cluster_descri, index_col = "index")

        with open(self.cluster_clf_storage,'rb') as clf_f:
            cluster_clf = pickle.load(clf_f)

            cluster_clf.names = descri_df.loc[:,"name"].to_dict()
            cluster_clf.warnings = descri_df.loc[:,"warning"].to_dict()

        with open(self.cluster_clf_storage,'wb') as clf_f:
            pickle.dump(cluster_clf,clf_f)

    def predict_clustering(self, smiles_to_predict):

        """
        Predicts the cluster for one smiles
        """

        #simpler not to set up bulk smiles processing
        smiles_df_to_predict = pd.DataFrame()
        smiles_df_to_predict.loc[0,"smiles"] = smiles_to_predict

        fps_df_to_predict = Smiles_DF_To_Folded_FP_DF(smiles_df_to_predict, None, store_data = False)

        with open(self.cluster_clf_storage,'rb') as clf_f:
            cluster_clf = pickle.load(clf_f)

        cluster_info = cluster_clf.predict(fps_df_to_predict)
        return(cluster_info)

    def fit_classification(self, cell_type = "b_cell"):
        """
        Fits a classifier for each cluster of the cell type and stores the 
        clf as pickle, the file name contains the cell_type and cluster,
        which can be used in the prediction method to select the correct clf.

        The clf is a custom classifier class, which can use any sci-kit learn clf,
        but also has some specific attributes which allows to store additional
        associated information.
        """

        #creates the folder to dump the pickled classifier (do not overwrite a pre existing one)
        os.makedirs(self.classifier_clf_storage, exist_ok= True)

        target = pd.read_csv(self.target_storage, index_col = "index")

        for fp_file in os.listdir(self.fp_unfolded_storage):

            #switch to ML terms here X (data, FPs), y (target, epitope or not)

            cluster_id = fp_file.strip(".csv")

            print("Fitting cluster {0}, cell type {1}".format(cluster_id, cell_type))

            #get the data
            path = os.path.join(self.fp_unfolded_storage, fp_file)
            X = pd.read_csv(path, index_col = "index")

            #get the target
            selected_target = target.loc[X.index,:]
            selected_target.loc[(selected_target[cell_type] == 1), "binary_target"] = 1
            selected_target.loc[(selected_target[cell_type] == 0), "binary_target"] = 0

            #target as binary series
            y = selected_target.loc[:,"binary_target"]
            y_smiles = selected_target.loc[:,"smiles"]

            #assign clf
            clf = RandomForestClassifier(n_estimators = 100, n_jobs = -1)

            #initiate custom epitope class
            e_clf = Epitope_Predictor(clf, X, y, y_smiles)

            e_clf.fit()

            #store as pickle
            out_path = os.path.join(self.classifier_clf_storage, \
                "cluster_{0}_type_{1}_clf.pickle".format(cluster_id, cell_type))

            with open(out_path,'wb') as clf_f:
                pickle.dump(e_clf,clf_f)

    def add_FP_imp_to_classification(self):
        """
        Adds the FP importance df to each classifier
        """

        for clf_f in os.listdir(self.classifier_clf_storage):

            print("Adding FP info to {0}".format(clf_f))

            #get the data
            path = os.path.join(self.classifier_clf_storage, clf_f)

            with open(path, "rb") as clf_f:
                e_clf = pickle.load(clf_f)
                e_clf.fit_FP_importance()

            with open(path,'wb') as clf_f:
                pickle.dump(e_clf,clf_f)

    def predict_classification(
        self, 
        smiles_to_predict, 
        cluster_id, 
        cell_type, 
        **kwargs,
        ):
        """
        Predict the probability of the smiles to be an epitope, based on the provided
        cluster_id and cell type
        """

        #search for the matching pickled classifier
        clf_to_use = None
        for clf_file in os.listdir(self.classifier_clf_storage):
            if str(cluster_id) in clf_file and cell_type in clf_file:
                clf_to_use = clf_file

        if not clf_to_use:
            raise Exception("Clf with cluster_id {0} and cell type {1} could not be found.".format(cluster_id, cell_type))
        else:

            #use clf to predict epitope probability
            path = os.path.join(self.classifier_clf_storage,clf_to_use)
            with open(path,'rb') as clf_f:
                e_clf = pickle.load(clf_f)

                classify_pred = {}
                classify_pred["proba"] = e_clf.predict(smiles_to_predict)
                classify_pred["fp_imp"] = e_clf.get_FP_importance_with_sim(
                    smiles_to_predict, 
                    **kwargs,
                    )

                #merge the predicted data with the default info of the clf
                classify_info = {**classify_pred, **e_clf.clf_information}

                return(classify_info)

    def prediction_chain(self, smiles, **kwargs):
        """
        Performs all steps for the epitope prediction:

            - predicts the cluster
            - predicts probability for each cell type

        """

        results = {}

        results["cluster_info"] = self.predict_clustering(smiles)
        results["classify_info_b"] = self.predict_classification(smiles, results["cluster_info"]["Cluster"], "b_cell", **kwargs)
        results["classify_info_t"] = self.predict_classification(smiles, results["cluster_info"]["Cluster"], "t_cell", **kwargs)

        return(results)

    def fit_chain(self, smiles_df):
        """
        Performs all fitting steps (FP generation, clf fitting)
        needed to set up the class for prediction
        each step below can also be run separately.

        Last run took 7262 seconds (2h)
        """

        print("Copy target df")
        self.load_target(smiles_df)

        print("Compile FPs (bit,folded) for clustering")
        self.compile_pre_cluster_FPs()

        print("Set correct order for FP and target")
        self.sanitize_pre_cluster_data()

        print("Fit the clusters")
        self.fit_clustering()

        print("Compile FPs (count, unfolded) for each cluster for classification")
        self.compile_cluster_FPs()

        print("Fit clf for b_cells")
        self.fit_classification(cell_type = "b_cell")

        print("Fit clf for t_cells")
        self.fit_classification(cell_type = "t_cell")

        print("Add FP importance to classification")
        self.add_FP_imp_to_classification()

        print("!!! DONE !!! \n Ready to predict.")


##########################
#Tests
##########################

#DATA_PATH = "ML_data"
# SMILES_PATH = os.path.join(DATA_PATH, "chebi_san_assigned.csv")
# smiles_df = pd.read_csv(SMILES_PATH, index_col = "index")
# smiles_df = smiles_df.iloc[:200,:]

# smiles_df.to_csv("test_data/chebi_san_assigned_for_tests.csv")

# exit()

# print("initiate class")
# predictor = NP_Epitope_Prediction(data_storage = DATA_PATH)
#predictor.fit_chain(smiles_df)

#predictor.get_IDs_for_onto_mapping()

#predictor.update_onto_mapping()
#exit()
#print("Data infos")
#predictor.get_data_infos()
#print("load target")
#predictor.load_target(smiles_df)
#print("precluster sanitizer")
#predictor.sanitize_pre_cluster_data()
# print("Compile pre cluster FPs")
# predictor.compile_pre_cluster_FPs()
#print("Fit clusterig")
#predictor.fit_clustering()
# print("Cluster infos")

# exit()

# predictor.get_cluster_stats()
# exit()
# #cluster = predictor.predict_clustering("CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC")
# #print(cluster)
# print("Compile Cluster FPs")
# predictor.compile_cluster_FPs()

# print("Fit classification b_cell")
# predictor.fit_classification(cell_type = "b_cell")

# print("Fit classification t_cell")
# predictor.fit_classification(cell_type = "t_cell")

# print("add_FP_imp")
# predictor.add_FP_imp_to_classification()

# print("Updating the cluster description")
# predictor.update_cluster_descri()

# print("add_ontolog")
# predictor.add_ontolog_to_classification()

# print("Predict classification")

##########################
#Test a bunch of smiles
##########################

# test_df = pd.read_csv(os.path.join(DATA_PATH, "target.csv"),index_col = "index")
# test_df = test_df.loc[test_df["b_cell"] == 1,:]

# res = []

# for counter, (index, row) in enumerate(test_df.iterrows()):
#     # print(counter)
#     results = predictor.prediction_chain(row["smiles"])
#     print(results)
#     # results["b_cell"] = row["b_cell"]
#     # results["t_cell"] = row["t_cell"]
#     # results["clusters"] = row["clusters"]
#     # results["ID"] = index
#     # res.append(results)

#     if counter == 5:
#         break

# res_df = pd.DataFrame(res)
# print(res_df)

#################################
#Check the fp alog
#################################

# test_df = pd.read_csv(os.path.join(DATA_PATH, "target.csv"),index_col = "index").loc["CHEBI:103210",:]
# test_fp = pd.read_csv(os.path.join(DATA_PATH, "fp_folded.csv"), index_col = "index").loc["CHEBI:103210",:]

# print(test_fp)

# smiles_df_to_predict = pd.DataFrame()
# smiles_df_to_predict.loc[0,"smiles"] = test_df["smiles"]

# fp_df = Smiles_DF_To_Folded_FP_DF(smiles_df_to_predict, None, store_data = False)
# print(fp_df)

#results = predictor.prediction_chain(test_df["smiles"])
#print(test_df)
#print(results)

##########################
#Test with any smiles
##########################

#smiles = "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O" #glucose
# smiles = "C[C@@H](C(=O)N[C@H](CCC(=O)N)C(=O)N[C@@H](CCC[C@@H](C(=O)O)N)C(=O)N[C@H](C)C(=O)N[C@H](C)C(=O)O)NC(=O)[C@@H](C)O[C@@H]1[C@H]([C@H](O[C@@H]([C@H]1O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)NC(=O)C)CO)O)NC(=O)C" #Peptidoglycan Pentapeptide-mdap3

# smiles = "CC(=O)SCCNC(=O)CCNC(=O)C(C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C(N=CN=C32)N)O)OP(=O)(O)O)O" #CoA

# smiles = "CCCCCCCCCCCCCCCCCC(=O)NCCCC[C@@H](NC(=O)CC[C@@H](NC(=O)[C@H](C)NC(=O)[C@@H](C)O[C@@H]1[C@@H](NC(C)=O)[C@H](O)O[C@H](CO)[C@H]1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1NC(C)=O)C(N)=O)C(O)=O"

# smiles = "CC1=C[C@H]2O[C@@H]3[C@H](O)C[C@@](C)([C@@]34CO4)[C@@]2(CO)[C@H](O)C1=O"

# smiles = "CC(=O)O[C@H]1C[C@@H]2CC[C@@H]3[C@H](CC[C@@]4(C)[C@H]3C[C@H]([N+]3(C)CCCCC3)[C@@H]4OC(C)=O)[C@@]2(C)C[C@@H]1N1CCCCC1"

# smiles = "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"

# smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@H](Cc1ccccc1)C(=O)N(C)[C@H](C(=O)N[C@H](C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](C)C(=O)OC)C(C)CC)C(C)C"

# smiles = "CCc1nn(C)c2c(=O)[nH]c(nc12)c3cc(ccc3OCC)S(=O)(=O)N4CCN(C)CC4" #viagra

# smiles = "C/C(=C\CO)/C=C/C=C(/C)\C=C\C1=C(C)CCCC1(C)C" #vit A

# smiles = "CO[C@H]1[C@@H](O)[C@H](NC(=O)[C@@H](O)CCO)[C@@H](C)O[C@@H]1O[C@H]1[C@@H](O)[C@H](NC(=O)[C@@H](O)CCO)[C@@H](C)O[C@@H]1O[C@H]1[C@@H](O)[C@H](NC(=O)[C@@H](O)CCO)[C@@H](C)O[C@@H]1O[C@H]1[C@@H](O)[C@H](NC(=O)[C@@H](O)CCO)[C@@H](C)O[C@@H]1O[C@@H]1[C@@H](O)O[C@H](C)[C@@H](NC(=O)[C@@H](O)CCO)[C@@H]1O"

# smiles = "Cc1onc(-c2c(F)cccc2Cl)c1C(=O)N[C@H](C(=O)NCCCC[C@H](N)C(=O)O)[C@@H]1N[C@@H](C(=O)O)C(C)(C)S1"

# results = predictor.prediction_chain(smiles, only_epitopes = False, compute_k_best = 5, show_k_best = None, sort_order = "E")

# json_text = Results_To_Json(results)

# with open("example_output.json", "w") as o_file:
#     o_file.write(json_text)
# print(json_text)

# html_text = Results_To_Html(results)

# with open("output-text.html", "w") as html_out:
#     html_out.write(html_text)

