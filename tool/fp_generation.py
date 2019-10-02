import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from io import StringIO
import sys
import os

def convert_to_count(val):
    if val == 0:
        return(0)
    else:
        return(1)

def Smiles_DF_To_Folded_FP_DF(smiles_df, 
                                OUTPATH = None, 
                                count_vector = False, 
                                nBits = 1024, 
                                radius = 3, 
                                useFeatures = False, 
                                useChirality = True, 
                                verbose = False,
                                **kwargs):

    result_df = pd.DataFrame(columns = list(range(nBits)))

    for index, (mol_id, row) in enumerate(smiles_df.iterrows()):

        #show process 
        if verbose:
            if index % 1000 == 0:
                print(index)

        mol = Chem.MolFromSmiles(row["smiles"])

        fp = AllChem.GetHashedMorganFingerprint(
                             mol,
                             radius, 
                             nBits = nBits,
                             useFeatures = useFeatures, 
                             useChirality = useChirality,
                             ).GetNonzeroElements()

        result_df.loc[mol_id,:] = 0
        result_df.loc[mol_id, fp.keys()] = fp.values()

    #simple way to convert count to bit
    if not count_vector:
        result_df = result_df.applymap(convert_to_count)

    #convert to int (saves memory)
    result_df = result_df.applymap(int)

    result_df.index.name = "index"

    if OUTPATH:
        result_df.to_csv(OUTPATH)

    return(result_df)


#################
#smiles_df_to_fp_df_bit (test)
#################


# DATA_PATH = "ML_data"
# SMILES_PATH = os.path.join(DATA_PATH, "chebi_san_assigned.csv")
# test_df = pd.read_csv(SMILES_PATH, index_col = "index")
# test_df = test_df.iloc[:100,:]

# result_df = smiles_df_to_fp_df_folded(test_df)

# exit()


def get_fp_molcount(smiles_df, radius = 3, useFeatures = False, useChirality = True, **kwargs):

    """
    For each Fingerprint (FP) get the count of mols which have this FP,
    the output CSV can be used to determine the FP cuf-off (the minimal
    FPs count needed to take this FP for the FP dataset)
    """

    fp_count = {}

    #create df which counts the occurences of each fingerprint
    for index, (mol_id, row) in enumerate(smiles_df.iterrows()):

        if index % 1000 == 0:
            print(""" *********** {0} *********** """.format(index))

        mol = Chem.MolFromSmiles(row["smiles"])

        fp = AllChem.GetMorganFingerprint(
                                        mol,
                                        radius,
                                        useFeatures = useFeatures, 
                                        useChirality = useChirality,
                                        ).GetNonzeroElements()

        #crate dict which counts the fps for all mols
        for idx, count in fp.items():
            if idx in fp_count.keys():
                fp_count[idx] += 1
            else:
                fp_count[idx] = 1

    FP_count_df = pd.DataFrame(pd.Series(fp_count), columns = ["count"])
    FP_count_df.index.name = "index"

    return(FP_count_df)


def get_FP_filtered(FP_count_df, cut_off = 10):
    """
    Returns the FPs which are occure in ore then "cut_off" mols
    """

    FP_to_use = FP_count_df.loc[FP_count_df["count"] >= cut_off].index.to_list()
    return(FP_to_use)

def create_all_fps(smiles_df, 
                    FP_to_use, 
                    OUTPATH, 
                    overwrite = True, 
                    count_vector = True, 
                    radius = 3,
                    useFeatures = False,
                    useChirality = True,
                    **kwargs):

    #otherwise other files might be appended
    if os.path.isfile(OUTPATH):

        if overwrite:
            os.remove(OUTPATH)
        else:
            raise Exception('{0} exists and overwrite kwarg is False'.format(OUTPATH))

    result_df_temp = pd.DataFrame(columns = FP_to_use)

    first_run = True

    #create df which counts the occurences of each fingerprint
    for index, (mol_id, row) in enumerate(smiles_df.iterrows()):

        if first_run:
            result_df = result_df_temp.copy()
            first_run = False

        mol = Chem.MolFromSmiles(row["smiles"])

        fp = AllChem.GetMorganFingerprint(
                                        mol,
                                        radius,
                                        useFeatures = useFeatures, 
                                        useChirality = useChirality,
                                        ).GetNonzeroElements()


        new_fp = {k: fp[k] for k in FP_to_use if k in fp}

        result_df.loc[mol_id,:] = 0
        result_df.loc[mol_id, new_fp.keys()] = new_fp.values()

        #store every 1000 mols
        if index % 1000 == 0 and index != 0:
            print("Dump mols up to index {0}".format(index))

            with open(OUTPATH, 'a') as f:
                do_header = (f.tell()==0)
                #do_header = False
                result_df.index.name = "index"
                result_df.to_csv(f, header=do_header)
                result_df = result_df_temp.copy()


    #store the rest
    print("Dump mols up to index {0}".format(index))

    with open(OUTPATH, 'a') as f:
        do_header = (f.tell()==0)
        #do_header = False
        result_df.index.name = "index"
        result_df.to_csv(f, header=do_header)
        result_df = result_df_temp.copy()


    final_df = pd.read_csv(OUTPATH, index_col = "index")

    final_df.index.name = "index"
    final_df = final_df.applymap(int)

    if not count_vector:
        final_df = final_df.applymap(convert_to_count)
        final_df.to_csv(OUTPATH)

    final_df.to_csv(OUTPATH)

    return(final_df)

def Smiles_DF_To_Unfolded_FP_DF(smiles_df, OUTPATH, cut_off = 10, df_memmory_warning = 25000000, **kwargs):

    """
    Chain the unfolded FP funs, pass kwargs to each.
    """

    result_df =  get_fp_molcount(smiles_df, **kwargs)
    FP_to_use = get_FP_filtered(result_df, cut_off = cut_off)

    print("FPs with cutoff {0} = {1}".format(cut_off, len(FP_to_use)))

    #warn if df gets to big (current warning 5000 rows * 5000 columns)
    value_count = len(FP_to_use) * smiles_df.shape[0]
    if value_count > df_memmory_warning:
        print("""With this FP cutoff, there will be {0} values in the created df, 
            this might lead to a memmory error on many machines (tested 16 GB ram),
            if this is the case incerase the cutoff or create smalles smiles df (clustering !)
            """.format(value_count))

    create_all_fps(smiles_df, FP_to_use, OUTPATH, **kwargs)

#################
#get_fp_molcount (test)
#################

# test_df = pd.DataFrame()
# test_df.loc[0,"smiles"] = "CCCCCCCCCCCC"
# #test_df.loc[1,"smiles"] = "CCCCCCCCCCCCX"
# test_df.loc[2,"smiles"] = "CCC"

# INPUT_PATH = "ML_data"
# TEST_PATH = os.path.join(INPUT_PATH, "chebi_san.csv")

# test_df = pd.read_csv(TEST_PATH)
# test_df = test_df.iloc[:2000,:]

# DATA_PATH = "FPs"
# OUTPATH = os.path.join(DATA_PATH, "FPs_count_all.csv")

# smiles_df_to_fp_df_unfolded(test_df, OUTPATH)

def Single_Smiles_To_Unfolded_FP_DF(smiles, 
                    FP_to_use, 
                    count_vector = True, 
                    radius = 3,
                    useFeatures = False,
                    useChirality = True,
                    verbose = False,
                    **kwargs):


    FP_to_use = [int(x) for x in FP_to_use]
    result_df = pd.DataFrame(columns = FP_to_use)

    mol = Chem.MolFromSmiles(smiles)

    fp = AllChem.GetMorganFingerprint(
                                    mol,
                                    radius,
                                    useFeatures = useFeatures, 
                                    useChirality = useChirality,
                                    ).GetNonzeroElements()

    #check if the fingerprints should be used, create new trimmed fp dict 
    #an count the not used fps
    new_fp = {}
    FP_not_used_count = 0
    for fp_id, fp_count in fp.items():
        if fp_id in FP_to_use:
            new_fp[fp_id] = fp_count
        else:
            FP_not_used_count += 1

    if verbose:
        print("Not used FPs: {0}".format(FP_not_used_count))

    # new_fp = {k: fp[k] for k in FP_to_use if k in fp}

    result_df.loc[0,:] = 0
    result_df.loc[0, new_fp.keys()] = new_fp.values()

    result_df.index.name = "index"
    result_df = result_df.applymap(int)

    if not count_vector:
        result_df = result_df.applymap(convert_to_count)

    return(result_df)

#################
#create all fps single (test)
#################

# test_fp_to_use = [26234434,266675433,368531313,864662311,1824088295,1842114921,2283643469,2549196227,2697642734,2905660137,2968968094,2976033787,3026394695,3083149185,3189457552,3192617127,3217380708,3218693969,3561006593,3698257053,10565946,516041383,864942730,1861965050,2117068077,2119439498,2246728737,2435602000,2976816164,515008442,533204632,864674487,1298690312,1338473746,1510328189,1535166686,2042559479,2245273601,2245384272,2246699815,2342113506,2807496773,3927890045,4003049590,4022716898,372101879,517457164,817588605,820930172,990985259,1276993226,1461258694,1842898132,2558947516,2629723425,3051231161,3117292559,3334525955,3577933307,4171172814,98513984,847433064,1016841875,2245900962,2424973678,2551483158,3452535345,3692055567,1506563592,1533864325,1542633699,1614704303,3542456614,3855312692,88831783,176403689,191340984,530227067,530250207,616155759,893111846,1016093624,1167268241,1167322652,1216095805,1237111057,1239802128,1705954535,2018127776,2246703798,2381125308,2410697484,2627605718,2648831588,2853515621,2856018280,3025629386,3202305582,3440734560,3545365497,3912994563,4089138501,3999906991,404368768,412434102,670706593,844214084,1206991059,1510461303,1542631284,1583009052,1583799011,1796099497,1823304278,1891071766,2067772118,2070976245,2245277810,2303564182,2423543607,2558947507,3051231160,3327625720,3375313142,3537119515,3537123720,3854544826,4025178337,4171172815,1710205153,2449355948,368713488,1313967653,1352399629,2076190208,2353112200,2944555726,2991151069,3145023872,3162837314,3631761933,3975275337,3982076256,3983062349,4194776273,4216335232,422715066,1429883190,2849741637,4121755354,584893129,847957139,1070821685,1083852209,1100037548,2041434490,2132511834,2297887526,2591432844,135162653,165937895,586457026,982949046,1324206100,1447784308,1759589175,1785683843,2412926691,2440473325,2440473326,2448950615,2456262944,2456262945,2667063169,3129492592,3194612514,4166791073,4166909777,396311007,994494548,1762382665,1772871781,1911793089,2163077160,2843136905,2955470977,3482813767,3516924900,3833608434,4096713990,369500293,379527158,742000539,1186671090,2246997334,2763854213,2984966880,3552309302,3696402029,3981762496,4197146180,265892994,737501821,1910766508,2157360080,97356054,847961216,920763749,1068491328,1396485351,1609756202,1858577693,1898829889,2181784098,2446961558,2460574971,2728180859,3153477100,3666284028,3721576064,3730747580,3792697566,4166778911,4207387074,140691206,328936174,352832185,400825689,540046244,553412256,1101907775,2092489639,2220038478,2237001212,2287269006,2296565457,2525418735,2750665496,2917424199,3176806076,3284564601,4008337421,4023654873,47607049,531116637,725322217,984189120,1171809476,1445852595,2192318254,2235918822,2257970297,2300762321,2805248962,2849870857,3345491492,3750444818,856674059,859799282,885218005,1275864840,1373661802,1688805770,1919750645,1956053223,1967344608,1978727030,2036328569,2181191047,2486452475,2834620058,3003632372,3153912292,3479423841,3544803785,3576922584,3628883864,3980805843,849275503,1054767590,2438720939,2669055056,639556223,690456263,909362583,1188701919,1214264900,1693331843,1961848637,2599973650,3351556771,3579962709,3925172229,3968746622,4222851645,404368769,926215114,1323888424,2181784125,2558119258,2566289208,2651067438,3083149214,2752034647,3118255683,3776905034,108183414,530219946,717512901,787069595,787132183,1038181128,1053116798,1173125914,1189967413,1211792670,1685248591,1849814881,1986202029,1998475861,2049407613,2065547089,2098603236,2127585437,2528877984,2967798987,3315826729,3745584548,4080822452,982949041,1732628033,3730747581,3824063894,505767871,541399149,1292826808,1353626381,1796099496,2636442481,3152170373,3344181102,3791047555,74537039,348315680,899522707,1182622762,1362518133,1506993418,1868697328,2154935424,2863170677,2939120473,4036774035,4278941385,161963127,589930590,800085213,932712697,1135286194,2222715027,2513505929,3599139274,4194366826,848128881,1067305901,2014255590,2636442480,2715536700,2811394787,3641942698,3657471097,3818546315,3921625505,416356657,642367021,831463242,997097697,2803848648,999334238,1583930800,2228063684,3510196525,199163361,3624155,1369588494,1674219510,2115476908,3120784405,847698334,2455552319,3824347764,2982494297,3375316278,3696389118,171200514,2142032900,2447748155,64070612,784020300,963553313,1792562719,1979311206,2516197204,2654043257,3120642300,3135357859,3328145258,4021328119,4078658161,1791271074,603510687,1136130381,1868602760,2784506312,2994748777,4041573576,369315086,3705139132,2185384887,2360741695,3976623167,571978829,1026928756,2246340824,2258843522,2604440622,2896319167,3848534832,2238101436,2479916646,3464861272,617798156,1010142783,1287474980,1790668568,1922435322,1978844073,2308348490,2592252298,3650394863,107748362,569172881,3337170256,951239203,3796841743,3850856377,4003054055,1073159142,847336149,977565923,2157595070,2806018737,2835079913,530226997,135162652,1043996772,2460574906,2473092670,2514004367,3639749812,3849240367,3881365835,3928910824,576060395,616374597,864942795,1510323402,2092554703,4037090848,185347322,979999009,2393016670,2782453342,3389297871,4086265842,570546812,2693580355,3194239760,2309124039,2827868305,3182824521,3731090,280868368,605229659,1098808928,1229209235,1506026072,1862632092,2006288109,2211525615,2432593989,2462011738,2529809651,2821901331,3062669668,3063495240,3087085494,3495750718,4239050898,95089064,1446633534,2245541727]

# df = create_all_fps_single("CCCCCCC", test_fp_to_use)

