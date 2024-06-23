import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

def load_double_mutants_data(file_path):
    """
    Loads the double mutants data from a file.

    Parameters:
    - file_path (str): Path to the file containing double mutants data.

    Returns:
    - pandas.DataFrame: DataFrame containing the loaded double mutants data.
    """
    full_path = os.path.join('data', file_path)
    df_doubles = pd.read_csv(full_path, sep="\t")
    return df_doubles


def filter_og_tsg(df_doubles):
    """
    Filters the double mutants DataFrame to select entries where Function1 is OG and Function2 is TSG,
    or Function1 is TSG and Function2 is OG.

    Parameters:
    - df_doubles (pandas.DataFrame): DataFrame containing double mutants data.

    Returns:
    - pandas.DataFrame: Filtered DataFrame containing OG-TSG couples.
    """
    df_og = df_doubles[((df_doubles["Function1"] == "OG") & (df_doubles["Function2"] == "TSG")) |
                       ((df_doubles["Function1"] == "TSG") & (df_doubles["Function2"] == "OG"))]
    df_og.reset_index(inplace=True, drop=True)
    return df_og


def create_og_couples(df_og):
    """
    Creates a set of OG-TSG couples from the filtered DataFrame.

    Parameters:
    - df_og (pandas.DataFrame): Filtered DataFrame containing OG-TSG couples.

    Returns:
    - set: Set of OG-TSG couples.
    """
    og_couples = set()
    for ind in df_og.index:
        mut1 = df_og.loc[ind, "Mut1"]
        mut2 = df_og.loc[ind, "Mut2"]
        g1 = mut1.split("_")[0]
        g2 = mut2.split("_")[0]
        og_couples.add((g1, g2))
    return og_couples


def load_gene_mutation_tumor_dict(file_path):
    """
    Loads the GeneMutation_tumor_dictionary from a pickle file.

    Parameters:
    - file_path (str): Path to the pickle file containing GeneMutation_tumor_dictionary.

    Returns:
    - dict: Loaded GeneMutation_tumor_dictionary.
    """
    full_path = os.path.join('data', file_path)
    with open(full_path, 'rb') as file:
        gene_mutation_tumor_dictionary = pickle.load(file)
    return gene_mutation_tumor_dictionary


def create_og_couples_tumor_dict(df_og, gene_mutation_tumor_dictionary, og_couples):
    """
    Creates a dictionary mapping OG-TSG couples to tumor sets based on GeneMutation_tumor_dictionary.

    Parameters:
    - df_og (pandas.DataFrame): Filtered DataFrame containing OG-TSG couples.
    - gene_mutation_tumor_dictionary (dict): Dictionary mapping mutations to tumor sets.
    - og_couples (set): Set of OG-TSG couples.

    Returns:
    - dict: Dictionary mapping OG-TSG couples to tumor sets.
    """
    og_couples_tumor_dict = {couple: set() for couple in og_couples}
    for ind in df_og.index:
        mut1 = df_og.loc[ind, "Mut1"]
        mut2 = df_og.loc[ind, "Mut2"]
        g1 = mut1.split("_")[0]
        g2 = mut2.split("_")[0]
        inters = gene_mutation_tumor_dictionary[mut1].intersection(gene_mutation_tumor_dictionary[mut2])
        og_couples_tumor_dict[(g1, g2)] = og_couples_tumor_dict[(g1, g2)].union(inters)
    return og_couples_tumor_dict


def load_tissue_tumor_dict(file_path):
    """
    Loads the Tissue_Tumor_Dict from a text file containing tissue information.

    Parameters:
    - file_path (str): Path to the text file containing tissue information.

    Returns:
    - dict: Dictionary mapping tissue to tumors.
    """
    full_path = os.path.join('data', file_path)
    tissue_tumor_dict = {}
    with open(full_path, "r") as infile:
        for line in infile:
            splitted = line.rstrip("\n").split("\t")
            tissue = splitted[2]
            if tissue not in tissue_tumor_dict.keys():
                tissue_tumor_dict[tissue] = set()
                tissue_tumor_dict[tissue].add(splitted[0])
            elif tissue in tissue_tumor_dict.keys():
                tissue_tumor_dict[tissue].add(splitted[0])
    return tissue_tumor_dict


def prepare_plotting_environment():
    """
    Prepares the plotting environment with specific font and SVG settings.

    No parameters or returns, modifies global matplotlib settings.
    """
    mpl.use('svg')
    new_rc_params = {
        "font.family": 'Arial',
        "font.size": 12,
        "font.serif": [],
        "svg.fonttype": 'none'
    }
    mpl.rcParams.update(new_rc_params)


def main():
    # File paths relative to 'data' folder
    double_mutants_file = "DiffGeneAnnotatedDoubles_Genie_vol6_Tcga.txt"
    gene_mutation_tumor_dict_file = "GeneMutation_tumor_dictionary.p"
    tissue_info_file = "NonhyperMutated_TCGA_GENIE_Vol6_VAF_0125_Subtype_Tissue_Primary_Metastasis_Info_ccBioPortal_API.txt"
    
    # Load data
    df_doubles = load_double_mutants_data(double_mutants_file)
    df_og = filter_og_tsg(df_doubles)
    og_couples = create_og_couples(df_og)
    gene_mutation_tumor_dictionary = load_gene_mutation_tumor_dict(gene_mutation_tumor_dict_file)
    og_couples_tumor_dict = create_og_couples_tumor_dict(df_og, gene_mutation_tumor_dictionary, og_couples)
    tissue_tumor_dict = load_tissue_tumor_dict(tissue_info_file)
    
    # Prepare plotting environment
    prepare_plotting_environment()
    
    # Further processing and plotting code...
    # Add your remaining code here for processing and plotting based on loaded data.


if __name__ == "__main__":
    main()
