
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 21:28:55 2023
@author: bengi ruken yavuz
"""
import os
import pandas as pd
from scipy.stats import fisher_exact
import numpy as np
import statsmodels.stats.multitest

# Constants and File Paths
DATA_DIR = "Data"
OGTSG_FILE = os.path.join(DATA_DIR, "OncoKB_CancerGeneList_OG_TSG_Label.txt")
HYPER_NOT_HYPER_FILE = os.path.join(DATA_DIR, "Rev#2__Hyper_NotHyperLabellingPatient_Mutation_Counts_for_HyperMutated_parameter=8_IQR.txt")
MUTATION_DATA_FILE = os.path.join(DATA_DIR, "Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt")
POTENTIAL_DOUBLETS_FILE = os.path.join(DATA_DIR, "DifferentGene_Potential_Doublets_Tcga_GenieVol6_Nonhyper_GEQ3.txt")
DOUBLETS_COUNT_FILE = os.path.join(DATA_DIR, "DifferentGeneDoubletsTumorCounts_Tcga_GenieVol6_Nonhyper_GEQ3.txt")
STATISTICS_FILE = os.path.join(DATA_DIR, "DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6.txt")
FDR_CORRECTED_FILE = os.path.join(DATA_DIR, "DifferentGeneStatisticsAllDoublets_AmongAll_NonHyper_Mutants__Tcga_GenieVol6_FDR_corrected_Q_0_3.txt")
KRAS_OUTPUT_FILE = os.path.join(DATA_DIR, "kras_g12_tcga_genieVol6.txt")
PIK3CA_ESR1_OUTPUT_FILE = os.path.join(DATA_DIR, "PIK3CA_ESR1_tcga_genieVol6.txt")

def load_og_tsg_data(file_path):
    """
    Load oncogene and tumor suppressor gene list from a file.
    
    Args:
    file_path (str): Path to the OGTSG file.
    
    Returns:
    dict: Dictionary with gene symbols as keys and OG/TSG labels as values.
    list: List of unique gene symbols.
    """
    df = pd.read_csv(file_path, sep="\t")
    df = df[df['OG_TSG'].isin(["OG", "TSG"])]
    og_tsg_dict = dict(zip(df['Hugo Symbol'], df['OG_TSG']))
    og_tsg_list = list(set(df['Hugo Symbol'].to_list()))
    return og_tsg_dict, og_tsg_list

def load_hypermutation_data(file_path):
    """
    Load hypermutated and non-hypermutated sample data.
    
    Args:
    file_path (str): Path to the file containing hypermutation information.
    
    Returns:
    set: Set of hypermutated samples.
    set: Set of non-hypermutated samples.
    """
    HyperMut = set()
    NotHyperMut = set()
    with open(file_path, "r") as infile:
        for line in infile:
            splitted = line.rstrip("\n").split("\t")
            if splitted[-1] == "hypermutated":
                HyperMut.add(splitted[0])
            elif splitted[-1] == "not_hypermutated":
                NotHyperMut.add(splitted[0])
    return HyperMut, NotHyperMut

def load_all_tumor_data(file_path):
    """
    Load all tumor data from a file.
    
    Args:
    file_path (str): Path to the tumor data file.
    
    Returns:
    set: Set of all tumor data.
    """
    AllTumorData = set()
    with open(file_path, "r") as infile:
        for line in infile:
            if "Tumor" not in line:
                splitted = line.rstrip("\n").split("\t")
                AllTumorData.add(splitted[2])
    return AllTumorData

def process_mutation_data(file_path, og_tsg_list, NotHyperMut):
    """
    Process mutation data to form potential doublets from mutations observed in >=3 tumors.
    
    Args:
    file_path (str): Path to the mutation data file.
    og_tsg_list (list): List of oncogenes and tumor suppressor genes.
    NotHyperMut (set): Set of non-hypermutated samples.
    
    Returns:
    dict: Dictionary with gene mutations as keys and sets of tumors as values.
    """
    Gene_Mutation_List = {}
    GeneMut_TumorDict = {}
    with open(file_path, "r") as infile:
        for line in infile:
            if "Unkonwn" not in line:
                splitted = line.rstrip("\n").split("\t")
                gene = splitted[0]
                mutation = splitted[7]
                genemut = gene + "_" + mutation
                patient = splitted[2]
                
                if gene not in Gene_Mutation_List.keys():
                    Gene_Mutation_List[gene] = []
                    Gene_Mutation_List[gene].append(mutation)
                else:
                    Gene_Mutation_List[gene].append(mutation)
                    
                if genemut not in GeneMut_TumorDict.keys():
                    GeneMut_TumorDict[genemut] = set()
                    GeneMut_TumorDict[genemut].add(patient)
                else:
                    GeneMut_TumorDict[genemut].add(patient)
    return GeneMut_TumorDict

def filter_gene_mutations(GeneMut_TumorDict, NotHyperMut, min_samples=3):
    """
    Filter gene mutations to include only those present in at least a specified number of non-hypermutated samples.
    
    Args:
    GeneMut_TumorDict (dict): Dictionary with gene mutations as keys and sets of tumors as values.
    NotHyperMut (set): Set of non-hypermutated samples.
    min_samples (int): Minimum number of samples for filtering. Default is 3.
    
    Returns:
    dict: Filtered dictionary with gene mutations as keys and sets of tumors as values.
    """
    GeneMut_TumorDict_GEQ3 = {}
    for genemut in GeneMut_TumorDict.keys():
        if len(GeneMut_TumorDict[genemut].intersection(NotHyperMut)) >= min_samples:
            GeneMut_TumorDict_GEQ3[genemut] = GeneMut_TumorDict[genemut].intersection(NotHyperMut)
    return GeneMut_TumorDict_GEQ3

def find_potential_doublets(Mutation_List, GeneMut_TumorDict_GEQ3, output_file):
    """
    Find potential doublets from the list of mutations.
    
    Args:
    Mutation_List (list): List of mutations.
    GeneMut_TumorDict_GEQ3 (dict): Dictionary with filtered gene mutations.
    output_file (str): Path to the output file for potential doublets.
    
    Returns:
    int: Count of potential doublets.
    """
    DoubleCount = 0
    with open(output_file, "a") as outfile:
        for mut1 in Mutation_List:
            for mut2 in Mutation_List:
                if (mut1 < mut2) and (mut1.split("_")[0] != mut2.split("_")[0]):
                    if GeneMut_TumorDict_GEQ3[mut1].intersection(GeneMut_TumorDict_GEQ3[mut2]) != set():
                        DoubleCount += 1
                        outfile.write(mut1 + "\t" + mut2 + "\n")
    return DoubleCount

def calculate_statistics(GeneMut_TumorDict_GEQ3, NotHyperMut, input_file, output_file):
    """
    Calculate statistical significance for potential doublets using Fisher's exact test.
    
    Args:
    GeneMut_TumorDict_GEQ3 (dict): Dictionary with filtered gene mutations.
    NotHyperMut (set): Set of non-hypermutated samples.
    input_file (str): Path to the input file containing potential doublets.
    output_file (str): Path to the output file for statistical results.
    """
    with open(output_file, "a") as outfile:
        outfile.write("Mut1\tMut2\tDoubleMut#\tOnlyMut1\tOnlyMut2\tRemainingTumor#\tAllTumor#VAF>\tp4\tOddsR4\n")
    
    df = pd.read_csv(input_file, sep="\t")
    
    for i in df.index:
        a = df.iloc[i]['DoubleMut#']
        b = df.iloc[i]["OnlyMut1"]
        c = df.iloc[i]["OnlyMut2"]
        d4 = len(NotHyperMut) - (a + b + c)
        
        table4 = np.array([[a, b], [c, d4]])  # Values independent, evaluated among all tumors with VAF > 0.25
        oddsr4, p4 = fisher_exact(table4, alternative='two-sided')
        
        with open(output_file, "a") as outfile:
            outfile.write(df.iloc[i]['Mut1'] + "\t" + df.iloc[i]['Mut2'] + "\t" + str(a) + "\t" + str(b) + "\t" + str(c) + "\t" + str(d4)
                          + "\t" + str(len(NotHyperMut)) + "\t" + str(p4) + "\t" + str(oddsr4) + "\n")

def fdr_correction(input_file, output_file, alpha=0.5, q_threshold=0.3):
    """
    Apply FDR correction to p-values and filter results.
    
    Args:
    input_file (str): Path to the input file containing statistical results.
    output_file (str): Path to the output file for FDR corrected results.
    alpha (float): Significance level for FDR correction. Default is 0.5.
    q_threshold (float): Q-value threshold for filtering. Default is 0.3.
    """
    df = pd.read_csv(input_file, sep="\t")
    p_values = df['p4'].values
    reject, q_values = statsmodels.stats.multitest.fdrcorrection(p_values, alpha=alpha)
    
    df['Q_value'] = q_values
    df['Reject'] = reject
    filtered_df = df[df['Q_value'] <= q_threshold]
    
    filtered_df.to_csv(output_file, sep="\t", index=False)

def main():
    # Load data
    og_tsg_dict, og_tsg_list = load_og_tsg_data(OGTSG_FILE)
    HyperMut, NotHyperMut = load_hypermutation_data(HYPER_NOT_HYPER_FILE)
    AllTumorData = load_all_tumor_data(MUTATION_DATA_FILE)
    
    # Process mutation data
    GeneMut_TumorDict = process_mutation_data(MUTATION_DATA_FILE, og_tsg_list, NotHyperMut)
    GeneMut_TumorDict_GEQ3 = filter_gene_mutations(GeneMut_TumorDict, NotHyperMut)
    
    # Find potential doublets
    Mutation_List = list(GeneMut_TumorDict_GEQ3.keys())
    find_potential_doublets(Mutation_List, GeneMut_TumorDict_GEQ3, POTENTIAL_DOUBLETS_FILE)
    
    # Calculate statistics
    calculate_statistics(GeneMut_TumorDict_GEQ3, NotHyperMut, DOUBLETS_COUNT_FILE, STATISTICS_FILE)
    
    # Apply FDR correction
    fdr_correction(STATISTICS_FILE, FDR_CORRECTED_FILE)

if __name__ == "__main__":
    main()
