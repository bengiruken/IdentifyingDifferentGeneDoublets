#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 18:31:13 2024

@author: yavuzb2
"""


import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import os


def load_rnaseq_data(folder = "Data"):
    """
    Load RNAseq data from the 'Data' folder. Returns data frame with RNAseq values for TCGA PAAD samples
    """
    file_path = os.path.join(folder, "PAAD_data_mrna_seq_v2_rsem.txt") #PAAD RNAseq values
    df = pd.read_csv(file_path, sep="\t")
    df.dropna(axis=0, how='any', inplace=True)
    df.set_index("Hugo_Symbol", drop=True, inplace=True)
    df.drop(columns=["Entrez_Gene_Id"], inplace=True)
    df = df.loc[~df.index.duplicated(keep='first')]
    return df
folder = 'data'
file_path = os.path.join(folder, "PAAD_data_mrna_seq_v2_rsem.txt")
file_path

df = pd.read_csv(file_path, sep="\t")

def get_kras_tp53_mutations(folder):
    """
    Get the list of KRASG12 and TP53 double mutations 
    """
    file_path = os.path.join(folder, "DiffGeneAnnotatedDoubles_Genie_vol6_Tcga.txt")
    df_doubles = pd.read_csv(file_path, sep="\t")
    AllDoubles = set()
    for ind in df_doubles.index:
        mut1 = df_doubles.loc[ind, "Mut1"]
        gene1 = mut1.split("_")[0]
        mut2 = df_doubles.loc[ind, "Mut2"]
        gene2 = mut2.split("_")[0]
        if (mut1 == "KRAS_G12") or (mut2 == "KRAS_G12"):
            if (gene1 == "TP53") or (gene2 == "TP53"):
                AllDoubles.add((mut1, mut2))
    return AllDoubles

get_kras_tp53_mutations(folder = 'Data')


def load_mutation_data(folder):
    """
    Load mutation data and return dictionaries tracking mutations and their associated tumors.
    """
    file_path = os.path.join(folder, "Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt")
    Gene_Mutation_List = {} #{Gene:[Mutated Residue List]}
    GeneMut_TumorDict = {} #{GeneMut: set(mutated tumors)}
    with open(file_path, "r") as infile:
        for line in infile:
            if ("Unknown" not in line) and ("Hugo" not in line):
                splitted = line.rstrip("\n").split("\t")
                gene = splitted[0]
                mutation = splitted[7]
                genemut = gene + "_" + mutation
                patient = splitted[2]
                if patient not in ['TCGA-FZ-5920-01', 'TCGA-FZ-5921-01','TCGA-FZ-5919-01','TCGA-FZ-5922-01', 'TCGA-FZ-5924-01']:
                    
                    if gene not in Gene_Mutation_List.keys():
                        Gene_Mutation_List[gene] = []
                    Gene_Mutation_List[gene].append(mutation)
    
                    if genemut not in GeneMut_TumorDict.keys():
                        GeneMut_TumorDict[genemut] = set()
                    GeneMut_TumorDict[genemut].add(patient)

    return {"Gene_Mutation_List":Gene_Mutation_List, "GeneMut_TumorDict":GeneMut_TumorDict}

x=load_mutation_data(folder = 'Data')
Gene_MutantTumors_Dict = load_mutation_data(folder = 'Data')["GeneMut_TumorDict"]


def get_gene_mutant_tumors(folder):
    """
    Form dictionary of genes and their associated mutant tumors.
    """
    GeneMut_TumorDict = load_mutation_data(folder = 'Data')["GeneMut_TumorDict"]
    Gene_Mutation_Tumor_GEQ3_Dict = {genemut:GeneMut_TumorDict[genemut] for genemut in GeneMut_TumorDict.keys() if len(GeneMut_TumorDict[genemut]) >= 1}
    #we only evaluate mutations observed on at least three tumors
    Gene_MutantTumors_Dict = {}
    for genemut in Gene_Mutation_Tumor_GEQ3_Dict.keys():
        gene = genemut.split("_")[0]
        patient = Gene_Mutation_Tumor_GEQ3_Dict[genemut]

        if gene not in Gene_MutantTumors_Dict.keys():
            Gene_MutantTumors_Dict[gene] = set()
        Gene_MutantTumors_Dict[gene] = Gene_MutantTumors_Dict[gene].union(patient)

    return Gene_MutantTumors_Dict
len(get_gene_mutant_tumors(folder = 'Data')["KRAS"])
len(get_gene_mutant_tumors(folder = 'Data')["PIK3CA"])


def load_hypermutated_data(folder= 'Data'):
    """
    Load data to distinguish between hypermutated and non-hypermutated tumors from the specified folder.
    """
    file_path = os.path.join(folder, "Rev#2__Hyper_NotHyperLabellingPatient_Mutation_Counts_for_HyperMutated_parameter=8_IQR.txt")
    HyperMut = set()
    NotHyperMut = set()
    with open(file_path, "r") as infile:
        for line in infile:
            splitted = line.rstrip("\n").split("\t")
            if splitted[-1] == "hypermutated":
                HyperMut.add(splitted[0])
            elif splitted[-1] == "not_hypermutated":
                NotHyperMut.add(splitted[0])
    return {'HyperMut':HyperMut, 'NotHyperMut':NotHyperMut}
y = load_hypermutated_data()
len(y['HyperMut'])


def load_subtype_tumors(folder= 'Data'):
    """
    Load tumor subtype information as a dictionary {Subtype:set(Tumors)}
    """
    file_path = os.path.join(folder, "TCGA_Primary_Metastasis_Info.txt")
    Subtype_Tumor_Dict = {}
    with open(file_path, "r") as infile:
        for line in infile:
            if 'ONCOTREE' not in line:
                splitted = line.rstrip("\n").split("\t")
                if splitted[2] not in Subtype_Tumor_Dict.keys():
                    Subtype_Tumor_Dict[splitted[2]] = set()
                Subtype_Tumor_Dict[splitted[2]].add(splitted[1])
    return Subtype_Tumor_Dict
z=load_subtype_tumors()
AllDoubles = get_kras_tp53_mutations(folder = 'Data')
paad = load_subtype_tumors()['PAAD']
len(paad)
def get_kras_tp53_double_mutants(folder= 'Data'):
    """
    Identify PAAD tumors with both KRAS_G12 and TP53 mutations.
    """
    GeneMut_TumorDict = load_mutation_data(folder = 'Data')["GeneMut_TumorDict"]
    Gene_Mutation_Tumor_GEQ3_Dict = {genemut:GeneMut_TumorDict[genemut] for genemut in GeneMut_TumorDict.keys() if len(GeneMut_TumorDict[genemut]) >= 1}
    #we only evaluate mutations observed on at least three tumors
    
    AllDoubles = get_kras_tp53_mutations(folder = 'Data')
    paad = load_subtype_tumors()['PAAD']
    KRAS_TP53_DoubleMutants = set()
    for double in AllDoubles:
        Doublet = Gene_Mutation_Tumor_GEQ3_Dict[double[0]].intersection(Gene_Mutation_Tumor_GEQ3_Dict[double[1]])
        KRAS_TP53_DoubleMutants = (KRAS_TP53_DoubleMutants.union(Doublet)).intersection(paad)
    return KRAS_TP53_DoubleMutants

q=get_kras_tp53_double_mutants()

def get_kras_single_mutants(folder = 'Data'):
    """
    Identify PAAD tumors with only KRAS_G12 mutations and no TP53 mutations.
    """
    GeneMut_TumorDict = load_mutation_data(folder = 'Data')["GeneMut_TumorDict"]
    Gene_Mutation_Tumor_GEQ3_Dict = {genemut:GeneMut_TumorDict[genemut] for genemut in GeneMut_TumorDict.keys() if len(GeneMut_TumorDict[genemut]) >= 1}
    #we only evaluate mutations observed on at least three tumors
    
    Gene_All_Mutants = get_gene_mutant_tumors(folder)
    KRAS_TP53_DoubleMutants = get_kras_tp53_double_mutants()
    paad = load_subtype_tumors()['PAAD']
    KRAS_Single_Mutant = list(((Gene_Mutation_Tumor_GEQ3_Dict["KRAS_G12"].union(Gene_All_Mutants["TP53"])).difference(KRAS_TP53_DoubleMutants)).intersection(paad))
    return KRAS_Single_Mutant

k = get_kras_single_mutants()

df = load_rnaseq_data(folder = "Data")
column_single1 = list(get_kras_single_mutants())
column_double1 = list(get_kras_tp53_double_mutants())
##get only 'TCGA' samples
column_single = [tumor for tumor in column_single1 if tumor.startswith('TCGA-')]
column_double = [tumor for tumor in column_double1 if tumor.startswith('TCGA-')]
#
set(column_single).difference({'TCGA-FZ-5922-01', 'TCGA-FZ-5924-01'})
df_single=df[list(set(column_single).difference({'TCGA-FZ-5922-01', 'TCGA-FZ-5924-01'}))]
df_double=df[list(set(column_double).difference({'TCGA-FZ-5920-01', 'TCGA-FZ-5921-01','TCGA-FZ-5919-01'}))]

import pandas as pd
import numpy as np
from scipy import stats


df_significants = pd.DataFrame(columns=["gene", "pval", "log2FC"])
significants = []
problem = []

for gene in df_single.index:
    S = df_single.loc[gene].values.tolist()
    D = df_double.loc[gene].values.tolist()

    np_S = np.array(S)
    np_D = np.array(D)

    # Avoid division by zero and log of zero issues
    if np_S.mean() == 0 or np_D.mean() == 0:
        continue

    logFC = np.log2(np_D.mean() / np_S.mean())

    try:
        statis = stats.mannwhitneyu(D, S)
        if (statis.pvalue <= 0.05) and (np.abs(logFC) > 0.5):
            to_append = pd.DataFrame([[gene, statis.pvalue, logFC]], columns=df_significants.columns)
            df_significants = pd.concat([df_significants, to_append], ignore_index=True)
            significants.append(gene)
    except ValueError:
        problem.append(gene)




def find_significant_genes(folder = "Data"):
    """
    Identify significantly differentially expressed genes between single and double mutants. 
    return TPM values as a data frame
    """
    df = load_rnaseq_data(folder = "Data")
    column_single1 = list(get_kras_single_mutants())
    column_double1 = list(get_kras_tp53_double_mutants())
    ##get only 'TCGA' samples
    column_single = [tumor for tumor in column_single1 if tumor.startswith('TCGA-')]
    column_double = [tumor for tumor in column_double1 if tumor.startswith('TCGA-')]
    #
    df_single=df[list(set(column_single).difference({'TCGA-FZ-5922-01', 'TCGA-FZ-5924-01'}))]
    df_double=df[list(set(column_double).difference({'TCGA-FZ-5920-01', 'TCGA-FZ-5921-01','TCGA-FZ-5919-01'}))]
    
    
    df_significants = pd.DataFrame(columns=["gene", "pval", "log2FC"])
    significants = []
    problem = []

    for gene in df_single.index:
        S = df_single.loc[gene].values.tolist()
        D = df_double.loc[gene].values.tolist()

        np_S = np.array(S)
        np_D = np.array(D)

        # Avoid division by zero and log of zero issues
        if np_S.mean() == 0 or np_D.mean() == 0:
            continue

        logFC = np.log2(np_D.mean() / np_S.mean())

        try:
            statis = stats.mannwhitneyu(D, S)
            if (statis.pvalue <= 0.05) and (np.abs(logFC) > 0.5):
                to_append = pd.DataFrame([[gene, statis.pvalue, logFC]], columns=df_significants.columns)
                df_significants = pd.concat([df_significants, to_append], ignore_index=True)
                significants.append(gene)
        except ValueError:
            problem.append(gene)

    return df_significants

p=find_significant_genes()

####denemeE
df = load_rnaseq_data(folder = "Data")
column_single1 = list(get_kras_single_mutants())
column_double1 = list(get_kras_tp53_double_mutants())
##get only 'TCGA' samples
column_single = [tumor for tumor in column_single1 if tumor.startswith('TCGA-')]
column_double = [tumor for tumor in column_double1 if tumor.startswith('TCGA-')]
#
df_single=df[list(set(column_single).difference({'TCGA-FZ-5922-01', 'TCGA-FZ-5924-01'}))]
df_double=df[list(set(column_double).difference({'TCGA-FZ-5920-01', 'TCGA-FZ-5921-01','TCGA-FZ-5919-01'}))]

All_Columns = list(column_single) + list(column_double)
df1 = df[All_Columns]
genes = find_significant_genes()['gene'].to_list()
df1 = df1[df1.index.isin(genes)]
df1T = df1.T
df_new1 = pd.DataFrame()
for col in df1T.columns:
    try:
        df_new1[col] = (df1T[col] - df1T[col].mean()) / df1T[col].std(ddof=0)
    except KeyError:
        print(f"KeyError in column: {col}")
df1_zscore = df_new1.T


def create_zscore_dataframe(folder = "Data"):
    """
    Create a z-score normalized dataframe for specified genes and columns.
    """
    df = load_rnaseq_data(folder = "Data")
    column_single1 = list(get_kras_single_mutants())
    column_double1 = list(get_kras_tp53_double_mutants())
    ##get only 'TCGA' samples
    column_single = [tumor for tumor in column_single1 if tumor.startswith('TCGA-')]
    column_double = [tumor for tumor in column_double1 if tumor.startswith('TCGA-')]
    #
    df_single=df[list(set(column_single).difference({'TCGA-FZ-5922-01', 'TCGA-FZ-5924-01'}))]
    df_double=df[list(set(column_double).difference({'TCGA-FZ-5920-01', 'TCGA-FZ-5921-01','TCGA-FZ-5919-01'}))]
    
    All_Columns = list(column_single) + list(column_double)
    df1 = df[All_Columns]
    genes = find_significant_genes()['gene'].to_list()
    df1 = df1[df1.index.isin(genes)]
    df1T = df1.T
    df_new1 = pd.DataFrame()
    for col in df1T.columns:
        try:
            df_new1[col] = (df1T[col] - df1T[col].mean()) / df1T[col].std(ddof=0)
        except ValueError:
            print(f"ValueError in column: {col}")
    df1_zscore = df_new1.T
    return df1_zscore

x=create_zscore_dataframe()
########
#df_zscore, genes, column_double, column_single,
def generate_heatmap( output_file):
    """
    Generate a heatmap from the z-score normalized dataframe.
    """
    df = load_rnaseq_data(folder = "Data")
    column_single1 = list(get_kras_single_mutants())
    column_double1 = list(get_kras_tp53_double_mutants())
    ##get only 'TCGA' samples
    column_single = [tumor for tumor in column_single1 if tumor.startswith('TCGA-')]
    column_double = [tumor for tumor in column_double1 if tumor.startswith('TCGA-')]
    #
    genes = find_significant_genes()['gene'].to_list()
    df_zscore = create_zscore_dataframe()
    
   
    
    df_gene = df_zscore[column_double]
    df_gene['mean'] = df_gene.mean(axis=1)
    positives = [ind for ind in df_gene.index if df_gene.loc[ind]["mean"] > 0]
    negatives = [ind for ind in df_gene.index if df_gene.loc[ind]["mean"] <= 0]
   
    
    ORDER = positives + negatives

    df_zscore = df_zscore.reindex(ORDER)
    df2 = pd.DataFrame(df_zscore.T.to_dict()).reset_index().melt(id_vars=['index'], var_name='Gene', value_name='Score')
    df2.columns = ['Patient', 'Gene', 'Score']

    heatmap1_data = pd.pivot_table(df2, values='Score', index=['Gene'], columns='Patient')
    heatmap1_data = heatmap1_data.reindex(ORDER)

    plt.figure(figsize=(0.50,0.8))
    sns.heatmap(heatmap1_data, center=0, cmap='coolwarm', vmin=-2, vmax=2)
    plt.savefig(output_file, bbox_inches="tight", dpi=1000)




if __name__ == "__main__":
    # Define file paths and folders
    data_folder = "Data"
    rna_seq_file = os.path.join(data_folder, "PAAD_data_mrna_seq_v2_rsem.txt")
    mutation_file = os.path.join(data_folder, "Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt")
    hypermutation_file = os.path.join(data_folder, "Rev#2__Hyper_NotHyperLabellingPatient_Mutation_Counts_for_HyperMutated_parameter=8_IQR.txt")
    
   
    # Load data
    print("Loading RNAseq data...")
    df = load_rnaseq_data(data_folder)
    
    print("Loading KRASG12 and TP53 mutations...")
    AllDoubles = get_kras_tp53_mutations(data_folder)
    
    print("Loading mutation data...")
    Gene_Mutation_List, GeneMut_TumorDict = load_mutation_data(data_folder)
    
    print("Loading hypermutated data...")
    HyperMut, NotHyperMut = load_hypermutated_data(data_folder)
    
    print("Loading subtype tumors...")
    Subtype_Tumor_Dict = load_subtype_tumors(data_folder)
    
    # Define genes of interest
    paad = set(df.columns)
    KRAS_TP53_DoubleMutants = get_kras_tp53_double_mutants(data_folder)
    KRAS_Single_Mutant = get_kras_single_mutants(data_folder)
    
    print("Finding significant genes...")
    df_single = df[list(KRAS_Single_Mutant)]
    df_double = df[list(KRAS_TP53_DoubleMutants)]
    df_significants = find_significant_genes(data_folder)
    genes = df_significants['gene'].tolist()
    
    print("Creating z-score dataframe...")
    df_zscore = create_zscore_dataframe(data_folder)
    
    # Generate heatmap
    output_file = "heatmap.png"
    print("Generating heatmap...")
    generate_heatmap(output_file)
    
    print(f"Heatmap saved to {output_file}")

