Certainly! Below is the script with added comments explaining each function and their arguments.

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 21:10:32 2023
Updated on Mon Jun 23 18:00:00 2024

@author: bengi
"""

import pandas as pd
import numpy as np
import networkx as nx

def load_data():
    """
    Load all necessary data files.

    Returns:
    df_ppi (DataFrame): Protein-protein interaction data.
    df_double_comps (DataFrame): Double mutant genes network data.
    df_uniprot (DataFrame): UniProt to gene symbol mapping data.
    df_uni (DataFrame): UniProt gene name conversion data.
    df_trrust (DataFrame): TRRUST transcription factor and target genes data.
    df_double_couples (DataFrame): Double mutant genes couples data.
    df_og_tsg (DataFrame): Oncogene and tumor suppressor gene list data.
    """
    df_ppi = pd.read_csv("data/OmniPath_Directed_Interactions.txt", sep="\t", low_memory=False)
    df_double_comps = pd.read_csv("data/DiffGene_doublets_component_genes_Network.txt", sep="\t")
    df_uniprot = pd.read_excel("data/GeneSymbolToUniprotDoubleComponents.xlsx")
    df_uni = pd.read_csv("data/UniProtGeneNameConversion.txt", sep="\t", dtype=str)
    df_trrust = pd.read_csv("data/Subnetwork_GeneSymbol_TRRUST_input_THR0_0_001key_regulators.tsv", sep="\t")
    df_double_couples = pd.read_excel("data/DiffGene_doublets_Couples_DoubleMutantGenes_Network.xlsx")
    df_og_tsg = pd.read_csv("data/OncoKB_CancerGeneList_OG_TSG_Label.txt", sep="\t")
    return df_ppi, df_double_comps, df_uniprot, df_uni, df_trrust, df_double_couples, df_og_tsg

def create_uniprot_gene_dict(df_uni):
    """
    Create a dictionary mapping UniProt IDs to gene symbols.

    Arguments:
    df_uni (DataFrame): DataFrame containing UniProt gene name conversion data.

    Returns:
    Uniprot_Gene (dict): Dictionary with UniProt IDs as keys and gene symbols as values.
    """
    Uniprot_Gene = {}
    for ind in df_uni.index:
        uni = df_uni.loc[ind, 'UniProt ID(supplied by UniProt)']
        if uni != np.nan:
            Uniprot_Gene[uni] = df_uni.loc[ind, 'Approved symbol']
    return Uniprot_Gene

def map_genes(df_ppi, Uniprot_Gene):
    """
    Map gene symbols to the protein-protein interaction data.

    Arguments:
    df_ppi (DataFrame): DataFrame containing PPI data.
    Uniprot_Gene (dict): Dictionary with UniProt IDs as keys and gene symbols as values.

    Returns:
    df_filtered (DataFrame): Filtered DataFrame with mapped gene symbols.
    """
    df_ppi["source_gene"] = df_ppi["source"].map(Uniprot_Gene).fillna("None")
    df_ppi["target_gene"] = df_ppi["target"].map(Uniprot_Gene).fillna("None")
    df_filtered = df_ppi[(df_ppi["source_gene"] != "None") & (df_ppi["target_gene"] != "None")]
    return df_filtered

def create_pagerank_network(df_ppi, df_uniprot):
    """
    Create a PageRank network and calculate PageRank scores.

    Arguments:
    df_ppi (DataFrame): Filtered DataFrame containing PPI data with mapped gene symbols.
    df_uniprot (DataFrame): DataFrame containing UniProt to gene symbol mapping data.

    Returns:
    G (networkx.DiGraph): Directed graph of the PPI network.
    pagerank_scores (dict): Dictionary of PageRank scores for each node in the graph.
    """
    G = nx.from_pandas_edgelist(df_ppi, source='source', target='target', create_using=nx.DiGraph())
    seed_genes = df_uniprot["Entry"].to_list()
    personalization = {gene: 1 for gene in seed_genes}
    pagerank_scores = nx.pagerank(G, alpha=0.85, personalization=personalization)
    return G, pagerank_scores

def create_subnetwork(G, pagerank_scores, threshold=0.001):
    """
    Create a subnetwork based on PageRank scores.

    Arguments:
    G (networkx.DiGraph): Directed graph of the PPI network.
    pagerank_scores (dict): Dictionary of PageRank scores for each node in the graph.
    threshold (float): PageRank score threshold for including nodes in the subnetwork.

    Returns:
    subnetwork (networkx.DiGraph): Subgraph of nodes with PageRank scores above the threshold.
    """
    subnetwork_genes = [gene for gene, score in pagerank_scores.items() if score > threshold]
    subnetwork = G.subgraph(subnetwork_genes)
    return subnetwork

def save_subnetwork(subnetwork, Uniprot_Gene, sif_filename):
    """
    Save the subnetwork in SIF format.

    Arguments:
    subnetwork (networkx.DiGraph): Subgraph of the main PPI network.
    Uniprot_Gene (dict): Dictionary with UniProt IDs as keys and gene symbols as values.
    sif_filename (str): Filename to save the subnetwork in SIF format.
    """
    with open(sif_filename, 'w') as sif_file:
        for edge in subnetwork.edges():
            source, target = edge
            source_gene = Uniprot_Gene.get(source, source)
            target_gene = Uniprot_Gene.get(target, target)
            sif_file.write(f"{source_gene} interactsWith {target_gene}\n")

def create_style_file(df_trrust, df_double_couples, subnetwork, style_filename):
    """
    Create a style file for the subnetwork.

    Arguments:
    df_trrust (DataFrame): DataFrame containing TRRUST transcription factor and target genes data.
    df_double_couples (DataFrame): DataFrame containing double mutant genes couples data.
    subnetwork (networkx.DiGraph): Subgraph of the main PPI network.
    style_filename (str): Filename to save the style file.
    """
    DoubleComponents = list(df_double_couples.itertuples(index=False, name=None))
    TF_target = [(row['Key TF'], target) for _, row in df_trrust.iterrows() for target in row['List of overlapped genes'].split(",")]
    TfList = set(df_trrust['Key TF'])
    
    with open(style_filename, 'w') as outfile:
        outfile.write("Gene1\tGene2\tEdgeDoubleMut\tEdgeTF-target\tTFGene1\tTFGene2\n")
        for edge in subnetwork.edges():
            source, target = edge
            couple1 = (source, target)
            couple2 = (target, source)
            tag1 = "DoubleMutation" if (couple1 in DoubleComponents) or (couple2 in DoubleComponents) else "NotDoubleMutation"
            tag2 = "TF-target" if (couple1 in TF_target) or (couple2 in TF_target) else "NotTF-target"
            tagGene1 = "TF" if source in TfList else "notTF"
            tagGene2 = "TF" if target in TfList else "notTF"
            outfile.write(f"{source}\t{target}\t{tag1}\t{tag2}\t{tagGene1}\t{tagGene2}\n")

def save_og_tsg_info(df_og_tsg, subnetwork, og_tsg_filename):
    """
    Save oncogene and tumor suppressor gene information for the subnetwork.

    Arguments:
    df_og_tsg (DataFrame): DataFrame containing oncogene and tumor suppressor gene list data.
    subnetwork (networkx.DiGraph): Subgraph of the main PPI network.
    og_tsg_filename (str): Filename to save the oncogene and tumor suppressor gene information.
    """
    og_tsg_dict = dict(zip(df_og_tsg['Hugo Symbol'], df_og_tsg['OG_TSG']))
    DONE = set()
    with open(og_tsg_filename, 'w') as outfile:
        for edge in subnetwork.edges():
            source, target = edge
            if source not in DONE and source in og_tsg_dict:
                outfile.write(f"{source}\t{og_tsg_dict[source]}\n")
                DONE.add(source)
            if target not in DONE and target in og_tsg_dict:
                outfile.write(f"{target}\t{og_tsg_dict[target]}\n")
                DONE.add(target)

def main():
    """
    Main function to orchestrate the workflow.
    """
    df_ppi, df_double_comps, df_uniprot, df_uni, df_trrust, df_double_couples, df_og_tsg = load_data()
    Uniprot_Gene = create_uniprot_gene_dict(df_uni)
    df_filtered = map_genes(df_ppi, Uniprot_Gene)
    G, pagerank_scores = create_pagerank_network(df_filtered, df_uniprot)
    subnetwork = create_subnetwork(G, pagerank_scores, threshold=0.001)
