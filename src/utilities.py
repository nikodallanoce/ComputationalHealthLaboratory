import requests
import numpy as np
import gseapy as gp
import pandas as pd
import networkx as nx
from config import BASE_URL, ACCESS_KEY


def __retrieve_interactions_from_biogrid(proteins_list: list) -> dict:
    request_url = BASE_URL + "/interactions"
    data = {}

    step = 5
    for i in range(0, len(proteins_list), step):
        end = i + step
        if end >= len(proteins_list):
            end = len(proteins_list)

        # List of genes to search for
        gene_list = proteins_list[i:end]

        params = {
            "accesskey": ACCESS_KEY,
            "format": "json",  # Return results in TAB2 format
            "geneList": "|".join(gene_list),  # Must be | separated
            "searchNames": "true",  # Search against official names
            "includeInteractors": "true",
            # Set to true to get any interaction involving EITHER gene, set to false to get interactions between genes
            "includeInteractorInteractions": "false",
            # Set to true to get interactions between the geneList’s first order interactors
            "includeEvidence": "false",
            # If false "evidenceList" is evidence to exclude, if true "evidenceList" is evidence to show
            "selfInteractionsExcluded": "true",  # If true no self-interactions will be included
        }

        r = requests.get(request_url, params=params)
        interactions = r.json()

        # Check if the interactions are more than the allowed number
        if len(interactions) == 10000:
            assert False

        # Create a hash of results by interaction identifier
        for interaction_id, interaction in interactions.items():
            data[interaction_id] = interaction

    return data


def __remove_useless_interactions(dataset: pd.DataFrame) -> pd.DataFrame:
    # Look for duplicated interactions
    duplicated_interactions = pd.DataFrame(np.sort(dataset[["InteractorA", "InteractorB"]].values, 1)).duplicated()
    print("Duplicated interactions:\n{0}".format(duplicated_interactions.value_counts()))

    # Delete such interactions from the dataset
    dataset = dataset[~duplicated_interactions.values]

    # Look for interactions where both proteins are the same
    same_proteins_interactions = pd.Series(dataset[["InteractorA", "InteractorB"]].nunique(axis=1) == 1)
    print("Useless interactions:\n{0}".format(same_proteins_interactions.value_counts()))

    # Delete such interactions from the dataset
    dataset = dataset[~same_proteins_interactions.values]
    return dataset


def retrieve_interactions(proteins_list: list, gene_interactions: pd.DataFrame) -> pd.DataFrame:
    # Load the data into a pandas dataframe
    data = __retrieve_interactions_from_biogrid(proteins_list)
    dataset = pd.DataFrame.from_dict(data, orient="index")

    # Re-order the columns and select only the columns we want to see
    columns = ["OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B"]
    dataset = dataset[columns]

    # Rename the columns and make all the values uppercase
    dataset = dataset.rename(columns={"OFFICIAL_SYMBOL_A": "InteractorA", "OFFICIAL_SYMBOL_B": "InteractorB"})
    dataset["InteractorA"] = dataset["InteractorA"].str.upper()
    dataset["InteractorB"] = dataset["InteractorB"].str.upper()

    # Remove duplicated and self-interactions
    dataset = __remove_useless_interactions(dataset)

    # Concatenate the found interactions with the ones involving the starting gene
    dataset = pd.concat([dataset, gene_interactions])
    return dataset


def retrieve_diseases(nodes: list, threshold: float = 0.1) -> pd.DataFrame:
    enr = gp.enrichr(gene_list=pd.DataFrame(nodes),
                     gene_sets=['DisGeNET'],  # Datasets from the gp.get_library_name() method
                     organism='Human',
                     description='DEGs_up_1d',
                     outdir='test'
                     )

    # Keep those pathways with an adjusted p-value < threshold
    df_diseases = enr.results[enr.results["Adjusted P-value"] < threshold][
        ["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]]
    return df_diseases


def intersection(lst1: list, lst2: list) -> list:
    inters = list()
    if not (len(lst1) == 0 or len(lst2) == 0):
        set1 = set(lst1)
        inters = [elem for elem in lst2 if elem in set1]
    return inters


def build_protein_graph(nodes: list, diseases: dict, interactions: pd.DataFrame) -> nx.Graph:
    protein_graph = nx.Graph(name='Protein Interactions Graph')

    # Build the nodes
    for node in nodes:
        protein_graph.add_node(node, diseases=[])  # Each node will have a list with the disease pathways it belongs to

    # Insert into the nodes their respective diseases
    for i, disease in diseases.items():
        disease_genes = disease['genes']
        for gene in disease_genes:
            protein_graph.nodes[gene]["diseases"].append(i)

    # Build the edges
    for _, interaction in interactions.iterrows():
        first_protein, second_protein = interaction[0], interaction[1]  # Proteins involved in the interaction

        # Retrieve the proteins' diseases
        prot1_dis = protein_graph.nodes()[interaction[0]]['diseases']
        prot2_dis = protein_graph.nodes()[interaction[1]]['diseases']

        # Build the edge
        protein_graph.add_edge(first_protein, second_protein, weight=len(intersection(prot1_dis, prot2_dis)))

    return protein_graph


def nodes_no_diseases(protein_graph: nx.Graph, out: str = "size"):
    nodes_no_disease = list()
    for node in protein_graph.nodes:
        if len(protein_graph.nodes[node]["diseases"]) == 0:
            nodes_no_disease.append(str(node))

    if out == "size":
        return len(nodes_no_disease)
    elif out == "nodes":
        return nodes_no_disease
    else:
        return nodes_no_disease, len(nodes_no_disease)


def edges_no_diseases(protein_graph: nx.Graph, out: str = "size"):
    edges_no_disease = list()
    for edge in protein_graph.edges:
        if protein_graph.edges[edge]["weight"] == 0:
            edges_no_disease.append(str(edge))

    if out == "size":
        return len(edges_no_disease)
    elif out == "nodes":
        return edges_no_disease
    else:
        return edges_no_disease, len(edges_no_disease)
