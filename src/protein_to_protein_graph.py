import pandas as pd
import networkx as nx
import networkx.algorithms.centrality as nx_c
from utilities import intersection


def build_protein_graph(nodes: list, diseases: dict, interactions: pd.DataFrame) -> nx.Graph:
    """
    Build the protein-to-protein graph by coloring rach node with the disease pathways it belongs to and weight
    the edges based on the shared diseases between the nodes
    :param nodes: list of nodes which were retrieved from the BioGRID dataset
    :param diseases: dict of diseases
    :param interactions: dataframe of interactions between nodes
    :return: networkx graph of the protein-to-protein network
    """
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
        first_prot_dis = protein_graph.nodes()[interaction[0]]['diseases']
        second_prot_dis = protein_graph.nodes()[interaction[1]]['diseases']

        # Build the edge
        protein_graph.add_edge(first_protein, second_protein, weight=len(intersection(first_prot_dis, second_prot_dis)))

    return protein_graph


def nodes_no_diseases(protein_graph: nx.Graph, out: str = "size"):
    """
    Nodes with node color/disease pathway, the method can also return their number or both values
    :param protein_graph: networkx graph of the protein-to-protein graph
    :param out: string that tells which output is requested
    :return: the number of nodes with no diseases, the nodes or both values depending on "out"
    """
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
    """
    Edges with no weight, the method can also return their number or both values
    :param protein_graph: networkx graph of the protein-to-protein graph
    :param out: string that tells which output is requested
    :return: the number of edges with no weights, the edges or both values depending on "out"
    """
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


def identify_biomarkers(protein_graph: nx.Graph, n: int = 10, protein: str = None) -> pd.DataFrame:
    """
    Identifies the biomarkers of the protein-to-protein network, which are the nodes with the highest degree
    :param protein_graph: networkx graph of the protein-to-protein graph
    :param n: number of biomarkers to consider
    :param protein: protein that needs to be considered a biomakers no matter its centrality
    :return: dataframe with the n nodes with the highest centrality
    """
    # Compute the centrality of each node based on their degree
    nodes_degree = pd.DataFrame.from_dict(nx_c.degree_centrality(protein_graph), orient="index", columns=["centrality"])
    nodes_degree = nodes_degree.sort_values(by='centrality', ascending=False)

    # Choose the first n nodes
    biomarkers = nodes_degree.iloc[:n]

    # Insert the starting gene
    if protein is not None:
        biomarkers = pd.concat([biomarkers, nodes_degree[nodes_degree.index == protein]])

    return biomarkers
