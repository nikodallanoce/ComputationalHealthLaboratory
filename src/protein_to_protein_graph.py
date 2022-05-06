import pandas as pd
import networkx as nx
import networkx.algorithms.centrality as nx_c
from utilities import intersection


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
        first_prot_dis = protein_graph.nodes()[interaction[0]]['diseases']
        second_prot_dis = protein_graph.nodes()[interaction[1]]['diseases']

        # Build the edge
        protein_graph.add_edge(first_protein, second_protein, weight=len(intersection(first_prot_dis, second_prot_dis)))

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


def identify_biomarkers(protein_graph: nx.Graph, n: int = 10, protein: str = None) -> pd.DataFrame:
    nodes_degree = pd.DataFrame.from_dict(nx_c.degree_centrality(protein_graph), orient="index", columns=["centrality"])
    nodes_degree = nodes_degree.sort_values(by='centrality', ascending=False)
    biomarkers = nodes_degree.iloc[:n]

    # Insert the starting gene
    if protein is not None:
        biomarkers = pd.concat([biomarkers, nodes_degree[nodes_degree.index == protein]])

    return biomarkers
