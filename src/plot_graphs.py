import networkx as nx
from pyvis.network import Network
import numpy as np
from communities import look_for_gene_community


def __build_network_disease__(protein_graph: nx.Graph, disease_genes: list = None) -> Network:
    net = Network(width=1080, height=720)
    node_index = dict()
    for i, node in enumerate(protein_graph.nodes()):
        node_index[node] = i
        if node not in disease_genes:
            net.add_node(i, label=node, size=8)
        else:
            net.add_node(i, label=node, size=16, color="red")

    for edge_from, edge_to in protein_graph.edges():
        if edge_from in disease_genes or edge_to in disease_genes:
            net.add_edge(node_index[edge_from], node_index[edge_to], color="red", value=1)
        else:
            net.add_edge(node_index[edge_from], node_index[edge_to])

    return net


def __build_network_protein__(protein_graph: nx.Graph, protein: str = None) -> Network:
    net = Network(width=1080, height=720)
    node_index = dict()
    for i, node in enumerate(protein_graph.nodes()):
        node_index[node] = i
        if node != protein:
            net.add_node(i, label=node, size=8)
        else:
            net.add_node(i, label=node, size=16, color="red")

    for edge_from, edge_to in protein_graph.edges():
        if edge_from == protein or edge_to == protein:
            net.add_edge(node_index[edge_from], node_index[edge_to], color="red", value=1)
        else:
            net.add_edge(node_index[edge_from], node_index[edge_to])

    return net


def plot_protein_network(protein_graph: nx.Graph, disease_genes: list = None, biomarkers: list = None) -> None:
    if biomarkers is not None:
        plot_graph = protein_graph.subgraph(biomarkers)
    else:
        plot_graph = protein_graph.copy()

    net = __build_network_disease__(plot_graph, disease_genes)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("protein_graph.html")


def plot_disease_network(protein_graph: nx.Graph, disease_genes: list, protein: str = None) -> None:
    sub_graph = protein_graph.subgraph(disease_genes)
    net = __build_network_protein__(sub_graph, protein)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("disease_graph.html")


def plot_community_protein(protein_graph: nx.Graph, communities: list, protein: str = None) -> None:
    if protein is not None:
        community = communities[look_for_gene_community(protein, communities)]
    else:
        community = np.random.randint(0, len(communities))

    sub_graph = protein_graph.subgraph(community)
    net = __build_network_protein__(sub_graph, protein)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("community_protein_graph.html")


def plot_community_disease(protein_graph: nx.Graph, disease_genes: list, community: set) -> None:
    sub_graph = protein_graph.subgraph(community)
    net = __build_network_disease__(sub_graph, disease_genes)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("community_disease_graph.html")
