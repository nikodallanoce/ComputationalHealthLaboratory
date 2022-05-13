import networkx as nx
import seaborn as sn
from pyvis.network import Network
import numpy as np
from communities import look_for_gene_community


def __build_network_disease__(protein_graph: nx.Graph,
                              genes_disease: list = None,
                              protein: str = None,
                              plot_other_edges: bool = True) -> Network:
    net = Network(width=1080, height=720)
    node_index = dict()
    for i, node in enumerate(protein_graph.nodes()):
        shape = "dot"
        size = 8
        color = "blue"
        node_index[node] = i
        if genes_disease is not None and node in genes_disease:
            color = "orange"
            size = 16
            shape = "diamond"

        if node == protein:
            color = "yellowgreen"
            size = 16

        net.add_node(i, label=node, size=size, color=color, shape=shape)

    for edge_from, edge_to in protein_graph.edges():
        if genes_disease is not None and ((edge_from == protein and edge_to in genes_disease) or
                                          (edge_from in genes_disease and edge_to == protein)):
            net.add_edge(node_index[edge_from], node_index[edge_to], color="green", value=1)
        elif genes_disease is not None and edge_from in genes_disease and edge_to in genes_disease:
            net.add_edge(node_index[edge_from], node_index[edge_to], color="orangered", value=1)
        elif plot_other_edges:
            net.add_edge(node_index[edge_from], node_index[edge_to], color="blue")

    return net


def __build_network_protein__(protein_graph: nx.Graph,
                              protein: str = None,
                              plot_other_edges: bool = True) -> Network:
    net = Network(width=1080, height=720)
    node_index = dict()
    for i, node in enumerate(protein_graph.nodes()):
        node_index[node] = i
        if node != protein:
            net.add_node(i, label=node, size=8)
        else:
            net.add_node(i, label=node, size=16, color="yellowgreen")

    for edge_from, edge_to in protein_graph.edges():
        if edge_from == protein or edge_to == protein:
            net.add_edge(node_index[edge_from], node_index[edge_to], color="green", value=1)
        elif plot_other_edges:
            net.add_edge(node_index[edge_from], node_index[edge_to])

    return net


def plot_protein_network(protein_graph: nx.Graph,
                         genes_disease: list = None,
                         biomarkers: list = None,
                         protein: str = None,
                         plot_other_edges: bool = True) -> None:
    if biomarkers is not None:
        plot_graph = protein_graph.subgraph(biomarkers)
    else:
        plot_graph = protein_graph.copy()

    net = __build_network_disease__(plot_graph, genes_disease, protein, plot_other_edges)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("protein_graph.html")


def plot_disease_network(protein_graph: nx.Graph, genes_disease: list, protein: str = None) -> None:
    sub_graph = protein_graph.subgraph(genes_disease)
    net = __build_network_protein__(sub_graph, protein)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("disease_graph.html")


def plot_community_protein(protein_graph: nx.Graph,
                           communities: list,
                           protein: str = None,
                           plot_other_edges: bool = True) -> None:
    if protein is not None:
        community = communities[look_for_gene_community(protein, communities)]
    else:
        community = np.random.randint(0, len(communities))

    sub_graph = protein_graph.subgraph(community)
    net = __build_network_protein__(sub_graph, protein, plot_other_edges)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("community_protein_graph.html")


def plot_community_disease(protein_graph: nx.Graph,
                           genes_disease: list,
                           community: set,
                           protein: str = None,
                           plot_other_edges: bool = True) -> None:
    sub_graph = protein_graph.subgraph(community)
    net = __build_network_disease__(sub_graph, genes_disease, protein, plot_other_edges)
    net.toggle_drag_nodes(False)
    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("community_disease_graph.html")


def plot_communities(protein_graph: nx.Graph, communities: list, protein: str = None) -> None:
    palette = sn.color_palette("tab10", len(communities)).as_hex()
    net = Network(width=1080, height=720)
    node_index = dict()
    index = 0
    gene_community = look_for_gene_community(protein, communities)
    if gene_community == -1:
        raise Warning("Protein {0} not found, no community will be highlighted with gold".format(protein))

    for i in range(len(communities)):
        community = communities[i]
        color = palette[i]
        sub_graph = protein_graph.subgraph(community)
        for _, node in enumerate(sub_graph.nodes()):
            node_index[node] = index
            if node == protein:
                net.add_node(index, label=node, size=32)
            else:
                net.add_node(index, label=node, size=8)

            index += 1

        for edge_from, edge_to in sub_graph.edges():
            net.add_edge(node_index[edge_from], node_index[edge_to], color=color)

    net.show_buttons(['physics'])
    net.force_atlas_2based(spring_strength=0.02)
    net.show("communities.html")
