from tqdm import tqdm
from numpy.ma.core import mean
import networkx as nx


def largest_conn_comp(protein_graph: nx.Graph, diseases_dict: dict) -> list:
    lcc_score = list()
    for _, disease_dict in tqdm(diseases_dict.items()):
        sub_graph = protein_graph.subgraph(disease_dict['genes'])  # Subgraph of the current disease
        largest_cc = max(nx.connected_components(sub_graph), key=len)
        lcc_score.append(len(largest_cc) / len(sub_graph.nodes()))

    return lcc_score


def distance_pathway_comps(protein_graph: nx.Graph, diseases_dict: dict) -> list:
    dpc_score = list()
    for _, disease_dict in tqdm(diseases_dict.items()):
        sub_graph = protein_graph.subgraph(disease_dict['genes'])
        conn_comps = list(nx.connected_components(sub_graph))
        distances = list()
        for i, comp in enumerate(conn_comps):
            for j in range(i + 1, len(conn_comps)):
                dist = 0
                for first_comp_protein in comp:
                    for second_comp_protein in conn_comps[j]:
                        dist += nx.shortest_path_length(protein_graph, source=first_comp_protein,
                                                        target=second_comp_protein)

                distances.append(dist / (len(comp) * len(conn_comps[j])))

        dpc_score.append(mean(distances))

    return dpc_score
