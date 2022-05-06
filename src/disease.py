from tqdm import tqdm
import gseapy as gp
import pandas as pd
from numpy.ma.core import mean
import networkx as nx


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


def mean_genes_diseases(diseases: dict, verbose: bool = False) -> float:
    mean_size = 0
    for _, disease in diseases.items():
        mean_size += len(disease['genes'])

    mean_size /= len(diseases.keys())
    if verbose:
        print("Mean number of genes for disease: {0}".format(mean_size))

    return mean_size


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
