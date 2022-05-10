from tqdm import tqdm
import gseapy as gp
import pandas as pd
from numpy.ma.core import mean
import networkx as nx


def retrieve_diseases(nodes: list, threshold: float = 0.1) -> pd.DataFrame:
    """
    Performs pathway enrichment on the passed nodes and keeps all the disease pathways that have an adjusted
    p-value less than the threshold
    :param nodes: list of all the nodes in the graph
    :param threshold: prune the disease pathways with an adjusted p-value bigger than this value
    :return: a dataframe with the disease pathways from DisGeNET using the gseapy package
    """
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


def build_diseases_dict(df_diseases: pd.DataFrame) -> dict:
    """
    Builds a dict of all diseases retrieved during pathway enrichment
    :param df_diseases: dataframe of all the disease pathways obtained by the pathway enrichment
    :return: a dict of the dataframe
    """
    diseases = dict()
    for i, disease in df_diseases.iterrows():
        disease_genes = disease['Genes'].split(";")
        diseases[i] = {"name": disease["Term"], "genes": disease_genes}

    return diseases


def mean_genes_diseases(diseases: dict) -> float:
    """
    Computes the mean number of genes of the diseases
    :param diseases: dict of the diseases coming from the build_diseases_dict method
    :return: the mean number of genes for disease
    """
    mean_size = 0
    for _, disease in diseases.items():
        mean_size += len(disease['genes'])

    mean_size /= len(diseases.keys())
    return mean_size


def largest_conn_comp(protein_graph: nx.Graph, diseases_dict: dict) -> list:
    """
    Fraction of disease proteins that lie in the disease’s largest pathway component (i.e., the relative size of the
    largest connected component (LCC) of the disease)
    :param protein_graph: networkx graph of the protein-to-protein network
    :param diseases_dict: dict of the diseases coming from the build_diseases_dict method
    :return: the largest connected component metric for each disease
    """
    lcc_score = list()
    for _, disease_dict in tqdm(diseases_dict.items()):
        sub_graph = protein_graph.subgraph(disease_dict['genes'])  # Subgraph of the current disease
        largest_cc = max(nx.connected_components(sub_graph), key=len)
        lcc_score.append(len(largest_cc) / len(sub_graph.nodes()))

    return lcc_score


def distance_pathway_comps(protein_graph: nx.Graph, diseases_dict: dict) -> list:
    """
    For each pair of pathway components, we calculate the average shortest path length between each set of proteins,
    and then, the average of this is taken over all pairs of the components
    :param protein_graph: networkx graph of the protein-to-protein network
    :param diseases_dict: dict of the diseases coming from the build_diseases_dict method
    :return: the distance of pathway components metric for each disease
    """
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
