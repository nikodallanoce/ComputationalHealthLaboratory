import networkx as nx
import numpy as np
import pandas as pd
from utilities import intersection


def create_communities(protein_graph: nx.Graph, threshold: int = 1) -> tuple:
    louvain_communities = list(nx.algorithms.community.louvain_communities(protein_graph))
    communities = [community for community in louvain_communities if len(community) > threshold]
    return louvain_communities, communities


def mean_size_communities(communities: list) -> float:
    mean_size = 0
    for community in communities:
        mean_size += len(community)

    mean_size /= len(communities)
    return mean_size


def diseases_in_community(protein_graph: nx.Graph, community: set) -> set:
    diseases_community = set()
    for protein in list(community):
        diseases_protein = protein_graph.nodes[protein]["diseases"]
        diseases_community.update(diseases_protein)

    return diseases_community


def mean_diseases_communities_size_n(communities: list, protein_graph: nx.Graph, n: int = 1) -> float:
    mean_diseases = 0
    n_size_commmunities = 0
    for community in communities:
        community = list(community)
        if len(community) == n:
            protein = community[0]
            protein_diseases = protein_graph.nodes[protein]['diseases']
            n_size_commmunities += 1
            mean_diseases += len(protein_diseases)

    if n_size_commmunities == 0:
        print("There are no communities with {0} nodes".format(n))
    else:
        mean_diseases /= n_size_commmunities

    return mean_diseases


def are_communities_distinct(communities: list) -> bool:
    for i, first_community in enumerate(communities):
        for j in range(i+1, len(communities)):
            second_community = communities[j]
            if len(intersection(first_community, second_community)) > 0:
                return False

    return True


def look_for_gene_community(protein: str, communities: list) -> int:
    for i, community in enumerate(communities):
        if protein in community:
            return i

    return -1


def communities_metrics(communities: list, diseases: dict) -> pd.DataFrame:
    df_ranks = list()  # pd.DataFrame(columns=["com", "disease", "rel_val", "common_genes"])
    for i, community in enumerate(communities):
        tot_genes = dict()
        shared_genes = dict()
        for k, disease in diseases.items():
            genes = disease['genes']
            shared_genes_community = intersection(genes, community)
            tot_genes[k] = len(genes)
            shared_genes[k] = len(shared_genes_community)

        for j in range(len(tot_genes)):
            n_genes, n_shared_genes = tot_genes[j], shared_genes[j]
            if n_shared_genes > 1:
                df_ranks.append({"Community": i, "Disease": diseases[j]['name'], "Shared genes": n_shared_genes,
                                 "Disease genes": n_genes, "Community size": len(community)})

    df_ranks = pd.DataFrame(df_ranks)

    # Ratio of the shared genes (between community and disease pathway) and the number of genes in the disease
    df_ranks["Ratio disease"] = df_ranks['Shared genes'] / df_ranks['Disease genes']

    # Ratio of the shared genes (between community and disease pathway) and the size of the community
    df_ranks['Ratio community'] = df_ranks["Shared genes"] / df_ranks["Community size"]

    # Relevance of the  based on the previous computed metrics
    df_ranks["Relevance"] = df_ranks["Ratio disease"] * df_ranks["Ratio community"]
    return df_ranks


def communities_distance(communities_rank: pd.DataFrame,
                         first_index: int,
                         second_index: int,
                         metric: str = "Relevance") -> float:
    # Retrieve all the diseases linked to the communities
    tmp_first = communities_rank[communities_rank["Community"] == first_index][["Disease", metric]]
    tmp_second = communities_rank[communities_rank["Community"] == second_index][["Disease", metric]]

    # Shared diseases between the two communities
    shared_diseases = tmp_first.merge(tmp_second, left_on="Disease", right_on="Disease")

    # Diseases that are not shared between the two communities
    tmp_first = tmp_first[~tmp_first['Disease'].isin(shared_diseases["Disease"])]
    tmp_second = tmp_second[~tmp_second['Disease'].isin(shared_diseases["Disease"])]

    distance = 0
    for _, disease in shared_diseases.iterrows():
        distance += np.power(disease[metric + "_x"] - disease[metric + "_y"], 2)

    for _, first_disease in tmp_first.iterrows():
        distance += np.power(first_disease[metric], 2)

    for _, second_disease in tmp_second.iterrows():
        distance += np.power(second_disease[metric], 2)

    return distance
