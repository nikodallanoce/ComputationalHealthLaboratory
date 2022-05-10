import networkx as nx
import numpy as np
import pandas as pd
from utilities import intersection
from disease import build_diseases_dict


def create_communities(protein_graph: nx.Graph, threshold: int = 1) -> tuple:
    """
    Create the communities of the protein-to-protein graph and prunes those with a size less than a threshold
    :param protein_graph: networkx graph of the protein-to-protein network
    :param threshold: prune those communities with a size less than the threshold
    :return: a list of all communities and a list of the communities after the pruning
    """
    louvain_communities = list(nx.algorithms.community.louvain_communities(protein_graph))

    # Prune the communities with size < threshold
    communities = [community for community in louvain_communities if len(community) > threshold]
    return louvain_communities, communities


def mean_size_communities(communities: list) -> float:
    """
    Computes the mean size of the communities
    :param communities: a list of all the communities
    :return: the mean size of the passed communities
    """
    mean_size = 0
    for community in communities:
        mean_size += len(community)

    mean_size /= len(communities)
    return mean_size


def diseases_in_community(protein_graph: nx.Graph, community: set) -> set:
    """
    Lists all the diseases that are found in a community, a disease is inside the community if at least one node
    is colored with its index
    :param protein_graph: networkx graph of the protein-to-protein network
    :param community: set of the nodes that belong to a community
    :return: a set of the diseases in a community
    """
    diseases_community = set()
    for protein in list(community):
        diseases_protein = protein_graph.nodes[protein]["diseases"]
        diseases_community.update(diseases_protein)

    return diseases_community


def mean_diseases_communities_size_n(communities: list, protein_graph: nx.Graph, n: int = 1) -> float:
    """
    Computes the mean number of diseases for those communities with a certain size
    :param communities: list of all communities
    :param protein_graph: networkx graph of the protein-to-protein network
    :param n: the size of the community to consider
    :return: the mean number of diseases for the communities
    """
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
    """
    By looking at the intersection between nodes for each pair of communities in order to tell if they are distinct
    which means that such intersection contains no nodes for all the pairs
    :param communities: list of all communities
    :return: a boolean indicating if the communities are distinct or not
    """
    for i, first_community in enumerate(communities):
        for j in range(i+1, len(communities)):
            second_community = communities[j]
            if len(intersection(first_community, second_community)) > 0:
                return False

    return True


def look_for_gene_community(protein: str, communities: list) -> int:
    """
    Look for the index of the protein's community, if it exists
    :param protein: a string that represents the protein's term
    :param communities: a list of all communities
    :return: the index of the community the protein belongs to, if it exists
    """
    for i, community in enumerate(communities):
        if protein in community:
            return i

    return -1


def communities_metrics(communities: list, diseases: dict) -> pd.DataFrame:
    """
    Computes three metrics:
    - Ratio of the shared genes (between community and disease pathway) and the number of genes in the disease
    - Ratio of the shared genes (between community and disease pathway) and the size of the community
    - Relevance of the disease in the community, based on the previous computed metrics
    And puts them inside a community-disease pandas dataframe
    :param communities: a list of all communities
    :param diseases: a dict of all the disease retrieved by the gseapy package
    :return: a pandas community-disease dataframe with the computed metrics as attributes
    """
    df_ranks = list()
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

    # Relevance of the disease in the community, based on the previous computed metrics
    df_ranks["Relevance"] = df_ranks["Ratio disease"] * df_ranks["Ratio community"]
    return df_ranks


def communities_distance(communities_rank: pd.DataFrame,
                         first_index: int,
                         second_index: int,
                         metric: str = "Relevance") -> float:
    """
    Computes the distance between two communities given their index and metric, you should pass the dataframe
    coming the communities_metrics method to let this one work properly
    :param communities_rank: a community-disease dataframe with some relevant disease metrics
    :param first_index: index of the first community
    :param second_index: index of the second commmunity
    :param metric: which metric to base the computation on ("Ratio disease", "Ratio community" or "Relevance")
    :return: the euclidean distance of the two communities based on a metric
    """
    # Retrieve all the diseases linked to the communities
    tmp_first = communities_rank[communities_rank["Community"] == first_index][["Disease", metric]]
    tmp_second = communities_rank[communities_rank["Community"] == second_index][["Disease", metric]]

    # Shared diseases between the two communities
    shared_diseases = tmp_first.merge(tmp_second, left_on="Disease", right_on="Disease")

    # Diseases that are not shared between the two communities
    tmp_first = tmp_first[~tmp_first['Disease'].isin(shared_diseases["Disease"])]
    tmp_second = tmp_second[~tmp_second['Disease'].isin(shared_diseases["Disease"])]

    distance = 0

    # Compute the euclidean distance for the shared diseases
    for _, disease in shared_diseases.iterrows():
        distance += np.power(disease[metric + "_x"] - disease[metric + "_y"], 2)

    # Compute the distance for the diseases that belong only to the first community
    for _, first_disease in tmp_first.iterrows():
        distance += np.power(first_disease[metric], 2)

    # Compute the distance for those diseases only in the second community
    for _, second_disease in tmp_second.iterrows():
        distance += np.power(second_disease[metric], 2)

    return distance


def retrieve_gene_community_diseases(communities: list,
                                     protein: str,
                                     diseases: dict,
                                     protein_graph: nx.Graph,
                                     df_diseases: pd.DataFrame,
                                     communities_rank: pd.DataFrame = None,
                                     threshold: int = 10) -> tuple:
    """
    Retrieves two dataframes (in both case the disease shares some genes in the protein community):
    - one with the diseases that do not have the passed protein into their pathways
    - the other with the disease that have the passed protein into their pathways
    :param communities: a list of all communities
    :param protein: a string that represents the protein's term
    :param diseases: a dict of all diseases
    :param protein_graph: networkx graph of the protein-to-protein network
    :param df_diseases: dataframe of the diseases pathways retrieved by using the gseapy package
    :param communities_rank: the community-disease dataframe coming from the community_metrics method
    :param threshold: the number of genes that the community must at least share with a disease pathway
    :return: two dataframes, one containing the diseases inside the community that do not have the
    passed protein and the other with all the diseases that have the passed protein
    """
    gene_community = look_for_gene_community(protein, communities)
    if gene_community == -1:
        raise Exception("The protein {0} is not in one of the communities".format(protein))

    if diseases is None:
        diseases = build_diseases_dict(df_diseases)

    if communities_rank is None:
        communities_rank = communities_metrics(communities, diseases)

    # Retrieve all the diseases in the gene community
    disease_rank = communities_rank[(communities_rank["Community"] == gene_community) &
                                    (communities_rank["Shared genes"] > threshold)].\
        sort_values(by="Relevance", ascending=False).\
        drop(["Community", "Shared genes", "Disease genes", "Community size"], axis=1)

    # Look for all the diseases pathways the gene belongs to
    gene_index_diseases = protein_graph.nodes()[protein]['diseases']
    gene_diseases = df_diseases.loc[gene_index_diseases, 'Term']

    # Retrieve those diseases in the same community of the gene, but the latter is not in their pathway
    diseases_no_gene = disease_rank[~disease_rank['Disease'].isin(gene_diseases)].sort_values(by='Relevance',
                                                                                              ascending=False)

    # Retrieve those diseases in the same community of the gene, in which the latter appears
    diseases_gene = disease_rank[disease_rank['Disease'].isin(gene_diseases)].sort_values(by="Relevance",
                                                                                          ascending=False)
    return diseases_no_gene, diseases_gene
