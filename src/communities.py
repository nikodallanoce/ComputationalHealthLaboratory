import networkx as nx
from utilities import intersection


def create_communities(protein_graph: nx.Graph, threshold: int = 1) -> list:
    louvain_communities = list(nx.algorithms.community.louvain_communities(protein_graph))
    communities = [community for community in louvain_communities if len(community) > threshold]
    return communities


def mean_size_communities(communities: list, verbose: bool = False) -> float:
    mean_size = 0
    for community in communities:
        mean_size += len(community)

    mean_size /= len(communities)
    if verbose:
        print("Mean size of communities: {0}".format(str(mean_size)))

    return mean_size


def mean_diseases_communities(communities: list, protein_graph: nx.Graph, size: int = 1) -> float:
    mean_diseases = 0
    n_size_commmunities = 0
    for community in communities:
        community = list(community)
        if len(community) == size:
            protein = community[0]
            protein_diseases = protein_graph.nodes[protein]['diseases']
            n_size_commmunities += 1
            mean_diseases += len(protein_diseases)

    if n_size_commmunities == 0:
        print("There are no communities with only one node")
    else:
        mean_diseases /= n_size_commmunities
        print("Mean diseases for those communities with one node: {0}".format(str(mean_diseases)))

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
