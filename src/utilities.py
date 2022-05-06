import requests
import numpy as np
import pandas as pd
from config import BASE_URL, ACCESS_KEY


def __retrieve_interactions_from_biogrid__(proteins_list: list) -> dict:
    request_url = BASE_URL + "/interactions"
    data = {}

    step = 5
    for i in range(0, len(proteins_list), step):
        end = i + step
        if end >= len(proteins_list):
            end = len(proteins_list)

        # List of genes to search for
        gene_list = proteins_list[i:end]

        params = {
            "accesskey": ACCESS_KEY,
            "format": "json",  # Return results in TAB2 format
            "geneList": "|".join(gene_list),  # Must be | separated
            "searchNames": "true",  # Search against official names
            "includeInteractors": "true",
            # Set to true to get any interaction involving EITHER gene, set to false to get interactions between genes
            "includeInteractorInteractions": "false",
            # Set to true to get interactions between the geneList’s first order interactors
            "includeEvidence": "false",
            # If false "evidenceList" is evidence to exclude, if true "evidenceList" is evidence to show
            "selfInteractionsExcluded": "true",  # If true no self-interactions will be included
        }

        r = requests.get(request_url, params=params)
        interactions = r.json()

        # Check if the interactions are more than the allowed number
        if len(interactions) == 10000:
            assert False

        # Create a hash of results by interaction identifier
        for interaction_id, interaction in interactions.items():
            data[interaction_id] = interaction

    return data


def __remove_useless_interactions__(dataset: pd.DataFrame) -> pd.DataFrame:
    # Look for duplicated interactions
    duplicated_interactions = pd.DataFrame(np.sort(dataset[["InteractorA", "InteractorB"]].values, 1)).duplicated()
    print("Duplicated interactions:\n{0}".format(duplicated_interactions.value_counts()))

    # Delete such interactions from the dataset
    dataset = dataset[~duplicated_interactions.values]

    # Look for interactions where both proteins are the same
    same_proteins_interactions = pd.Series(dataset[["InteractorA", "InteractorB"]].nunique(axis=1) == 1)
    print("Useless interactions:\n{0}".format(same_proteins_interactions.value_counts()))

    # Delete such interactions from the dataset
    dataset = dataset[~same_proteins_interactions.values]
    return dataset


def retrieve_interactions(proteins_list: list, gene_interactions: pd.DataFrame) -> pd.DataFrame:
    # Load the data into a pandas dataframe
    data = __retrieve_interactions_from_biogrid__(proteins_list)
    dataset = pd.DataFrame.from_dict(data, orient="index")

    # Re-order the columns and select only the columns we want to see
    columns = ["OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B"]
    dataset = dataset[columns]

    # Rename the columns and make all the values uppercase
    dataset = dataset.rename(columns={"OFFICIAL_SYMBOL_A": "InteractorA", "OFFICIAL_SYMBOL_B": "InteractorB"})
    dataset["InteractorA"] = dataset["InteractorA"].str.upper()
    dataset["InteractorB"] = dataset["InteractorB"].str.upper()

    # Remove duplicated and self-interactions
    dataset = __remove_useless_interactions__(dataset)

    # Concatenate the found interactions with the ones involving the starting gene
    dataset = pd.concat([dataset, gene_interactions])
    return dataset


def intersection(lst1: list, lst2: list) -> list:
    inters = list()
    if not (len(lst1) == 0 or len(lst2) == 0):
        set1 = set(lst1)
        inters = [elem for elem in lst2 if elem in set1]
    return inters
