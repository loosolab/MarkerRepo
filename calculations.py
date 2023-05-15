import numpy as np
import marker_repo as mr
import pandas as pd
import yaml
import urllib.request
import os

def compare_marker_lists(REPO_LISTS_PATH, keywords, case_sensitive=False, exact=False):
    """
    This function compares and scores markers from selected marker lists using Inverse Document Frequency (IDF).

    Parameters
    ----------
    keywords : str, dict
        Keywords for selecting marker lists. If a string, the function will check if the string is contained anywhere in the marker lists. If a dictionary, the keys are the column names and the values are the keywords to search for in those columns.
    REPO_LISTS_PATH : str
        The path to the directory containing the marker lists.
    case_sensitive : bool, default: False
        If True, the search will be case-sensitive. If False, the search will be case-insensitive.
    exact : bool, default: False
        If True, the search will look for exact matches. If False, the search will look for substrings.

    Returns
    --------
    pd.DataFrame :
        A DataFrame containing the name ("Info"), the marker ("Marker") and the score ("Score") for each marker in the selected marker lists.
    """

    df = mr.search_db(mr.get_db(REPO_LISTS_PATH), keywords, case_sensitive=case_sensitive, exact=exact)
    uids = [int(idx) for idx in df.index]
    files = mr.get_uid_paths(REPO_LISTS_PATH, uids)

    # Initialize dictionary
    marker_dict = {"Info": [], "Marker": []}

    # Load marker lists from selected files
    for file in files:
        with open(file, 'r') as f:
            data = yaml.safe_load(f)
            marker_list_section = data.get('marker_list', [])
            
            for marker_list in marker_list_section:
                name = marker_list.get('name', '')
                markers = marker_list.get('markers', [])
                marker_dict["Info"].extend([name]*len(markers))
                marker_dict["Marker"].extend(markers)

    df = pd.DataFrame(marker_dict)
    
    # Calculate scores
    total_lists = len(df['Info'].unique())
    marker_counts = df['Marker'].value_counts()
    df['Score'] = df['Marker'].apply(lambda x: np.log(total_lists / marker_counts[x]))

    # Sort by Score in descending order
    df.sort_values('Score', ascending=False, inplace=True)

    return df


def invert_score(df):
    """
    Inverts and scales the score values in the DataFrame (marker list).
    
    The function calculates the inverse of the scores in the DataFrame, 
    in addition the highest score becomes 0 and the lowest score becomes 1. 
    All other scores are scaled accordingly. Thus, higher values correspond to lower uniqueness.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing the Info, Marker and Score columns of the marker list.
        
    Returns
    --------
    pandas.DataFrame :
        The DataFrame containing the modified scores.
    """

    # Get min and max score
    min_score = df['Score'].min()
    max_score = df['Score'].max()

    # Invert and scale scores
    df['Score'] = (max_score - df['Score']) / (max_score - min_score)

    return df


def download_homologene_data(url="ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"):
    """
    Download and parse the HomoloGene data.

    Parameters
    ----------
    url : str
        URL to the HomoloGene data file.

    Returns
    --------
    pandas.DataFrame : 
        DataFrame with the HomoloGene data.
    """

    # Download the HomoloGene data
    homologene_data = pd.read_csv(url, sep='\t', header=None, index_col=0)

    # Rename the columns
    homologene_data.columns = ["Taxonomy ID", "Gene ID", "Gene Symbol", "Protein gi", "Protein accession"]
    homologene_data.index.names = ["HID"]

    return homologene_data


def transfer_markers(df, source_organism, target_organism):
    """
    Transfer markers between organisms based on homology.
    
    Parameters
    ----------
    df : DataFrame
        Input DataFrame containing 'Marker' and 'Info' columns.
    source_organism : int
        Taxonomy ID of the source organism.
    target_organism : int
        Taxonomy ID of the target organism.
        
    Returns
    --------
    DataFrame:
        Transferred markers DataFrame with columns corresponding to source_organism and target_organism.
    """
    
    pass
