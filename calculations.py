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


def download_homologene_data():
    """
    Downloads the latest HomoloGene data from NCBI and stores it in a Pandas DataFrame.

    Returns
    -------
    pandas.DataFrame :
        A DataFrame containing the HomoloGene data. The columns are: 
        'HID' (HomoloGene group ID), 'Taxonomy ID', 'Gene ID', 'Gene Symbol', and 'Protein ID'.
    """

    url = 'ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data'
    filename = 'homologene.data'

    # Download DB
    urllib.request.urlretrieve(url, filename)

    # Save it into df
    col_names = ['HID', 'Taxonomy ID', 'Gene ID', 'Gene Symbol', 'Protein ID']
    df = pd.read_csv(filename, sep='\t', header=None, names=col_names)

    os.remove(filename)

    return df
