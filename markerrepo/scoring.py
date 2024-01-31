import numpy as np
from .homology import get_panglao_ui
from .homology import transfer_ui_to_homologs
from .marker_repo import combine_dfs, combine_lists, search_df
from sklearn.preprocessing import MinMaxScaler


def compare_marker_lists(repo_path=".", lists_path=None, keywords=None, marker_df=None, case_sensitive=False, exact=False):
    """
    This function compares and scores markers from selected marker lists using Ubiquitousness Index.
    A score of '0' signifies that the marker is the most specific within this selection, 
    while a score of '1' indicates that the marker is the most prevalent.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    keywords : str, dict
        Keywords for selecting marker lists. If a string, the function will check if the string is contained anywhere in the marker lists. 
        If a dictionary, the keys are the column names and the values are the keywords to search for in those columns.
    marker_df : pd.DataFrame
        Input DataFrame with columns "Info", "Marker"
    case_sensitive : bool, default: False
        If True, the search will be case-sensitive. If False, the search will be case-insensitive.
    exact : bool, default: False
        If True, the search will look for exact matches. If False, the search will look for substrings.

    Returns
    --------
    pd.DataFrame :
        A DataFrame containing the name ("Info"), the marker ("Marker") and the score ("Score") for each marker in the selected marker lists.
    """

    if marker_df is not None:
        df = marker_df
    else:
        df = search_df(combine_dfs(repo_path=repo_path, lists_path=lists_path), keywords, case_sensitive=case_sensitive, exact=exact)
        uids = [int(idx) for idx in df.index]
        df = combine_lists(uids, repo_path=repo_path, lists_path=lists_path)

    df = df.drop_duplicates()

    # Calculate scores
    total_lists = len(df['Info'].unique())
    marker_counts = df['Marker'].value_counts()
    df['Score'] = df['Marker'].apply(lambda x: marker_counts[x] / total_lists)

    # Scaling scores to be between 0 and 1
    scaler = MinMaxScaler()
    df['Score'] = scaler.fit_transform(df[['Score']])
    df.sort_values('Score', ascending=True, inplace=True)

    return df


def update_scores(df, organism="human", repo_path=".", biomart_target=None):
    """
    Update the scores in the dataframe using the ubiquitousness index from the panglao database.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns "Marker", "Info", and "Score". 
    organism : str, default human
        Organism to consider when retrieving the ubiquitousness index.
    repo_path : str, default "."
        The path of the Marker Repo.
    biomart_target : str, default None
        Name of the target organism in the BioMart database. If None, the name is inferred from the organism parameter.

    Returns
    -------
    pd.DataFrame :
        Updated dataframe with new scores.
    """

    # Retrieve the ubiquitousness index dictionary for the given organism
    panglao_organisms = ["human", "mouse"]

    if organism in panglao_organisms:
        ui_dict = get_panglao_ui(repo_path=repo_path, organism=organism)
    else:
        print(f"Transferring Panglao ubiquitousness index to {organism}...")
        if biomart_target is None:
            ui_dict = transfer_ui_to_homologs(target_organism=organism, repo_path=repo_path)
        else:
            ui_dict = transfer_ui_to_homologs(biomart_target=biomart_target, repo_path=repo_path, target_organism=organism)

    # Split the "Marker" column and take the first part
    df['MainMarker'] = df['Marker'].str.split().str[0].str.upper()

    # Check if 'Score' column is in df
    if 'Score' not in df.columns:
        df['Score'] = np.nan

    # Update scores where the main marker is in ui_dict
    df['Score'] = df['MainMarker'].map(ui_dict).fillna(df['Score'])


    df = df.drop(columns='MainMarker')
    df = df.dropna()
    df.sort_values('Score', ascending=True, inplace=True)

    return df