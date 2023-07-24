import numpy as np
import pandas as pd
from annotate_by_marker_and_features.markerrepo.homology import transfer_ui_to_homologs
from annotate_by_marker_and_features.markerrepo.marker_repo import combine_dfs, combine_lists, get_gene_dict, search_df
from sklearn.preprocessing import MinMaxScaler


def compare_marker_lists(repo_path=".", keywords=None, marker_df=None, case_sensitive=False, exact=False):
    """
    This function compares and scores markers from selected marker lists using Ubiquitousness Index.
    A score of '0' signifies that the marker is the most specific within this selection, 
    while a score of '1' indicates that the marker is the most prevalent.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
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
        df = search_df(combine_dfs(repo_path=repo_path), keywords, case_sensitive=case_sensitive, exact=exact)
        uids = [int(idx) for idx in df.index]
        df = combine_lists(uids, repo_path=repo_path)

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


def get_panglao_ui(panglao_file="panglao_markers", organism="human", id_type='symbol', repo_path="."):
    """
    Create a dictionary with gene symbols or Ensembl IDs (based on id_type) and average ubiquitousness index as values.
    Only considers rows with the given organism.

    Parameters
    ----------
    panglao_file : str
        Path to the panglao markers.
    organism : str, default human
        Organism to consider (e.g. "human" or "mouse").
    id_type : str, default 'symbol'
        Type of gene identifier to use as dictionary keys. 'symbol' for gene symbols, 'ensembl' for Ensembl IDs.
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    -------
    dict :
        Dictionary with gene identifiers and average ubiquitousness index as values.
    """
    df = pd.read_csv(f"{repo_path}/{panglao_file}", sep="\t")

    # Convert full organism name to short code for filtering dataframe
    organism_dict = {'human': 'Hs', 'mouse': 'Mm'}
    if organism not in organism_dict:
        raise ValueError(f'Invalid organism {organism}. Expected "human" or "mouse".')
    species_code = organism_dict[organism]

    # Filter dataframe by organism
    df = df[df['species'].str.contains(species_code, na=False)]

    panglao_ui_dict = {}
    gene_dict = get_gene_dict(organism=organism, repo_path=repo_path)

    for _, row in df.iterrows():
        # Get gene symbols and nicknames
        symbols = [row['official gene symbol']]
        if pd.notna(row['nicknames']):
            symbols.extend(row['nicknames'].split('|'))

        # Convert symbols to uppercase and to Ensembl IDs if id_type is 'ensembl'
        symbols = [symbol.upper() for symbol in symbols]
        if id_type == 'ensembl':
            symbols = [gene_dict.get(symbol) for symbol in symbols]
            # Remove None values (symbols that couldn't be converted to Ensembl IDs)
            symbols = [symbol for symbol in symbols if symbol is not None]


        # Fill dictionary
        for symbol in symbols:
            symbol = symbol.upper()
            if symbol in panglao_ui_dict:
                # If the symbol is already in the dictionary, update the value to the average
                panglao_ui_dict[symbol] = round((panglao_ui_dict[symbol] + row['ubiquitousness index']) / 2, 3)
            else:
                panglao_ui_dict[symbol] = row['ubiquitousness index']

    return panglao_ui_dict


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
            ui_dict = transfer_ui_to_homologs(biomart_target=biomart_target, repo_path=repo_path)

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