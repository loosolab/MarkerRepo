from .marker_repo import search_db, combine_lists, get_db
import pandas as pd
import os
import urllib.request
from .utils import read_whitelist
from sklearn.preprocessing import MinMaxScaler
from pybiomart import Server


def compare_marker_lists(repo_lists_path="./lists", keywords=None, marker_df=None, case_sensitive=False, exact=False):
    """
    This function compares and scores markers from selected marker lists using Ubiquitousness Index.
    A score of '0' signifies that the marker is the most specific within this selection, 
    while a score of '1' indicates that the marker is the most prevalent.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path to the directory containing the marker lists.
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
        df = search_db(get_db(repo_lists_path=repo_lists_path), keywords, case_sensitive=case_sensitive, exact=exact)
        uids = [int(idx) for idx in df.index]
        df = combine_lists(uids, repo_lists_path=repo_lists_path)

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


def download_homologene_data(file_name="homologene.data", url="ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"):
    """
    Download and parse the HomoloGene data.

    Parameters
    ----------
    file_name : str, default "homologene.data"
        Name of the local HomoloGene data file.
    url : str, default "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
        URL to the HomoloGene data file.

    Returns
    --------
    pandas.DataFrame : 
        DataFrame with the HomoloGene data.
    """

    # Check if file already exists
    if os.path.exists(file_name):
        overwrite = input(f"'{file_name}' already exists. Do you want to overwrite it? (yes/no): ").lower()
        
        if overwrite == 'no':
            # Load existing data
            homologene_data = pd.read_csv(file_name, sep='\t', header=None, index_col=0)
        else:
            # Download new data and overwrite existing file
            urllib.request.urlretrieve(url, file_name)
            homologene_data = pd.read_csv(file_name, sep='\t', header=None, index_col=0)
    else:
        # Download data
        urllib.request.urlretrieve(url, file_name)
        homologene_data = pd.read_csv(file_name, sep='\t', header=None, index_col=0)

    # Rename columns
    homologene_data.columns = ["Taxonomy ID", "Gene ID", "Gene Symbol", "Protein GI", "Protein accession"]
    homologene_data.index.names = ["HID"]

    return homologene_data


def get_supported_taxonomy_ids():
    """
    Returns all supported taxonomy IDs in the downloaded HomoloGene data.
    
    Returns
    -------
    list of str:
        List of supported taxonomy IDs.
    """
    
    organisms = []

    if os.path.exists("homologene.data"):
        homologene_data = pd.read_csv("homologene.data", sep='\t', header=None, index_col=0)
        homologene_data.columns = ["Taxonomy ID", "Gene ID", "Gene Symbol", "Protein GI", "Protein accession"]
        homologene_data.index.names = ["HID"]
    else:
        homologene_data = download_homologene_data()

    # Get unique taxonomy IDs from HomoloGene db and convert them to strings
    unique_taxonomy_ids = homologene_data['Taxonomy ID'].unique().astype(str).tolist()
    # Get support organisms from whitelist repository
    supported_organisms = read_whitelist("organism")['whitelist']

    for so in supported_organisms:
        name, tax = so.split(" ")
        if tax in unique_taxonomy_ids:
            organisms.append(f"{name} {tax}")

    return organisms


def transfer_markers(df, source_organism, target_organism, hg_db):
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
    hg_db : pd.DataFrame
        DataFrame containing the HomoloGene db.
        
    Returns
    --------
    DataFrame :
        Transferred markers DataFrame with columns corresponding to source_organism and target_organism.
    """

    # Create an explicit copy of df (to avoid SettingWithCopyWarning)
    df_copy = df.copy()

    # Adjust the Marker column in df_copy to contain only the first marker identifier (gene symbol)
    df_copy.loc[:, 'Marker'] = df_copy['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x).str.upper()    
    # Filter homologene_data for the source and target organisms
    source_data = hg_db[hg_db['Taxonomy ID'] == int(source_organism)]
    target_data = hg_db[hg_db['Taxonomy ID'] == int(target_organism)]
    
    # Merge source and target data on HID
    merged_data = pd.merge(source_data, target_data, left_index=True, right_index=True, suffixes=('_source', '_target'))
    merged_data.rename(columns={'Gene Symbol_source': 'Marker', 'Gene Symbol_target': 'Transferred Marker'}, inplace=True)
    merged_data['Marker'] = merged_data['Marker'].str.upper()

    # Merge input df_copy with merged_data on Marker
    result_df = pd.merge(df_copy, merged_data, on='Marker')
    result_df = result_df[['Transferred Marker', 'Info']]
    result_df['Transferred Marker'] = result_df['Transferred Marker'].str.upper()
    result_df.drop_duplicates(inplace=True)
    
    return result_df


def get_panglao_ui(panglao_file="panglao_markers", organism="Hs"):
    """
    Create a dictionary with gene symbols and nicknames as keys and average ubiquitousness index as values.
    Only considers rows with the given organism.

    Parameters
    ----------
    panglao_file : str
        Path to the panglao markers.
    organism : str
        Organism to consider (e.g. "Hs" or "Mm").

    Returns
    -------
    dict :
        Dictionary with gene symbols and nicknames as keys and average ubiquitousness index as values.
    """

    df = pd.read_csv(panglao_file, sep="\t")
    
    # Filter dataframe by organism
    df = df[df['species'].str.contains(organism, na=False)]
    
    panglao_ui_dict = {}

    for _, row in df.iterrows():
        # Get gene symbols and nicknames
        symbols = [row['official gene symbol']]
        if pd.notna(row['nicknames']):
            symbols.extend(row['nicknames'].split('|'))

        # Fill dictionary
        for symbol in symbols:
            symbol = symbol.upper()
            if symbol in panglao_ui_dict:
                # If the symbol is already in the dictionary, update the value to the average
                panglao_ui_dict[symbol] = round((panglao_ui_dict[symbol] + row['ubiquitousness index']) / 2, 3)
            else:
                panglao_ui_dict[symbol] = row['ubiquitousness index']

    return panglao_ui_dict


def update_scores(df, organism="Hs", panglao_file="panglao_markers"):
    """
    Update the scores in the dataframe using the ubiquitousness index from the panglao database.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns "Marker", "Info", and "Score". 
    organism : str
        Organism to consider when retrieving the ubiquitousness index (e.g. "Hs" or "Mm").
    panglao_file : str
        Path to the panglao markers.

    Returns
    -------
    pd.DataFrame :
        Updated dataframe with new scores.
    """

    # Retrieve the ubiquitousness index dictionary for the given organism
    ui_dict = get_panglao_ui(panglao_file, organism)

    # Split the "Marker" column and take the first part
    df['MainMarker'] = df['Marker'].str.split().str[0].str.upper()

    # Update scores where the main marker is in ui_dict
    df['Score'] = df['MainMarker'].map(ui_dict).fillna(df['Score'])


    df = df.drop(columns='MainMarker')
    df.sort_values('Score', ascending=True, inplace=True)

    return df


def fetch_homologs(source_organism, target_organism):
    """
    Fetch homologous genes using BioMart.
    
    Parameters
    ----------
    source_organism : str
        Name of the source organism.
    target_organism : str
        Name of the target organism.
    
    Returns
    -------
    pd.DataFrame :
        DataFrame containing homologous genes.
    """
    # Initialize BioMart server
    server = Server(host='http://www.ensembl.org')

    # Define source and target datasets
    source_dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets[source_organism + '_gene_ensembl']
    target_homolog_attribute = target_organism + '_homolog_ensembl_gene'

    # Query BioMart database
    attributes = ['ensembl_gene_id', 'external_gene_name', target_homolog_attribute]
    data = source_dataset.query(attributes=attributes)
    
    return data


def create_dataset_dict():
    """
    Creates a dictionary mapping the display names of the datasets to their actual names.
    
    Parameters
    ----------
    datasets : dict
        A dictionary of available datasets from the Biomart server.

    Returns
    --------
    dict :
        A dictionary with display names as keys and actual dataset names as values.
    """

    # Get the available datasets from the Biomart server
    server = Server(host='http://www.ensembl.org')
    datasets = server.marts['ENSEMBL_MART_ENSEMBL'].datasets

    dataset_dict = {}

    for name, dataset in datasets.items():
        dataset_dict[dataset.display_name] = name.split("_")[0]

    return dataset_dict


def get_dataset_names(organism_name):
    """
    Retrieves the names of datasets corresponding to a specific organism.
    
    Parameters
    ----------
    organism_name : str
        The name of the organism to search for.
    dataset_dict : dict
        A dictionary with display names as keys and actual dataset names as values.

    Returns
    --------
    list :
        A list of dataset names that correspond to the input organism_name.
    """

    dataset_dict = create_dataset_dict()
    matching_names = [dataset_name for display_name, dataset_name in dataset_dict.items() if organism_name.lower() in display_name.lower()]

    return matching_names


def transfer_markers_biomart(biomart_df, source_df):
    """
    This function merges two dataframes based on a common column.

    Parameters
    ----------
    biomart_df : pd.DataFrame
        DataFrame obtained from the BioMart database, with columns corresponding to 'Gene stable ID', 'Gene name', and a column containing the homologs of interest.
    source_df : pd.DataFrame
        Source DataFrame, with columns 'Marker', 'Info'. The 'Marker' column contains two gene names separated by a space.

    Returns
    --------
    pd.DataFrame :
        Target DataFrame containing 'Marker' and 'Info' columns. The 'Marker' column contains the homologs of interest from the BioMart DataFrame.
    """

    # Split the 'Marker' column and keep the ensembl ID
    source_df['Marker'] = source_df['Marker'].apply(lambda x: x.split(' ')[1] if len(x.split(' ')) > 1 else x)

    # Merge the two dataframes on the common column ('Marker' from source_df and 'Gene stable ID' from biomart_df)
    merged_df = pd.merge(biomart_df, source_df, left_on='Gene stable ID', right_on='Marker', how='inner')

    # Create the target dataframe
    target_df = merged_df[[biomart_df.columns[2], 'Info']]
    target_df.rename(columns={biomart_df.columns[2]: 'Marker'}, inplace=True)
    target_df.drop_duplicates(inplace=True)

    return target_df


def get_supported_biomart_organisms():
    """
    Returns all supported BioMart organisms in the downloaded Ensembl db.

    Returns
    -------
    list of str:
        List of supported organisms.
    """
    
    organisms = []

    # Get organisms from Ensembl db and convert them to strings
    dataset_list = list(create_dataset_dict().keys())
    ensembl_organisms = [s.split(" genes")[0].lower() for s in dataset_list]

    # Get support organisms from whitelist repository
    supported_organisms = read_whitelist("organism")['whitelist']

    for so in supported_organisms:
        name, tax = so.split(" ")

        if name.lower() in ensembl_organisms:
            organisms.append(f"{name} {tax}")

    return organisms


def select_db(biomart, homologene):
    """
    Ask the user to select a database from the provided list of supported organisms in each database.

    Parameters
    ----------
    biomart : list
        List of supported organisms in the Biomart database.
    homologene : list
        List of supported organisms in the HomoloGene database.

    Returns
    -------
    str :
        The chosen database, either "biomart" or "homologene".
    """

    while True:
        print("Supported organisms in the Biomart database:")
        for organism in biomart:
            print(organism)

        print("\nSupported organisms in the HomoloGene database:")
        for organism in homologene:
            print(organism)

        db_choice = input("\nPlease choose a database (biomart/homologene): ")

        if db_choice.lower() in ['biomart', 'homologene']:
            return db_choice
        else:
            print("\nInvalid choice. Please choose either 'biomart' or 'homologene'.")
