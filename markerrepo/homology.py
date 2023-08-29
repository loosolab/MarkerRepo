import pandas as pd
import os
import urllib.request
from pybiomart import Server
from .marker_repo import get_gene_dict
from .marker_repo import combine_dfs, get_gene_dict, get_whitelists, guided_search, search_df, select, update_markers
from .plotting import plot_gene_counts
from .utils import read_whitelist
from IPython.display import display


def check_organisms(biomart_orgs, homologene_orgs, source_orgs, target_org):
    """
    Checks if the source organisms are supported by the BioMart or HomoloGene approach.

    Parameters
    ----------
    biomart_orgs : list
        List of organisms supported by the BioMart approach.
    homologene_orgs : list
        List of organisms supported by the HomoloGene approach.
    source_orgs : list
        List of source organisms that the user wants to use for gene transfer.
    target_org : str
        Name and tax ID of the target organism.

    Returns
    -------
    lists of str :
        - Organisms available in both approaches
        - Organisms available in neither approach
        - Organisms available only in the BioMart approach
        - Organisms available only in the HomoloGene approach
    """

    if target_org in source_orgs:
        source_orgs.remove(target_org)
        print(f"\nRemoved {target_org} from source organisms as it matches the target organism.")

    source_orgs = set(source_orgs)

    in_both = list(source_orgs.intersection(biomart_orgs).intersection(homologene_orgs))
    in_neither = list(source_orgs.difference(biomart_orgs).difference(homologene_orgs))
    only_in_biomart = list(source_orgs.intersection(biomart_orgs).difference(homologene_orgs))
    only_in_homologene = list(source_orgs.intersection(homologene_orgs).difference(biomart_orgs))

    return in_both, in_neither, only_in_biomart, only_in_homologene


def download_homologene_data(file_name="homologene.data", url="ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data", check=False):
    """
    Download and parse the HomoloGene data.

    Parameters
    ----------
    file_name : str, default "homologene.data"
        Name of the local HomoloGene data file.
    url : str, default "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
        URL to the HomoloGene data file.
    check : bool, default False
        If True, ask wether you want to overwrite the HomoloGene file.

    Returns
    --------
    pandas.DataFrame : 
        DataFrame with the HomoloGene data.
    """

    # Check if file already exists
    if os.path.exists(file_name):
        if check:
            overwrite = input(f"'{file_name}' already exists. Do you want to overwrite it? (yes/no): ").lower()

            if overwrite == 'yes':
                # Download new data and overwrite existing file
                urllib.request.urlretrieve(url, file_name)
    else:
        # Download data
        urllib.request.urlretrieve(url, file_name)

    # Load HomoloGene db
    homologene_data = pd.read_csv(file_name, sep='\t', header=None, index_col=0)

    # Rename columns
    homologene_data.columns = ["Taxonomy ID", "Gene ID", "Gene Symbol", "Protein GI", "Protein accession"]
    homologene_data.index.names = ["HID"]

    return homologene_data


def get_supported_taxonomy_ids(repo_path="."):
    """
    Returns all supported taxonomy IDs in the downloaded HomoloGene data.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.

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
    supported_organisms = read_whitelist("organism", repo_path=repo_path)['whitelist']

    for so in supported_organisms:
        name, tax = so.split(" ")
        if tax in unique_taxonomy_ids:
            organisms.append(f"{name} {tax}")

    return organisms


def get_target_genes_counts(df, gene_dict, source_column='Marker', target_column='Transferred Marker'):
    """
    Show plot and DataFrame of count of target genes per source gene.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame resulting from a gene transfer, containing columns for source genes and target genes.
    gene_dict : dict
        Dictionary containing Gene Symbols and Ensembl IDs.
    source_column : str, default 'Marker'
        The name of the column in df that contains the source genes.
    target_column : str, default 'Transferred Marker'
        The name of the column in df that contains the target genes.
    """

    gene_counts_df = get_transfer_counts(df, source_column=source_column, target_column=target_column)
    gene_counts_df = update_markers(gene_counts_df, gene_dict, column='Source Gene')

    return gene_counts_df


def process_and_filter_genes(transfer_counts_df, source_df, source_whitelist, target_counts=None, plots=False, id_type='symbol'):
    """
    Process and filter genes based on given conditions.

    Parameters
    ----------
    transfer_counts_df : pd.DataFrame
        DataFrame containing the transfer counts data.
    source_df : pd.DataFrame
        DataFrame containing the source genes data.
    source_whitelist : list
        List of source genes to be whitelisted.
    target_counts : int, optional
        Maximum count of target genes allowed. If not specified, all genes are included.
    plots : bool, default False
        Whether to generate plots or not.
    id_type : str, default 'symbol'
        The type of identifiers in the transfer_counts_df. 'symbol' for gene symbols and 'ensembl' for Ensembl IDs.

    Returns
    -------
    pd.DataFrame :
        DataFrame containing the filtered source genes data.
    """

    id_index = 0 if id_type == 'symbol' else 1

    # Calculate target genes count per source genes
    gene_counts_df = get_target_genes_counts(transfer_counts_df, get_gene_dict(w_markers=source_whitelist))
    print("\nCount of target genes per source gene:")
    display(gene_counts_df)

    # Plot count of target genes per source gene
    if plots:
        plot_gene_counts(gene_counts_df)

    if target_counts:
        filter_df = gene_counts_df.copy()
        filter_df = filter_df.loc[filter_df['Count target genes'] <= target_counts]
        filter_df['Source Gene'] = filter_df['Source Gene'].apply(lambda x: x.split(' ')[id_index] if len(x.split(' ')) > 1 else x)
        source_df_copy = source_df.copy()
        source_df_copy['Marker'] = source_df_copy['Marker'].apply(lambda x: x.split(' ')[id_index] if len(x.split(' ')) > 1 else x)
        filtered_source_df = pd.merge(source_df_copy, filter_df, left_on='Marker', right_on='Source Gene', how='inner')
        filtered_source_df = filtered_source_df[['Marker', 'Info']]

        print(f"Filtered source DataFrame with target_counts <= {target_counts}:")
        display(update_markers(filtered_source_df, get_gene_dict(w_markers=source_whitelist)))
        num_filtered_genes = source_df_copy.shape[0] - filtered_source_df.shape[0]

        print(f"Filtered: {num_filtered_genes} source genes, {100 - filtered_source_df.shape[0] / source_df_copy.shape[0] * 100:.2f}%\n")

        return filtered_source_df


def calculate_homology_proportions(df, source_genes, target_genes, id_type='symbol'):
    """
    Calculate the percentage of source and target genes present in a given DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame that contains gene identifiers for source and target organism.
        First column: source organism, second column: target organism.
    source_genes : list of str
        Whole set of genes from the source organism. Format: 'GeneName Ensembl ID'.
    target_genes : list of str
        Whole set of genes from the target organism. Format: 'GeneName Ensembl ID'.
    id_type : str, default 'symbol'
        The type of identifiers in the whitelist and df. 'symbol' for gene symbols and 'ensembl' for Ensembl IDs.

    Returns
    -------
    tuple of float
        Tuple containing the proportion of all possible transferred genes and of these genes in the target organism.
    """

    id_index = 0 if id_type == 'symbol' else 1

    # Split each string in the lists by space
    source_genes_ids = [gene.split()[id_index].upper() for gene in source_genes]
    target_genes_ids = [gene.split()[id_index].upper() for gene in target_genes]

    # Calculate the number of source and target genes
    num_source_genes = len(source_genes_ids)
    num_target_genes = len(target_genes_ids)

    # Find the intersection of the source/target genes and the DataFrame
    source_intersection = df.iloc[:, 0].isin(source_genes_ids)
    target_intersection = df.iloc[:, 1].isin(target_genes_ids)

    # Calculate the percentage of genes that are in the DataFrame
    source_percentage = (source_intersection.sum() / num_source_genes) * 100
    target_percentage = (target_intersection.sum() / num_target_genes) * 100

    return source_percentage, target_percentage


def transfer_markers_homologene(df, source_organism, target_organism, hg_db, target_whitelist, source_whitelist, calc_proportions=False, plots=False, target_counts=None):
    """
    Transfer markers between organisms based on homology and calculate the proportion of transferred markers.

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
    target_whitelist : list of str
        Whitelist containing all gene names of the target organism.
    source_whitelist : list of str
        Whitelist containing all gene names of the source organism. 
    calc_proportions : bool, default False
        If true, the proportions of all source genes and target genes are calculated.
    plots : bool, default False
        If true, show plots of transfer statistics.
    target_counts : int, default None
        If not None, filter those target genes whose number of target genes per source gene is <= target_counts.

    Returns
    --------
    DataFrame :
        Transferred markers DataFrame with columns corresponding to source_organism and target_organism.
    """

    df_copy = df.copy()

    # Adjust the Marker column in df_copy to contain only the first marker identifier (gene symbol)
    df_copy.loc[:, 'Marker'] = df_copy['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x).str.upper()

    # Filter homologene_data for the source and target organisms
    source_data = hg_db[hg_db['Taxonomy ID'] == int(source_organism)]
    target_data = hg_db[hg_db['Taxonomy ID'] == int(target_organism)]

    # Merge source and target data on HID
    merged_data = pd.merge(source_data, target_data, left_index=True, right_index=True, suffixes=('_source', '_target'))
    merged_data.rename(columns={'Gene Symbol_source': 'Marker', 'Gene Symbol_target': 'Transferred Marker'}, inplace=True)
    merged_data = merged_data[['Marker', 'Transferred Marker']]
    merged_data['Marker'] = merged_data['Marker'].str.upper()
    merged_data['Transferred Marker'] = merged_data['Transferred Marker'].str.upper()
    merged_data.drop_duplicates(inplace=True)

    # Merge input df_copy with merged_data on Marker
    target_df = pd.merge(df_copy, merged_data, on='Marker')
    target_df = target_df[['Transferred Marker', 'Info']]
    target_df['Transferred Marker'] = target_df['Transferred Marker'].str.upper()
    target_df.drop_duplicates(inplace=True)

    # Calculate the percentage of transferred markers
    marker_transfer_rate = target_df.shape[0] / df.shape[0] * 100
    print(f"The Marker Transfer Rate is {marker_transfer_rate:.2f}%, indicating the percentage of successfully transferred markers from the source to the target organism.")

    target_genes_list = [entry.split(' ')[0].upper() for entry in target_whitelist]
    filtered_df = target_df[target_df['Transferred Marker'].isin(target_genes_list)]

    # # Calculate the percentage of transferred markers after filtering
    # marker_filtered_transfer_rate = filtered_df.shape[0] / df_copy.shape[0] * 100
    # print(f"Marker transfer rate after filtering: {marker_filtered_transfer_rate:.2f}%")

    target_df.rename(columns={'Transferred Marker': 'Marker'}, inplace=True)

    if calc_proportions:
        if source_whitelist is not None and target_whitelist is not None:
            # Calculate the proportions
            possible_transfer_rate, transferred_genes_in_target = calculate_homology_proportions(merged_data, source_whitelist, target_whitelist, id_type='symbol')

            print(f"\nGeneral information:")
            print(f"Percentage of genes that can potentially be transferred from the source organism: {possible_transfer_rate:.2f}%")
            print(f"If all genes from the source organism were transferred, they would cover {transferred_genes_in_target:.2f}% of all genes of the target organism.")

    if target_counts:
        return process_and_filter_genes(merged_data, df, source_whitelist, target_counts=target_counts, plots=plots)

    markers_extended = update_markers(target_df, get_gene_dict(w_markers=target_whitelist))

    return markers_extended


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
    attributes = ['ensembl_gene_id', target_homolog_attribute]
    data = source_dataset.query(attributes=attributes)

    return data


def get_dataset_names(organism_name):
    """
    Retrieves the names of datasets corresponding to a specific organism.

    Parameters
    ----------
    organism_name : str
        The name of the organism to search for.

    Returns
    --------
    list :
        A list of dataset names that correspond to the input organism_name.
    """

    dataset_dict = create_dataset_dict()
    matching_names = [dataset_name for display_name, dataset_name in dataset_dict.items() if organism_name.lower() in display_name.lower()]

    return matching_names


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


def transfer_ui_to_homologs(target_organism="human", repo_path=".", biomart_target=None):
    """
    Transfer ubiquitousness index (ui) from source organism to target organism using homologous genes.

    Parameters
    ----------
    target_organism : str, default "human"
        Name of the target organism.
    repo_path : str, default "."
        The path of the Marker Repo.
    biomart_target : str, default None
        Name of the target organism in the BioMart database. If None, the name is inferred from the target_organism parameter.

    Returns
    -------
    pd.DataFrame :
        DataFrame with homologous genes and their ui.
    """

    if biomart_target is None:
        biomart_target = select(whitelist=get_dataset_names(target_organism), heading="BioMart target organism", repo_path=repo_path)

    biomart_dict = {'human': 'hsapiens', 'mouse': 'mmusculus'}
    panglao_organisms = ["human", "mouse"]

    transferred_panglao_dfs = []

    for organism in panglao_organisms:
        organism_ui_dict = get_panglao_ui(repo_path=repo_path, id_type='ensembl', organism=organism)
        biomart_source = biomart_dict[organism]

        # Fetch homologous genes
        organism_homologs_df = fetch_homologs(biomart_source, biomart_target)

        # Map "ui" values from the organism_ui_dict to the DataFrame
        organism_homologs_df['ui'] = organism_homologs_df['Gene stable ID'].map(organism_ui_dict)

        # Append DataFrame of the current organism to the list of DataFrames
        transferred_panglao_dfs.append(organism_homologs_df)

    # Concatenate all DataFrames in the list
    transferred_scores = pd.concat(transferred_panglao_dfs, ignore_index=True)

    # Remove rows with NaNs, 'Gene stable ID' column and duplicates
    transferred_scores = transferred_scores.dropna()
    transferred_scores = transferred_scores.drop(columns=['Gene stable ID'])
    transferred_scores = transferred_scores.drop_duplicates()

    # Extend genes
    gene_dict = get_gene_dict(organism=target_organism, repo_path=repo_path)
    transferred_scores = update_markers(transferred_scores, gene_dict, column=transferred_scores.columns[0])

    # Convert to dict with keys (gene symbols) = first column, values (UI) = second column
    ui_dict = pd.Series(transferred_scores[transferred_scores.columns[1]].values,
                            index=transferred_scores[transferred_scores.columns[0]].str.split(' ').str[0]).to_dict()
    return ui_dict


def create_dataset_dict():
    """
    Creates a dictionary mapping the display names of the datasets to their actual names.

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


def transfer_markers_biomart(biomart_df, source_df, target_whitelist, source_whitelist, calc_proportions=False, plots=False, target_counts=None):
    """
    This function merges two dataframes based on a common column and calculates the proportion of transferred markers.

    Parameters
    ----------
    biomart_df : pd.DataFrame
        DataFrame obtained from the BioMart database, with columns corresponding to 'Gene stable ID', 'Gene name', and a column containing the homologs of interest.
    source_df : pd.DataFrame
        Source DataFrame, with columns 'Marker', 'Info'. The 'Marker' column contains two gene names separated by a space.
    target_whitelist : list of str
        Whitelist containing all gene names and IDs of the target organism, separated by a space.
    source_whitelist : list of str
        Whitelist containing all gene names and IDs of the source organism, separated by a space.
    calc_proportions : bool, default False
        If true, the proportions of all source genes and target genes are calculated.
    plots : bool, default False
        If true, show plots of transfer statistics.
    target_counts : int, default None
        If not None, filter those target genes whose number of target genes per source gene is <= target_counts.

    Returns
    --------
    pd.DataFrame :
        Target DataFrame containing 'Marker' and 'Info' columns. The 'Marker' column contains the homologs of interest from the BioMart DataFrame.
    """

    source_df = source_df.copy()

    # Split the 'Marker' column and keep the ensembl ID
    source_df['Marker'] = source_df['Marker'].apply(lambda x: x.split(' ')[1] if len(x.split(' ')) > 1 else x)

    # Merge the two dataframes on the common column ('Marker' from source_df and 'Gene stable ID' from biomart_df)
    merged_df = pd.merge(biomart_df, source_df, left_on='Gene stable ID', right_on='Marker', how='inner')
    merged_df.drop_duplicates(inplace=True)
    transfer_counts_df = merged_df.rename(columns={biomart_df.columns[1]: 'Transferred Marker'}).copy()

    # Calculate the percentage of transferred markers
    marker_transfer_rate = merged_df.shape[0] / source_df.shape[0] * 100
    print(f"The Marker Transfer Rate is {marker_transfer_rate:.2f}%, indicating the percentage of successfully transferred markers from the source to the target organism.")

    target_genes_list = [entry.split(' ')[1] for entry in target_whitelist]
    filtered_df = merged_df[merged_df[biomart_df.columns[1]].isin(target_genes_list)]

    # # Calculate the percentage of transferred markers after filtering
    # marker_filtered_transfer_rate = filtered_df.shape[0] / source_df.shape[0] * 100
    # print(f"Marker transfer rate after filtering: {marker_filtered_transfer_rate:.2f}%")

    # Create the target dataframe
    target_df = merged_df[[biomart_df.columns[1], 'Info']].copy()
    target_df.rename(columns={biomart_df.columns[1]: 'Marker'}, inplace=True)
    target_df.drop_duplicates(inplace=True)

    if calc_proportions:
        if source_whitelist is not None:
            # Calculate the proportions
            possible_transfer_rate, transferred_genes_in_target = calculate_homology_proportions(biomart_df, source_whitelist, target_whitelist, id_type="ensembl")

            print(f"\nGeneral information:")
            print(f"Percentage of genes that can potentially be transferred from the source organism: {possible_transfer_rate:.2f}%")
            print(f"If all genes from the source organism were transferred, they would cover {transferred_genes_in_target:.2f}% of all genes of the target organism.")

    if target_counts:
        return process_and_filter_genes(transfer_counts_df, source_df, source_whitelist, target_counts=target_counts, plots=plots, id_type='ensembl')

    markers_extended = update_markers(target_df, get_gene_dict(w_markers=target_whitelist))

    return markers_extended


def get_supported_biomart_organisms(repo_path="."):
    """
    Returns all supported BioMart organisms in the downloaded Ensembl db.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.

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
    supported_organisms = read_whitelist("organism", repo_path=repo_path)['whitelist']

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


def get_transfer_counts(df, source_column='Marker', target_column='Transferred Marker'):
    """
    Counts the number of target genes for each source gene after a gene transfer.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame resulting from a gene transfer, containing columns for source genes and target genes.
    source_column : str, default 'Marker'
        The name of the column in df that contains the source genes.
    target_column : str, default 'Transferred Marker'
        The name of the column in df that contains the target genes.

    Returns
    -------
    pd.DataFrame :
        DataFrame with two columns: 'Source Gene' and 'Count target genes'.
    """

    # Count target genes per source gene, create DF, sort descending, return DF.
    df = df[[source_column, target_column]].drop_duplicates()
    gene_counts = df.groupby(source_column)[target_column].nunique()
    gene_counts_df = gene_counts.reset_index()
    gene_counts_df.columns = ['Source Gene', 'Count target genes']
    gene_counts_df.sort_values('Count target genes', ascending=False, inplace=True)

    return gene_counts_df


def prepare_gene_transfer(search_terms=None, case_sensitive=False, exact=False, repo_path="."):
    """
    Prepares the transfer of marker genes from a source organism to a target organism 
    by querying and fetching all necessary data.

    Parameters
    ----------
    search_terms : list of str
        Search terms to use for the search. Terms can be prefixed with '+' to denote that they must be included,
        or with '-' to denote that they must not be included. Terms without a prefix will include rows that contain them,
        but will not exclude rows that do not.
    case_sensitive : bool, default False
        If True, the search will be case-sensitive. If False, the search will be case-insensitive.
    exact : bool, default False
        If True, the search will look for exact matches. If False, the search will look for substrings.
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    -------
    misc

    if BioMart:
        biomart df, source df, target genes whitelist, source genes whitelist
    else:
        source df, source tax ID, target tax ID, HomoloGene db, target genes whitelist, source genes whitelist
    """

    get_whitelists()

    biomart_orgs = get_supported_biomart_organisms(repo_path=repo_path)
    homologene_orgs = get_supported_taxonomy_ids(repo_path=repo_path)

    db_choice = select_db(biomart_orgs, homologene_orgs)

    if db_choice == "biomart":
        organisms = biomart_orgs
    else:
        organisms = homologene_orgs

    source_organism, source_tax = select(whitelist=organisms, heading="source organism", repo_path=repo_path).split(" ")
    target_organism, target_tax = select(whitelist=organisms, heading="target organism", repo_path=repo_path).split(" ")
    print(f"Loading genes of {source_organism}...")
    source_genes = read_whitelist(f"genes/{source_organism}", repo_path=repo_path)['whitelist']
    print("Done!\n")
    print(f"Loading genes of {target_organism}...")
    target_genes = read_whitelist(f"genes/{target_organism}", repo_path=repo_path)['whitelist']
    print("Done!\n")

    if db_choice == "biomart":
        print("Specify BioMart organism selection:")
        source_organism_bm = select(whitelist=get_dataset_names(source_organism), heading="BioMart source organism", repo_path=repo_path)
        target_organism_bm = select(whitelist=get_dataset_names(target_organism), heading="BioMart target organism", repo_path=repo_path)

        print("Fetch necessary data from BioMart...")
        biomart_db = fetch_homologs(source_organism_bm, target_organism_bm).dropna()
    else:
        print("Get HomoloGene db...")
        hg_db = download_homologene_data()

    if search_terms:
        print("\nGenerate marker DataFrame based on the given search terms: ")
        for search_term in search_terms:
            print(search_term)
        source_df = search_df(combine_dfs(repo_path=repo_path), search_terms, case_sensitive=case_sensitive, exact=exact, out="marker_list", repo_path=repo_path)
    else:
        print("\nSelect the marker lists to be transferred to the target organism:")
        source_df = guided_search(repo_path=repo_path, out="marker_list")

    print(f"\nDataFrame of the markers of the source organism ({source_organism}) to be transferred to the target organism ({target_organism}):")
    display(source_df)

    if db_choice == "biomart":
        print(f'The preparations have been completed successfully.\nPlease continue by following the steps in the section "Transfer markers from one organism to another using BioMart".')
        return biomart_db, source_df, target_genes, source_genes
    else:
        print(f'The preparations have been completed successfully.\nPlease continue by following the steps in the section "Transfer markers from one organism to another using HomoloGene db".')
        return source_df, source_tax, target_tax, hg_db, target_genes, source_genes


def set_biomart_defaults(repo_path="."):
    """
    Processes the BioMart dataset dictionary to handle organisms with more than one value. In such cases, prompts the user to select a value.
    The processed dictionary is then saved to a file, with one value per line.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    -------
    str :
        The path of the saved file.
    """

    default_orgs = []

    orgs = get_supported_biomart_organisms(repo_path=repo_path)
    for org in orgs:
        biomart_orgs = (get_dataset_names(org.split(" ")[0]))
        if len(biomart_orgs) > 1:
            print(f"Please select the BioMart organism for {org}:")
            default_orgs.append(select(whitelist=biomart_orgs, heading="BioMart organism", repo_path=repo_path))

    defaults_path = os.path.join(repo_path, "biomart_defaults")
    with open(defaults_path, "w") as f:
        for org in default_orgs:
            f.write(f"{org}\n")

    return defaults_path


def get_biomart_defaults(repo_path="."):
    """
    Retrieves the default BioMart organisms from the biomart_defaults file.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    -------
    list of str :
        The default BioMart organisms.
    """

    defaults_path = os.path.join(repo_path, "biomart_defaults")
    if os.path.exists(defaults_path):
        with open(defaults_path, "r") as f:
            defaults = f.readlines()
        return [default.strip() for default in defaults]
    else:
        return []
    

def get_supported_organisms(repo_path="."):
    """
    Retrieves all supported organisms from the Biomart and HomoloGene databases.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    -------
    list of str :
        The supported organisms.
    """
    
    biomart_organisms = get_supported_biomart_organisms(repo_path=repo_path)
    homologene_organisms = get_supported_taxonomy_ids(repo_path=repo_path)
    supported_organisms = list(set(biomart_organisms + homologene_organisms))

    return supported_organisms