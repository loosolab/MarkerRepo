import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from .parsing import get_marker_lists, process_file, split_marker_elements, get_list, dataframe_to_dict
from .utils import read_whitelist
import yaml
from git import Repo, GitCommandError
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from IPython.display import display
import math
import re
import time
import random
import warnings


def search_df(df, search_terms, col_to_search=None, case_sensitive=False, exact=False, out="metadata", repo_path=".", lists_path=None, column_specific_terms=None, suffix=None):
    """
    This function filters a given DataFrame based on the provided keywords. Depending on the 'out' parameter,
    the function either returns the filtered DataFrame or a combined list of markers.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to be filtered.
    search_terms : list of str
        Search terms to use for the search. Terms can be prefixed with '+' to denote that they must be included,
        or with '-' to denote that they must not be included. Terms without a prefix will include rows that contain them,
        but will not exclude rows that do not.
    col_to_search : str, default None
        Column to perform the search in. If None, the search will be performed in all columns.
    case_sensitive : bool, default False
        If True, the search will be case-sensitive. If False, the search will be case-insensitive.
    exact : bool, default False
        If True, the search will look for exact matches. If False, the search will look for substrings.
    out : str, default "metadata"
        Determines the output of the function. If 'metadata', the function returns the filtered DataFrame. If
        'marker_list', the function returns a combined list of markers.
    repo_path : str, default "."
        The path of the Marker Repo. Required if out = 'marker_list'.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If none, the lists folder of the repo_path is used.
    column_specific_terms : dict, default None
        A dictionary with column names as keys and lists of search terms as values. If provided, 'col_to_search' and 'search_terms' are ignored.
    suffix : str, default None
        The key of the metadata section whose value should be appended to the marker names.

    Returns
    -------
    pd.DataFrame :
        Either the filtered search results as metadata or as a combined list of markers.
    """
    
    if not os.path.exists(repo_path):
        raise FileNotFoundError(f"The specified repository path '{repo_path}' does not exist.")
    
    if out not in ["metadata", "marker_list"]:
        raise ValueError("The parameter 'out' must be either 'metadata' or 'marker_list'.")
    
    df  = df.reset_index()
    
    # Convert col_to_search and search_terms to column_specific_terms if None
    if column_specific_terms is None:
        column_specific_terms = {col_to_search: search_terms} if col_to_search is not None else {col: search_terms for col in df.columns}

    flags = 0 if case_sensitive else re.IGNORECASE
    
    for col, terms in column_specific_terms.items():
        try:
            if col not in df.columns:
                print(f"Warning: The specified column '{col}' does not exist in the DataFrame. Skipping this column.")
                continue

            if isinstance(terms, str):
                terms = [terms]

            must_include_terms = [term.lstrip('+') for term in terms if term.startswith('+')]
            positive_terms = [term for term in terms if not term.startswith('-') and not term.startswith('+')]
            negative_terms = [term.lstrip('-') for term in terms if term.startswith('-')]

            if exact:
                must_include_terms = [f"^{term}$" for term in must_include_terms]
                positive_terms = [f"^{term}$" for term in positive_terms]
                negative_terms = [f"^{term}$" for term in negative_terms]

            # Apply filtering logic for each column
            for term in must_include_terms:
                df = df[df[col].astype(str).apply(lambda x: bool(re.search(term, x, flags=flags)))]

            if positive_terms:
                positive_mask = pd.concat([df[col].astype(str).apply(lambda x: bool(re.search(term, x, flags=flags))) for term in positive_terms], axis=1).any(axis=1)
                df = df[positive_mask]

            for term in negative_terms:
                df = df[~df[col].astype(str).apply(lambda x: bool(re.search(term, x, flags=flags)))]

        except Exception as e:
            print(f"An error occurred for column '{col}': {str(e)}")

    if out == "marker_list":
        if repo_path is None:
            raise ValueError("repo_path must be provided when out='marker_list'")
        uids = [int(idx) for idx in df.index]
        return combine_lists(uids, repo_path=repo_path, lists_path=lists_path, suffix=suffix)
    
    df.set_index("ID", inplace=True)

    return df


def guided_search(repo_path=".", lists_path=None, df=None, out="metadata", exact=False, case_sensitive=False, suffix=None):
    """
    An interactive function that guides the user through the process of searching the DataFrame.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    df : pd.DataFrame, default None
        The DataFrame to search in. If not provided, the function will create one from the repo_path.
    out : str, default "metadata"
        Determines the output of the function. If 'metadata', the function returns the filtered DataFrame. If
        'marker_list', the function returns a combined list of markers.
    exact : bool, default False
        If True, the search will look for exact matches. If False, the search will look for substrings.
    case_sensitive : bool, default False
        If True, the search will be case-sensitive. If False, the search will be case-insensitive.
    suffix : str, default None
        The key of the metadata section whose value should be appended to the marker names.

    Returns
    -------
    pd.DataFrame
        Either the filtered search results as metadata or as a combined list of markers.
    """

    if not os.path.exists(repo_path):
        raise FileNotFoundError(f"The specified repository path '{repo_path}' does not exist.")
    
    if out not in ["metadata", "marker_list"]:
        raise ValueError("The parameter 'out' must be either 'metadata' or 'marker_list'.")

    # Get the DataFrame if not provided
    if df is None:
        df = combine_dfs(repo_path=repo_path, lists_path=lists_path)

    df_copy = df.copy()
    df_copy = df_copy.reset_index()
    columns = df_copy.columns.tolist()
    page = 1
    per_page = 10
    num_pages = math.ceil(len(columns) / per_page)

    while True:
        start_index = (page - 1) * per_page
        end_index = start_index + per_page

        print("Available columns for search:")
        for i, col in enumerate(columns[start_index:end_index], start=start_index):
            print(f"{i+1}: {col}")

        if page < num_pages:
            print("n: Next page")
        if page > 1:
            print("p: Previous page")

        column = input("Enter identifier of column to search in (leave blank to search in all columns): ")
        if column.lower() == 'n' and page < num_pages:
            page += 1
            continue
        elif column.lower() == 'p' and page > 1:
            page -= 1
            continue

        col_to_search = None
        if column:
            try:
                col_to_search = columns[int(column) - 1]
            except (IndexError, ValueError):
                print("Invalid column identifier. Please try again.")
                continue

            while True:
                unique_entries = input("Do you want to see all unique entries in this column? (yes/no) ")
                if unique_entries.lower() in ['yes', 'no']:
                    if unique_entries.lower() == 'yes':
                        if df_copy[col_to_search].dtype == 'object':
                            unique_values = df_copy[col_to_search].explode().unique()
                            print("Unique entries:")
                            for val in unique_values:
                                print(val)
                        else:
                            print(df_copy[col_to_search].unique())
                        break
                    else:
                        break
                else:
                    print("Please answer with 'yes' or 'no'.")

        search_terms = input("Enter search terms (separated by commas, '-' for negative search): ")
        search_terms = [term.strip() for term in search_terms.split(",")]

        df_copy = search_df(df_copy, search_terms, col_to_search=col_to_search, case_sensitive=case_sensitive, exact=exact, suffix=suffix)
        print(f"Number of results: {len(df_copy)}")

        while True:
            see_results = input("Do you want to see the results? (yes/no): ")
            if see_results.lower() == 'yes':
                display(df_copy)
                break
            elif see_results.lower() == 'no':
                break
            else:
                print("Please answer with 'yes' or 'no'.")

        while True:
            continue_search = input("Do you want to continue searching? (yes/no): ")
            if continue_search.lower() in ['yes', 'no']:
                if continue_search.lower() == 'no':
                    if out == "marker_list":
                        if repo_path is None:
                            raise ValueError("repo_path must be provided when out='marker_list'")
                        uids = [int(idx) for idx in df_copy.index]
                        return combine_lists(uids, repo_path=repo_path, lists_path=lists_path, suffix=suffix)

                    return df_copy
                break
            else:
                print("Please answer with 'yes' or 'no'.")


def get_db(repo_path=".", lists_path=None, parallel=True):
    """
    Get the database of the Marker Repo as DataFrame, containing metadata information.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    parallel : bool, default True
        If True, uses parallel processing to improve performance.

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing metadata information of all lists.
    """

    if not lists_path:
        lists_path = f"{repo_path}/lists"

    # Check if the provided repository path exists
    if not os.path.exists(lists_path):
        raise FileNotFoundError(f"The specified path '{lists_path}' does not exist.")

    # Get all YAML file paths in the directory
    file_paths = [os.path.join(root, file) for root, dirs, files in os.walk(lists_path) for file in files if file.endswith(".yaml")]

    # Check if there are any YAML files in the path
    if not file_paths:
        raise FileNotFoundError(f"No YAML files found in the path '{lists_path}'.")
    
    if parallel:
        # Use a ProcessPoolExecutor to read and parse files in parallel
        with ProcessPoolExecutor() as executor:
            data = list(executor.map(process_file, file_paths))
    else:
        data = []

        for file_path in file_paths:
            # Read YAML file and extract all leaf values from "metadata"
            data.append(process_file(file_path))

    # Create DataFrame and set "ID" as index
    df = pd.DataFrame(data)
    df.rename(columns=get_display_names(repo_path=repo_path), inplace=True)
    if "ID" in df.columns:
        df["ID"] = pd.to_numeric(df["ID"])
        df.set_index("ID", inplace=True)
        df.sort_values("ID", inplace=True)

    return df


def combine_dfs(repo_path=".", lists_path=None, parallel=True, preprocessed=False, meta_lists="meta_lists", marker_lists="marker_lists"):
    """
    Combine the outputs of 'get_db' and 'get_marker_lists' based on the given columns.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
    parallel : bool, default True
        If True, uses parallel processing to improve performance.
    preprocessed : bool, default True
        If True, uses preprocessed data from .tsv files.
    meta_lists : str, default "meta_lists"
        Path to the preprocessed metadata .tsv file.
    marker_lists : str, default "marker_lists"
        Path to the preprocessed marker .tsv file.

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing combined information.
    """

    if preprocessed:
        if not (os.path.isfile(meta_lists) and os.path.isfile(marker_lists)):
            raise FileNotFoundError("Preprocessed .tsv files not found. Please check the file paths.")
        df_meta = pd.read_csv(meta_lists, sep='\t')
        df_marker = pd.read_csv(marker_lists, sep='\t')
    else:
        df_meta = get_db(repo_path=repo_path, lists_path=lists_path, parallel=parallel)
        df_marker = get_marker_lists(repo_path=repo_path, lists_path=lists_path, parallel=parallel)
    
    df_marker = df_marker.groupby('ID').agg({
        'Marker': lambda x: list(set(x)),
        'Info': lambda x: list(set(x))
    }).reset_index()
    
    df_combined = df_meta.merge(df_marker, on='ID', how='left')
    df_combined.set_index('ID', inplace=True)

    # Split marker elements
    df_combined['Marker'] = df_combined['Marker'].apply(split_marker_elements)

    return df_combined


def find_leaf_keys(dictionary, prefix=''):
    """
    Find all leaf keys in a nested dictionary.
    
    Parameters
    ----------
    dictionary : dict
        The dictionary to search through.
    prefix : str
        The prefix for the keys, used for nested dictionaries.
        
    Returns
    -------
    dict
        A dictionary with leaf keys (including their path if nested) and their values.
    """

    leaves = {}
    for key, value in dictionary.items():
        full_key = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            leaves.update(find_leaf_keys(value, full_key))
        else:
            leaves[full_key] = value

    return leaves


def get_marker_list(file_path, suffix=None):
    """
    Reads a YAML file containing a "marker_list" section and optionally appends a leaf key value from the metadata section to the "Info" column. 
    Warns if the suffix is provided but not a leaf key in the metadata.

    Parameters
    ----------
    file_path : str
        The path to the input YAML file.
    suffix : str, optional
        The leaf key of the metadata section whose value should be appended to the "name".

    Returns
    -------
    pd.DataFrame :
        The resulting DataFrame containing "Marker" and "Info" columns.
    """

    with open(file_path, 'r') as file:
        yaml_data = yaml.safe_load(file)
    
    
    metadata_suffix = ""
    if suffix:
        # Find all leaf keys in metadata
        leaf_keys = find_leaf_keys(yaml_data.get('metadata', {}))
        if suffix in leaf_keys:
            metadata_suffix = "_" + str(leaf_keys[suffix])
        else:
            warnings.warn(f"Suffix '{suffix}' not found as a leaf key in metadata.", UserWarning)
    
    marker_list = yaml_data['marker_list']
    
    marker_data = []
    for item in marker_list:
        name = item['name'] + metadata_suffix
        markers = item['markers']
        for marker in markers:
            marker_data.append({"Marker": marker, "Info": name})

    df = pd.DataFrame(marker_data)

    return df

def export_marker_list(df, path="./exported_lists", file_name=None, header=False, marker_id=None):
    """
    Exports a marker list (df) to path/file_name. If a file with this name already exists,
    a timestamp suffix is added to the filename.

    Parameters
    ----------
    df : pd.DataFrame
        The marker list to be exported.
    path : str, default "."
        The path where the marker list should be saved.
    file_name : str, default None
        The filename of the marker list.
    header : bool, default False
        If True, the header will also be exported
    marker_id : str, default None
        If "ensembl", only the second marker in the 'Marker' column is kept before exporting the df.
        If "symbol", only the first marker in the 'Marker' column is kept before exporting the df.

    Returns
    -------
    export_path : str
        The full path where the marker list was saved.
    """

    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Folder {os.path.abspath(path)} created.")

    if marker_id == "ensembl":
        df['Marker'] = df['Marker'].apply(lambda x: x.split(' ')[1] if len(x.split(' ')) > 1 else x.split(' ')[0])
    elif marker_id == "symbol":
        df['Marker'] = df['Marker'].apply(lambda x: x.split(' ')[0])

    if not path:
        path = "."

    if not file_name:
        file_name = get_valid_filename()

    # Generate the full file path
    export_path = os.path.join(path, f"{file_name}")
    
    # Check if a file with this name already exists
    if os.path.exists(export_path):
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        export_path = os.path.join(path, f"{file_name}_{timestamp}")

    # Export the marker list
    df.to_csv(export_path, sep="\t", index=False, header=header)
    print(f"Marker list saved: {os.path.abspath(export_path)}")
    
    return os.path.abspath(export_path)


def get_uid_paths(uids, repo_path=".", lists_path=None):
    """
    Searches for files in the specified folder and its subfolders with names in the format "name_UID.yaml",
    where UID is an integer. Returns the paths of the files that contain the UIDs from the given list.

    Parameters
    ----------
    uids : list of int
        A list of integers representing the UIDs to search for.
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If none, the lists folder of the repo_path is used.

    Returns
    -------
    list of str :
        A list of file paths containing the specified UIDs.
    """

    matching_files = []

    if not lists_path:
        lists_path = f"{repo_path}/lists"

    for root, _, files in os.walk(f"{lists_path}"):
        for file in files:
            if file.endswith('.yaml'):
                uid = int(file.split('_')[-1].split('.')[0])
                if uid in uids:
                    matching_files.append(os.path.join(root, file))

    return matching_files


def combine_lists(uids, repo_path=".", lists_path=None, suffix=None):
    """
    Combine multiple lists to one custom list.

    Parameters
    ----------
    uids : list of str
        The uids of the lists which will be combined.
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    suffix : str, default None
        The key of the metadata section whose value should be appended to the marker names.

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing the combinend list
    """

    # Read lists which are going to be combined
    dfs = []
    for file in get_uid_paths(uids, repo_path=repo_path, lists_path=lists_path):
        dfs.append(get_marker_list(file, suffix=suffix))

    if not dfs:
        print("No lists found.")
        return pd.DataFrame()
        
    # Perform outer join
    combined_df = pd.concat(dfs).reset_index(drop=True)
    combined_df.drop_duplicates(inplace=True)

    return combined_df


def update_markers(df, marker_dict, column='Marker'):
    """
    Updates markers by extending gene names with ensembl IDs and the other way round.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing the marker list.
    marker_dict : dict
        Dictionary containing the names and IDs as keys and values.
    column : str, default 'Marker'
        The column which is going to be extended.

    Returns
    --------
    pandas.DataFrame :
        The DataFrame containing updated markers and corresponding information.
    """

    df = df.copy()
    df[column] = df[column].apply(lambda marker: marker + ' ' + marker_dict[marker] if marker in marker_dict else marker)

    # Make sure that the Ensembl ID is always in second position
    df[column] = df[column].apply(lambda x: ' '.join(x.split(' ')[::-1]) if x.split(' ')[0].startswith('ENS') else x)

    return df


def select(whitelist=None, key=None, heading=None, repo_path="."):
    """
    Shows selection of whitelist and returns selected value.
    If only a key is passed, the corresponding whitelist is used as a selection.
    If a whitelist (list of strings) is passed, the key is the heading only.

    Parameters
    ----------
    whitelist : list of str, default None
        The selection to choose from.
    key : str, default None
        The key of the whitelist. For example "organism".
    heading: str, default None
        The heading (description) of the whitelist.
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    --------
    str :
        The selected string of the whitelist.
    """

    if not whitelist:
        if key:
            whitelist = read_whitelist(key, repo_path=repo_path)['whitelist']
        elif heading:
            raise Exception(f"No values for '{heading}' available. Please try again using other parameters.")
        else:
            raise Exception("No values available. Please try again using other parameters.")

    if len(whitelist) > 1:
        if not heading:
            print(f"Select {key}")
        else:
            print(f"Select {heading}")

        for i, value in enumerate(whitelist):
            print(str(i+1) + ":\t" + value)

    if len(whitelist) == 1:
        selection = whitelist[0]
    else:
        while True:
            try:
                user_input = int(input())
                if user_input < 1 or user_input > len(whitelist):
                    raise ValueError
                else:
                    selection = whitelist[user_input-1]
                    break
            except ValueError:
                print("Invalid input. Please enter a number corresponding to the options above.")
            
    print(f"Selection: {selection}\n")

    return selection


def get_gene_dict(organism=None, w_markers=None, repo_path="."):
    """
    Creates dictionary of whitelist of genes of specific organism.

    Parameters
    ----------
    organism : str, default None
        The organism that owns the corresponding genes.
    w_markers : list of str, default None
        A list of Gene Symbols and Ensembl IDs separated by space.
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    --------
    dict :
        Dictionary which contains the gene names and ensembl IDs.
    """

    gene_dict = {}

    if organism:
        if len(organism.split(' ')) > 1:
            w_markers = read_whitelist(f"genes/{organism.split(' ')[0]}", repo_path=repo_path)['whitelist']
        else:
            w_markers = read_whitelist(f"genes/{organism}", repo_path=repo_path)['whitelist']
    elif not w_markers:
        raise ValueError("Provide organism or whitelist of genes (w_markers).")

    for marker in w_markers:
        name, ensg = marker.split(" ")[0].upper(), marker.split(" ")[1].upper()
        gene_dict[name] = ensg
        gene_dict[ensg] = name

    return gene_dict


def get_uid():
    """
    Creates a new ID by combining a timestamp and a random number.

    Returns
    --------
    str :
        A string containing the new ID.
    """

    timestamp = int(time.time())
    random_number = random.randint(1, 1000000)

    # Combine the timestamp and random number to create a (almost) unique ID
    uid = f"{timestamp}{random_number}"

    return uid


def extract_display_names(d, parent_key='', sep='_'):
    """
    Extract the display names from a nested dictionary.

    Parameters
    ----------
    d : dict
        The input dictionary.
    parent_key : str, default ''
        The parent key used during recursion.
    sep : str, default '_'
        The separator used to concatenate keys.

    Returns
    -------
    dict :
        A dictionary containing the display names as values and concatenated names as keys.
    """

    display_names = {}
    for k, v in d.items():
        # Determine the new key. If the current key is "value", keep the parent key.
        new_key = parent_key if k == "value" else parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            # Extract display name of the current key if exists
            display_name = v.get('display_name')
            if display_name:
                display_names[new_key] = display_name
            # Extract display names from the nested dictionary
            display_names.update(extract_display_names(v, new_key, sep=sep))

    return display_names


def get_display_names(repo_path="."):
    """
    Extracts the display names of the keys.yaml.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns:
    -------
    dict :
        A dictionary containing the display names as values and concatenated names as keys.
    """

    # Load the YAML file
    with open(f"{repo_path}/keys.yaml", 'r', encoding='utf-8') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)

    # Extract the display names
    display_names = extract_display_names(yaml_data['metadata'])

    return display_names


def push_marker_list(list_path, repo_path="."):
    """
    Creates a new branch with the given list name, adds the new list,
    commits the changes and pushes the new branch to the remote repository.

    Parameters
    ----------
    list_path : str
        The path of the new list that is to be added. The list name and branch name
        will be extracted from this path.
    repo_path : str, default "."
        The path of the Marker Repo.
    """

    # Extract the list name from the list path
    list_name = os.path.splitext(os.path.basename(list_path))[0]  # Removes the .yaml extension

    try:
        print("Opening repository")
        repo = Repo(repo_path)
        assert not repo.bare

        # Pull the latest changes
        print("Pulling latest changes from the repository")
        repo.remotes['origin'].pull()

        # Check out new branch
        print(f"Creating and checking out new branch '{list_name}'")
        repo.git.checkout('HEAD', b=list_name)

        # Add and commit changes
        print("Adding new list and committing changes")
        repo.git.add(list_path)
        repo.git.commit('-m', f'Add new list: {list_name}')

        # Push changes
        print("Pushing changes to the repository")
        repo.git.push('--set-upstream', 'origin', list_name)

        print("Changes pushed successfully.")

    except GitCommandError as e:
        print(f"An error occurred while pushing to the repository: {e}")
    except AssertionError:
        print(f"The provided path '{repo_path}' does not appear to be a valid Git repository.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def preprocess_lists_to_tsv(repo_path=".", output_path_meta="meta_lists.tsv", output_path_markers="marker_lists.tsv"):
    """
    Generates preprocessed .tsv files of the database and marker lists.
    This function can be used to speed up subsequent reads of the data.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    output_path_meta : str, default "meta_lists.tsv"
        The file path for output metadata .tsv file.
    output_path_markers : str, default "marker_lists.tsv"
        The file path for output marker .tsv file.

    Returns
    -------
    str, str :
        The absolute paths to the generated .tsv files (metadata and markers respectively).
    """

    # Get metadata and marker lists
    df_meta = get_db(repo_path=repo_path)
    df_markers = get_marker_lists(repo_path=repo_path)
    
    # Get absolute file paths
    output_path_meta = os.path.abspath(output_path_meta)
    output_path_markers = os.path.abspath(output_path_markers)

    # Write to tsv
    df_meta.to_csv(output_path_meta, sep='\t', index=True)
    df_markers.to_csv(output_path_markers, sep='\t', index=True)

    return output_path_meta, output_path_markers


def get_valid_filename(prompt="Enter file name or path: "):
    """
    Allows the user to enter a path or file name, which will be checked for validity.

    Parameters
    ----------
    prompt : str, default "Enter file name: "
        The text of the prompt.

    Returns
    -------
    str :
        A valid path or file name.
    """

    while True:
        file_name = input(prompt)
        
        # Check if file name is not empty
        if not file_name.strip():
            print("File name or path cannot be empty.")
            continue

        # Check if file name ends with a slash
        if file_name.endswith("/"):
            print("File name or path cannot end with a slash.")
            continue


        return file_name


def get_selected_lists(keywords=None, metadata_df=None, repo_path=".", lists_path=None, case_sensitive=False, exact=False, order=["Info", "Marker"], suffix=None):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.

    Parameters
    ----------
    keywords : dict or str, default None
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    metadata_df : pd.DataFrame, default None
        A DataFrame containing the metadata of a selection of marker lists.
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    case_sensitive : bool, default False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.
    order : list of str, default ["Info", "Marker"]
        The order of the columns in the resulting DataFrame.
    suffix : str, default None
        The key of the metadata section whose value should be appended to the marker names.
        
    Returns
    --------
    pd.DataFrame :
        The DataFrame conaining the combined lists.
    """

    if keywords:
        db = combine_dfs(repo_path=repo_path, lists_path=lists_path)
        df = search_df(db, keywords, case_sensitive=case_sensitive, exact=exact, suffix=suffix)
        if df.empty:
            raise Exception(
                        f"No search results available!")
    elif metadata_df is not None:
        df = metadata_df
    else:
        raise Exception(
                    f"You need to specify keywords or a metadata DataFrame!")


    # Get UIDs and combine lists
    uids = [int(idx) for idx in df.index]
    combined_df = combine_lists(uids, repo_path=repo_path, lists_path=lists_path, suffix=suffix)

    # Drop duplicates, keep one marker only, rearrange column order
    markers_filtered = combined_df.drop_duplicates()
    markers_filtered = markers_filtered[order]

    return markers_filtered


def update_organism(organism, repo_path):
    """
    Updates the organism input to include both the organism name and the taxon ID, based on the whitelist.

    Parameters
    ----------
    organism : str or int
        The input organism as either the name, taxon ID, or both.
    repo_path : str
        Path to the marker repository.

    Returns
    -------
    str :
        The updated organism name and taxon ID as a string, or None if not found.
    """
    
    valid_organisms = read_whitelist("organism", repo_path=repo_path)['whitelist']
    organism_str = str(organism).lower()
    
    for entry in valid_organisms:
        entry_lower = entry.lower()
        if organism_str in entry_lower:
            return entry
            
    return None


def check_ensembl(adata, verbose=False):
    """
    Checks if the majority of gene identifiers in the .var index of an Anndata object are Ensembl IDs.
    Returns True if Ensembl IDs are the majority.

    Parameters
    ----------
    adata : anndata.AnnData
        An Anndata object containing the gene expression data.
    verbose : bool, default False
        If True, prints the number of Ensembl IDs and other IDs.

    Returns
    -------
    bool :
        True if Ensembl IDs are the majority in the .var index of the Anndata object. 
    """

    ensembl_pattern = re.compile(r'ENS[A-Z]{0,1}G\d+')

    ensembl_count = 0
    other_count = 0

    for gene_id in adata.var.index:
        if ensembl_pattern.match(gene_id):
            ensembl_count += 1
        else:
            other_count += 1

    if verbose:
        print(f"Ensembl IDs: {ensembl_count}")
        print(f"Other IDs: {other_count}")

    if ensembl_count > other_count:
        return True
    else:
        return False


def delete_files(file_paths):
    """
    Delete files from a list of file paths.

    Parameters
    ----------
    file_paths : list of str
        A list containing the file paths of the files to be deleted.
    """

    for file_path in file_paths:
        try:
            os.remove(file_path)
            print(f"File deleted: {file_path}")
        except FileNotFoundError:
            print(f"File not found: {file_path}")
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")


def transform_marker_list(LIST_PATH, info_col, marker_col, MARKER_TYPE, ORGANISM):
    """
    Transforms a marker list by filtering and extending markers based on their type and gene dictionary.

    Parameters
    ----------
    LIST_PATH : str
        The path to the marker list.
    info_col : int
        The column number in the list that contains additional information.
    marker_col : int
        The column number in the list that contains markers.
    MARKER_TYPE : str
        The type of markers, e.g., "Genes" or "Genomic regions".
    ORGANISM : str
        The organism name to be used for fetching the gene dictionary.

    Returns
    -------
    list of dict :
        A list of dictionaries, each containing a marker name and its corresponding markers.
    """

    markers = get_list(LIST_PATH, info_col=info_col, marker_col=marker_col, marker_type=MARKER_TYPE).drop_duplicates()
    print("All markers of provided list:")
    display(markers)

    if MARKER_TYPE == "Genes":
        markers['Marker'] = markers['Marker'].str.upper()
        gene_dict = get_gene_dict(ORGANISM)
        markers_removed = markers[~markers['Marker'].isin(gene_dict.keys())]
        print("Removed markers:")
        display(markers_removed)

        markers_filtered = markers[markers['Marker'].isin(gene_dict.keys())]
        markers_extended = update_markers(markers_filtered, gene_dict)
        print("Filtered and extended markers:")
        display(markers_extended)

        marker_dict = dataframe_to_dict(markers_extended)
    elif MARKER_TYPE == "Genomic regions":
        # genomic regions support
        marker_dict = dataframe_to_dict(markers)
    else:
        raise ValueError(f"Unsupported marker type: {MARKER_TYPE}")
    
    marker_list = []
    for name in marker_dict.keys():
        marker_list.append({'name': name, 'markers': marker_dict[name]})

    return marker_list
