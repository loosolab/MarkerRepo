import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import git
from .utils import read_whitelist
import yaml
from git import Repo
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from IPython.display import display
import math
import re

def search_df(df, search_terms, col_to_search=None, case_sensitive=False, exact=False, out="metadata", repo_lists_path="./lists"):
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
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
        Required if out = 'marker_list'.

    Returns
    -------
    pd.DataFrame :
        Either the filtered search results as metadata or as a combined list of markers.
    """

    must_include_terms = [term.lstrip('+') for term in search_terms if term.startswith('+')]
    positive_terms = [term for term in search_terms if not term.startswith('-') and not term.startswith('+')]
    negative_terms = [term.lstrip('-') for term in search_terms if term.startswith('-')]

    if exact:
        must_include_terms = [f"^{term}$" for term in must_include_terms]
        positive_terms = [f"^{term}$" for term in positive_terms]
        negative_terms = [f"^{term}$" for term in negative_terms]

    if col_to_search:
        for term in must_include_terms:
            df = df[df[col_to_search].astype(str).apply(lambda x: bool(re.search(term, x, flags=0 if case_sensitive else re.IGNORECASE)))]
        for term in positive_terms:
            df = df[df[col_to_search].astype(str).apply(lambda x: bool(re.search(term, x, flags=0 if case_sensitive else re.IGNORECASE)))]
        for term in negative_terms:
            df = df[~df[col_to_search].astype(str).apply(lambda x: bool(re.search(term, x, flags=0 if case_sensitive else re.IGNORECASE)))]
    else:
        for term in must_include_terms:
            df = df[df.apply(lambda x: x.astype(str).str.contains(term, flags=0 if case_sensitive else re.IGNORECASE, regex=True).any(), axis=1)]
        for term in positive_terms:
            df = df[df.apply(lambda x: x.astype(str).str.contains(term, flags=0 if case_sensitive else re.IGNORECASE, regex=True).any(), axis=1)]
        for term in negative_terms:
            df = df[~df.apply(lambda x: x.astype(str).str.contains(term, flags=0 if case_sensitive else re.IGNORECASE, regex=True).any(), axis=1)]

    if out == "marker_list":
        if repo_lists_path is None:
            raise ValueError("repo_lists_path must be provided when out='marker_list'")
        uids = [int(idx) for idx in df.index]
        return combine_lists(uids, repo_lists_path=repo_lists_path)

    return df


def guided_search(repo_lists_path="./lists", df=None, out="metadata"):
    """
    An interactive function that guides the user through the process of searching the DataFrame.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored.
    df : pd.DataFrame, default None
        The DataFrame to search in. If not provided, the function will create one from the repo_lists_path.
    out : str, default "metadata"
        Determines the output of the function. If 'metadata', the function returns the filtered DataFrame. If
        'marker_list', the function returns a combined list of markers.

    Returns
    -------
    pd.DataFrame
        Either the filtered search results as metadata or as a combined list of markers.
    """

    # Get the DataFrame if not provided
    if df is None:
        df = combine_dfs(repo_lists_path=repo_lists_path)

    df_copy = df.copy()
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
            col_to_search = columns[int(column) - 1]

            unique_entries = input("Do you want to see all unique entries in this column? (yes/no) ")
            if unique_entries.lower() == 'yes':
                if df_copy[col_to_search].dtype == 'object':
                    unique_values = df_copy[col_to_search].explode().unique()
                    print("Unique entries:")
                    for val in unique_values:
                        print(val)
                else:
                    print(df_copy[col_to_search].unique())

        search_terms = input("Enter search terms (separated by commas, '-' for negative search): ")
        search_terms = [term.strip() for term in search_terms.split(",")]

        exact = input("Perform an exact search? (yes/no): ")
        exact = exact.lower() == 'yes'

        case_sensitive = input("Consider case sensitivity? (yes/no): ")
        case_sensitive = case_sensitive.lower() == 'yes'

        df_copy = search_df(df_copy, search_terms, col_to_search=col_to_search, case_sensitive=case_sensitive, exact=exact)
        print(f"Number of results: {len(df_copy)}")

        see_results = input("Do you want to see the results? (yes/no): ")
        if see_results.lower() == 'yes':
            display(df_copy)

        continue_search = input("Do you want to continue searching? (yes/no): ")
        if continue_search.lower() != 'yes':
            break
    
    if out == "marker_list":
        if repo_lists_path is None:
            raise ValueError("repo_lists_path must be provided when out='marker_list'")
        uids = [int(idx) for idx in df_copy.index]
        return combine_lists(uids, repo_lists_path=repo_lists_path)

    return df_copy


def get_db(repo_lists_path="./lists", parallel=True):
    """
    Get the database of the Marker Repo as DataFrame, containing metadata information.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    parallel : bool, default True
        If True, uses parallel processing to improve performance.

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing metadata information of all lists.
    """
    
    file_paths = [os.path.join(root, file) for root, dirs, files in os.walk(repo_lists_path) for file in files if file.endswith(".yaml")]

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
    df.rename(columns=get_display_names(repo_lists_path.split("/lists")[0]), inplace=True)
    if "ID" in df.columns:
        df["ID"] = pd.to_numeric(df["ID"])
        df.set_index("ID", inplace=True)
        df.sort_values("ID", inplace=True)

    return df


def process_file(file_path):
    """
    Reads and parses a marker list file.
    
    Parameters
    ----------
    file_path : str
        Path of the marker list file (yaml-file).
    
    Returns
    -------
    dict :
        A dictionary containing flattened metadata information from the file.
    """

    with open(file_path, 'r', encoding='utf-8') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)
        metadata = yaml_data.get("metadata", {})
        flattened_metadata = flatten_dict(metadata)

        return flattened_metadata
    

def get_marker_lists(repo_lists_path="./lists", parallel=True):
    """
    Get the marker list from the Marker Repo as DataFrame.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    parallel : bool, default True
        If True, uses parallel processing to improve performance.

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing all markers and their designations.
    """
    
    file_paths = [os.path.join(root, file) for root, dirs, files in os.walk(repo_lists_path) for file in files if file.endswith(".yaml")]

    if parallel:
        # Use a ProcessPoolExecutor to read and parse files in parallel
        with ProcessPoolExecutor() as executor:
            data = list(executor.map(process_file_markers, file_paths))
    else:
        data = []
        for file_path in file_paths:
            data.append(process_file_markers(file_path))

    # Flatten the list of lists into a single list and create DataFrame
    data = [item for sublist in data for item in sublist]
    df = pd.DataFrame(data)
    df["ID"] = pd.to_numeric(df["ID"])
    df.set_index("ID", inplace=True)

    return df


def process_file_markers(file_path):
    """
    Reads and parses a marker list file to extract the marker list.
    
    Parameters
    ----------
    file_path : str
        Path of the marker list file (yaml-file).
    
    Returns
    -------
    list of dict :
        A list of dictionaries containing marker information for each marker, with 'Marker', 'Info', and 'ID' as keys.
    """
    
    with open(file_path, 'r', encoding='utf-8') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)
        marker_list_data = yaml_data.get("marker_list", [])

        # Extract ID from filename
        filename = os.path.basename(file_path)
        id = filename.split("_")[-1].replace(".yaml", "")

        # Extract markers and their names
        markers_data = []
        for item in marker_list_data:
            markers = item.get("markers", [])
            name = item.get("name", "")
            for marker in markers:
                markers_data.append({"Marker": marker, "Info": name, "ID": id})

        return markers_data


def split_marker_elements(marker_list):
    """
    This function takes in a list of marker elements and splits each element into two separate elements if 
    it contains two strings separated by a space. 
    
    Parameters
    ----------
    marker_list : list of str
        The list of markers that should be split into separate elements.

    Returns
    -------
    list of str :
        The list containing split elements - one element -> one marker
    """

    new_marker_list = []
    for marker in marker_list:
        # Split the marker into two elements if it contains a space
        new_marker_list.extend(marker.split())
    return new_marker_list


def combine_dfs(repo_lists_path="./lists", parallel=True, preprocessed=False, meta_lists="meta_lists", marker_lists="marker_lists"):
    """
    Combine the outputs of 'get_db' and 'get_marker_lists' based on the given columns.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the marker lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
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
        df_meta = get_db(repo_lists_path=repo_lists_path, parallel=parallel)
        df_marker = get_marker_lists(repo_lists_path=repo_lists_path, parallel=parallel)
    
    df_marker = df_marker.groupby('ID').agg({
        'Marker': lambda x: list(set(x)),
        'Info': lambda x: list(set(x))
    }).reset_index()
    
    df_combined = df_meta.merge(df_marker, on='ID', how='left')
    df_combined.set_index('ID', inplace=True)

    # Split marker elements
    df_combined['Marker'] = df_combined['Marker'].apply(split_marker_elements)

    return df_combined


def get_marker_list(file_path):
    """
    Reads a YAML file containing a section named "marker_list". The "marker_list" section consists of a list,
    where each element contains the keys "name" and "markers". The key "name" contains a string, and the key
    "markers" contains a list of strings. The function returns a DataFrame with two columns: "Marker" and "Info".
    The "Marker" column contains all elements of the "markers" key values, and the "Info" column contains the
    corresponding "name" key values.

    Parameters
    ----------
    file_path : str
        The path to the input YAML file.

    Returns
    -------
    pd.DataFrame :
        The resulting DataFrame containing "Marker" and "Info" columns.
    """

    with open(file_path, 'r') as file:
        yaml_data = yaml.safe_load(file)

    marker_list = yaml_data['marker_list']

    marker_data = []
    for item in marker_list:
        name = item['name']
        markers = item['markers']
        for marker in markers:
            marker_data.append({"Marker": marker, "Info": name})

    df = pd.DataFrame(marker_data)

    return df


def flatten_dict(d, parent_key='', sep='_', list_sep='\n'):
    """
    Flatten a nested dictionary, concatenating keys with a separator.

    Parameters
    ----------
    d : dict
        The input dictionary to be flattened.
    parent_key : str, default ''
        The parent key used during recursion.
    sep : str, default '_'
        The separator used to concatenate keys.
    list_sep : str, default '\n'
        The separator used to join list elements in a single cell.

    Returns
    dict :
        The flattened dictionary with concatenated keys.
    """

    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, dict):
            # If the value is a dictionary, flatten it
            items.extend(flatten_dict(v, new_key, sep=sep, list_sep=list_sep).items())
        elif isinstance(v, list):
            # If the value is a list, process its elements
            list_items = []
            for i, elem in enumerate(v):
                if isinstance(elem, dict):
                    # If the element is a dictionary, join its key-value pairs with the list separator
                    list_items.append(list_sep.join(f"{key}: {value}" for key, value in elem.items()))
                else:
                    # Otherwise, convert the element to a string
                    list_items.append(str(elem))
            # Join the list items with the list separator and store them in a single cell
            items.append((new_key, list_sep.join(list_items)))
        else:
            # If the value is not a dictionary or a list, store it directly
            items.append((new_key, v))
            
    return dict(items)


def get_list(path, info_col=1, marker_col=0):
    """
    Reads the marker lists and converts it to a DataFrame using the information
    of info_col and marker_col.

    Parameters
    ----------
    path : str
        The path where the list is stored.
    info_col : int, default 1
        The column which contains additional information like cell type or phase.
    marker_col: int, default 0
        The column which contains the marker (gene or genomic region).

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing the list
    """

    headers = ["Marker", "Info"]
    if type(info_col) == int:
        header = [headers[marker_col], headers[info_col]]
        df = pd.read_csv(path, sep='\t', names=header)
    else:
        df = pd.read_csv(path, sep='\t', names=[headers[0]])
        df[headers[1]] = info_col

    return df


def dataframe_to_dict(df):
    """
    Converts DataFrame of marker list to dictionary,
    using info_col as keys and marker_col as values.

    # TODO only marker column available

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing the marker list.

    Returns
    --------
    dict :
        The dictionary containing markers and corresponding information
    """

    result = {}
    for index, row in df.iterrows():
        key = row["Info"]
        value = row["Marker"]
        
        if key in result:
            result[key].append(value)
        else:
            result[key] = [value]

    return result


def export_marker_list(df, path=".", file_name=None, header=False, marker_id=None):
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

    if marker_id == "ensembl":
        df['Marker'] = df['Marker'].apply(lambda x: x.split(' ')[1])
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


def get_uid_paths(uids, repo_lists_path="./lists"):
    """
    Searches for files in the specified folder and its subfolders with names in the format "name_UID.yaml",
    where UID is an integer. Returns the paths of the files that contain the UIDs from the given list.

    Parameters
    ----------
    uids : list of int
        A list of integers representing the UIDs to search for.
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.

    Returns
    -------
    list of str :
        A list of file paths containing the specified UIDs.
    """

    matching_files = []

    for root, _, files in os.walk(repo_lists_path):
        for file in files:
            if file.endswith('.yaml'):
                uid = int(file.split('_')[-1].split('.')[0])
                if uid in uids:
                    matching_files.append(os.path.join(root, file))

    return matching_files


def combine_lists(uids, repo_lists_path="./lists"):
    """
    Combine multiple lists to one custom list.

    Parameters
    ----------
    uids : list of str
        The uids of the lists which will be combined.
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing the combinend list
    """

    # Read lists which are going to be combined
    dfs = []
    for file in get_uid_paths(uids, repo_lists_path=repo_lists_path):
        dfs.append(get_marker_list(file))
        
    # Perform outer join
    combined_df = pd.concat(dfs).reset_index(drop=True)
    combined_df.drop_duplicates(inplace=True)
    
    # TODO: inner join, etc ...

    return combined_df


def show_statistics(metadata, repo_lists_path="./lists", dpi=120):
    """
    Shows content of whole Marker Repo.

    Parameters
    ----------
    metadata : dict
        The dictionary containing the metadata information.
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    dpi : int, default 120
    """

    # TODO tissue plot - show count only

    sns.set_style("darkgrid")
    sns.set(rc={"figure.dpi": dpi, "savefig.dpi": dpi})
    fig, axes = plt.subplots(2, 3)
    axes_arr = [(0,0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    # Load all lists
    df = getDB(repo_lists_path)

    # Keep values which are not None
    filters = {}
    for key in metadata:
        if metadata[key]:
            filters[key] = metadata[key]

    # Plot statistics
    for count, key in enumerate(filters):
        stat_df = df[key].value_counts()
        stat_df.plot(kind='bar', title=key, ax=axes[axes_arr[count]])
    
    plt.show()


def get_whitelists(repo_path="."):
    """
    Fetches whitelists of the metadata_whitelist repository.
    
    Parameters
    ----------
    repo_path : str, default "."
        The path of the marker repository.
    """

    # Based on https://gitlab.gwdg.de/loosolab/software/metadata-organizer/-/blob/main/metaTools.py
    print('Fetching whitelists...\n')
    if not os.path.exists(f"{repo_path}/metadata_whitelists"):
        repo = git.Repo.clone_from('https://gitlab.gwdg.de/loosolab/software/metadata_whitelists.git/', f"{repo_path}/metadata_whitelists")
    else:
        repo = git.Repo('metadata_whitelists')
        o = repo.remotes.origin
        o.pull()
    print("Done!")


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


def select(whitelist=None, key=None, heading=None):
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

    Returns
    --------
    str :
        The selected string of the whitelist.
    """

    if not whitelist:
        if key:
            whitelist = read_whitelist(key)['whitelist']
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
        selection = whitelist[int(input())-1]
        
    print(f"Selection: {selection}\n")

    return selection


def get_gene_dict(organism=None, w_markers=None):
    """
    Creates dictionary of whitelist of genes of specific organism.

    Parameters
    ----------
    organism : str, default None
        The organism that owns the corresponding genes.
    w_markers : list of str, default None
        A list of Gene Symbols and Ensembl IDs separated by space.

    Returns
    --------
    dict :
        Dictionary which contains the gene names and ensembl IDs.
    """

    gene_dict = {}

    if organism:
        if len(organism.split(' ')) > 1:
            w_markers = read_whitelist(f"genes/{organism.split(' ')[0]}")['whitelist']
        else:
            w_markers = read_whitelist(f"genes/{organism}")['whitelist']
    elif not w_markers:
        raise ValueError("Provide organism or whitelist of genes (w_markers).")

    for marker in w_markers:
        name, ensg = marker.split(" ")[0].upper(), marker.split(" ")[1].upper()
        gene_dict[name] = ensg
        gene_dict[ensg] = name

    return gene_dict


def get_uid(path="./lists"):
    """
    Creates a new UID by iterating through all files in path.

    Parameters
    ----------
    path : str, default "./lists"
        The root of the files containing UIDs.

    Returns
    --------
    str :
        A string containing an unused UID
    """

    existing_uids = []
    
    # Iterate through all files and subdirectories in the given path
    for root, _, files in os.walk(path):
        for file in files:
            # Extract the UID from the filename
            file_parts = file.split('_')
            uid = file_parts[-1].split('.yaml')[0]
            existing_uids.append(uid)

    # Generate a new UID and make sure it is unique
    new_uid = 1
    while str(new_uid) in existing_uids:
        new_uid += 1

    return str(new_uid)


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
        The path of the marker repository.

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
        The path of the repository. Defaults to the current directory.
    """

    # Extract the list name from the list path
    list_name = os.path.splitext(os.path.basename(list_path))[0]  # Removes the .yaml extension

    repo = Repo(repo_path)
    assert not repo.bare

    # Pull the latest changes
    repo.remotes['origin'].pull()

    # Preprocess lists
    meta_path, markers_path = preprocess_lists_to_tsv(repo_lists_path=f"{repo_path}/lists")

    # Check out new branch
    repo.git.checkout('HEAD', b=list_name)
    repo.git.add(list_path)
    repo.git.add(meta_path)
    repo.git.add(markers_path)
    repo.git.commit('-m', f'Add new list: {list_name}')
    repo.git.push('--set-upstream', 'origin', list_name)


def preprocess_lists_to_tsv(repo_lists_path="./lists", output_path_meta="meta_lists.tsv", output_path_markers="marker_lists.tsv"):
    """
    Generates preprocessed .tsv files of the database and marker lists.
    This function can be used to speed up subsequent reads of the data.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
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
    df_meta = get_db(repo_lists_path)
    df_markers = get_marker_lists(repo_lists_path)
    
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