import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import git
from .utils import read_whitelist
import yaml
import string
from git import Repo
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor


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

    return df


def process_file_markers(file_path):
    """
    Reads and parses a marker list file to extract the marker list.
    
    Parameters
    ----------
    file_path : str
        Path of the marker list file (yaml-file).
    """
    
    with open(file_path, 'r', encoding='utf-8') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)
        marker_list_data = yaml_data.get("marker_list", [])

        # Extract markers and their names
        markers_data = []
        for item in marker_list_data:
            markers = item.get("markers", [])
            name = item.get("name", "")
            for marker in markers:
                markers_data.append({"marker": marker, "name": name})

        return markers_data


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


def search_db(df, keywords, case_sensitive=False, exact=False, out="metadata", repo_lists_path="./lists"):
    """
    This function filters a given DataFrame based on the provided keywords. Depending on the 'out' parameter,
    the function either returns the filtered DataFrame or a combined list of markers.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to be filtered.
    keywords : dict or str
        The keywords to filter the DataFrame.
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
    pd.DataFrame
        Either the filtered search results as metadata or as a combined list of markers.
    """

    if not case_sensitive:
        df = df.applymap(lambda x: str(x).lower() if isinstance(x, str) else x)
        if isinstance(keywords, dict):
            keywords = {k: v.lower() for k, v in keywords.items()}
        else:
            keywords = keywords.lower()

    if isinstance(keywords, dict):
        filtered_df = df.copy()
        for column, value in keywords.items():
            if exact:
                filtered_df = filtered_df[filtered_df[column] == value]
            else:
                filtered_df = filtered_df[filtered_df[column].str.contains(value, na=False, regex=False)]
    else:
        if exact:
            mask = df.applymap(lambda x: keywords == str(x)).any(axis=1)
        else:
            mask = df.applymap(lambda x: keywords in str(x)).any(axis=1)
        filtered_df = df[mask]

    if out == "marker_list":
        if repo_lists_path is None:
            raise ValueError("mr_path must be provided when out='marker_list'")
        uids = [int(idx) for idx in filtered_df.index]
        return combine_lists(uids, repo_lists_path=repo_lists_path)

    return filtered_df


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
        df = get_db(repo_lists_path=repo_lists_path)

    columns = df.columns.tolist()
    identifiers = [str(i) for i in range(1, 10)] + list(string.ascii_lowercase)[:len(columns)-9]

    # Split columns into groups of 10 for pagination
    page = 0
    pages = [columns[i:i+10] for i in range(0, len(columns), 10)]

    while True:
        # Print identifiers and column names for the current page
        print("Available columns for search:")
        for identifier, column in zip(identifiers, pages[page]):
            print(f"{identifier}: {column}")
        
        # Ask for column to search in
        command_options = []
        if page > 0: command_options.append("'p' for previous page")
        if page < len(pages) - 1: command_options.append("'n' for next page")
        command_options = ", ".join(command_options)
        col_to_search_identifier = input(f"Enter identifier of column to search in (leave blank to search in all columns)\nEnter {command_options}: ")
        
        if col_to_search_identifier == 'n':
            page = (page + 1) % len(pages)
            continue
        elif col_to_search_identifier == 'p':
            page = (page - 1) % len(pages)
            continue
        elif col_to_search_identifier == '':
            col_to_search = None  # Search in all columns
            break
        elif col_to_search_identifier in identifiers:
            col_to_search = pages[page][identifiers.index(col_to_search_identifier)]
            break

    # Ask for value to search for
    if col_to_search:
        show_possible_values = input("Do you want to see all possible values for this column? (yes/no): ").lower() == "yes"
        if show_possible_values:
            unique_values = df[col_to_search].unique()
            for value in unique_values:
                print(value)
    
    search_terms = input("Enter search terms (separate multiple terms with a comma): ").split(',')
    exact = input("Perform an exact search? (yes/no): ").lower() == "yes"
    case_sensitive = input("Consider case sensitivity? (yes/no): ").lower() == "yes"

    # Perform the search for each term and combine the results
    results = pd.DataFrame()
    for search_term in search_terms:
        if col_to_search:
            keywords = {col_to_search: search_term.strip()}
        else:
            keywords = search_term.strip()
        result = search_db(df, keywords, exact=exact, case_sensitive=case_sensitive)
        results = pd.concat([results, result])

    print(f"Number of results: {len(results)}")
    see_results = input("Do you want to see the results? (yes/no): ").lower() == "yes"

    if see_results:
        display(results)

    # Further filtering?
    further_filter = input("Do you want to filter the results further? (yes/no): ").lower() == "yes"

    if further_filter:
        return guided_search(repo_lists_path=repo_lists_path, df=results, out=out)
    
    if out == "marker_list":
        uids = [int(idx) for idx in results.index]
        return combine_lists(uids, repo_lists_path=repo_lists_path)
    
    return results


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


def export_marker_list(df, path=".", file_name="marker_list", header=False, marker_id=None):
    """
    Exports a marker list (df) to path/file_name. If a file with this name already exists,
    a timestamp suffix is added to the filename.

    Parameters
    ----------
    df : pd.DataFrame
        The marker list to be exported.
    path : str, default "."
        The path where the marker list should be saved.
    file_name : str, default "marker_list"
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

    # Check out new branch
    repo.git.checkout('HEAD', b=list_name)
    repo.git.add(list_path)
    repo.git.commit('-m', f'Add new list: {list_name}')
    repo.git.push('--set-upstream', 'origin', list_name)
