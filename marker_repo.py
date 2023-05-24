import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import git
import src.utils as utils
import yaml
import string

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


def search_db(df, keywords, case_sensitive=False, exact=False):
    """
    This function filters a given DataFrame based on the provided keywords. The keywords can be either
    a dictionary or a string. If the keywords are provided as a dictionary, the DataFrame will be filtered
    using the dictionary keys as column names and the corresponding values as the keywords to search within
    those columns. If the keywords are provided as a string, the function will check if the DataFrame contains
    the string anywhere and retains only the rows that fulfill the search criteria. The search can be made
    either exact or partial (substring) and case-sensitive or case-insensitive.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to be filtered.
    keywords : dict or str
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    case_sensitive : bool, default: False
        If True, the search will be case-sensitive. If False, the search will be case-insensitive.
    exact : bool, default: False
        If True, the search will look for exact matches. If False, the search will look for substrings.

    Returns
    -------
    pd.DataFrame :
        The filtered DataFrame containing only the rows that meet the search criteria.
    """

    if not case_sensitive:
        df = df.applymap(lambda x: str(x).lower() if isinstance(x, str) else x)
        if isinstance(keywords, dict):
            keywords = {k: v.lower() for k, v in keywords.items()}
        else:
            keywords = keywords.lower()

    if isinstance(keywords, dict):
        # If keywords is a dictionary, filter by matching column names and values
        filtered_df = df.copy()
        for column, value in keywords.items():
            if exact:
                filtered_df = filtered_df[filtered_df[column] == value]
            else:
                filtered_df = filtered_df[filtered_df[column].str.contains(value, na=False, regex=False)]
    else:
        # If keywords is a string, filter by checking if the string exists anywhere in the DataFrame
        if exact:
            mask = df.applymap(lambda x: keywords == str(x)).any(axis=1)
        else:
            mask = df.applymap(lambda x: keywords in str(x)).any(axis=1)
        filtered_df = df[mask]

    return filtered_df


def guided_search(REPO_LISTS_PATH, df=None):
    """
    An interactive function that guides the user through the process of searching the DataFrame.

    Parameters
    ----------
    REPO_LISTS_PATH : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    df : pd.DataFrame, default: None
        The DataFrame to search in. If not provided, the function will create one from the REPO_LISTS_PATH.

    Returns
    -------
    pd.DataFrame :
        The DataFrame that contains the search results.
    """

    # Get the DataFrame if not provided
    if df is None:
        df = get_db(REPO_LISTS_PATH)

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
        col_to_search_identifier = input("Enter identifier of column to search in (leave blank to search in all columns)\nEnter 'n' for next page, 'p' for previous page: ")
        
        if col_to_search_identifier == 'n':
            page = (page + 1) % len(pages)
            continue
        elif col_to_search_identifier == 'p':
            page = (page - 1) % len(pages)
            continue

        if col_to_search_identifier in identifiers:
            col_to_search = pages[page][identifiers.index(col_to_search_identifier)]
            break
    
    # Ask for value to search for
    search_term = input("Enter search term: ")
    exact = input("Perform an exact search? (yes/no): ").lower() == "yes"
    case_sensitive = input("Consider case sensitivity? (yes/no): ").lower() == "yes"
    
    # Build the keywords
    if col_to_search:
        keywords = {col_to_search: search_term}
    else:
        keywords = search_term

    # Perform the search
    results = search_db(df, keywords, exact=exact, case_sensitive=case_sensitive)

    print(f"Number of results: {len(results)}")
    see_results = input("Do you want to see the results? (yes/no): ").lower() == "yes"

    if see_results:
        display(results)

    # Further filtering?
    further_filter = input("Do you want to filter the results further? (yes/no): ").lower() == "yes"

    if further_filter:
        return guided_search(REPO_LISTS_PATH, results)
    
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
            # if the value is a dictionary, recursively flatten it
            items.extend(flatten_dict(v, new_key, sep=sep, list_sep=list_sep).items())
        elif isinstance(v, list):
            # if the value is a list, process its elements
            list_items = []
            for i, elem in enumerate(v):
                if isinstance(elem, dict):
                    # if the element is a dictionary, join its key-value pairs with the list separator
                    list_items.append(list_sep.join(f"{key}: {value}" for key, value in elem.items()))
                else:
                    # otherwise, convert the element to a string
                    list_items.append(str(elem))
            # join the list items with the list separator and store them in a single cell
            items.append((new_key, list_sep.join(list_items)))
        else:
            # if the value is not a dictionary or a list, store it directly
            items.append((new_key, v))
            
    return dict(items)


def get_db(REPO_LISTS_PATH):
    """
    Get the database of the Marker Repo as DataFrame.

    Parameters
    ----------
    REPO_LISTS_PATH : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.

    Returns
    --------
    pandas.DataFrame :
        Dataframe containing all lists
    """

    data = []

    # iterate through all files in the folder and subfolders
    for root, dirs, files in os.walk(REPO_LISTS_PATH):
        for file in files:
            if file.endswith(".yaml"):
                file_path = os.path.join(root, file)

                # read YAML file and extract all leaf values from "metadata"
                with open(file_path, 'r', encoding='utf-8') as yaml_file:
                    yaml_data = yaml.safe_load(yaml_file)
                    metadata = yaml_data.get("metadata", {})

                    # flatten the dictionary and save the results
                    flattened_metadata = flatten_dict(metadata)
                    data.append(flattened_metadata)

    # create DataFrame and set "id" as index
    df = pd.DataFrame(data)
    df.rename(columns=get_display_names(), inplace=True)
    if "ID" in df.columns:
        df.set_index("ID", inplace=True)

    return df


def get_list(path, info_col=1, marker_col=0):
    """
    Reads the marker lists and converts it to a DataFrame using the information
    of info_col and marker_col.

    Parameters
    ----------
    path : str
        The path where the list is stored.
    info_col : integer, default 1
        The column which contains additional information like cell type or phase.
    marker_col: integer, default 0
        The column which contains the marker (gene or genomic region).

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing the list
    """

    # ltype_dict = {"celltype": ["Cell type", "Marker"], "cellcycle": ["Marker", "Phase"], "mito": "Marker", "gender": "Marker", "blacklist": ["Chr", "Start", "Stop"]}
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


def export_marker_list(path, file_name, df):
    """
    Export marker list (df) to path/file_name

    Parameters
    ----------
    path : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    file_name : str
        The file name of the combined list.
    df : pd.DataFrame
        The marker list which will be exported.
    """

    export_path = f"{path}/{file_name}"

    df.to_csv(export_path, sep="\t", index=False)
    print(f"Combined list saved: {export_path}")


def get_uid_paths(REPO_LISTS_PATH, uids):
    """
    Searches for files in the specified folder and its subfolders with names in the format "name_UID.yaml",
    where UID is an integer. Returns the paths of the files that contain the UIDs from the given list.

    Parameters
    ----------
    REPO_LISTS_PATH : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    uids : list of int
        A list of integers representing the UIDs to search for.

    Returns
    -------
    list of str :
        A list of file paths containing the specified UIDs.
    """

    matching_files = []

    for root, _, files in os.walk(REPO_LISTS_PATH):
        for file in files:
            if file.endswith('.yaml'):
                uid = int(file.split('_')[-1].split('.')[0])
                if uid in uids:
                    matching_files.append(os.path.join(root, file))

    return matching_files


def combine_lists(REPO_LISTS_PATH, uids):
    """
    Combine multiple lists to one custom list.

    Parameters
    ----------
    REPO_LISTS_PATH : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    uids : list of str
        The uids of the lists which will be combined.

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing the combinend list
    """

    # read lists which are going to be combined
    dfs = []
    for file in get_uid_paths(REPO_LISTS_PATH, uids):
        dfs.append(get_marker_list(file))
        
    # perform outer join
    combined_df = pd.concat(dfs).reset_index(drop=True)

    # TODO: inner join, etc ...

    return combined_df


def show_statistics(REPO_LISTS_PATH, metadata, dpi=120):
    """
    Shows content of whole Marker Repo.

    Parameters
    ----------
    REPO_LISTS_PATH : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    metadata : dict
        The dictionary containing the metadata information.
    dpi : integer, default 120
    """

    # TODO tissue plot - show count only

    sns.set_style("darkgrid")
    sns.set(rc={"figure.dpi": dpi, "savefig.dpi": dpi})
    fig, axes = plt.subplots(2, 3)
    axes_arr = [(0,0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]

    # Load all lists
    df = getDB(REPO_LISTS_PATH)

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


def get_whitelists():
    """
    Fetches whitelists of the metadata_whitelist repository.
    """

    # Based on https://gitlab.gwdg.de/loosolab/software/metadata-organizer/-/blob/main/metaTools.py
    print('Fetching whitelists...\n')
    if not os.path.exists('metadata_whitelists'):
        repo = git.Repo.clone_from('https://gitlab.gwdg.de/loosolab/software/metadata_whitelists.git/', 'metadata_whitelists')
    else:
        repo = git.Repo('metadata_whitelists')
        o = repo.remotes.origin
        o.pull()
    print("Done!")


def update_markers(df, marker_dict):
    """
    Updates markers by extending gene names withi ensembl IDs and the other way round.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame containing the marker list.
    marker_dict : dict
        Dictionary containing the names and IDs as keys and values.

    Returns
    --------
    dict :
        The dictionary containing markers and corresponding information
    """

    def apply_update(marker):
        return marker + ' ' + marker_dict[marker] if marker in marker_dict else marker

    df['Marker'] = df['Marker'].apply(apply_update)

    return df


def select(whitelist=None, key=None, heading=None):
    """
    Shows selection of whitelist and returns selected value.
    If only a key is passed, the corresponding whitelist is used as a selection.
    If a whitelist (list of strings) is passed, the key is the heading only.

    Parameters
    ----------
    whitelist : list of str, default: None
        The selection to choose from.
    key : str, default: None
        The key of the whitelist. For example "organism".
    heading: str, default: None
        The heading (description) of the whitelist.

    Returns
    --------
    str :
        The selected string of the whitelist.
    """
    if not whitelist:
        whitelist = utils.read_whitelist(key)['whitelist']
    
    if not heading:
        print(f"Select {key}")
    else:
        print(f"Select {heading}")

    for i, value in enumerate(whitelist):
        print(str(i+1) + ":\t" + value)
    
    selection = whitelist[int(input())-1]
    print(f"Selection: {selection}\n")

    return selection


def get_gene_dict(organism):
    """
    Creates dictionary of whitelist of genes of specific organism.

    Parameters
    ----------
    organism : str
        The organism that owns the corresponding genes.

    Returns
    --------
    dict :
        Dictionary which contains the gene names and ensembl IDs.
    """
    gene_dict = {}
    w_markers = utils.read_whitelist(f"genes/{organism.split(' ')[0]}")['whitelist']
    for marker in w_markers:
        name, ensg = marker.split(" ")[0].upper(), marker.split(" ")[1].upper()
        gene_dict[name] = ensg
        gene_dict[ensg] = name

    return gene_dict


def get_uid(path):
    """
    Creates a new UID by iterating through all files in path.

    Parameters
    ----------
    path : str
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
            # extract display name of the current key if exists
            display_name = v.get('display_name')
            if display_name:
                display_names[new_key] = display_name
            # recursively extract display names from the nested dictionary
            display_names.update(extract_display_names(v, new_key, sep=sep))

    return display_names


def get_display_names():
    """
    Extracts the display names of the keys.yaml.

    Returns:
    -------
    dict :
        A dictionary containing the display names as values and concatenated names as keys.
    """

    # load the YAML file
    with open('keys.yaml', 'r', encoding='utf-8') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)

    # extract the display names
    display_names = extract_display_names(yaml_data['metadata'])

    return display_names
