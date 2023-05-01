import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import git
import src.utils as utils
import yaml


def read_marker_list(file_path):
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


def search_db(df, keywords, exact=False, case_sensitive=False):
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

    exact : bool, default: False
        If True, the search will look for exact matches. If False, the search will look for substrings.

    case_sensitive : bool, default: False
        If True, the search will be case-sensitive. If False, the search will be case-insensitive.

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
    
    #return df.reset_index(drop=True)


def flatten_dict(d, parent_key='', sep='_', list_sep='\n'):
    """
    Flatten a nested dictionary, concatenating keys with a separator.

    Parameters
    ----------
    d : dict
        The input dictionary to be flattened.
    parent_key : string, default ''
        The parent key used during recursion (default is an empty string).
    sep : string, default '_'
        The separator used to concatenate keys (default is an underscore).
    list_sep : string, default '\n'
        The separator used to join list elements in a single cell (default is a newline character).

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
    Get the database of the Marker Repo as dataframe.

    Parameters
    ----------
    REPO_LISTS_PATH : string
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
    if "id" in df.columns:
        df.set_index("id", inplace=True)
    return df


def get_list(path, info_col=1, marker_col=0):
    """
    Reads the marker lists and converts it to a dataframe using the information
    of info_col and marker_col.

    Parameters
    ----------
    path : string
        The path where the list is stored.
    info_col : integer
        The column which contains additional information like cell type or phase.
    marker_col: integer
        The column which contains the marker (gene or genomic region).

    Returns
    --------
    pandas.DataFrame :
        Dataframe containing the list
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


def combine_lists(paths, list_type, file_name="custom_list"):
    """
    Combine multiple lists to one custom list.

    Parameters
    ----------
    paths : array of strings
        The paths of the lists which will be combined.
    file_name : string, default "custom_list"
        The file name of the combined list.
    Returns
    --------
    pandas.DataFrame :
        Dataframe containing the combinend list
    """

    # read lists which are going to be combined
    dfs = []
    for path in paths:
        dfs.append(getList(path, list_type))
        
    # perform outer join
    combined_df = pd.concat(dfs).reset_index(drop=True)
    display(combined_df)

    # TODO: inner join, etc ...

    # save custom list
    combined_df.to_csv(file_name, sep="\t", index=False)
    print(f"Combined list saved: {os.path.abspath(file_name)}")

    return combined_df


def show_statistics(REPO_LISTS_PATH, metadata, dpi=120):
    """
    Shows content of whole Marker Repo.

    Parameters
    ----------
    REPO_LISTS_PATH : string
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


def dataframe_to_dict(df, info_col=0, marker_col=1):
    """
    Converts dataframe of marker list to dictionary,
    using info_col as keys and marker_col as values.

    # TODO only marker column available

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe containing the marker list.
    info_col : integer
        The column which contains additional information like cell type or phase.
    marker_col: integer
        The column which contains the marker (gene or genomic region).
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


def update_markers(df, marker_dict):
    """
    Updates markers by extending gene names withi ensembl IDs and the other way round.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe containing the marker list.
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


def select(key):
    """
    Shows selection of whitelist and returns selected value.

    Parameters
    ----------
    key : string
        The key of the whitelist. For example "organism".

    Returns
    --------
    string :
        The selection of the whitelist.
    """
    whitelist = utils.read_whitelist(key)['whitelist']
    print(f"Select {key}")
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
    organism : string
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