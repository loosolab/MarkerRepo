import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
import yaml


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


def get_marker_lists(repo_path=".", lists_path=None, parallel=True):
    """
    Get the marker list from the Marker Repo as DataFrame.

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
        DataFrame containing all markers and their designations.
    """

    if not lists_path:
        lists_path = f"{repo_path}/lists"

    # Check if the provided repository path exists
    if not os.path.exists(lists_path):
        raise FileNotFoundError(f"The specified path '{lists_path}' does not exist.")

    file_paths = [os.path.join(root, file) for root, dirs, files in os.walk(lists_path) for file in files if file.endswith(".yaml")]

    # Check if there are any YAML files in the path
    if not file_paths:
        raise FileNotFoundError(f"No YAML files found in the path '{lists_path}'.")

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


def get_list(path, info_col=1, marker_col=0, marker_type="Genes"):
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
    marker_type : str, default "Genes"
        The type of the marker ("Genes" or "Genomic regions").

    Returns
    --------
    pandas.DataFrame :
        DataFrame containing the list
    """

    if marker_type == "Genes":
        # Detect the number of columns
        with open(path, 'r') as f:
            line = f.readline()
            num_columns = len(line.split('\t'))

        headers = ["Marker", "Info"]
        if type(info_col) == int and num_columns > 2:
            df = pd.read_csv(path, sep='\t', usecols=[marker_col, info_col], names=headers)
        elif type(info_col) == int and num_columns <= 2:
            df = pd.read_csv(path, sep='\t', names=headers)
        else:
            df = pd.read_csv(path, sep='\t', names=[headers[0]])
            df[headers[1]] = info_col

    elif marker_type == "Genomic regions":
        rows = []

        with open(path, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                # Assuming the format is CHR START STOP INFO
                if info_col == 1 and marker_col == 0:
                    # marker = f"{parts[0]}:{parts[1]}-{parts[2]}"
                    marker = parts[0]
                    info = parts[-1]
                # Assuming the format is INFO CHR START STOP
                elif info_col == 0 and marker_col == 1:
                    # marker = f"{parts[1]}:{parts[2]}-{parts[3]}"
                    marker = parts[1]
                    info = parts[0]
                else:
                    raise ValueError("Invalid column configuration. Please set info_col and marker_col correctly.")
                
                rows.append({"Marker": marker, "Info": info}) 
        
        df = pd.DataFrame(rows, columns=["Marker", "Info"]) 

    else:
        raise ValueError(f"Invalid marker type '{marker_type}'. Please choose 'Genes' or 'Genomic regions'.")

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


def filter_strings_by_prefix(strings, valid_prefixes = ["leiden", "louvain", "lovain", "kmeans"]):
    """
    Filters a list of strings by a list of prefixes.

    Parameters
    ----------
    strings : list of str
        The list of strings that should be filtered.
    valid_prefixes : list of str, default ["leiden", "louvain", "lovain", "kmeans"]
        The list of prefixes that should be used for filtering.

    Returns
    -------
    list of str :
        The list of strings that start with one of the valid prefixes.
    """

    filtered_strings = []

    for string in strings:
        if any(string.lower().startswith(prefix) for prefix in valid_prefixes):
            filtered_strings.append(string)

    return filtered_strings

