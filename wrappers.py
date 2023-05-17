import marker_repo as mr
import calculations as calc
import os

def get_selected_lists(REPO_LISTS_PATH, keywords, case_sensitive=False, exact=False):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.

    Parameters
    ----------
    REPO_LISTS_PATH : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    keywords : dict or str
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    case_sensitive : bool, default: False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default: False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.

    Returns
    --------
    pandas.DataFrame :
        The DataFrame conaining the combined lists.
    """

    db = mr.get_db(REPO_LISTS_PATH)
    df = mr.search_db(db, keywords, case_sensitive=case_sensitive, exact=exact)
    if df.empty:
        raise Exception(
                    f"No search results available!")
    print(df)

    # Get UIDs and combine lists
    uids = [int(idx) for idx in df.index]
    combined_df = mr.combine_lists(REPO_LISTS_PATH, uids)

    # Drop duplicates, keep one marker only, rearrange column order
    markers_filtered = combined_df.drop_duplicates()
    markers_filtered['Marker'] = markers_filtered['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x)
    markers_filtered = markers_filtered[['Info', 'Marker']]


def get_two_column_markers(REPO_LISTS_PATH, keywords, path=None, case_sensitive=False, exact=False):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.
    Optionally, it can export the DataFrame to a file.

    Parameters
    ----------
    REPO_LISTS_PATH : str
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    keywords : dict or str
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    path : str, default: None
        The path to the file where the new marker list will be saved. If not specified, the function will not save the DataFrame to a file.
    case_sensitive : bool, default: False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default: False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.

    Returns
    --------
    pandas.DataFrame or str :
        If path is specified, the function returns the absolute path to the file where the marker list was saved.
        If path is not specified, the function returns the DataFrame.
    """

    df = get_selected_lists(REPO_LISTS_PATH, keywords, case_sensitive=case_sensitive, exact=exact)

    if path:
        # Export marker list
        mr.export_marker_list(REPO_LISTS_PATH, path, df)
        return os.path.abspath(path)
    else:
        return df


def transform_list_to_panglao(df, organism="Hs", tissue="all"):
    """
    Adds additional columns to the given DataFrame and rearranges column order.
    This converts the format of the marker list to the format of the Panglao db.
    
    Parameters
    ----------
    df : DataFrame
        Input DataFrame with columns "Info", "Marker", and "Score".
        
    Returns
    -------
    pd.DataFrame :
        Updated DataFrame with new columns, specified column order, and updated Marker column.
    """

    # Update Marker column to retain only the first marker
    df['Marker'] = df['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x)
    
    # Add new columns
    df['Organism'] = organism
    df['Nicknames'] = df['Marker']
    df['Tissue'] = tissue
    
    # Rearrange column order
    df = df[['Organism', 'Marker', 'Info', 'Nicknames', 'Score', 'Tissue']]
    
    return df


def get_panglao_style_markers():
    #TODO
    pass