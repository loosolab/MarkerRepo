from .marker_repo import get_db, search_db, combine_lists, export_marker_list
from .calculations import compare_marker_lists, update_scores

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
    case_sensitive : bool, default False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.

    Returns
    --------
    pandas.DataFrame :
        The DataFrame conaining the combined lists.
    """

    db = get_db(REPO_LISTS_PATH)
    df = search_db(db, keywords, case_sensitive=case_sensitive, exact=exact)
    if df.empty:
        raise Exception(
                    f"No search results available!")
    # print(df)

    # Get UIDs and combine lists
    uids = [int(idx) for idx in df.index]
    combined_df = combine_lists(REPO_LISTS_PATH, uids)

    # Drop duplicates, keep one marker only, rearrange column order
    markers_filtered = combined_df.drop_duplicates()
    markers_filtered['Marker'] = markers_filtered['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x)
    markers_filtered = markers_filtered[['Info', 'Marker']]

    return markers_filtered


def convert_markers(REPO_LISTS_PATH, keywords=None, df=None, path=None, file_name="marker_list", case_sensitive=False, exact=False, style="two_column", organism="Hs", tissue="all"):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.
    Optionally, it can export the DataFrame to a file.

    Parameters
    ----------
    REPO_LISTS_PATH : str, default None
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    keywords : dict or str, default None
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    df : pd.DataFrame, default None
        A DataFrame containing a two column marker list.
    path : str, default None
        The path to the file where the new marker list will be saved. If not specified, the function will not save the DataFrame to a file.
    file_name : str, default "marker_list"
        The filename of the marker list.
    case_sensitive : bool, default False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.
    style : str, default "two_column"
        The format style of which the exported marker list should look like.
        Currently there are three options available: "two_column", "score" and "panglao"
    organism : str, default "Hs"
        Organism of panglao style markers.
    tissue : str, default "all"
        Tissue of panglao style markers.

    Returns
    --------
    pd.DataFrame or str :
        If path is specified, the function returns the absolute path to the file where the marker list was saved.
        If path is not specified, the function returns the DataFrame.
    """

    if REPO_LISTS_PATH and keywords:
        marker_list = get_selected_lists(REPO_LISTS_PATH, keywords, case_sensitive=case_sensitive, exact=exact)
    else:
        marker_list = df
    
    match style:
        case "two_column":
            print("Preparing two column style marker list...")
            
        case "score":
            print("Preparing score style marker list...")
            marker_list = compare_marker_lists(marker_df=marker_list)
            
        case "panglao":
            print("Preparing panglao style marker list...")
            marker_list = compare_marker_lists(marker_df=marker_list)
            marker_list = update_scores(marker_list)
            marker_list = transform_list_to_panglao(df=marker_list, organism=organism, tissue=tissue)
        case _:
            print("Style not recognized. Try 'two_column', 'score' or 'panglao'")


    if path or file_name:
        # Export marker list
        return export_marker_list(marker_list, path=path, file_name=file_name)
    else:
        return marker_list


def transform_list_to_panglao(df, organism="Hs", tissue="all"):
    """
    Adds additional columns to the given DataFrame and rearranges column order.
    This converts the format of the marker list to the format of the Panglao db.
    
    Parameters
    ----------
    df : DataFrame
        Input DataFrame with columns "Info", "Marker", and "Score".
    organism : str, default "Hs"
        Organism of panglao style markers.
    tissue : str, default "all"
        Tissue of panglao style markers.
        
    Returns
    -------
    pd.DataFrame :
        Updated DataFrame with new columns, specified column order, and updated Marker column.
    """

    # Update Marker column to retain only the first marker
    df['Marker'] = df['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x)
    
    # Add new columns
    df['Organism'] = organism
    df['Aliases'] = df['Marker']
    df['Tissue'] = tissue
    
    # Rearrange column order
    df = df[['Organism', 'Marker', 'Info', 'Aliases', 'Score', 'Tissue']]
    
    return df
