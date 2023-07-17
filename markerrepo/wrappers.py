from .marker_repo import search_df, combine_lists, export_marker_list, guided_search, combine_dfs
from .calculations import compare_marker_lists, update_scores

def get_selected_lists(keywords=None, metadata_df=None, repo_path=".", case_sensitive=False, exact=False):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.

    Parameters
    ----------
    keywords : dict or str, default None
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    repo_path : str, default "."
        The path of the Marker Repo.
    metadata_df : pd.DataFrame, default None
        A DataFrame containing the metadata of a selection of marker lists.
    case_sensitive : bool, default False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.

    Returns
    --------
    pd.DataFrame :
        The DataFrame conaining the combined lists.
    """

    if keywords:
        db = combine_dfs(repo_path=repo_path)
        df = search_df(db, keywords, case_sensitive=case_sensitive, exact=exact)
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
    combined_df = combine_lists(uids, repo_path=repo_path)

    # Drop duplicates, keep one marker only, rearrange column order
    markers_filtered = combined_df.drop_duplicates()
    markers_filtered = markers_filtered[['Info', 'Marker']]

    return markers_filtered


def convert_markers(repo_path=".", keywords=None, df=None, path="exported_lists", file_name="marker_list", case_sensitive=False, exact=False, style="two_column", organism="Hs", tissue="all", gs=False, ensembl=False):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.
    Optionally, it can export the DataFrame to a file.

    Parameters
    ----------
    repo_lists_path : str, default "."
        The path of the Marker Repo.
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
    gs : bool, default False
        If true, guided search is enabled and the resulting DataFrame will be used.
    ensembl : bool, default False
        If True, the Ensembl IDs will be used instead of the gene symbols.

    Returns
    --------
    pd.DataFrame or str :
        If path is specified, the function returns the absolute path to the file where the marker list was saved.
        If path is not specified, the function returns the DataFrame.
    """

    if gs:
        marker_list = guided_search(repo_path=repo_path, out="marker_list")
    elif repo_path and keywords:
        marker_list = get_selected_lists(keywords, repo_path=repo_path, case_sensitive=case_sensitive, exact=exact)
    elif df is not None:
        marker_list = df
    else:
        raise Exception(
            f"You need to specify keywords, a marker list DataFrame or use the guided search ('gs=True') !")

    if ensembl:
        marker_list['Marker'] = marker_list['Marker'].apply(lambda x: x.split(' ')[1] if len(x.split(' ')) > 1 else x)
    else:
        marker_list['Marker'] = marker_list['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x)
    
    match style:
        case "two_column":
            print("Preparing two column style marker list...")
            
        case "score":
            print("Preparing score style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list)
            
        case "panglao":
            print("Preparing panglao style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list)
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
