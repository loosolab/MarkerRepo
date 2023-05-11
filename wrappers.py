import marker_repo as mr
import os

def get_two_column_markers(keywords, REPO_LISTS_PATH, path=None, case_sensitive=False, exact=False):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.
    Optionally, it can export the DataFrame to a file.

    Parameters
    ----------
    keywords : dict or str
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    REPO_LISTS_PATH : string
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
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
    db = mr.get_db(REPO_LISTS_PATH)
    df = mr.search_db(db, keywords, case_sensitive=case_sensitive, exact=exact)
    if df.empty:
        raise Exception(
                    f"No search results available!")
    print(df)

    # Get UIDs and combine lists
    uids = [int(idx) for idx in df.index]
    combined_df = mr.combine_lists(REPO_LISTS_PATH, uids)

    # Drop duplicates
    markers_filtered = combined_df.drop_duplicates()

    if path:
        # Export marker list
        mr.export_marker_list(REPO_LISTS_PATH, path, markers_filtered)
        return os.path.abspath(path)
    else:
        return markers_filtered
