import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def checkFiles(LIST_PATH, METADATA_PATH, list_type):
    """
    Checks the input files: does the list and - if not None - the metadata file exist? 
    Checks the format of the list: correct amount of columns? Separation correct?

    Parameters
    ----------
    LIST_PATH : string
        The path where the list is stored.
    METADATA_PATH : string
        The path where the metadata file of the list is stored.
    list_type : string
        The type of the list (gene/region).

    Returns
    --------
    boolean :
        True if the file(s) exist and the input format seems to be correct, False else 
    """

    # Checking path(s) of input files
    files = [LIST_PATH]
    if METADATA_PATH:
        files.append(METADATA_PATH)
    for file in files:
        if os.path.isfile(file):
            print(f"{file} exists.")
        else:
            print(f"Please make sure that your input is correct. {file} does not exist.")
            return False
     
    # Checking format of input list
    correct = True
    with open(f"{LIST_PATH}", "r") as list_file:
        lines = list_file.readlines()
        if list_type == "celltype" or "cellcycle":
            for line in lines:
                if len(line.split("\t")) != 2:
                    print("Please make sure that your input file consists of two columns, separated by tabs.")
                    if list_type == "celltype":
                        print(f"First column: cell type\nSecond column: gene")
                    if list_type == "cellcycle":
                        print(f"First column: cellcycle gene\nSecond column: phase")
                    correct = False
                    break
        elif list_type == "mito" or "gender":
            for line in lines:
                if len(line.split("\t")) != 1:
                    print("Please make sure that your input file consists of one column.")
                    print("This column should contain gene names.")
                    correct = False
                    break
        # TODO: blacklist
        if correct:
            print("The format of the list seems correct.")
            return True
        else:
            return False


def addList(REPO_LISTS_PATH, LIST_PATH, metadata):
    """
    Adds a new list to the Marker Gene Repo.

    Parameters
    ----------
    REPO_LISTS_PATH : string
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    LIST_PATH : string
        The path where the list is stored.
    metadata : dictionary
        The dictionary containing the metadata information.

    Returns
    --------
    string :
        The path of the newly added list
    """

    # create folder(s) and get path of new list
    folder = f"{REPO_LISTS_PATH}/{metadata['Kind']}/{metadata['Organism']}/{metadata['Tissue']}/{metadata['Year']}/{metadata['List type']}"
    new_file_path = f"{folder}/{metadata['Title']}"

    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder {folder} created.")

    # check whether list already exists
    if os.path.isfile(new_file_path):
        go_on = input(f"The file {new_file_path} already exists.\nDo you want to override the existing file?\
        Enter yes or no: ")
        override = True if go_on == "yes" else False
        # override existing file
        if override:
            shutil.copyfile(LIST_PATH, f"{folder}/{metadata['Title']}")
            print(f"Replaced list in folder {folder}.")
    else:
        # copy list to Marker Repo
        shutil.copyfile(LIST_PATH, f"{folder}/{metadata['Title']}")
        print(f"Copied list to {folder}.")

    return f"{folder}/{metadata['Title']}"


def searchDB(df, keywords, exact=False):
    """
    Search for specific lists of the Marker Gene Repo.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe which contains the data to be searched.
    keywords : dictionary
        The dictionary containing the keywords for filtering.
    exact : boolean, default False
        If True perform an exact search, otherwise use "contains".
    
    Returns
    --------
    pandas.DataFrame :
        Dataframe containing all the hits
    """

    # keep values which are not None
    filters = {}
    for key in keywords:
        if keywords[key]:
            filters[key] = keywords[key]

    # filter dataframe
    for key in filters:
        if exact:
            df = df[df[key] == filters[key]]
        else:
            df = df[df[key].str.contains(filters[key])]
    
    return df.reset_index(drop=True)


def getDB(REPO_LISTS_PATH):
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

    # get all paths of lists of the Marker Repo
    files = [os.path.join(root, name) for root, dirs, files in os.walk(REPO_LISTS_PATH) for name in files]
    kinds, organisms, tissues, years, ltypes, titles = ([] for i in range(6))

    # create dataframe of paths
    for file in files:
        file = file.split("lists/")[1]
        kind, organism, tissue, year, ltype, title = file.split("/")
        kinds.append(kind)
        organisms.append(organism)
        tissues.append(tissue)
        years.append(year)
        ltypes.append(ltype)
        titles.append(title)

    list_dict = {"Kind": kinds, "Organism": organisms, "Tissue": tissues, 
                "Year": years, "List type": ltypes, "Title": titles}

    df = pd.DataFrame(list_dict)

    return df


def getPaths(REPO_LISTS_PATH, df):
    """
    Converts the dataframe of the Marker Repo DB to paths.

    Parameters
    ----------
    REPO_LISTS_PATH : string
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    df : pandas.DataFrame
        The dataframe containing lists of the Marker Repo.
    Returns
    --------
    array of strings :
        The paths of the lists inside the dataframe.
    """

    # read dataframe, convert rows to strings
    paths = []
    rows = df.to_string(header=False, index=False, index_names=False).split('\n')

    # convert strings to paths
    files = ['/'.join(row.split()) for row in rows]
    for file in files:
        path = f"{REPO_LISTS_PATH}/{file}"
        paths.append(path)
        print(path)
    
    return paths


def getList(path, list_type):
    """
    Checks the input files: does the list and - if not None - the metadata file exist? 
    Checks the format of the list: correct amount of columns? Separation correct?

    Parameters
    ----------
    path : string
        The path where the list is stored.
    list_type : string
        The type of the list (gene/region).

    Returns
    --------
    pandas.DataFrame :
        Dataframe containing the list
    """

    # create dictionary containing headers of different list types
    # TODO use whitelist instead
    ltype_dict = {"celltype": ["Cell type", "Marker"], "cellcycle": ["Marker", "Phase"], "mito": "Marker", "gender": "Marker", "blacklist": ["Chr", "Start", "Stop"]}
    header = ltype_dict[list_type]

    df = pd.read_csv(path, sep='\t', names=header)

    return df


def combineLists(paths, list_type, file_name="custom_list"):
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


def showStatistics(REPO_LISTS_PATH, metadata, dpi=120):
    """
    Shows content of whole Marker Repo.

    Parameters
    ----------
    REPO_LISTS_PATH : string
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    metadata : dictionary
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


def convertList():
    pass
