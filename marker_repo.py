import os
import shutil
import pandas as pd


def check_input_files(LIST_PATH, METADATA_PATH, list_type):
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

    folder = f"{REPO_LISTS_PATH}/{metadata['Kind']}/{metadata['Organism']}/{metadata['Tissue']}/{metadata['Year']}/{metadata['List type']}"
    new_file_path = f"{folder}/{metadata['Title']}"

    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Folder {folder} created.")

    if os.path.isfile(new_file_path):
        go_on = input(f"The file {new_file_path} already exists. Do you want to override the existing file?\
        Enter yes or no: ")
        override = True if go_on == "yes" else False
        if override:
            shutil.copyfile(LIST_PATH, f"{folder}/{metadata['Title']}")
            print(f"Replaced list in folder {folder}.")
    else:
        shutil.copyfile(LIST_PATH, f"{folder}/{metadata['Title']}")
        print(f"Copied list to {folder}.")

    return(f"{folder}/{metadata['Title']}")

def searchDB(REPO_LISTS_PATH, keywords):
    """
    Search for specific lists of the Marker Gene Repo.

    Parameters
    ----------
    REPO_LISTS_PATH : string
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    keywords : dictionary
        The dictionary containing the keywords for filtering.

    Returns
    --------
    pandas.DataFrame :
        Dataframe containing all the hits
    """
    hits = []
    filtered_files = [os.path.join(root, name) for root, dirs, files in os.walk(REPO_LISTS_PATH) for name in files]

    # Filter title
    if keywords['Title']:
        filtered_files = [file for file in filtered_files if keywords['Title'] in file.split("/")[-1]]
        
    # Filter list_type
    if keywords['List type']:
        filtered_files = [file for file in filtered_files if keywords['List type'] in file.split("/")[-2]]
        
    # Filter year
    if keywords['Year']:
        filtered_files = [file for file in filtered_files if keywords['Year'] in file.split("/")[-3]]
        
    # Filter tissue
    if keywords['Tissue']:
        filtered_files = [file for file in filtered_files if keywords['Tissue'] in file.split("/")[-4]]
        
    # Filter organism
    if keywords['Organism']:
        filtered_files = [file for file in filtered_files if keywords['Organism'] in file.split("/")[-5]]
        
    # Filter kind
    if keywords['Kind']:
        filtered_files = [file for file in filtered_files if keywords['Kind'] in file.split("/")[-6]]

    count = len(filtered_files)
    if count > 0:
        print(f"{str(count)} lists have been found:")
    else:
        print("Unfortunately there was no hit. Maybe try fewer keywords.")
        
    for file in filtered_files:
        hit = file.split("lists/")[1]
        hits.append(hit)
        print(hit)

    kinds, organisms, tissues, years, ltypes, titles = ([] for i in range(6))
    for hit in hits:
        kind, organism, tissue, year, ltype, title = hit.split("/")
        kinds.append(kind)
        organisms.append(organism)
        tissues.append(tissue)
        years.append(year)
        ltypes.append(ltype)
        titles.append(title)
        
    list_dict = {"Kind": kinds, "Organism": organisms, "Tissue": tissues, 
                "Year": years, "List type": ltypes, "Title": titles}

    df = pd.DataFrame(list_dict)
    
    return(df)


def convertList():
    pass


def mergeLists():
    pass

