import os
import sys

def get_all_uids(repo_lists_path="./lists"):
    """ 
    Traverse the specified directory and its subdirectories to extract UIDs from all marker lists.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.

    Returns
    --------
    list :
        A list of all UIDs.
    """

    uids = []
    for _, _, files in os.walk(repo_lists_path):
        for file in files:
            if '.yaml' in file:
                uid = file.split('_')[-1].split('.yaml')[0]
                uids.append(uid)

    return uids


def check_uid(repo_lists_path="./lists"):
    """
    Check whether the new UID is unique.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.

    Returns
    --------
    bool :
        True if the new UID is unique, False otherwise.
    """

    # Get the new UID from the environment variable
    new_uid = os.getenv('NEW_UID')
    uids = get_all_uids(repo_lists_path=repo_lists_path)
    
    if new_uid in uids:
        print(f"Duplicate UID found: {new_uid}")
        return False
    else:
        print(f"No duplicate UID found.")
        return True

# Run the check and exit with error if the UID is not unique
if not check_uid():
    raise ValueError(f'UID {os.getenv("NEW_UID")} already exists.')
