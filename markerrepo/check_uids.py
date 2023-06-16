import os
import sys

def get_all_uids(repo_lists_path="./lists"):
    """ 
    Traverse the specified directory and its subdirectories to extract UIDs from filenames.
    The filenames are expected to be in the format '<name>_UID.<extension>', where '<name>' can contain underscores.

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

    print(uids)

    return uids

def check_uid(new_uid, repo_lists_path="./lists"):
    """
    Check whether the new UID is unique.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    new_uid : str
        The UID to be checked for uniqueness.
    """

    uids = get_all_uids(repo_lists_path=repo_lists_path)
    if new_uid in uids:
        print(f"Duplicate UID found: {new_uid}")
        sys.exit(1)
    else:
        print(f"No duplicate UID found.")
