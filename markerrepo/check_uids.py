import argparse
import os
from collections import Counter

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("uid", help="The UID to be checked for uniqueness.")
args = parser.parse_args()
new_uid = args.uid

def get_all_uids(repo_path="."):
    """ 
    Traverse the specified directory and its subdirectories to extract UIDs from filenames.
    The filenames are expected to be in the format '<name>_UID.<extension>', where '<name>' can contain underscores.
    You can specify a UID to exclude from the collected UIDs.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    exclude_uid : str, default None
        A UID to exclude from the collected UIDs.

    Returns
    --------
    list :
        A list of all UIDs.
    """

    uids = []
    for _, _, files in os.walk(repo_path):
        for file in files:
            if '.yaml' in file:
                uid = file.split('_')[-1].split('.yaml')[0]
                uids.append(uid)

    return uids


def check_uid(new_uid, repo_path="."):
    """
    Check whether the new UID is unique.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    new_uid : str
        The UID to be checked for uniqueness.

    Returns
    --------
    bool :
        True if the new UID is unique, False otherwise.
    """

    uids = get_all_uids(repo_path=repo_path)
    
    print(f"Your UID: {new_uid}")
    print(f"All UIDs: {uids}")

    uid_counts = Counter(uids)

    if uid_counts[new_uid] > 1:
        print(f"Duplicate UID found: {new_uid}")
        return False
    else:
        print(f"No duplicate UID found.")
        return True

# Run the check and exit with error if the UID is not unique
if not check_uid(new_uid, repo_path="."):
    raise ValueError(f'UID {new_uid} already exists.')
