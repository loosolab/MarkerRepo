import argparse
import os

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("uid", help="The UID to be checked for uniqueness.")
args = parser.parse_args()
new_uid = args.uid

def get_all_uids(repo_lists_path="./lists", exclude_uid=None):
    """ 
    Traverse the specified directory and its subdirectories to extract UIDs from filenames.
    The filenames are expected to be in the format '<name>_UID.<extension>', where '<name>' can contain underscores.
    You can specify a UID to exclude from the collected UIDs.

    Parameters
    ----------
    repo_lists_path : str, default "./lists"
        The path where the lists of the Marker Repo are stored - probable 'REPO_PATH/lists'.
    exclude_uid : str, default None
        A UID to exclude from the collected UIDs.

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
                if uid != exclude_uid:
                    uids.append(uid)
    
    return uids


def check_uid(new_uid, repo_lists_path="./lists"):
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

    uids = get_all_uids(repo_lists_path=repo_lists_path, exclude_uid=new_uid)
    
    print(f"Your UID: {new_uid}")
    print(f"All UIDs: {uids.sort()}")

    if new_uid in uids:
        print(f"Duplicate UID found: {new_uid}")
        return False
    else:
        print(f"No duplicate UID found.")
        return True

# Run the check and exit with error if the UID is not unique
if not check_uid(new_uid, repo_lists_path="./lists"):
    raise ValueError(f'UID {os.getenv("NEW_UID")} already exists.')
