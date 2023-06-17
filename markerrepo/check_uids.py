import argparse
import os

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("uid", help="The UID to be checked for uniqueness.")
args = parser.parse_args()
new_uid = args.uid

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

    uids = get_all_uids(repo_lists_path=repo_lists_path)
    print(new_uid, type(new_uid))
    print(uids)
    print(uids[0], type(uids[0]))
    if new_uid in uids:
        print(f"Duplicate UID found: {new_uid}")
        return False
    else:
        print(f"No duplicate UID found.")
        return True

# Run the check and exit with error if the UID is not unique
if not check_uid(new_uid, repo_lists_path="./lists"):
    raise ValueError(f'UID {os.getenv("NEW_UID")} already exists.')
