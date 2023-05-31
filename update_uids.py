import os
import re
import yaml


def get_all_paths_and_uids(path):
    """
    Get all files and UIDs in the given path, including subdirectories.

    Parameters
    ----------
    path : str
        The path to the folder containing the marker lists.

    Returns
    -------
    dict :
        A dictionary with the UIDs as keys and lists of file paths as values.
    """
    files_with_uids = {}
    pattern = re.compile(r"_([0-9]+)\.yaml$")
    for root, _, files in os.walk(path):
        for file in files:
            if file.endswith(".yaml"):
                match = pattern.search(file)
                if match:
                    uid = int(match.group(1))
                    file_path = os.path.join(root, file)
                    if uid in files_with_uids:
                        files_with_uids[uid].append(file_path)
                    else:
                        files_with_uids[uid] = [file_path]

    return files_with_uids


def update_uids(repo_lists_path):
    """
    Update the UIDs of the files in the given path to ensure uniqueness.

    Parameters
    ----------
    repo_lists_path : str
        The path to the folder containing the marker lists.

    """
    files_with_uids = get_all_paths_and_uids(repo_lists_path)
    uids_in_use = set(files_with_uids.keys())

    for uid, file_paths in files_with_uids.items():
        if len(file_paths) > 1:
            file_paths.sort(key=os.path.getmtime)
            for file_path in file_paths[1:]:
                new_uid = uid
                while new_uid in uids_in_use:
                    new_uid += 1
                uids_in_use.add(new_uid)
                new_file_path = re.sub(r"_([0-9]+)\.yaml$", f"_{new_uid}.yaml", file_path)
                
                # Update the value of the "id" key in the "metadata" section
                with open(file_path, 'r') as file:
                    content = yaml.safe_load(file)
                content['metadata']['id'] = new_uid
                with open(new_file_path, 'w') as file:
                    yaml.safe_dump(content, file)

                os.remove(file_path)
                print(f"File '{file_path}' was renamed to '{new_file_path}' and its ID was updated to {new_uid}")


def update():
    lists_path = "lists"
    update_uids(lists_path)
