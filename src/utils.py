import yaml
import os
import copy

# The following functions were copied from Mampok
# https://gitlab.gwdg.de/loosolab/software/mampok/-/blob/master/mampok/utils.py


def save_as_yaml(dictionary, file_path):
    """
    save dictionary as YAML file
    :param dictionary: a dictionary that should be saved
    :param file_path: the path of the yaml file to be created
    """
    with open(file_path, 'w') as file:
        documents = yaml.dump(dictionary, file, sort_keys=False)


def read_in_yaml(yaml_file):
    """
    read yaml, auto lower all keys
    :param yaml_file: the path to the yaml file to be read in
    :return: low_output: a dictionary containing the information of the yaml
    """
    with open(yaml_file) as file:
        output = yaml.load(file, Loader=yaml.FullLoader)
    low_output = {k.lower(): v for k, v in output.items()}
    return low_output


# The following function was found on Stackoverflow
# https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-
# nested-dictionaries-and-lists
# original function name: findkeys
# submitted by user 'arainchi' on Nov 9, 2013 at 3:14,
# edited on Sep 12, 2019 at 21:50

def find_keys(node, kv):
    """
    generator to return all values of given key in dictionary
    :param node: a dictionary to be searched
    :param kv: the key to search for
    """
    if isinstance(node, list):
        for i in node:
            for x in find_keys(i, kv):
                yield x
    elif isinstance(node, dict):
        if kv in node:
            yield node[kv]
        for j in node.values():
            for x in find_keys(j, kv):
                yield x


# The following function is based on 'findkeys' and customized to solve
# related problems

def find_values(node, kv):
    """
    generator to return all keys in dictionary that fit a given key value
    :param node: a dictionary to be searched
    :param kv: the key value to search for
    """
    if isinstance(node, list):
        for i in node:
            for x in find_values(i, kv):
                yield x
    elif isinstance(node, dict):
        for val in node.values():
            if isinstance(val, dict) or isinstance(val, list):
                for x in find_values(val, kv):
                    yield x
            else:
                if ((type(kv) is int or type(
                        kv) is float) and (
                            type(val) is int or type(val) is float)) or (
                        type(kv) is bool and type(val) is bool):
                    if kv == val:
                        yield kv
                else:
                    if all(elem in str(val).lower() for elem in str(kv).lower().split(' ')):
                        yield val


def read_whitelist(key):
    """
    This function reads in a whitelist and returns it.
    :param key: the key that contains a whitelist
    :return: whitelist: the read in whitelist
    """
    try:
        whitelist = read_in_yaml(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                         'metadata_whitelists', 'whitelists', key))
    except (AttributeError, FileNotFoundError):
        whitelist = None
    return whitelist


def read_grouped_whitelist(whitelist, filled_object):
    """
    This function parses a whitelist of type 'group'. If there are more than 30
    values it is formed into a plain whitelist.
    :param whitelist: the read in whitelist
    :param filled_object: a dictionary containing filled information
    :return: whitelist: the read in whitelist
    """
    headers = {}
    for key in whitelist['whitelist']:
        if not isinstance(whitelist['whitelist'][key], list) and \
                os.path.isfile(
                os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             '..', 'metadata_whitelists', 'whitelists',
                             whitelist['whitelist'][key])):
            whitelist['whitelist'][key] = \
                get_whitelist(whitelist['whitelist'][key], filled_object)
            if isinstance(whitelist['whitelist'][key], dict):
                if whitelist['whitelist'][key]['whitelist_type'] == 'depend':
                    if whitelist['whitelist'] and 'whitelist' in \
                            whitelist['whitelist']:
                        whitelist['whitelist'][key] = whitelist['whitelist']
                    else:
                        whitelist['whitelist'][key] = None
                elif whitelist['whitelist'][key]['whitelist_type'] == 'group':
                    whitelist['whitelist'][key] = \
                        [x for xs in list(whitelist['whitelist'][key].values())
                         for x in xs]
                elif whitelist['whitelist'][key]['whitelist_type'] == 'plain':
                    if 'headers' in whitelist['whitelist'][key]:
                        headers[key] = whitelist['whitelist'][key]['headers']
                    whitelist['whitelist'][key] = \
                        whitelist['whitelist'][key]['whitelist']
    w = [f'{x} ({xs})' for xs in list(whitelist['whitelist'].keys()) if whitelist['whitelist'][xs] is
    not None for x in whitelist['whitelist'][xs] if x is not None]
    #w = [f'{x}' for xs in list(whitelist['whitelist'].keys()) if
    #     whitelist['whitelist'][xs] is not None for x in
    #     whitelist['whitelist'][xs] if x is not None]

    if len(w) > 30:
        new_whitelist = copy.deepcopy(whitelist)
        new_whitelist['whitelist'] = w
        new_whitelist['whitelist_keys'] = list(whitelist['whitelist'].keys())
        whitelist = new_whitelist
        whitelist['whitelist_type'] = 'plain_group'
    if len(list(headers.keys())) > 0:
        whitelist['headers'] = headers
    return whitelist


def read_depend_whitelist(whitelist, depend):
    """
    This function parses a whitelist of type 'depend' in order to get the
    values fitting the dependency.
    :param whitelist: the read in whitelist
    :param depend: the key whose values the whitelist depends on
    :return: whitelist: the read in whitelist
    """
    if depend in whitelist:
        whitelist = whitelist[depend]
    elif os.path.isfile(os.path.join(os.path.dirname(
            os.path.abspath(__file__)), '..', 'metadata_whitelists',
            'whitelists', depend)):
        whitelist = read_whitelist(depend)
    if not isinstance(whitelist, list) and not isinstance(whitelist, dict) \
            and os.path.isfile(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '..',
            'metadata_whitelists', 'whitelists',
            whitelist)):
        whitelist = read_whitelist(whitelist)
    return whitelist


def get_whitelist(key, filled_object):
    """
    This function reads in a whitelist and parses it depending on its type.
    :param key: the key that contains a whitelist
    :param filled_object: a dictionary containing filled information
    :return: whitelist: the parsed whitelist
    """
    group = False
    stay_depend = False
    plain = False
    abbrev = False
    whitelist = read_whitelist(key)

    while isinstance(whitelist,
                     dict) and not group and not stay_depend and not \
            plain and not abbrev:
        if whitelist['whitelist_type'] == 'group':
            whitelist = read_grouped_whitelist(whitelist, filled_object)
            group = True
        elif whitelist['whitelist_type'] == 'plain':
            plain = True
        elif whitelist['whitelist_type'] == 'abbrev':
            abbrev = True
        elif whitelist['whitelist_type'] == 'depend':
            depend = list(find_keys(filled_object, whitelist['ident_key']))
            if len(depend) == 0:
                if whitelist['ident_key'] == 'organism_name':
                    depend = list(find_keys(filled_object, 'organism'))
            if len(depend) > 0:
                whitelist = read_depend_whitelist(whitelist['whitelist'],
                                                  depend[0].split(' ')[0])
            else:
                stay_depend = True

    if group:
        if whitelist['whitelist_type'] != 'plain_group':
            for key in whitelist['whitelist']:
                if whitelist['whitelist'][key] is not None and key != 'headers' \
                        and key != 'whitelist_type' and key != 'whitelist_keys':
                    whitelist['whitelist'][key] = sorted(
                        whitelist['whitelist'][key])

    elif whitelist and not stay_depend and not plain and not abbrev:
        whitelist = sorted(whitelist)
    return whitelist


def find_list_key(item, l):
    """
    This function finds an item in a list of dictionaries.
    :param item: the key that should be found
    :param l: a list of dictionaries
    :return: item: a list containing the values of all found keys
    """
    for k in l.split(':'):
        item = list(find_keys(item, k))
    return item
