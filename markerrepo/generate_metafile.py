import sys
from tabulate import tabulate
from .utils import read_in_yaml, save_as_yaml, find_keys, get_whitelist
import datetime
import os
import readline
import re
import copy

size = 100 # size of screen
factor = []
not_editable = ['list_name', 'id', 'marker_list', 'organism', 'marker_type']
id = ''
exp_fac = {}
list_def = []

class WhitelistCompleter:
    def __init__(self, whitelist):
        self.whitelist = whitelist

    def complete(self, text, state):
        results = [x for x in self.whitelist if
                   re.search(text, x, re.IGNORECASE) is not None] + [None]
        if len(results) > 30:
            results = results[:30]
        return results[state]


# ---------------------------------GENERATE-------------------------------------

def generate_file(path, input_id, name, mandatory_mode, marker_list, organism, marker_type):
    """
    This function is used to generate metadata by calling functions to compute
    user input. It writes the metadata into a yaml file after validating it.
    :param path: the path to the folder the metadata file should be saved to
    :param input_id: the ID of the experiment
    :param mandatory_mode: if True only mandatory files are filled out
    """

    global id
    id = input_id

    file_name = f'{name}_{input_id}.yaml'

    # test if list already exists
    if os.path.exists(
            os.path.join(path, file_name)) or os.path.exists(
            os.path.join(path, file_name)):
        print(f'The marker list {input_id} already exists.')
        overwrite = parse_list_choose_one(['True ', 'False '],
                              f'\nDo you want to overwrite the file?')
        if not overwrite:
            sys.exit(f'Program terminated.')

    # read in structure file
    key_yaml = read_in_yaml(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                     'keys.yaml'))

    # create metadata dictionary and fill it with the given organism and marker type
    org_dict = {'organism_name': organism.split()[0], 'taxonomy_id': organism.split()[1]}
    result_dict = {'metadata': {'id': input_id, 'list_name': name, 'organism': org_dict, 'marker_type': marker_type}}

    # parse through metadata structure and fill it for every key
    for item in key_yaml:
        # print(f"ITEM: {key_yaml[item]}")
        if item == 'marker_list':
            result_dict['marker_list'] = marker_list

        elif item in result_dict:
            result_dict[item] = {**result_dict[item],
                                 **get_redo_value(key_yaml[item], item, False,
                                                  mandatory_mode, result_dict,
                                                  True, False, True)}
        else:
            result_dict[item] = get_redo_value(key_yaml[item], item, False,
                                               mandatory_mode, result_dict,
                                               True, False, True)
        print(f'{"".center(size, "-")}\n'
              f'{"SUMMARY".center(size, " ")}\n'
              f'{"".center(size, "-")}\n')
        sum = print_summary(result_dict[item], 1, False)
        print(sum)
        print(f'\n\n')
        print(f'{"".center(size, "-")}\n')

        correct = True
        #correct = parse_list_choose_one(['True ', 'False '],
        #                                f'\nIs the input correct? You can redo it by selecting \'False\'')
        while not correct:
            result_dict[item] = edit_item(item, result_dict[item], key_yaml[item], result_dict, mandatory_mode)

            print(f'{"".center(size, "-")}\n'
                  f'{"SUMMARY".center(size, " ")}\n'
                  f'{"".center(size, "-")}\n')
            sum = print_summary(result_dict[item], 1, False)
            print(sum)
            print(f'\n\n')
            print(f'{"".center(size, "-")}\n')

            correct = parse_list_choose_one(['True ', 'False '],
                                            f'\nIs the input correct? You can redo it by selecting \'False\'')


    # TODO: extra function for summary
    # print summary
    print(f'{"".center(size, "-")}\n'
          f'{"SUMMARY".center(size, " ")}\n'
          f'{"".center(size, "-")}\n')
    sum = print_summary(result_dict, 1, False)
    print(sum)
    print(f'\n\n')
    print(f'{"".center(size, "-")}\n')

    # print(f'{"".center(size, "-")}\n'
    #       f'{"FILE VALIDATION".center(size, " ")}\n'
    #       f'{"".center(size, "-")}\n')
    # valid, missing_mandatory_keys, invalid_keys, invalid_entries, \
    #     invalid_values, pool_warn, ref_genome_warn = validate_yaml.\
    #     validate_file(result_dict)
    # if not valid:
    #     validate_yaml.print_validation_report(
    #         result_dict, missing_mandatory_keys, invalid_keys,
    #         invalid_entries, invalid_values)
    # elif len(pool_warn) > 0 or len(ref_genome_warn) > 0:
    #     validate_yaml.print_warning(result_dict, pool_warn, ref_genome_warn)
    # else:
    #     print(f'Validation complete. No errors found.\n')

    list_path = os.path.join(path, file_name)
    save_as_yaml(result_dict, list_path)
    print(f'Marker list saved:\n{list_path}')
    
    # not needed for marker repo
    # print_sample_names(result_dict, input_id, path)


def edit_item(key_name, item, key_yaml, result_dict, mandatory_mode):
    # TODO: add + remove list elements
    
    if isinstance(item, list):
        if all(not isinstance(x, dict) for x in item):
            item = get_redo_value(key_yaml, key_name, False, mandatory_mode, result_dict, True, False, True)
        else:
            all_options = []
            for i in range(len(item)):
                if isinstance(item[i], dict):
                    all_options.append('\n'.join(f'{x}: {item[i][x]}' for x in item[i]))
                else:
                    all_options.append(item[i])
            if len(all_options) > 1:
                print(f'Please choose the list elements you want to edit (1-{len(item)}) divided by comma.')
                print_option_list(all_options, False)
                chosen_options = parse_input_list(all_options, False)
                options = []
                for i in range(len(all_options)):
                    if all_options[i] in chosen_options:
                        options.append(i)
                for i in options:
                    item[i] = edit_item(key_name, item[i], key_yaml, result_dict, mandatory_mode)
            else:
                item = edit_item(key_name, item[0], key_yaml, result_dict, mandatory_mode)
    elif isinstance(item, dict):
        print(f'Please choose the keys (1-{len(item.keys())}) you want to edit divided by comma.')
        options = [key for key in key_yaml['value'] if key not in not_editable]
        options.insert(0, 'all')
        print_option_list(options, False)
        options = parse_input_list(options, False)
        if 'all' in options:
            new_item = get_redo_value(key_yaml, key_name, False,
                    mandatory_mode, result_dict, True, False, False)
            if isinstance(new_item, list):
                item = new_item
            else:
                item = {**result_dict[key_name],
                        **new_item}
        else:
            for key in options:
                if key == 'organism':
                    item = get_redo_value(key_yaml, 'experimental_setting', False, mandatory_mode, result_dict, True, False, False)
                elif key == 'experimental_factors' and 'organism' not in options:
                    new_yaml = copy.deepcopy(key_yaml)
                    new_yaml['value'].pop('organism')
                    item = get_redo_value(new_yaml, 'experimental_setting',
                                          False, mandatory_mode, copy.deepcopy(item),
                                          True, False, False)
                elif key == 'conditions' and 'experimental_factors' not in options and 'organism' not in options:
                    for i in range(len(item['experimental_factors'])):
                        for fac in list_def:
                            if item['experimental_factors'][i]['factor'] == fac['factor']:
                                item['experimental_factors'][i] = fac
                                break
                    new_yaml = copy.deepcopy(key_yaml)
                    new_yaml['value'].pop('organism')
                    new_yaml['value'].pop('experimental_factors')
                    item = get_redo_value(new_yaml, 'experimental_setting',
                                          False, mandatory_mode,
                                          copy.deepcopy(item),
                                          True, False, False)
                else:
                    if key in item:
                        item[key] = edit_item(key, item[key], key_yaml['value'][key], result_dict, mandatory_mode)
                    else:
                        item[key] = get_redo_value(key_yaml['value'][key], key, False, mandatory_mode, result_dict, True, False, True)
    else:
        item = parse_input_value(key_name, key_yaml['desc'], key_yaml['whitelist'], key_yaml['input_type'], item)
    return item


def get_redo_value(node, item, optional, mandatory_mode, result_dict,
                   first_node, is_factor, do_redo):
    """
    This function tests whether a list must be specified for a value and,
    depending on this, calls a function to enter the value.
    :param node: a part of the keys.yaml that is being iterated over
    :param item: the name of the key that is being filled
    :param optional: a bool that states if a key is optional
    :param mandatory_mode: a bool that states if the mandatory mode is active
    :param result_dict: a dictionary that contains the filled metadata
    :param first_node: a bool that states how a header should be printed
    :param is_factor: a bool that states if the current value is an
                      experimental factor
    :return: value: the value that was entered for the key
    """
    # test if the input value is of type list
    if node['list']:
        # test if one list element contains a dictionary
        if isinstance(node['value'], dict) and not \
                set(['mandatory', 'list', 'desc', 'display_name', 'value']) \
                <= set(node['value'].keys()):
            # test if the input value is an experimental factor
            if is_factor:

                value = fill_metadata_structure(node['value'], item, {},
                                                optional,
                                                mandatory_mode, result_dict,
                                                first_node, is_factor)
            else:
                # repeat the input prompt for the keys of the list element
                # until the user specifies it is complete
                redo = True
                value = []
                while redo:
                    value.append(
                        fill_metadata_structure(node['value'], item, {}, optional,
                                                mandatory_mode, result_dict,
                                                first_node, is_factor))
                    if do_redo:
                        redo = parse_list_choose_one(['True ', 'False '],
                                                 f'\nDo you want to add '
                                                 f'another {item}?')
                    else:
                        redo = False
        else:
            # enable the input of multiple elements of the whitelist
            value = get_input_list(node, item, result_dict)

    else:

        # test if the input value is an experimental factor and is either not
        # a dictionary or a dictionary that contains the key 'merge' as special
        # case (means that the input is treated like a single value and then
        # split into a dictionary, e.g. gene -> gene_name, ensembl_id)
        if is_factor and \
                (not isinstance(node['value'], dict) or
                    (isinstance(node['value'], dict) and
                        set(
                            ['mandatory', 'list', 'desc', 'display_name',
                             'value']) <= set(node['value'].keys()) or
                        ('special_case' in node and 'merge' in
                            node['special_case']))):

            # ask the user to put in a list of experimental factors
            value = get_input_list(node, item, result_dict)

        else:
            # call function to fill in value
            value = fill_metadata_structure(node, item, {}, optional,
                                            mandatory_mode, result_dict,
                                            first_node, is_factor)
    return value


def fill_metadata_structure(node, key, return_dict, optional, mandatory_mode,
                            result_dict, first_node, is_factor, custom=None):
    """
    This function calls other functions to fill in metadata information for a
    key depending on its type.
    :param node: a part of the keys.yaml that contains the key to be filled
    :param key: the name of the key to be filled
    :param return_dict: a dictionary that contains all filled in information
    :param optional: a bool that states if a key is optional
    :param mandatory_mode: a bool that states if mandatory mode is active
    :param result_dict: a dictionary that contains all filled in information
    :param first_node: a bool that states how a header should be printed
    :param is_factor: a bool that states if the current value is an
                      experimental factor
    :return: return_dict: a dictionary containing the filled information
    """

    # test if the given part of the metadata structure is a dictionary
    if isinstance(node, dict) and not \
            set(['mandatory', 'list', 'desc', 'display_name', 'value']) <= \
            set(node.keys()):
        # test if the input is of type value_unit
        if len(node.keys()) == 2 and 'value' in node.keys() and 'unit' in \
                node.keys():

            if is_factor:

                # call function to get a list of value_units if the input value
                # is an experimental factor
                return get_list_value_unit(result_dict, key)

            else:

                # call function to fill a single value_unit
                return get_value_unit(result_dict)

        else:

            # create a lists to store optional keys and their description
            optionals = []
            desc = []

            # iterate through all keys in the given part of the metadata
            # structure
            for item in node:
                # parameter to declare if a key is optional, False per default
                optional = False

                # test if the key is editable
                if item not in not_editable:

                    # TODO: specify function
                    # if the key is mandatory, call the ... function on it

                    # Debugging:
                    # print(item)
                    # print(node)
                    # print(node[item])

                    if node[item]['mandatory']:
                        return_dict[item] = get_redo_value(node[item], item,
                                                           optional,
                                                           mandatory_mode,
                                                           return_dict, False,
                                                           is_factor, True)

                    else:

                        # TODO:
                        if node[item]['list'] or item not in factor or \
                                is_factor:
                            optionals.append(item)
                            desc.append(node[item]['desc'])

            if custom is not None:
                optionals.append(custom['display_name'])
                desc.append(custom['desc'])

            # if there are optional keys and mandatory mode is not active, ask
            # the user whether he wants to add optional information
            if len(optionals) > 0 and mandatory_mode == False:

                # set te optional parameter to True
                optional = True

                print(
                    f'\nDo you want to add any of the following optional keys?'
                    f' (1,...,{len(optionals)} or n)')

                # print a list of all optional keys
                print_option_list(optionals, desc)

                # parse the user input
                options = parse_input_list(optionals, True)

                if options:

                    # if the user chose optional keys, iterate through every
                    # chosen optional key
                    for option in options:

                        # Loop to enter Custom Keys
                        if option == 'Custom Tag':
                            redo = True
                            while redo:
                                print('Enter a key for your custom tag:')
                                custom_key = parse_input_value('Custom tag', '', False,
                                                    'str', result_dict)
                                custom_value = get_redo_value(custom, custom_key, optional, mandatory_mode, result_dict, first_node, is_factor, True)
                                return_dict[custom_key] = custom_value
                                redo = parse_list_choose_one(['True ', 'False '], 'Do you want to add another custom key?')

                        else:
                            # TODO: runter rücken?
                            new_element = True

                            # test if the value for the optional key is of type
                            # list
                            if node[option]['list']:

                                # TODO: ??
                                if option in return_dict and all(
                                        isinstance(x, dict) for x in
                                        return_dict[option]):

                                    # set the new_element parameter to False
                                    new_element = False

                                    # read in the structure file
                                    key_yaml = read_in_yaml(
                                        os.path.join(os.path.dirname(
                                            os.path.abspath(__file__)), '..',
                                            'keys.yaml'))

                                    # list all possible keys that can occur in one
                                    # list element of the optional key
                                    possible_keys = list(list(
                                        find_keys(
                                            key_yaml, option))[0]['value'].keys())

                                    # create lists for list elements that contain
                                    # optional keys that have no value yet and the
                                    # descriptions of the list elements
                                    elems = []
                                    desc = []

                                    # iterate through all list elements that
                                    # already exist for the optional key
                                    for i in range(len(return_dict[option])):

                                        # test if optional keys are missing in the
                                        # list element
                                        if not all(k in return_dict[option][i] for
                                                   k in possible_keys):

                                            # add the list element with missing
                                            # optional keys to the 'elems' list
                                            elems.append(', '.join(
                                                [f'{k}: '
                                                 f'{return_dict[option][i][k]}'
                                                 for k in return_dict[option][i]]))

                                            # add a description of the list element
                                            # with missing optional keys to the
                                            # 'desc' list
                                            keys = ', '.join(
                                                [x for x in possible_keys if x
                                                 not in return_dict[option][i]])
                                            desc.append(
                                                f'Possible information to add: '
                                                f'{keys}')

                                    # test if list with list elems containing
                                    # unfilled optional keys is not empty
                                    if len(elems) > 0:

                                        # add an option for adding a new list value
                                        # to the 'elems' list (with empty desc)
                                        elems.append(f'Add new {option}')
                                        desc.append('')

                                        # print information for user that there are
                                        # existing list elements for the optional
                                        # key and let the user select if he wants
                                        # to add information to an existing list
                                        # element or create a new one + parse the
                                        # user input
                                        print(
                                            f'\nThere are existing elements for '
                                            f'{option}. Please select the elements'
                                            f' for which you want to add '
                                            f'information.\n')
                                        print_option_list(elems, desc)
                                        list_elems = parse_input_list(
                                            range(len(return_dict[option]) + 1),
                                            False)

                                        # iterate through the list elements the
                                        # user chose to add information
                                        for indc in list_elems:

                                            # test if the user chose an existing
                                            # list element to add information
                                            if int(indc) - 1 < len(
                                                    return_dict[option]):

                                                # set a caption to show the user
                                                # which existing list element he is
                                                # about to edit and print it
                                                caption = elems[
                                                    int(indc) - 1].replace("\n",
                                                                           ", ").\
                                                    center(size, ' ')
                                                line = ''.center(size, '_')
                                                print(f'\n'
                                                      f'{line}\n\n'
                                                      f'List element: {caption}\n'
                                                      f'{line}\n')

                                                # save all unfilled optional keys
                                                # of the list element into a list
                                                possible_input = [
                                                    x for x in possible_keys
                                                    if x not in list(
                                                        return_dict[option]
                                                        [int(indc) - 1].keys())]

                                                # test if there are multiple
                                                # unfilled keys in the list element
                                                if len(possible_input) > 1:

                                                    # copy the metadata structure
                                                    # of the optional key
                                                    part_node = copy.deepcopy(
                                                        node[option])

                                                    # iterate through the keys
                                                    # contained as value of the
                                                    # optional key and those keys
                                                    # to a list if they are
                                                    # already filled
                                                    remove_keys = []
                                                    for k in part_node['value']:
                                                        if k not in possible_input:
                                                            remove_keys.append(k)

                                                    # remove the filled keys from
                                                    # the metadata structure of
                                                    # the optional key
                                                    for k in remove_keys:
                                                        part_node['value'].pop(k)

                                                    # call the
                                                    # fill_metadata_structure
                                                    # function for the structure
                                                    # of the optional key without
                                                    # the already filled keys
                                                    val = fill_metadata_structure(
                                                        part_node, option, {},
                                                        optional, mandatory_mode,
                                                        result_dict, False,
                                                        is_factor)

                                                    # merge the prefilled optional
                                                    # key with the new input
                                                    # information
                                                    return_dict[
                                                        option][int(indc) - 1] \
                                                        = merge_dicts(
                                                        return_dict[option][
                                                            int(indc) - 1], val)

                                                # if there is just one unfilled key
                                                else:

                                                    # find the structure of the
                                                    # unfilled key in the metadata
                                                    # structure and save it
                                                    part_node = list(
                                                        find_keys(
                                                            key_yaml,
                                                            possible_input[0]))[0]

                                                    # call function to fill the
                                                    # unfilled key
                                                    val = get_redo_value(
                                                        part_node,
                                                        possible_input[0],
                                                        optional, mandatory_mode,
                                                        result_dict, False,
                                                        is_factor, True)

                                                    # save the now filled key in
                                                    # the dictionary
                                                    return_dict[option][
                                                        int(indc) - 1][
                                                        possible_input[0]] = val

                                            else:

                                                # if the user chose to add a new
                                                # element to the list, print a
                                                # caption for the new element and
                                                # set new_element to True

                                                h_line = ''.center(size,
                                                                   '_')
                                                caption = f'New {option}'.center(
                                                    size, ' ')
                                                print(f'\n'
                                                      f'{h_line}\n\n'
                                                      f'{caption}\n'
                                                      f'{h_line}\n')
                                                new_element = True

                                    else:

                                        # set new_element to True if there are no
                                        # list elements with unfilled keys
                                        new_element = True

                            # create a new list element if new_element is set to
                            # True
                            if new_element:

                                if 'special_case' in node[option] and 'merge' in \
                                        node[option]['special_case']:
                                    value = parse_input_value(option,
                                                              node[option]['desc'],
                                                              True, 'str',
                                                              result_dict)
                                    return_dict[option] = value
                                else:
                                    val = get_redo_value(node[option],
                                                         option,
                                                         optional,
                                                         mandatory_mode,
                                                         result_dict, False,
                                                         is_factor, True)
                                    if node[option]['list']:
                                        if option in return_dict:
                                            return_dict[option] += val
                                        else:
                                            return_dict[option] = val
                                    else:
                                        return_dict[option] = val
    else:
        if node['mandatory'] or optional or is_factor:
            if 'special_case' in node and 'merge' in node['special_case']:
                value = parse_input_value(key, node['desc'], True, 'str',
                                          result_dict)
            else:
                value = enter_information(node, key, return_dict, optional,
                                          mandatory_mode, result_dict,
                                          first_node, is_factor)
            return value
    return return_dict


# ---------------------------------SUMMARY--------------------------------------


def print_summary(result, depth, is_list):
    """
    This function parses the dictionary into a string with the same structure
    as it will be saved to the yaml file.
    :param result: the filled dictionary
    :param depth: an integer that specifies the depth of indentation
    :param is_list: a bool that states if a key contains a list
    :return: summary: a string that contains all entered information
    """
    summary = ''
    if isinstance(result, dict):
        for key in result:
            if key == list(result.keys())[0] and is_list:
                summary = f'{summary}\n{"    " * (depth - 1)}{"  - "}{key}: ' \
                          f'{print_summary(result[key], depth + 1, is_list)}'
            else:
                summary = f'{summary}\n{"    " * depth}{key}: ' \
                          f'{print_summary(result[key], depth + 1, is_list)}'
    elif isinstance(result, list):
        for elem in result:
            if not isinstance(elem, list) and not isinstance(elem, dict):
                summary = f'{summary}\n{"    " * (depth - 1)}{"  - "}{elem}'
            else:
                summary = f'{summary}{print_summary(elem, depth, True)}'
    else:
        summary = f'{summary}{result}'
    return summary


# ---------------------------------UTILITIES------------------------------------


def enter_information(node, key, return_dict, optional, mandatory_mode,
                      result_dict, first_node, is_factor):
    """
    This function is used to create prompts for the user to enter information
    and parses the input.
    :param node: a part of the keys.yaml
    :param key: the name of the key that should be filled
    :param return_dict: a dictionary containing input information
    :param optional: a bool to state if the key is optional
    :param mandatory_mode: a bool to state if mandatory mode is active
    :param result_dict: a dictionary containing all filled information
    :param first_node: a bool to state how a header should be printed
    :param is_factor: a bool to state if the key is an experimental factor
    :return: the filled key
    """
    # test if the key contains a dictionary
    if isinstance(node['value'], dict) and \
            not set(['mandatory', 'list', 'desc', 'display_name', 'value']) \
            <= set(node['value'].keys()):
        display_name = node['display_name']
        if first_node:

            # if the key is on top level, print a bigger caption
            print(f'{"".center(size, "_")}\n\n'
                  f'{f"{display_name}".center(size, " ")}\n'
                  f'{"".center(size, "_")}\n')
        else:

            # if the key is on a lower level, print a smaller caption
            print(f'\n'
                  f'{"".center(size, "-")}\n'
                  f'{f"{display_name}".center(size, " ")}\n'
                  f'{"".center(size, "-")}\n')

        # print a description if one is stated in the metadata structure
        if node['desc'] != '':
            print(f'{node["desc"]}\n')

        if 'special_case' in node and 'custom_key' in node['special_case']:
            custom = node['special_case']['custom_key']
        else:
            custom = None
        # call fill_metadata_structure to fill in the dictionary
        return fill_metadata_structure(node['value'], key, return_dict,
                                       optional,
                                       mandatory_mode, result_dict, False,
                                       is_factor, custom)

    else:
        # call parse_input_value to fill in a single value
        return parse_input_value(key, node['desc'], node['whitelist'],
                                 node['input_type'], result_dict)


def parse_list_choose_one(whitelist, header):
    """
    This function prints an indexed whitelist and prompts the user to choose a
    value by specifying the index.
    :param whitelist: a list of values for the user to choose from
    :param header: a headline or description that should be printed
    :return: value: the whitelist value chosen by the user
    """

    # print the given header and indexed whitelist and parse the user input
    try:
        print(f'{header}')
        print_option_list(whitelist, False)
        value = whitelist[int(input()) - 1]

    # redo the input prompt if the user input is not an integer
    except (IndexError, ValueError) as e:
        print(f'Invalid entry. Please enter a number between 1 and '
              f'{len(whitelist)}')
        value = parse_list_choose_one(whitelist, header)

    if value == 'True ':
        value = True
    elif value == 'False ':
        value = False

    # return the user input
    return value


def parse_input_value(key, desc, has_whitelist, value_type, result_dict):
    """
    This function lets the user enter an input and tests if the input is valid.
    :param key: the key that should be filled
    :param desc: a description that should be printed
    :param has_whitelist: a bool to state if the key has a whitelist
    :param value_type: the type that is expected for the value
    :param result_dict: a dictionary containing filled information
    :return: input_value: the value that was input by the user
    """

    # read in whitelist if the key has one or set the whitelist to None
    if has_whitelist:
        whitelist = get_whitelist(key, result_dict)
    else:
        whitelist = None

    # test if the whitelist is not None
    if whitelist is not None:

        # test if the whitelist is a dictionary
        if isinstance(whitelist['whitelist'], dict):

            # create a list containing the values of every key in the whitelist
            w = [x for xs in list(whitelist['whitelist'].values()) if xs is
                 not None for x in xs]

            # use autocomplete if the whitelist contains more than 30 values
            if len(w) > 30:
                input_value = complete_input(w, key)

                # repeat if the input value does not match the whitelist
                while input_value not in w:
                    print(f'The value you entered does not match the '
                          f'whitelist. Try tab for autocomplete.')
                    input_value = complete_input(w, key)

            else:

                # parse grouped whitelist if the whitelist contains less than
                # 30 values
                input_value = parse_group_choose_one(whitelist['whitelist'],
                                                     w, f'{key}:')

        # use autocomplete if the whitelist is a list longer than 30
        elif len(whitelist['whitelist']) > 30:
            input_value = complete_input(whitelist['whitelist'], key)

            # repeat if the input value does not match the whitelist
            while input_value not in whitelist['whitelist']:
                print(f'The value you entered does not match the '
                      f'whitelist. Try tab for autocomplete.')
                input_value = complete_input(whitelist['whitelist'], key)

        # call parse_list_choose_one to print indexed whitelist and parse user
        # input
        else:
            input_value = parse_list_choose_one(whitelist['whitelist'],
                                                f'{key}:')
        w_key = None
        if whitelist['whitelist_type'] == 'plain_group':
            for k in whitelist['whitelist_keys']:
                if input_value.endswith(f' ({k})'):
                    input_value = input_value.replace(f' ({k})', '')
                    w_key = k

        # test if the whitelist contains the key 'headers'
        if 'headers' in whitelist:
            if whitelist['whitelist_type'] == 'group' or whitelist['whitelist_type'] == 'plain_group' and w_key is not None:
                if w_key in whitelist['headers']:
                    headers = whitelist['headers'][w_key].split(' ')
                    vals = input_value.split(' ')

                    # create a dictionary to store the new value
                    value = {}

                    # iterate through the headers and save the header and value of the
                    # same index into a dictionary with header as key
                    for i in range(len(headers)):
                        value[headers[i]] = vals[i]

                    # overwrite the input value with the dictionary
                    input_value = value

            else:

                # split the headers and the input value at ' ' and save each to
                # a list
                headers = whitelist['headers'].split(' ')
                vals = input_value.split(' ')

                # create a dictionary to store the new value
                value = {}

                # iterate through the headers and save the header and value of the
                # same index into a dictionary with header as key
                for i in range(len(headers)):
                    value[headers[i]] = vals[i]

                # overwrite the input value with the dictionary
                input_value = value

    # if there is no whitelist
    else:

        # print a description if one is stated in the metadata structure
        if desc != '':
            print(f'\n{desc}')

        # prompt the user to choose between True and False if the input value
        # is of type boolean
        if value_type == 'bool':
            input_value = parse_list_choose_one(['True ', 'False '],
                                                f'Is the sample {key}')

        else:

            # print the key, add a newline if a description was printed
            if desc == '':
                input_value = input(f'\n{key}: ')
            else:
                input_value = input(f'{key}: ')

            # repeat input if nothing is passed
            if input_value == '':
                print(f'Please enter something.')
                input_value = parse_input_value(key, desc, whitelist,
                                                value_type, result_dict)

            # test if input fits the input_type int and repeat the input if not
            if value_type == 'number':
                try:
                    input_value = int(input_value)
                except ValueError:
                    print(f'Input must be of type int. Try again.')
                    input_value = parse_input_value(key, desc, whitelist,
                                                    value_type,
                                                    result_dict)

            # if the input_type is date, test if the date fits the format
            # dd.mm.yyyy and repeat the input if not
            elif value_type == 'date':
                try:
                    input_date = input_value.split('.')
                    if len(input_date) != 3 or len(
                            input_date[0]) != 2 or len(
                            input_date[1]) != 2 or len(input_date[2]) != 4:
                        raise SyntaxError
                    input_value = datetime.date(int(input_date[2]),
                                                int(input_date[1]),
                                                int(input_date[0]))
                    input_value = input_value.strftime("%d.%m.%Y")
                except (IndexError, ValueError, SyntaxError) as e:
                    print(f'Input must be of type \'DD.MM.YYYY\'.')
                    input_value = parse_input_value(key, desc, whitelist,
                                                    value_type,
                                                    result_dict)

            else:

                # repeat input if it is of type string and contains an
                # invalid character
                if '\"' in input_value:
                    print(f'Ivalid symbol \'\"\'. Please try again.')
                    input_value = parse_input_value(key, desc, whitelist,
                                                    value_type, result_dict)
                elif '{' in input_value:
                    print(f'Ivalid symbol \'{"{"}\'. Please try again.')
                    input_value = parse_input_value(key, desc, whitelist,
                                                    value_type, result_dict)
                elif '}' in input_value:
                    print(f'Ivalid symbol \'{"}"}\'. Please try again.')
                    input_value = parse_input_value(key, desc, whitelist,
                                                    value_type, result_dict)
                elif '|' in input_value:
                    print(f'Ivalid symbol \'|\'. Please try again.')
                    input_value = parse_input_value(key, desc, whitelist,
                                                    value_type, result_dict)

    # return the user input
    return input_value


def parse_group_choose_one(whitelist, w, header):
    """
    This function prints a grouped whitelist with indexes and lets the user
    choose one value.
    :param whitelist: the grouped whitelist as a dictionary
    :param w: a list containing all whitelist values
    :param header: a header that should be printed
    :return:
    """

    try:

        # print a header or description
        print(f'{header}\n')

        # set index to 1
        i = 1

        # iterate through keys in whitelist
        for key in whitelist:

            # print the indexed whitelist values with the keys as captions
            if whitelist[key] is not None:
                print(f'\033[1m{key}\033[0m')
                for value in whitelist[key]:
                    print(f'{i}: {value}')
                    i += 1

        # select the ehitelist value that fits the input index
        value = w[int(input()) - 1]

    # redo the input if the passed value is not an integer >1
    except (IndexError, ValueError) as e:
        print(f'Invalid entry. Please enter a number between 1 and '
              f'{len(whitelist)}')
        value = parse_list_choose_one(whitelist, header)

    # return the user input
    return value


def merge_dicts(a, b):
    """
    This function merges two dictionaries with the same structure to create
    one.
    :param a: the first dictionary
    :param b: the second dictionary
    :return: res: the merged dictionary
    """
    if isinstance(a, list):
        res = []
        for i in range(len(a)):
            res.append(merge_dicts(a[i], b[i]))
    elif isinstance(a, dict):
        b_keys = list(b.keys())
        res = {}
        for key in a.keys():
            if key in b_keys:
                res[key] = merge_dicts(a[key], b[key])
                b_keys.remove(key)
            else:
                res[key] = a[key]
        for key in b_keys:
            res[key] = b[key]
    else:
        res = a
    return res


def get_input_list(node, item, filled_object):
    """
    This function asks the user to enter a list of values divided by comma and
    parses the input.
    :param node: a part of the keys.yaml
    :param item: the key that should be filled
    :param filled_object: a dictionary containing filled information
    :return: used_values: the filled in list of values
    """
    if 'whitelist' in node and node['whitelist'] or 'special_case' in node \
            and 'merge' in node['special_case']:
        whitelist = get_whitelist(item, filled_object)
        if whitelist:
            if len(whitelist['whitelist']) > 30:
                used_values = []
                redo = True
                print(f'\nPlease enter the values for tags '
                      f'{item}.')
                while redo:
                    input_value = complete_input(whitelist['whitelist'],
                                                 item)
                    if input_value in whitelist['whitelist']:
                        used_values.append(input_value)
                        redo = parse_list_choose_one(
                            ['True ', 'False '],
                            f'\nDo you want to add another {item}?')
                    else:
                        print(f'The value you entered does not match the '
                              f'whitelist. Try tab for autocomplete.')
            else:
                print('\nPlease enter the list')
                if whitelist['whitelist_type'] == 'plain' or whitelist['whitelist_type'] == 'plain_group':
                    print_option_list(whitelist['whitelist'], '')
                    used_values = parse_input_list(whitelist['whitelist'],
                                                   False)
                elif whitelist['whitelist_type'] == 'group':
                    w = [x for xs in list(whitelist['whitelist'].values()) for
                         x in xs]
                    i = 1
                    for w_key in whitelist['whitelist']:
                        print(f'\033[1m{w_key}\033[0m')
                        for value in whitelist['whitelist'][w_key]:
                            print(f'{i}: {value}')
                            i += 1
                    used_values = parse_input_list(w, False)

            w_keys = []
            if whitelist['whitelist_type'] == 'plain_group':
                for i in range(len(used_values)):
                    for k in whitelist['whitelist_keys']:
                        if used_values[i].endswith(f' ({k})'):
                            used_values[i] = used_values[i].replace(f' ({k})', '')
                            w_keys.append(k)

            if 'headers' in whitelist:
                if whitelist['whitelist_type'] == 'group' or whitelist['whitelist_type'] == 'plain_group' and len(w_keys)>0:
                    for i in range(len(used_values)):
                        if w_keys[i] in whitelist['headers']:
                            headers = whitelist['headers'][w_keys[i]].split(' ')
                            vals = used_values[i].split(' ')
                            used_values[i] = {}
                            for j in range(len(headers)):
                                used_values[i][headers[j]] = vals[j]
                else:
                    headers = whitelist['headers'].split(' ')
                    for i in range(len(used_values)):
                        vals = used_values[i].split(' ')
                        used_values[i] = {}
                        for j in range(len(headers)):
                            used_values[i][headers[j]] = vals[j]
        else:
            print('No whitelist')
            used_values = [0]
    else:
        value_type = node['input_type']
        print(f'\nPlease enter a list of'
              f' {item}, which describe the list, divided by comma:\n')
        used_values = parse_input_list(value_type, False)
    return used_values


def print_option_list(options, desc):
    """
    This function prints an indexed whitelist.
    :param options: the whitelist values
    :param desc: a description to be printed
    """
    if desc:
        data = [[f'{i + 1}:', f'{options[i]}', desc[i]] for i in
                range(len(options))]
        print(tabulate(data, tablefmt='plain',
                       maxcolwidths=[size * 1 / 8,
                                     size * 3 / 8,
                                     size * 4 / 8]))
    else:
        data = [[f'{i + 1}:', f'{options[i]}'] for i in range(len(options))]
        print(tabulate(data, tablefmt='plain',
                       maxcolwidths=[size * 1 / 8,
                                     size * 7 / 8]))


def parse_input_list(options, terminable):
    """
    This function parses the user input for a list.
    :param options: possible input options
    :param terminable: a bool to state in nothing can be input
    :return: input_list: a list containing the input values
    """
    input_list = input()
    if terminable and input_list.lower() == 'n':
        return None
    else:
        if isinstance(options, list):
            try:
                input_list = [options[int(i.strip()) - 1] for i in
                              input_list.split(',')]
            except (IndexError, ValueError) as e:
                print(f'Invalid entry, try again:')
                input_list = parse_input_list(options, terminable)
        else:
            try:
                input_list = [x.strip() for x in input_list.split(',')]
                if options == 'number':
                    for i in range(len(input_list)):
                        input_list[i] = int(input_list[i])
            except (ValueError, IndexError) as e:
                print(f'Invalid entry. Please enter {options} numbers divided '
                      f'by comma.')
                input_list = parse_input_list(options, terminable)
    return input_list


def complete_input(whitelist, key):
    """
    This function uses a completer to autofill user input and print matching
    values.
    :param whitelist: a whitelist with possible values
    :param key: the key that should be filled
    :return: inpu_value: the value that was input by the user
    """
    print(f'\nPress tab once for autofill if'
          f' possible or to get a list of up to'
          f' 30 possible input values.')
    readline.parse_and_bind("tab: complete")
    readline.parse_and_bind("set show-all-if-ambiguous On")
    readline.parse_and_bind("set show-all-if-unmodified On")
    completer = WhitelistCompleter(whitelist)
    readline.set_completer(completer.complete)
    readline.set_completer_delims('')
    input_value = input(f'{key}: ')
    readline.parse_and_bind('tab: self-insert')
    return input_value


def get_value_unit(result_dict):
    """
    This function prompts the user to input a value_unit.
    :param result_dict: a dictionary containing filled information
    :return: val_un: a dictionary containing the unit and value
    """
    val_un = {'unit': parse_input_value('unit', '', True, str, result_dict),
              'value': parse_input_value('value', '', False, int, result_dict)}
    return val_un


def get_list_value_unit(result_dict, ex_factor):
    """
    This function prompts the user to enter a list of value_units.
    :param result_dict: a dictionary containing filled information
    :param factor: the experimental factor that contains the value_unit
    :return: a dictionary containing a unit and a list of values
    """
    print(f'\nPlease enter the unit for factor {ex_factor}:')
    unit = parse_input_value('unit', '', True, 'select', result_dict)
    val_un = []
    print(f'\nPlease enter int values for factor {ex_factor} (in {unit}) '
          f'divided by comma:')
    value = parse_input_list('number', False)
    for val in value:
        val_un.append({'unit': unit, 'value': val})
    return val_un
