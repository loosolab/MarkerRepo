import os.path
import datetime
from src import utils


# This script includes functions for the validation of metadata yaml files

# ---------------------------------VALIDATION-----------------------------------

generated = ['condition_name', 'sample_name']
factor = None


def validate_file(metafile):
    """
    In this function all functions for the validation of a metadata file are
    called. The validation is based on the data in the file 'keys.yaml'. It is
    tested if all mandatory keys are included, if the included keys are valid
    and if the entered values correspond to the whitelist.
    :param metafile: the read in metadata yaml file
    :return:
    valid: bool, true if the file is valid, false if it is invalid
    missing_mandatory_keys: a list containing the missing mandatory keys
    invalid_keys: a list containing the invalid keys
    invalid_entries: a list containing the invalid entries -> (key, [values])
    """
    pool_warn = []
    ref_genome_warn = []
    valid = True
    key_yaml = utils.read_in_yaml(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'keys.yaml'))
    invalid_keys, invalid_entries, invalid_value = \
        new_test(metafile, key_yaml, [], '', [], [], [], None, [], None, metafile)
    missing_mandatory_keys = test_for_mandatory(metafile, key_yaml,
                                                [x.split(':')[-1] for x in
                                                 invalid_keys])
    if len(missing_mandatory_keys) > 0 or len(invalid_keys) > 0 or len(
            invalid_entries) > 0 or len(invalid_value) > 0:
        valid = False
    else:
        pool_warn, ref_genome_warn = validate_logic(metafile)
    return valid, missing_mandatory_keys, invalid_keys, invalid_entries, \
        invalid_value, pool_warn, ref_genome_warn


# -----------------------------------REPORT-------------------------------------


def print_validation_report(metafile, missing_mandatory_keys, invalid_keys,
                            invalid_values, invalid_value):
    """
    This function outputs a report on invalid files. The report contains the ID
     of the project, the path to the file, as well as the missing mandatory
     keys, invalid keys and invalid entries.
    :param invalid_value: a list containing invalid values
    :param metafile: the metafile that is validated
    :param missing_mandatory_keys: a list containing all missing mandatory keys
    :param invalid_keys: a list containing all invalid keys
    :param invalid_values: a list containing all invalid values
    """
    try:
        input_id = metafile['project']['id']
    except KeyError:
        input_id = 'missing'
    try:
        path = metafile['path']
    except KeyError:
        path = 'missing'
    invalid_entries = '\n- '.join(invalid_keys)
    missing = '\n- '.join(missing_mandatory_keys)
    whitelist_values = []
    for v in invalid_values:
        key = ':'.join(v.split(':')[:-1])
        entry = v.split(':')[-1]
        whitelist_values.append(entry + ' in ' + key + '\n')
    value = []
    for v in invalid_value:
        value.append(f'{v[0]}: {v[1]} -> {v[2]}')
    print(f'{"ERROR".center(80, "-")}\n'
          f'Project ID: {input_id}\n'
          f'Path: {path}\n\n'
          f'Report:\n')
    if len(invalid_keys) > 0:
        print(f'The following keys were invalid:\n'
              f'- {invalid_entries}\n')
    if len(missing_mandatory_keys) > 0:
        print(f'The following mandatory keys were missing:\n'
              f'- {missing}\n')
    if len(invalid_values) > 0:
        print(f'The following values do not match the whitelist:\n'
              f'- {"- ".join(whitelist_values)}')
    if len(invalid_value) > 0:
        print(f'The following values are invalid:\n'
              f'- {"- ".join(value)}')
    print(f'{"".center(80, "-")}')


def print_warning(metafile, pool_warn, ref_genome_warn):
    """
    This function prints a warning message.
    :param metafile: the metafile that contains the warning
    :param pool_warn: a list of warnings concerning pooled and donor_count
    :param ref_genome_warn: a list of warnings concerning the reference_genome
    """
    try:
        input_id = metafile['project']['id']
    except KeyError:
        input_id = 'missing'
    try:
        path = metafile['path']
    except KeyError:
        path = 'missing'
    print(f'{"WARNING".center(80, "-")}\n'
          f'Project ID: {input_id}\n'
          f'Path: {path}\n\n'
          f'Report:\n')
    if len(pool_warn) > 0:
        for elem in pool_warn:
            print(f'- Sample \'{elem[0]}\':\n{elem[1]}')
    if len(ref_genome_warn) > 0:
        for elem in ref_genome_warn:
            print(f'- Run from {elem[0]}:\n{elem[1]}')
    print(f'{"".center(80, "-")}')


# --------------------------------UTILITIES------------------------------------

def new_test(metafile, key_yaml, sub_lists, key_name, invalid_keys,
             invalid_entry, invalid_value, input_type, is_factor,
             local_factor, full_metadata):
    """
    This function test if all keys in the metadata file are valid.
    :param metafile: the metadata file
    :param key_yaml: the read in keys.yaml
    :param sub_lists: a list to save all items within a key if it has a list as
                      value
    :param key_name: the name of the key that is tested
    :param invalid_keys: a list of all invalid keys
    :param invalid_entry: a list of all invalid entries
    :param invalid_value: a list of all invalid values
    :param input_type: the input type that is expected for the value
    :param is_factor: a bool to state if the key is an experimental factor
    :param local_factor: a parameter to save the current experimental factor
    :return:
    invalid_keys: a list containing the invalid keys
    invalid_entries: a list containing the invalid entries
    invalid_value: a list containing the invalid values
    """
    if isinstance(metafile, dict) and not (
            'value' in metafile and 'unit' in metafile):
        for key in metafile:
            if not key_yaml and key_name.split(':')[-1] in is_factor or (key_name.split(':')[-1] == 'values' and local_factor is not None):
                new_yaml1 = utils.read_in_yaml(os.path.join(
                    os.path.dirname(os.path.abspath(__file__)), '..',
                    'keys.yaml'))
                if key_name.split(':')[-1] in is_factor:
                    new_yaml = list(utils.find_keys(new_yaml1, key_name.split(':')[-1]))
                else:
                    new_yaml = list(utils.find_keys(new_yaml1, local_factor))
                if len(new_yaml) > 0:
                    if 'whitelist' in new_yaml[0] and new_yaml[0]['whitelist']:
                        if key_name.split(':')[-1] in is_factor:
                            w = utils.get_whitelist(key_name.split(':')[-1],
                                                    full_metadata)
                        else:
                            w = utils.get_whitelist(local_factor, full_metadata)
                        if w and 'headers' in w:
                            if isinstance(w['headers'], dict):
                                if 'whitelist_keys' in w:
                                    headers = []
                                    for w_k in w['whitelist_keys']:
                                        if w_k in w['headers']:
                                            headers += w['headers'][w_k].split(' ')
                            else:
                                headers = w['headers'].split(' ')
                            new_yaml[0]['value'] = headers
                    if key not in new_yaml[0]['value']:
                        invalid_keys.append(key)
                    elif isinstance(metafile[key], list) != new_yaml[0]['list']:
                        invalid_keys.append(key)
                else:
                    invalid_keys.append(key)
            elif key not in key_yaml:
                invalid_keys.append(key)
            else:
                if key == 'factor':
                    global factor
                    factor = metafile[key]
                    local_factor = metafile[key]
                    is_factor.append(metafile[key])
                input_type = None
                if key == 'values' and factor is not None:
                    node = list(utils.find_keys(key_yaml, factor))
                    if len(node) > 0:
                        if 'input_type' in node:
                            input_type = node['input_type']
                elif isinstance(metafile[key], list) != key_yaml[key]['list']:
                    invalid_keys.append(key)
                elif 'input_type' in key_yaml[key]:
                    input_type = key_yaml[key]['input_type']
                res_keys, res_entries, res_values = new_test(
                    metafile[key], key_yaml[key]['value'], sub_lists,
                    f'{key_name}:{key}' if key_name != '' else key,
                    invalid_keys, invalid_entry, invalid_value, input_type,
                    is_factor, local_factor, full_metadata)
                invalid_keys = res_keys
    elif isinstance(metafile, list):
        for item in metafile:
            sub_lists.append(item)
            res_keys, res_entries, res_values = new_test(
                item, key_yaml, sub_lists, key_name, invalid_keys,
                invalid_entry, invalid_value, input_type, is_factor,
                local_factor, full_metadata)
            invalid_keys = res_keys
            sub_lists = sub_lists[:-1]
    else:
        invalid = new_test_for_whitelist(key_name.split(':')[-1], metafile,
                                         sub_lists)
        if invalid:
            invalid_entry.append(f'{key_name}:{metafile}')

        inv_value, message = validate_value(metafile, input_type,
                                            key_name.split(':')[-1])

        if not inv_value:
            invalid_value.append((key_name, metafile, message))

    return invalid_keys, invalid_entry, invalid_value


def new_test_for_whitelist(entry_key, entry_value, sublists):
    """
    This function tests if the value of a key matches the whitelist.
    :param entry_key: the key that is tested
    :param entry_value: the value that has to match the whitelist
    :param sublists: a list to save all items within a key if it has a list as
                      value
    :return: True if the entry does not match the whitelist else False
    """
    whitelist = utils.read_whitelist(entry_key)
    if whitelist and whitelist['whitelist_type'] == 'plain':
        whitelist = whitelist['whitelist']
    if isinstance(whitelist, dict):
        while isinstance(whitelist, dict) and whitelist[
                'whitelist_type'] == 'depend':
            whitelist_key = whitelist['ident_key']
            for i in reversed(range(len(sublists))):

                value = list(utils.find_keys(sublists[i], whitelist_key))

                if len(value) > 0:
                    if len(value) == 1:
                        break
                    else:
                        print("ERROR: multiple values")
                        break

            if value[0] in whitelist:
                whitelist = whitelist[value[0]]
            else:
                whitelist = utils.read_whitelist(value[0])
                if whitelist and whitelist['whitelist_type'] == 'plain':
                    whitelist = whitelist['whitelist']
        if isinstance(whitelist, dict) and whitelist[
                'whitelist_type'] == 'group':
            whitelist = utils.read_grouped_whitelist(whitelist, {})
            whitelist = [x for xs in list(whitelist['whitelist'].values())
                         if xs is not None for x in xs]
    if whitelist and not isinstance(whitelist, list) and not isinstance(
            whitelist, dict) \
            and os.path.isfile(os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '..',
            'metadata_whitelists', 'whitelists', whitelist)):
        whitelist = utils.read_whitelist(whitelist)
        if whitelist:
            whitelist = whitelist['whitelist']
    if whitelist and entry_value not in whitelist:
        return True
    return False


def test_for_mandatory(metafile, key_yaml, invalid_keys):
    """
    This function calls a function to get the missing mandatory keys for every 
    part of the metadata object.
    :param metafile: the read in metadata yaml file
    :param key_yaml: the read in structure file 'keys.yaml'
    :param invalid_keys: a list of keys that are invalid and should be ignored
    :return: missing_keys: a list containing the missing mandatory keys
    """
    missing_keys = []
    for key in key_yaml:
        missing_keys += get_missing_keys(key_yaml[key], metafile,
                                         invalid_keys, key, [], len(metafile[key]) if key in metafile and key_yaml[key]['list'] else 1)
    return missing_keys


def get_missing_keys(node, metafile, invalid_keys, pre, missing, list_len):
    """
    This function tests if all mandatory keys from the structure file
    'keys.yaml' are present in the metadata file.
    :param node: a key within the read in structure file 'keys.yaml'
    :param metafile: the read in metadata file
    :param invalid_keys: a list containing invalid keys that should be ignored
    :param pre: a string to save and chain keys in order to save their position
    :param missing: a list to save the missing mandatory keys
    """

    metafile = find_key(metafile, pre)

    if len(metafile) == 0:
        if node['mandatory']:
            missing.append(pre)
    else:
         if pre.split(':')[-1] not in invalid_keys:
            if isinstance(node['value'], dict) and not set(['mandatory', 'list',
                                                    'desc', 'display_name',
                                                    'value']) <= \
            set(node['value'].keys()):
                if isinstance(metafile[0], list):
                    for elem in metafile[0]:
                        for key in node['value']:
                            missing = get_missing_keys(node['value'][key], elem,
                                           invalid_keys, pre + ':' + key,
                                           missing, len(metafile[0]) if node['list'] else 1)
                else:
                    for key in node['value']:
                        missing = get_missing_keys(node['value'][key], metafile[0],
                                                   invalid_keys,
                                                   pre + ':' + key,
                                                   missing,
                                                   len(metafile[0]) if node[
                                                       'list'] else 1)
    return list(set(missing))


def find_key(metafile, key):
    """
    This function searches for a key of the keys.yaml in the metafile.
    :param metafile: the read in metadata file
    :param key: a string of chained keys (key1:key2...)
    :param is_list: bool, true if the instance in the structure is a list
    :param values: the default entry for the key within the structure file
    :param invalid_keys: a list containing invalid keys that should be ignored
    :return: missing_keys: a list containing the missing mandatory keys
    """
    for k in key.split(':'):
        new_metafile = list(utils.find_keys(metafile, k))
    return new_metafile


def validate_value(input_value, value_type, key):
    """
    This function tests if an entered value matches its type and contains
    invalid characters.
    :param input_value: the value to be valiated
    :param value_type: the type of which the value should be
    :param key: the key that contains the value
    :return:
    valid: a boolean that states if the value is valid
    message: a string that contains information about the error if tha value
             is invalid
    """
    valid = True
    message = None
    if input_value is not None:
        if value_type == 'bool':
            if input_value not in [True, False]:
                valid = False
                message = 'The value has to be of type bool (True or False).'
        elif value_type == 'number':
            if not isinstance(input_value, int):
                valid = False
                message = 'The value has to be an integer.'
        elif value_type == 'date':
            try:
                input_date = input_value.split('.')
                if len(input_date) != 3 or len(input_date[0]) != 2 or len(
                        input_date[1]) != 2 or len(input_date[2]) != 4:
                    raise SyntaxError
                input_value = datetime.date(int(input_date[2]),
                                            int(input_date[1]),
                                            int(input_date[0]))
            except (IndexError, ValueError, SyntaxError) as e:
                valid = False
                message = f'Input must be of type \'DD.MM.YYYY\'.'
        elif type(input_value) == str and ('\"' in input_value or '{' in input_value or '}' in
                 input_value or '|' in input_value) and key not in generated:
            valid = False
            message = 'The value contains an invalid character ' \
                      '(\", {, } or |).'
    return valid, message


def validate_logic(metafile):
    """
    This functions tests the logic of the input data.
    :param metafile: the metafile to be validated
    :return:
    pool_warn: a list containing warnings about the donor_count and pooled
    ref_genome_warn: a list containing warnings about the reference genome
    """
    pool_warn = []
    ref_genome_warn = []
    samples = list(utils.find_keys(metafile, 'samples'))
    for cond in samples:
        for sample in cond:
            warning, warn_message = validate_donor_count(sample['pooled'],
                                                         sample['donor_count'])
            if warning:
                pool_warn.append((sample['sample_name'], warn_message))
    organisms = list(utils.find_keys(metafile, 'organism_name'))
    runs = list(utils.find_keys(metafile, 'runs'))
    if len(runs) > 0:
        for run in runs[0]:
            if 'reference_genome' in run:
                warning, warn_message = validate_reference_genome(
                    organisms, run['reference_genome'])
                if warning:
                    ref_genome_warn.append((run['date'], warn_message))
    return pool_warn, ref_genome_warn


def validate_reference_genome(organisms, reference_genome):
    """
    This function tests if the reference genome matches the organism.
    :param organisms: a list of all organisms in the metadata file
    :param reference_genome: the reference genome that was specified
    :return:
    invalid: boolean to state if the reference genome is invalid
    message: a string explaining the logical error
    """
    invalid = False
    message = None
    ref_genome_whitelist = utils.get_whitelist(
        'reference_genome', None)['whitelist']
    if not any([reference_genome in ref_genome_whitelist[organism] for
                organism in organisms]):
        invalid = True
        organisms = [f'\'{organism}\'' for organism in organisms]
        message = (f'The reference genome \'{reference_genome}\' does not '
                   f'match the input organism ({", ".join(organisms)}).')
    return invalid, message


def validate_donor_count(pooled, donor_count):
    """
    This function tests if the donor_count matches the value stated for pooled.
    :param pooled: a boolean stating if a sample is pooled
    :param donor_count: the donor_count of a sample
    :return:
    invalid: boolean to state if the value is invalid
    message: a string explaining the logical error
    """
    invalid = False
    message = None
    if pooled and donor_count <= 1:
        invalid = True
        message = (f'Found donor count {donor_count} for pooled sample. '
                   f'The donor count should be greater than 1.')
    elif not pooled and donor_count > 1:
        invalid = True
        message = (f'Found donor count {donor_count} for sample that is not '
                   f'pooled. The donor count should be 1.')
    return invalid, message
