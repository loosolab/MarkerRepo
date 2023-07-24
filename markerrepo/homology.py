def check_organisms(biomart_orgs, homologene_orgs, source_organisms):
    """
    Checks if the source organisms are supported by the BioMart or HomoloGene approach.

    Parameters
    ----------
    biomart_orgs : list
        List of organisms supported by the BioMart approach.
    homologene_orgs : list
        List of organisms supported by the HomoloGene approach.
    source_organisms : list
        List of source organisms that the user wants to use for gene transfer.

    Returns
    -------
    lists of str :
        - Organisms available in both approaches
        - Organisms available in neither approach
        - Organisms available only in the BioMart approach
        - Organisms available only in the HomoloGene approach
    """

    biomart_set = set(biomart_orgs)
    homologene_set = set(homologene_orgs)
    source_set = set(source_organisms)

    in_both = list(source_set.intersection(biomart_set).intersection(homologene_set))
    in_neither = list(source_set.difference(biomart_set).difference(homologene_set))
    only_in_biomart = list(source_set.intersection(biomart_set).difference(homologene_set))
    only_in_homologene = list(source_set.intersection(homologene_set).difference(biomart_set))

    if in_both:
        print("The following organisms can be used in both approaches: " + ', '.join(in_both) + ".")
    if in_neither:
        print("The following organisms can't be used in either approach: " + ', '.join(in_neither) + ".")
    if only_in_biomart:
        print("The following organisms can only be used in the BioMart approach: " + ', '.join(only_in_biomart) + ".")
    if only_in_homologene:
        print("The following organisms can only be used in the HomoloGene approach: " + ', '.join(only_in_homologene) + ".")

    return in_both, in_neither, only_in_biomart, only_in_homologene