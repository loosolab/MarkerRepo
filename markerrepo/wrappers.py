from .marker_repo import search_df, combine_lists, export_marker_list, guided_search, combine_dfs, get_whitelists, select
from .calculations import compare_marker_lists, update_scores, get_supported_biomart_organisms, get_supported_taxonomy_ids, get_dataset_names, fetch_homologs, download_homologene_data, transfer_markers_biomart, transfer_markers_homologene
from .utils import read_whitelist

def get_selected_lists(keywords=None, metadata_df=None, repo_path=".", case_sensitive=False, exact=False):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.

    Parameters
    ----------
    keywords : dict or str, default None
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    repo_path : str, default "."
        The path of the Marker Repo.
    metadata_df : pd.DataFrame, default None
        A DataFrame containing the metadata of a selection of marker lists.
    case_sensitive : bool, default False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.

    Returns
    --------
    pd.DataFrame :
        The DataFrame conaining the combined lists.
    """

    if keywords:
        db = combine_dfs(repo_path=repo_path)
        df = search_df(db, keywords, case_sensitive=case_sensitive, exact=exact)
        if df.empty:
            raise Exception(
                        f"No search results available!")
    elif metadata_df is not None:
        df = metadata_df
    else:
        raise Exception(
                    f"You need to specify keywords or a metadata DataFrame!")


    # Get UIDs and combine lists
    uids = [int(idx) for idx in df.index]
    combined_df = combine_lists(uids, repo_path=repo_path)

    # Drop duplicates, keep one marker only, rearrange column order
    markers_filtered = combined_df.drop_duplicates()
    markers_filtered = markers_filtered[['Info', 'Marker']]

    return markers_filtered


def convert_markers(repo_path=".", keywords=None, df=None, path="exported_lists", file_name="marker_list", case_sensitive=False, exact=False, style="two_column", organism="Hs", tissue="all", gs=False, ensembl=False):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.
    Optionally, it can export the DataFrame to a file.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    keywords : dict or str, default None
        The keywords to filter the DataFrame. Can be either a dictionary with column names as keys and
        keywords as values, or a single string to search for in the entire DataFrame.
    df : pd.DataFrame, default None
        A DataFrame containing a two column marker list.
    path : str, default None
        The path to the file where the new marker list will be saved. If not specified, the function will not save the DataFrame to a file.
    file_name : str, default "marker_list"
        The filename of the marker list.
    case_sensitive : bool, default False
        If True, the function will consider the case of the keywords. If False, the function will ignore the case.
    exact : bool, default False
        If True, the function will search for exact matches of the keywords. If False, the function will search for the keywords as substrings.
    style : str, default "two_column"
        The format style of which the exported marker list should look like.
        Currently there are three options available: "two_column", "score" and "panglao"
    organism : str, default "Hs"
        Organism of panglao style markers.
    tissue : str, default "all"
        Tissue of panglao style markers.
    gs : bool, default False
        If true, guided search is enabled and the resulting DataFrame will be used.
    ensembl : bool, default False
        If True, the Ensembl IDs will be used instead of the gene symbols.

    Returns
    --------
    pd.DataFrame or str :
        If path is specified, the function returns the absolute path to the file where the marker list was saved.
        If path is not specified, the function returns the DataFrame.
    """

    if gs:
        marker_list = guided_search(repo_path=repo_path, out="marker_list")
    elif repo_path and keywords:
        marker_list = get_selected_lists(keywords, repo_path=repo_path, case_sensitive=case_sensitive, exact=exact)
    elif df is not None:
        marker_list = df
    else:
        raise Exception(
            f"You need to specify keywords, a marker list DataFrame or use the guided search ('gs=True') !")

    if ensembl:
        marker_list['Marker'] = marker_list['Marker'].apply(lambda x: x.split(' ')[1] if len(x.split(' ')) > 1 else x)
    else:
        marker_list['Marker'] = marker_list['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x)
    
    match style:
        case "two_column":
            print("Preparing two column style marker list...")
            
        case "score":
            print("Preparing score style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list)
            
        case "panglao":
            print("Preparing panglao style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list)
            marker_list = update_scores(marker_list)
            marker_list = transform_list_to_panglao(df=marker_list, organism=organism, tissue=tissue)
        case _:
            print("Style not recognized. Try 'two_column', 'score' or 'panglao'")


    if path or file_name:
        # Export marker list
        return export_marker_list(marker_list, path=path, file_name=file_name)
    else:
        return marker_list


def transform_list_to_panglao(df, organism="Hs", tissue="all"):
    """
    Adds additional columns to the given DataFrame and rearranges column order.
    This converts the format of the marker list to the format of the Panglao db.
    
    Parameters
    ----------
    df : DataFrame
        Input DataFrame with columns "Info", "Marker", and "Score".
    organism : str, default "Hs"
        Organism of panglao style markers.
    tissue : str, default "all"
        Tissue of panglao style markers.
        
    Returns
    -------
    pd.DataFrame :
        Updated DataFrame with new columns, specified column order, and updated Marker column.
    """

    # Update Marker column to retain only the first marker
    df['Marker'] = df['Marker'].apply(lambda x: x.split(' ')[0] if len(x.split(' ')) > 1 else x)
    
    # Add new columns
    df['Organism'] = organism
    df['Aliases'] = df['Marker']
    df['Tissue'] = tissue
    
    # Rearrange column order
    df = df[['Organism', 'Marker', 'Info', 'Aliases', 'Score', 'Tissue']]
    
    return df


def transfer_markers(target_org=None, source_df=None, repo_path=".", target_counts=1, weight_markers=True):
    """
    Performs all steps of transferring marker genes from source organism(s)
    to one target organism.
    
    Parameters
    ----------
    target_org : str, default None
        The target organism to which the genes of the source organism(s) are to be transferred. 
        The string must contain the organism name and the taxonomy id. If None, select a target organism.
        Example: "human 9606"
    source_df : DataFrame, default None
        Metadata DataFrame of source information.
    repo_path : str, default "."
        The path of the Marker Repo.
    target_counts : int, default 1
        If not None, filter those target genes whose number of target genes per source gene is <= target_counts.

    Returns
    -------
    pd.DataFrame :
        Updated DataFrame with new columns, specified column order, and updated Marker column.
    """

    get_whitelists()

    if not target_org:
        target_org = select(key="organism", heading="target organism:")
    target_organism, target_tax = target_org.split(" ")

    print(f"Loading genes of {target_organism}...")  
    target_genes = read_whitelist(f"genes/{target_organism}", repo_path=repo_path)['whitelist']
    print("Done!\n")

    if source_df is None:
        source_df = guided_search(out="metadata")
    
    unique_organisms = source_df[['Organism name', 'Taxonomy ID']].drop_duplicates()

    source_organisms = [' '.join(map(str, tup)) for tup in unique_organisms.values]
    biomart_organisms = get_supported_biomart_organisms(repo_path=repo_path)
    homologene_organisms = get_supported_taxonomy_ids(repo_path=repo_path)

    in_both, in_neither, only_in_biomart, only_in_homologene = check_organisms(biomart_organisms, homologene_organisms, source_organisms)
    biomart_organisms = in_both + only_in_biomart

    print("")
    if target_org in biomart_organisms:
        biomart_organisms.remove(target_org)
        print(f"Removed {target_org} from BioMart source organisms as it matches the target organism.")

    homologene_organisms = in_both + only_in_homologene
    if target_org in homologene_organisms:
        homologene_organisms.remove(target_org)
        print(f"Removed {target_org} from HomoloGene source organisms as it matches the target organism.")

    if len(homologene_organisms) > 0:
        print("\nStarting HomoloGene approach...")
        for source_organism in homologene_organisms:
            source_organism, source_tax = source_organism.split(" ")
            print(f"Loading genes of {source_organism}...")
            source_genes = read_whitelist(f"genes/{source_organism}", repo_path=repo_path)['whitelist']
            print("Done!\n")

            print("Get HomoloGene db...")
            hg_db = download_homologene_data()

            uids = source_df.loc[source_df['Organism name'] == source_organism].index.tolist()
            source_marker_list = combine_lists(uids, repo_path=repo_path)
            print(f"\nDataFrame of the markers of the source organism ({source_organism}) to be transferred to the target organism ({target_organism}):")
            display(source_marker_list)

            if target_counts:
                print(f"Filter source DataFrame by the number of target genes per source gene: remove all source genes that lead to more than {target_counts} target genes.")
                filtered_source_df = transfer_markers_homologene(source_marker_list, source_tax, target_tax, hg_db, target_genes, source_whitelist=source_genes, calc_proportions=True, 
                                                        plots=True, target_counts=target_counts)
            else:
                filtered_source_df = source_marker_list
            
            print(f"Create DataFrame containing the transferred genes based on the filter criteria.")
            transferred_list = transfer_markers_homologene(filtered_source_df, source_tax, target_tax, hg_db, target_genes,
                                                    source_whitelist=source_genes, calc_proportions=True, plots=True)
            print("Transferred markers:")
            display(transferred_list)

            if weight_markers:
                results_scored = compare_marker_lists(marker_df=transferred_list)
                transferred_list = update_scores(df=results_scored, repo_path=repo_path)

            export_marker_list(transferred_list, path="./transferred_markers", file_name=f"{source_organism}_{target_organism}_HomoloGene", marker_id="symbol")

    if len(biomart_organisms) > 0:
        print("\nStarting BioMart approach...")
        print("\nSpecify BioMart organism selection:")
        target_organism_bm = select(whitelist=get_dataset_names(target_organism), heading="BioMart target organism", repo_path=repo_path)
        for source_organism in biomart_organisms:
            source_organism = source_organism.split(" ")[0]
            source_organism_bm = select(whitelist=get_dataset_names(source_organism), heading="BioMart source organism", repo_path=repo_path)
            print(f"Loading genes of {source_organism}...")
            source_genes = read_whitelist(f"genes/{source_organism}", repo_path=repo_path)['whitelist']
            print("Done!\n")

            print("Fetch necessary data from BioMart...")
            biomart_db = fetch_homologs(source_organism_bm, target_organism_bm).dropna()

            uids = source_df.loc[source_df['Organism name'] == source_organism].index.tolist()
            source_marker_list = combine_lists(uids, repo_path=repo_path)
            print(f"\nDataFrame of the markers of the source organism ({source_organism}) to be transferred to the target organism ({target_organism}):")
            display(source_marker_list)

            if target_counts:
                print(f"Filter source DataFrame by the number of target genes per source gene: remove all source genes that lead to more than {target_counts} target genes.")
                filtered_source_df = transfer_markers_biomart(biomart_db, source_marker_list, target_genes, source_whitelist=source_genes,
                                                            calc_proportions=True, plots=True, target_counts=target_counts)
            else:
                filtered_source_df = source_marker_list
            
            print(f"Create DataFrame containing the transferred genes based on the filter criteria.")
            transferred_list = transfer_markers_biomart(biomart_db, filtered_source_df, target_genes, source_whitelist=source_genes,
                                                        calc_proportions=True, plots=True, target_counts=None)
            print("Transferred markers:")
            display(transferred_list)

            if weight_markers:
                results_scored = compare_marker_lists(marker_df=transferred_list)
                transferred_list = update_scores(df=results_scored, repo_path=repo_path)

            export_marker_list(transferred_list, path="./transferred_markers", file_name=f"{source_organism}_{target_organism}_BioMart", marker_id="symbol")
    

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
