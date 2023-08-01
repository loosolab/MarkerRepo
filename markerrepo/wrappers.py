from .scoring import compare_marker_lists, update_scores
from .homology import check_organisms, download_homologene_data, fetch_homologs, get_dataset_names, get_supported_biomart_organisms, get_supported_taxonomy_ids, transfer_markers_biomart, transfer_markers_homologene
from .marker_repo import combine_lists, export_marker_list, guided_search, get_whitelists, select, get_selected_lists, search_df, combine_dfs, get_valid_filename
from .homology import get_biomart_defaults
from .utils import read_whitelist
from IPython.display import display


def create_marker_lists(organism, repo_path=".", style="score", path=".", file_name=None):
    """
    Creates marker lists for a given organism.

    Parameters
    ----------
    organism : str
        The organism of the marker lists.
    repo_path : str, default "."
        The path of the Marker Repo.
    style : str, default "score"
        The style of the marker lists. Currently there are four options available: "two_column", "score", "ui" and "panglao"
    path : str, default "."
        The path of the exported marker lists.
    file_name : str, default None
        The name of the exported marker lists.

    Returns
    -------
    List of paths to the created marker lists.

    """

    paths = []

    weighted = True if style == "score" or style == "ui" else False
    ui = True if style == "ui" else False
    custom_file_name = False if file_name else True

    while True:  
        df = search_df(df=combine_dfs(repo_path=repo_path), col_to_search="Organism name", search_terms=[f"+{organism}"])

        if df.empty:
            print("No marker lists found for this organism.")
            print("Trying to create marker lists via homology...")

            paths.extend(transfer_markers(target_org=organism, source_df=None, repo_path=repo_path, target_counts=1, 
                          weight_markers=weighted, export_suffix="annotation", ui=ui, custom_file_name=custom_file_name))
        else:
            print(f"Found {len(df)} marker lists for the given organism {organism}.")
            display(df)
            print(f"Please specify the marker lists you want to use for the annotation.")
            paths.append(convert_markers(style=style, repo_path=repo_path, df=guided_search(repo_path=repo_path, df=df, out="marker_list"), path=path, file_name=file_name))

        user_input = input("Do you want to add another marker list? (yes/no): ")
        if user_input.lower() != "yes":
            break 

    return paths


def convert_markers(repo_path=".", keywords=None, df=None, path="exported_lists", file_name=None, case_sensitive=False, exact=False, style="two_column", organism="Hs", tissue="all", gs=False, ensembl=False):
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
        Currently there are four options available: "two_column", "score", "ui" and "panglao"
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
        case "ui":
            print("Preparing ui style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list)
            marker_list = update_scores(marker_list, repo_path=repo_path)
        case "panglao":
            print("Preparing panglao style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list)
            marker_list = update_scores(marker_list, repo_path=repo_path)
            marker_list = transform_list_to_panglao(df=marker_list, organism=organism, tissue=tissue)
        case _:
            print("Style not recognized. Try 'two_column', 'score', 'ui' or 'panglao'")

    if path or file_name:
        # Export marker list
        return export_marker_list(marker_list, path=path, file_name=file_name)
    elif path and not file_name:
        return export_marker_list(marker_list, path=path, file_name=get_valid_filename(prompt="Enter file name: "))
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


def transfer_markers(target_org=None, source_df=None, repo_path=".", target_counts=1, weight_markers=False, export_suffix=None, ui=False, custom_file_name=False):
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
    weight_markers : bool, default False
        If True, a third column containing scores is added to the transferred marker list.
    export_suffix : str, default None
        A suffix that will be added to the file name of the transferred marker list.
    ui : bool, default False
        If True, try to update scores using the Panglao ubiquitousness index.
    custom_file_name : bool, default False
        If True, the user can specify a custom file name for the exported marker list.

    Returns
    -------
    list of str : 
        The paths of the exported transferred marker lists.
    """

    paths = []

    get_whitelists()

    biomart_organisms = get_supported_biomart_organisms(repo_path=repo_path)
    homologene_organisms = get_supported_taxonomy_ids(repo_path=repo_path)
    supported_organisms = list(set(biomart_organisms + homologene_organisms))

    if not target_org:
        target_org = select(whitelist=supported_organisms, heading="target organism:")
    target_organism, target_tax = target_org.split(" ")

    if target_org not in biomart_organisms:
        biomart_organisms = []
    if target_org not in homologene_organisms:
        homologene_organisms = []

    print(f"Loading genes of {target_organism}...")  
    target_genes = read_whitelist(f"genes/{target_organism}", repo_path=repo_path)['whitelist']
    print("Done!\n")

    if source_df is None:
        print("Select lists of source markers.")
        source_df = guided_search(out="metadata", repo_path=repo_path)
    
    unique_organisms = source_df[['Organism name', 'Taxonomy ID']].drop_duplicates()
    source_organisms = [' '.join(map(str, tup)) for tup in unique_organisms.values]

    in_both, in_neither, only_in_biomart, only_in_homologene = check_organisms(biomart_organisms, homologene_organisms, source_organisms, target_org)
    print("")
    if in_both:
        print("The following source organisms can be used in both approaches: " + ', '.join(in_both) + ".")
    if in_neither:
        print("The following source organisms can't be used in either approach: " + ', '.join(in_neither) + ".")
    if only_in_biomart:
        print("The following source organisms can only be used in the BioMart approach: " + ', '.join(only_in_biomart) + ".")
    if only_in_homologene:
        print("The following source organisms can only be used in the HomoloGene approach: " + ', '.join(only_in_homologene) + ".")

    biomart_source_organisms = in_both + only_in_biomart
    homologene_source_organisms = in_both + only_in_homologene

    if len(homologene_source_organisms) > 0 and target_org in homologene_organisms:
        print("\nStarting HomoloGene approach...")
        for source_organism in homologene_source_organisms:
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
            
            if not filtered_source_df.empty:
                print(f"Create DataFrame containing the transferred genes based on the filter criteria.")
                transferred_list = transfer_markers_homologene(filtered_source_df, source_tax, target_tax, hg_db, target_genes,
                                                        source_whitelist=source_genes, calc_proportions=True, plots=True)
                print("Transferred markers:")
                display(transferred_list)

                if weight_markers:                   
                    if target_org in biomart_organisms and ui:
                        default_org = next((x for x in get_biomart_defaults(repo_path=repo_path) if x in get_dataset_names(target_organism)), None)
                        if default_org:
                            transferred_list = update_scores(df=transferred_list, repo_path=repo_path, organism=target_organism, biomart_target=default_org)
                        else:
                            transferred_list = update_scores(df=transferred_list, repo_path=repo_path, organism=target_organism)
                    transferred_list = compare_marker_lists(marker_df=transferred_list)
                    print("Weighted transferred markers:")
                    display(transferred_list)

                if not custom_file_name:
                    file_name=f"{source_organism}_{target_organism}_HomoloGene"
                    if export_suffix:
                        file_name = f"{file_name}_{export_suffix}"
                else:
                    file_name = get_valid_filename(prompt="Enter file name: ")

                paths.append(export_marker_list(transferred_list, path="./transferred_markers", file_name=file_name, marker_id="symbol"))
            else:
                print("Source DataFrame is empty.")

    if len(biomart_source_organisms) > 0 and target_org in biomart_organisms:
        print("\nStarting BioMart approach...")
        print("\nSpecify BioMart target organism selection:")
        default_org = next((x for x in get_biomart_defaults(repo_path=repo_path) if x in get_dataset_names(target_organism)), None)
        if default_org:
            print(f"Default organism: {default_org}")
            target_organism_bm = default_org
        else:
            target_organism_bm = select(whitelist=get_dataset_names(target_organism), heading="BioMart target organism", repo_path=repo_path)
        for source_organism in biomart_source_organisms:
            source_organism = source_organism.split(" ")[0]
            print("\nSpecify BioMart source organism selection:")
            default_org = next((x for x in get_biomart_defaults(repo_path=repo_path) if x in get_dataset_names(source_organism)), None)
            if default_org:
                print(f"Default organism: {default_org}")
                source_organism_bm = default_org
            else:
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
            
            if not filtered_source_df.empty:
                print(f"Create DataFrame containing the transferred genes based on the filter criteria.")
                transferred_list = transfer_markers_biomart(biomart_db, filtered_source_df, target_genes, source_whitelist=source_genes,
                                                            calc_proportions=True, plots=True, target_counts=None)
                print("Transferred markers:")
                display(transferred_list)

                if weight_markers:
                    if target_org in biomart_organisms and ui:
                        default_org = next((x for x in get_biomart_defaults(repo_path=repo_path) if x in get_dataset_names(target_organism)), None)
                        if default_org:
                            transferred_list = update_scores(df=transferred_list, repo_path=repo_path, organism=target_organism, biomart_target=default_org)
                        else:
                            transferred_list = update_scores(df=transferred_list, repo_path=repo_path, organism=target_organism, biomart_target=target_organism_bm)
                    transferred_list = compare_marker_lists(marker_df=transferred_list)
                    print("Weighted transferred markers:")
                    display(transferred_list)

                if not custom_file_name:
                    file_name=f"{source_organism}_{target_organism}_BioMart"
                    if export_suffix:
                        file_name = f"{file_name}_{export_suffix}"
                else:
                    file_name = get_valid_filename(prompt="Enter file name: ")

                paths.append(export_marker_list(transferred_list, path="./transferred_markers", file_name=file_name, marker_id="symbol"))
            else:
                print("Source DataFrame is empty.")
                
    print("\nFinished!")

    return paths

