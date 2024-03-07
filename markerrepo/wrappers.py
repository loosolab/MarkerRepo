from .scoring import compare_marker_lists, update_scores
from .homology import check_organisms, download_homologene_data, fetch_homologs, get_dataset_names, get_supported_biomart_organisms, get_supported_taxonomy_ids, transfer_markers_biomart, transfer_markers_homologene
from .marker_repo import combine_lists, export_marker_list, guided_search, select, get_selected_lists, search_df, combine_dfs, get_valid_filename, update_organism, transform_marker_list, get_uid
from .homology import get_biomart_defaults
from .utils import read_whitelist, get_whitelists
from .annotation import annot_ct, show_tables, reformat_marker_list, compare_cell_types, update_adata_with_markers
from .generate_metafile import generate_file
import scanpy as sc
from IPython.display import display
import os
import sys
import contextlib
import io
import logging
import warnings
import inspect
import numpy as np
import episcanpy as epi

try:
    from sctoolbox.tools import celltype_annotation
except ModuleNotFoundError:
    warnings.warn("Please install the latest sctoolbox version. Some functionality may not be available.", RuntimeWarning)


@contextlib.contextmanager
def suppress_logging(logger_name, level=logging.CRITICAL):
    logger = logging.getLogger(logger_name)
    old_level = logger.getEffectiveLevel()
    logger.setLevel(level)
    try:
        yield
    finally:
        logger.setLevel(old_level)


@contextlib.contextmanager
def suppress_output():
    new_stdout, new_stderr = io.StringIO(), io.StringIO()
    old_stdout, old_stderr = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_stdout, new_stderr
        yield
    finally:
        sys.stdout, sys.stderr = old_stdout, old_stderr


def create_marker_lists(organism=None, repo_path=".", lists_path=None, style="two_column", path=".", file_name=None, ensembl=False, 
                        col_to_search=None, search_terms=None, force_homology=False, show_lists=True,
                        column_specific_terms=None, adata=None, suffix=None):
    """
    Creates marker lists for a given organism.

    Parameters
    ----------
    organism : str, default None
        The organism of the marker lists.
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    style : str, default "two_column"
        The style of the marker lists. Currently there are four options available: "two_column", "score", "ui" and "panglao"
    path : str, default "."
        The path of the exported marker lists.
    file_name : str, default None
        The name of the exported marker lists.
    ensembl : bool, default False
        If True, the Ensembl IDs will be used instead of the gene symbols.
    col_to_search : str, default None
        The column of the DataFrame to search in.
    search_terms : list of str, default None
        The search terms to search for in the DataFrame.
    force_homology : bool, default False
        If True, the function will try to create marker lists via homology even if marker lists for the given organism already exist.
    show_lists : bool, default True
        If True, the function will show the marker lists of the query.
    column_specific_terms : dict, default None
        A dictionary with column names as keys and lists of search terms as values. If provided, 'col_to_search' and 'search_terms' are ignored.
    adata : AnnData, default None
        If provided, the function will add the marker list IDs to the .uns table of the AnnData object.
    suffix : str, default None
        The key of the metadata section whose value should be appended to the marker names.
        
    Note:
    ----
    It is recommended to use 'column_specific_terms' instead of 'col_to_search' and 'search_terms'.
    The parameters 'col_to_search' and 'search_terms' are planned to be deprecated in future releases.
        
    Returns
    -------
    List of str :
    Paths to the created marker lists.
    """

    if not os.path.exists(repo_path):
        raise FileNotFoundError(f"The specified repository path '{repo_path}' does not exist.")

    if style not in ["two_column", "score", "ui", "panglao"]:
        raise ValueError("The parameter 'style' must be one of 'two_column', 'score', 'ui', or 'panglao'.")

    if col_to_search is not None or search_terms is not None:
        print("Warning: 'col_to_search' and 'search_terms' are deprecated and will be removed in future versions.")
        print("It is recommended to use 'column_specific_terms' for more robust functionality.")
    
    # TODO: check also whether organism, col_to_search and keys (columns) of column_specific_terms are valid

    paths = []

    weighted = True if style == "score" or style == "ui" else False
    ui = True if style == "ui" else False
    custom_file_name = False if file_name else True

    while True:
        if organism:  
            df = search_df(df=combine_dfs(repo_path=repo_path, lists_path=lists_path), col_to_search="Organism name", search_terms=[f"+{organism.split(' ')[0]}"], suffix=suffix)
        else:
            df = combine_dfs(repo_path=repo_path, lists_path=lists_path)

        if df.empty or force_homology:
            if not force_homology:
                print("No marker lists found for this organism.")
            print("Trying to create marker lists via homology...")

            source_df = None
            if search_terms or column_specific_terms:
                source_df = search_df(df=combine_dfs(repo_path=repo_path, lists_path=lists_path), col_to_search=col_to_search, search_terms=search_terms, column_specific_terms=column_specific_terms, suffix=suffix)
            paths.extend(transfer_markers(target_org=organism, source_df=source_df, repo_path=repo_path, lists_path=lists_path, target_counts=1, 
                          weight_markers=weighted, export_suffix="annotation", ui=ui, custom_file_name=custom_file_name, ensemble=ensembl, suffix=suffix))
        else:
            if organism:
                print(f"Found {len(df)} marker lists for the given organism {organism.split(' ')[0]}.")
                
                if show_lists:
                    display(df)
            if search_terms or column_specific_terms:
                df = search_df(df=df, col_to_search=col_to_search, search_terms=search_terms, column_specific_terms=column_specific_terms, suffix=suffix)
                print(f"Found {len(df)} marker lists for the given search terms.")
                
                if show_lists:
                    display(df)
                df_combined = get_selected_lists(metadata_df=df, repo_path=repo_path, lists_path=lists_path, order=["Marker", "Info"], suffix=suffix)
                paths.append(convert_markers(style=style, repo_path=repo_path, lists_path=lists_path, df=df_combined, path=path, file_name=file_name, ensembl=ensembl, suffix=suffix))
            else:
                print(f"Please specify the marker lists you want to use for the annotation.")
                df = guided_search(repo_path=repo_path, lists_path=lists_path, df=df, out="metadata")
                df_combined = get_selected_lists(metadata_df=df, repo_path=repo_path, lists_path=lists_path, order=["Marker", "Info"], suffix=suffix)
                paths.append(convert_markers(style=style, repo_path=repo_path, lists_path=lists_path, df=df_combined, path=path, file_name=file_name, ensembl=ensembl, suffix=suffix))

        if adata and file_name:
            update_adata_with_markers(adata, file_name, df)

        if search_terms or column_specific_terms:
            return paths
        
        user_input = input("Do you want to add another marker list? (yes/no): ")
        if user_input.lower() != "yes":
            break 

    return paths


def create_multiple_marker_lists(cml_parameters=[{}], repo_path=".", lists_path=None, organism=None, style='two_column', path='.', file_name=None, ensembl=False, 
                                 col_to_search=None, search_terms=None, force_homology=False, show_lists=True, 
                                 column_specific_terms=None, adata=None, suffix=None):
    """
    Executes the 'create_marker_lists' function for multiple sets of parameters.

    Each set of parameters is provided as a dictionary within a list. This function
    iterates over each dictionary, using its contents to call 'create_marker_lists'.
    Default values are used for any missing parameters, which can be overridden by individual dictionaries.

    Parameters
    ----------
    cml_parameters : list of dict, default [{}]
        A list of dictionaries, where each dictionary contains parameters for a single
        call to 'create_marker_lists'. Keys in the dictionaries should match the parameter
        names of 'create_marker_lists', and values should be the desired values for those parameters.
    repo_path : str, default "."
        The path of the Marker Repo. This value is passed directly to 'create_marker_lists'.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    organism : str, default None
        Default organism to use for all 'create_marker_lists' calls.
    style : str, default 'two_column'
        Default style for all 'create_marker_lists' calls.
    path : str, default '.'
        Default path for all 'create_marker_lists' calls.
    file_name : str, default None
        Default file name for all 'create_marker_lists' calls.
    ensembl : bool, default False
        Default ensembl flag for all 'create_marker_lists' calls.
    col_to_search : str, default None
        Default column to search for all 'create_marker_lists' calls.
    search_terms : list of str, default None
        Default search terms for all 'create_marker_lists' calls.
    force_homology : bool, default False
        Default force homology flag for all 'create_marker_lists' calls.
    show_lists : bool, default True
        Default show lists flag for all 'create_marker_lists' calls.
    column_specific_terms : dict, default None
        Default column specific terms for all 'create_marker_lists' calls.
    adata : AnnData, default None
        Default AnnData object for all 'create_marker_lists' calls.
    suffix : str, default None
        Default suffix for all 'create_marker_lists' calls.

    Note:
    ----
    It is recommended to use 'column_specific_terms' instead of 'col_to_search' and 'search_terms'.
    The parameters 'col_to_search' and 'search_terms' are planned to be deprecated in future releases.

    Returns
    -------
    list of str :
        A combined list of all paths to the created marker lists from each call to 'create_marker_lists'.
    """

    if not cml_parameters:
        cml_parameters = [{}]

    if col_to_search is not None or search_terms is not None:
        print("Warning: 'col_to_search' and 'search_terms' are deprecated and will be removed in future versions.")
        print("It is recommended to use 'column_specific_terms' for more robust functionality.")

    all_paths = []

    for setting in cml_parameters:
        # Set default values for parameters of 'create_marker_lists'
        params = {
            'organism': organism,
            'style': style,
            'path': path,
            'file_name': file_name,
            'ensembl': ensembl,
            'col_to_search': col_to_search,
            'search_terms': search_terms,
            'force_homology': force_homology,
            'show_lists': show_lists,
            'column_specific_terms': column_specific_terms,
            'adata': adata,
            'repo_path': repo_path,
            'lists_path': lists_path,
            'suffix': suffix
        }

        # Update these defaults with values from the current dictionary
        params.update(setting)

        # Call 'create_marker_lists' with updated parameters
        paths = create_marker_lists(**params)
        all_paths.extend(paths)

    return all_paths


def run_annotation(adata, marker_repo=True, SCSA=True, marker_lists=None, mr_obs="mr", scsa_obs="scsa", 
                   rank_genes_column=None, clustering_column=None, reference_obs=None, keep_all=False, 
                   verbose=False, show_ct_tables=False, show_plots=False, show_comparison=False, ignore_overwrite=False,
                   celltype_column_name=None, omic="RNA"):
    """
    Performs annotations on single cell data and allows the user to choose between different annotation methods. 

    Parameters
    ----------
    adata : AnnData
        The anndata object to annotate.
    marker_repo : bool, default True
        Whether to use Marker Repo annotation.
    SCSA : bool, default True
        Whether to use SCSA annotation.
    marker_lists : list of str, default []
        Paths to marker list files.
    mr_obs : str, default "mr"
        .obs key for Marker Repo annotation.
    scsa_obs : str, default "scsa"
        .obs key for SCSA annotation.
    rank_genes_column : str, default None
        The column of the .uns table which contains the rank genes scores. E.g. "rank_genes_groups". 
        If None, the ranking will be performed on the clustering_column.
    clustering_column : str, default None
        The column of the .obs table which contains the clustering information. E.g. "louvain" or "leiden".
    reference_obs : str, default None
        A reference annotation already present in the .obs table that can be compared with the other annotations.
    keep_all : bool, default False
        If True, all annotation columns will be kept. If False, only the selected annotation column and the reference_obs will be kept.
    verbose : bool, default False
        If True, the function will print additional information.
    show_ct_tables : bool, default False
        If True, the function will show the tables of the annotation.
    show_plots : bool, default False
        If True, the function will show the plots of the annotation.
    show_comparison : bool, default False
        If True, the function will show the comparison of the annotations.
    ignore_overwrite : bool, default False
        If True, the function will not ask for confirmation before overwriting existing files.
    celltype_column_name : str, default None
        The name of the selected cell type annotation column. If None, all annotation columns will be kept.
    omic : str, default "RNA"
        The omic type of the AnnData object. E.g. "RNA", "ATAC"
    """

    if not marker_repo and not SCSA:
        raise ValueError("At least one of 'marker_repo' or 'SCSA' must be True.")

    if marker_lists is None or not marker_lists:
        raise ValueError("No marker lists provided. Please provide a list of marker list paths.")

    if isinstance(marker_lists, str):
        marker_lists = [marker_lists]
    
    for marker_list in marker_lists:
        if not os.path.exists(marker_list):
            raise FileNotFoundError(f"Marker list file not found: {marker_list}")
        
    if not clustering_column:
        clustering_column = select(whitelist=list(adata.obs.columns), heading="clustering column")

    if clustering_column not in adata.obs:
        raise ValueError(f"Clustering column '{clustering_column}' not found in adata.obs.")

    if rank_genes_column is not None and rank_genes_column not in adata.uns:
        raise ValueError(f"Rank genes column '{rank_genes_column}' not found in adata.uns.")

    if reference_obs is not None and reference_obs not in adata.obs:
        raise ValueError(f"Reference annotation column '{reference_obs}' not found in adata.obs.")
    
    if not rank_genes_column:
        rank_genes_column = rank_feature_groups(adata, clustering_column, show_plots=show_plots, verbose=verbose, omic=omic)

    annotation_columns = [] if reference_obs is None else [reference_obs]        

    if show_plots:
        if 'X_umap' in adata.obsm:
            sc.pl.umap(adata, color=clustering_column, wspace=0.5, cmap=None)

    for marker_list in marker_lists:
        name = marker_list.split('/')[-1]
        annotation_dir = f"./annotation/{name}"
        plot_columns = [] if reference_obs is None else [reference_obs]

        if marker_repo:
            ct_column = f"{mr_obs}_{name}_{clustering_column}"
            annotation_columns.append(ct_column)
            plot_columns.append(ct_column)
            
            # Execute Marker Repo annotation
            annot_ct(adata, output_path=annotation_dir, db_path=marker_list,
                           cluster_column=clustering_column, rank_genes_column=rank_genes_column, 
                           ct_column=ct_column, verbose=verbose, ignore_overwrite=ignore_overwrite)

            # Show tables and alternative cell types of each cluster
            if show_ct_tables:
                print(f"Tables of cell type annotation with clustering {clustering_column} and marker list {name}:")
                show_tables(annotation_dir=annotation_dir, n=5, clustering_column=clustering_column, show_diff=True)

        if SCSA:
            column_added = f"{scsa_obs}_{name}_{clustering_column}"
            annotation_columns.append(column_added)
            plot_columns.append(column_added)
            original_var_index = adata.var.index.copy()
            adata.var.index = adata.var.index.map(lambda x: x.upper())

            # Execute SCSA annotation
            if verbose:
                celltype_annotation.run_scsa(adata, 
                    gene_column=None, 
                    key=rank_genes_column, 
                    column_added=column_added,
                    inplace=True, 
                    species=None, 
                    fc=1.5, 
                    pvalue=0.05, 
                    user_db=reformat_marker_list(marker_list), 
                    celltype_column="cell_name")
            else:
                with suppress_logging(logger_name='sctoolbox'):
                    celltype_annotation.run_scsa(adata, 
                        gene_column=None, 
                        key=rank_genes_column, 
                        column_added=column_added,
                        inplace=True, 
                        species=None, 
                        fc=1.5, 
                        pvalue=0.05, 
                        user_db=reformat_marker_list(marker_list), 
                        celltype_column="cell_name")
                    
            adata.var.index = original_var_index

        # Show plots
        if show_plots:
            if 'X_umap' in adata.obsm:
                sc.pl.umap(adata, color=plot_columns, wspace=0.5, cmap=None)
            else:
                print("UMAP embedding not found in the AnnData object.")


    # Compare annotations
    if show_comparison:
        print("Comparison of cell type annotations:")
        display(compare_cell_types(adata, clustering_column, annotation_columns))

    # Select cell type annotation
    if celltype_column_name:
        annotation_column = select(whitelist=annotation_columns, heading="Select cell type annotation column:")
        adata.obs.rename(columns={annotation_column: celltype_column_name}, inplace=True)
    else:
        keep_all = True

    if not keep_all:
        # Keep only the selected annotation column and reference_obs if provided
        columns_to_keep = [annotation_column]
        if reference_obs is not None and reference_obs in adata.obs.columns:
            columns_to_keep.append(reference_obs)

        columns_to_remove = [col for col in annotation_columns if col not in columns_to_keep]
        adata.obs.drop(columns=columns_to_remove, inplace=True)


def rank_feature_groups(adata, clustering_column=None, show_plots=False, verbose=False, n=None, omic="RNA"):
    """
    Rank features for each cluster in the given AnnData object.

    Parameters
    ----------
    adata : AnnData
        The anndata object to rank.
    clustering_column : str, default None
        The column of the .obs table which contains the clustering information. E.g. "louvain" or "leiden".
    show_plots : bool, default False
        If True, the function will show the plots of the ranking.
    verbose : bool, default False
        If True, the function will print additional information.
    n : int, default None
        The number of features that appear in the returned table. If None, all features will be returned.

    Returns
    -------
    str :
        The name of the column in the .uns table which contains the features and scores.
    """

    if not n:
        # Use 10 percent of the features of the adata
        n = int(adata.n_vars * 0.1)

    adata.uns['omic'] = omic

    if not clustering_column:
        clustering_column = select(whitelist=list(adata.obs.columns), heading="clustering column")
    
    sc.tl.dendrogram(adata, groupby=clustering_column)

    if 'log1p' in adata.uns:
        adata.uns['log1p']['base'] = np.e
    else:
        adata.uns['log1p'] = {'base': np.e}
            
    rank_genes_column = f'rank_feat_groups_{clustering_column}'
    if verbose:
        print(f'Ranking feature groups for clusters using obs column {clustering_column}')
        epi.tl.rank_features(adata, groupby=f'{clustering_column}', use_raw=False, key_added=rank_genes_column, method='t-test', n_features=n, omic=omic)
        if show_plots:
            epi.pl.rank_feat_groups_matrixplot(adata, n_features=10, key=rank_genes_column, show=True)
    else:
        with suppress_output():
            epi.tl.rank_features(adata, groupby=f'{clustering_column}', use_raw=False, key_added=rank_genes_column, method='t-test', n_features=n, omic=omic)
            if show_plots:
                epi.pl.rank_feat_groups_matrixplot(adata, n_features=10, key=rank_genes_column, show=True)

    return rank_genes_column


def convert_markers(repo_path=".", lists_path=None, keywords=None, df=None, path="exported_lists", file_name=None, case_sensitive=False, exact=False, style="two_column", organism="Hs", tissue="all", gs=False, ensembl=False, suffix=None):
    """
    Searches the database for given keywords and combines the found marker lists into a new DataFrame.
    Optionally, it can export the DataFrame to a file.

    Parameters
    ----------
    repo_path : str, default "."
        The path of the Marker Repo.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
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
    suffix : str, default None
        The key of the metadata section whose value should be appended to the marker names.

    Returns
    --------
    pd.DataFrame or str :
        If path is specified, the function returns the absolute path to the file where the marker list was saved.
        If path is not specified, the function returns the DataFrame.
    """

    if not os.path.exists(repo_path):
        raise FileNotFoundError(f"The specified repository path '{repo_path}' does not exist.")

    if gs:
        marker_list = guided_search(repo_path=repo_path, lists_path=lists_path, out="marker_list", suffix=suffix)
    elif repo_path and keywords:
        marker_list = get_selected_lists(keywords, repo_path=repo_path, lists_path=lists_path, case_sensitive=case_sensitive, exact=exact, suffix=suffix)
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
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list, suffix=suffix)
        case "ui":
            print("Preparing ui style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list, suffix=suffix)
            marker_list = update_scores(marker_list, repo_path=repo_path)
        case "panglao":
            print("Preparing panglao style marker list...")
            marker_list = compare_marker_lists(repo_path=repo_path, marker_df=marker_list, suffix=suffix)
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


def transfer_markers(target_org=None, source_df=None, repo_path=".", lists_path=None, target_counts=1, weight_markers=False, export_suffix=None, ui=False, custom_file_name=False, ensemble=False, filter_transferred=True, verbose=False, suffix=None, biomart_target=None, homologene_target=None):
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
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
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
    ensemble : bool, default False
        If True, the Ensembl IDs will be used instead of the gene symbols.
    filter_transferred : bool, default True
        If True, filter all marker lists which have already been transferred.
    verbose : bool, default False
        If True, the function will print additional information.
    suffix : str, default None
        The key of the metadata section whose value should be appended to the marker names.

    Returns
    -------
    list of str : 
        The paths of the exported transferred marker lists.
    """

    #TODO: implement verbose

    if not os.path.exists(repo_path):
        raise FileNotFoundError(f"The specified repository path '{repo_path}' does not exist.")

    paths = []
    marker_id = "ensembl" if ensemble else "symbol"

    biomart_organisms = get_supported_biomart_organisms(repo_path=repo_path)
    homologene_organisms = get_supported_taxonomy_ids(repo_path=repo_path)
    supported_organisms = list(set(biomart_organisms + homologene_organisms))

    if not target_org:
        target_org = select(whitelist=supported_organisms, heading="target organism:")
    else:
        target_org = update_organism(target_org, repo_path)

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
        if filter_transferred:
            # Select only lists which were not transferred before
            df = search_df(df=combine_dfs(repo_path=repo_path, lists_path=lists_path), col_to_search="tags_transferred", search_terms=[f"nan"], suffix=suffix)
            source_df = guided_search(df=df, out="metadata", repo_path=repo_path, lists_path=lists_path)
    
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

            hg_db = download_homologene_data(repo_path=repo_path)

            uids = source_df.loc[source_df['Organism name'] == source_organism].index.tolist()
            source_marker_list = combine_lists(uids, repo_path=repo_path, lists_path=lists_path, suffix=suffix)
            print(f"\nDataFrame of the markers of the source organism ({source_organism}) to be transferred to the target organism ({target_organism}):")
            display(source_marker_list)

            if target_counts:
                print(f"Filter source DataFrame by the number of target genes per source gene: remove all source genes that lead to more than {target_counts} target genes.\n")
                filtered_source_df = transfer_markers_homologene(source_marker_list, source_tax, target_tax, hg_db, target_genes, source_whitelist=source_genes, calc_proportions=True, 
                                                        plots=True, target_counts=target_counts)
            else:
                filtered_source_df = source_marker_list
            
            if not filtered_source_df.empty:
                print(f"Create DataFrame containing the transferred genes based on the filter criteria.\n")
                transferred_list = transfer_markers_homologene(filtered_source_df, source_tax, target_tax, hg_db, target_genes,
                                                        source_whitelist=source_genes, calc_proportions=False, plots=True)
                print("\nTransferred markers:")
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

                paths.append(export_marker_list(transferred_list, path="./transferred_markers", file_name=file_name, marker_id=marker_id))
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
            source_marker_list = combine_lists(uids, repo_path=repo_path, lists_path=lists_path, suffix=suffix)
            print(f"\nDataFrame of the markers of the source organism ({source_organism}) to be transferred to the target organism ({target_organism}):")
            display(source_marker_list)

            if target_counts:
                print(f"Filter source DataFrame by the number of target genes per source gene: remove all source genes that lead to more than {target_counts} target genes.\n")
                filtered_source_df = transfer_markers_biomart(biomart_db, source_marker_list, target_genes, source_whitelist=source_genes,
                                                            calc_proportions=True, plots=True, target_counts=target_counts)
            else:
                filtered_source_df = source_marker_list
            
            if not filtered_source_df.empty:
                print(f"Create DataFrame containing the transferred genes based on the filter criteria.\n")
                transferred_list = transfer_markers_biomart(biomart_db, filtered_source_df, target_genes, source_whitelist=source_genes,
                                                            calc_proportions=False, plots=True, target_counts=None)
                print("\nTransferred markers:")
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

                paths.append(export_marker_list(transferred_list, path="./transferred_markers", file_name=file_name, marker_id=marker_id))
            else:
                print("Source DataFrame is empty.")
                
    print("\nFinished!")

    return paths


def validate_settings(cml_parameters=None, repo_path=None, lists_path=None, adata=None, organism=None, rank_genes_column=None, genes_column=None, 
                      clustering_column=None, ensembl=None, col_to_search=None, search_terms=None, column_specific_terms=None):
    """
    Validates user settings including file paths, anndata object columns, and specified organism.

    Parameters
    ----------
    cml_parameters : list of dict, default None
        A list of dictionaries, each containing parameters for 'create_marker_lists' function.
    repo_path : str, default None
        Path to the marker repository.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.
    adata : anndata.AnnData, default None
        The loaded AnnData object.
    organism : str or int, default None
        Organism name, taxon ID, or both. E.g., "mouse", 10090, or "mouse 10090".
    rank_genes_column : str, default None
        Column in .obs table where ranked genes are stored.
    genes_column : str, default None
        Column in .var table where gene symbols or IDs are stored.
    clustering_column : str, default None
        The column in .obs table of the clustering you want to annotate.
    ensembl : bool, default None
        Whether the genes in the genes_column are Ensembl IDs.
    col_to_search : str, default None
        Column to search in for marker list selection.
    search_terms : list of str, default None
        List of search terms to use for marker list selection.
    column_specific_terms : list of dicts, default None
        List of dictionaries with specific 'col_to_search' and 'search_terms' for each.

    Returns
    -------
    bool :
        Returns True if all settings are valid, otherwise prints the errors and returns False.
    """ 

    errors = []
    combined_df_columns = list(combine_dfs(repo_path=repo_path, lists_path=lists_path).columns)
    get_whitelists(repo_path=repo_path, silent_skip=True, update=False)

    # Validate settings dictionaries
    if cml_parameters:
        # Validate keys of dictionaries
        valid_params = set(param.name for param in inspect.signature(create_marker_lists).parameters.values())
        for setting in cml_parameters:
            invalid_params = set(setting) - valid_params
            if invalid_params:
                errors.append(f"Invalid parameters in settings: {', '.join(invalid_params)}")
        # Validate values of dictionaries
        if not invalid_params:
            for i, setting in enumerate(cml_parameters, 1):
                if "repo_path" in setting:
                    if not setting["repo_path"]:
                        errors.append(f"Settings for marker list {i}: No repo_path provided.")
                    else:
                        if not os.path.exists(setting["repo_path"]):
                            errors.append(f"Settings for marker list {i}: Repo path {setting['repo_path']} does not exist.")
                if "organism" in setting:
                    if not setting["organism"]:
                        errors.append(f"Settings for marker list {i}: No organism provided.")
                    else:
                        organism_str = str(setting["organism"]).strip()
                        is_valid_organism = any(organism_str == valid_entry.split(" ")[0] or organism_str == valid_entry.split(" ")[1] or organism_str == valid_entry for valid_entry in valid_organisms)
                        if not is_valid_organism:
                            formatted_valid_organisms = "\n  - " + "\n  - ".join(valid_organisms)
                            errors.append(f"Settings for marker list {i}: Invalid organism or taxon ID {setting['organism']}.\nAvailable options:{formatted_valid_organisms}\n")
                if "column_specific_terms" in setting:
                    if not setting["column_specific_terms"]:
                        errors.append(f"Settings for marker list {i}: No column specific terms provided.")
                    else:
                        for specific_column in setting["column_specific_terms"]:
                            if specific_column not in combined_df_columns:
                                formatted_combined_df_columns = "\n  - " + "\n  - ".join(combined_df_columns)
                                errors.append(f"Settings for marker list {i}: Invalid col_to_search {specific_column}.\nAvailable columns:{formatted_combined_df_columns}\n") 
                else:
                    if "col_to_search" in setting:
                        if not setting["col_to_search"]:
                            errors.append(f"Settings for marker list {i}: No col_to_search provided.")
                        else:
                            if setting["col_to_search"] not in combined_df_columns:
                                formatted_combined_df_columns = "\n  - " + "\n  - ".join(combined_df_columns)
                                errors.append(f"Settings for marker list {i}: Invalid col_to_search {setting['col_to_search']}.\nAvailable columns:{formatted_combined_df_columns}\n")
                    if "search_terms" in setting:
                        if not setting["search_terms"]:
                            errors.append(f"Settings for marker list {i}: No search_terms provided.")
                if "style" in setting:
                    if not setting["style"]:
                        errors.append(f"Settings for marker list {i}: No style provided.")
                    else:
                        if setting["style"] not in ["two_column", "score", "ui", "panglao"]:
                            errors.append(f"Settings for marker list {i}: Invalid style {setting['style']}. Available styles: two_column, score, ui, panglao")
                    
    # Validate individual parameters
    if repo_path and not os.path.exists(repo_path):
        errors.append(f"Repo path {repo_path} does not exist.")

    if adata is None:
        errors.append("No AnnData object provided.")
        
    # List of valid organisms and their tax IDs
    valid_organisms = read_whitelist("organism", repo_path=repo_path)['whitelist']

    # Validate the organism
    if organism:
        organism_str = str(organism).strip()
        is_valid_organism = any(organism_str == valid_entry.split(" ")[0] or organism_str == valid_entry.split(" ")[1] or organism_str == valid_entry for valid_entry in valid_organisms)
        if not is_valid_organism:
            formatted_valid_organisms = "\n  - " + "\n  - ".join(valid_organisms)
            errors.append(f"Invalid organism or taxon ID {organism}.\nAvailable options:{formatted_valid_organisms}\n")

    # Validate obs and var columns
    if rank_genes_column and rank_genes_column not in adata.obs.columns:
        formatted_obs_columns = "\n  - " + "\n  - ".join(adata.obs.columns)
        errors.append(f"Invalid rank_genes_column {rank_genes_column}.\nAvailable columns in adata.obs:{formatted_obs_columns}\n")

    if genes_column and genes_column not in adata.var.columns:
        formatted_var_columns = "\n  - " + "\n  - ".join(adata.var.columns)
        errors.append(f"Invalid genes_column {genes_column}.\nAvailable columns in adata.var:{formatted_var_columns}\n")

    if clustering_column and clustering_column not in adata.obs.columns:
        formatted_obs_columns = "\n  - " + "\n  - ".join(adata.obs.columns)
        errors.append(f"Invalid column {clustering_column}.\nAvailable columns in adata.obs:{formatted_obs_columns}\n")

    # Validate col_to_search and search_terms based on column_specific_terms if provided
    if column_specific_terms:
        for specific_column in column_specific_terms:
            if specific_column not in combined_df_columns:
                formatted_combined_df_columns = "\n  - " + "\n  - ".join(combined_df_columns)
                errors.append(f"Invalid col_to_search {specific_column}.\nAvailable columns:{formatted_combined_df_columns}\n")
    else:
        # Validate col_to_search using the columns from the combined DataFrames in the repo
        if col_to_search:
            if col_to_search not in combined_df_columns:
                formatted_combined_df_columns = "\n  - " + "\n  - ".join(combined_df_columns)
                errors.append(f"Invalid col_to_search {col_to_search}.\nAvailable columns:{formatted_combined_df_columns}\n")

    # Print errors or confirm validation
    if errors:
        print("Validation failed due to the following errors:")
        print("-" * 40)
        for error in errors:
            print(error)
        print("-" * 40)

        return False
    else:
        print("All settings are valid.")
        print(f"Summary of settings:")
        print("-" * 40)
        print("General parameters:")
        print(f"  Repo path: {repo_path}")
        print(f"  Organism: {organism}")
        print(f"  Rank genes column: {rank_genes_column}")
        print(f"  Genes column: {genes_column}")
        print(f"  Clustering column: {clustering_column}")
        print(f"  Ensembl IDs: {ensembl}")

        if column_specific_terms:
            print(f"  Column specific terms:")
            for i, specific_column in enumerate(column_specific_terms.keys()):
                print(f"    {i+1}. Column to search: {specific_column}, Search terms: {column_specific_terms[specific_column]}")
        else:
            print(f"  Column to search: {col_to_search}")
            print(f"  Search terms: {search_terms}")

        if cml_parameters:
            print("\nParameters from dictionary:")
            for i, setting in enumerate(cml_parameters, 1):
                print(f"  {i}. Marker list:")
                for key, value in setting.items():
                    print(f"    {key}: {value}")

        print("-" * 40)

    return True


def create_yaml_list(list_path, list_name=None, organism=None, marker_type=None, marker_col=0, info_col=1, output_path=None, repo_path="."):
    """
    Creates a YAML file containing the marker list.

    Parameters
    ----------
    list_path : str
        The path of the marker list (one or two tab-separated columns).
    list_name : str
        The name of the marker list.
    organism : str, default None
        The organism of the marker list.
    marker_type : str, default None
        The type of the marker list (Genes or Genomic regions).
    marker_col : str, default 0
        The column number where the markers are stored.
    info_col : str, default 1
        The column number where the additional information is stored (e.g. cell type).
    repo_path : str, default "."
        The path of the Marker Repo.

    Returns
    -------
    str :
        The absolute path of the created YAML file.
    """

    if not os.path.exists(list_path):
        raise FileNotFoundError(f"The specified marker list path '{list_path}' does not exist.")

    if not list_name:
        list_name = input("Please enter the name of the marker list: ")

    if not organism:
        organism = select(key="organism")
    else:
        organism = update_organism(organism, repo_path)

    if not marker_type:
        marker_type = select(key="marker_type")

    if not output_path:
        output_path = "."

    marker_list = transform_marker_list(list_path, info_col, marker_col, marker_type, organism)
    UID = get_uid()
    yaml_path = generate_file(UID, list_name, False, marker_list, organism, marker_type, repo_path=repo_path, output_path=output_path)

    return yaml_path


def export_markers_from_anndata(adata, n=50, rank_genes_column='rank_genes_groups', file_name='ranked_markers', omic="RNA"):
    """
    Export the top n marker genes per cluster from an anndata object's ranked genes groups.

    Parameters
    ----------
    adata : anndata.AnnData
        The anndata object containing the ranked genes.
    n : int, default: 50
        The number of top genes per cluster to export.
    rank_genes_column : str, default: 'rank_genes_groups'
        The key/column name where the ranked genes are stored in the anndata object.
    file_name : str, default: 'ranked_markers'
        The file name to save the exported marker list.
    omic : str, default: "RNA"
        The omic type of the anndata object.

    Returns
    -------
    str :
        The path where the marker list is saved.
    """

    if not rank_genes_column in adata.uns.keys() or not rank_genes_column:
        rank_genes_column = rank_feature_groups(adata, show_plots=True, omic=omic)

    df = sc.get.rank_genes_groups_df(adata, group=None, key=rank_genes_column)
    df_sorted = df.sort_values(by=['group', 'scores'], ascending=[True, False])
    df_sliced = df_sorted.groupby('group').head(n)

    marker_gene_df = df_sliced[['names', 'group']]

    # Display the top row of each unique group
    top_rows = df_sliced.groupby('group').first()
    top_rows.rename(columns={"names": "Marker gene", "scores": "Score"}, inplace=True)
    top_rows.index.name = "Cell type"
    top_rows_subset = top_rows[["Marker gene", "Score"]]
    display(top_rows_subset)

    path = export_marker_list(marker_gene_df, file_name=file_name)

    return path