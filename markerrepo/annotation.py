import math
import os
import statistics
import pandas as pd
import scanpy as sc
from IPython.display import display
from .marker_repo import read_whitelist, combine_dfs, export_marker_list


def annot_ct(genes_adata, adata=None, output_path=".", db_path=None, cluster_path=None, cluster_column=None, rank_genes_column=None, sample="sample", ct_column="cell_types", tissue="all", species="Hs", inplace=True, header=False, min_hits=4, verbose=False, ignore_overwrite=False):
    """
    This function calculates potential cell types per cluster and adds them to the obs table of the anndata object.

    Parameters
    ----------
    genes_adata : anndata.AnnData
        The anndata object which contains clustered data, gene IDs as index as well as rank genes groups.
    adata : anndata.AnnData, default None
        The anndata object to add the annotations to. If None, the annotations will be written to genes_adata.
    output_path : string, default "."
        The path to the folder where the annotation file will be written and where the ranks folder will be created.
    db_path : string, default None
        The path to the cell type marker gene database file.
    cluster_path : string, default None
        The path to the folder which contains the "cluster files": Tab-separated files containing the genes and
        the corresponding ranked scores. Use only if you already created your own cluster files.
    cluster_column : string, default None
        The column of the .obs table which contains the clustering information. E.g. "louvain" or "leiden".
    rank_genes_column : string, default None
        The column of the .uns table which contains the rank genes scores. E.g. "rank_genes_groups".
    sample : string, default "sample"
        The name of the sample. E.g. "sample1" or "zebrafish". This will be used for naming the output files.
    ct_column : string, default "cell_types"
        The column of the .obs table of the anndata object (adata or genes_adata) which will include the new cell type annotation.
    tissue : string, default "all"
        If tissue is not "all", only marker genes found in the entered tissue will be taken into account.
        This only works if you use the whole panglao database. If you use custom databases (such as combined lists of the marker repo), 
        you can ignore this parameter.
    species : string, default "hs"
        The species of the data. (Hs or Mm supported)
        This only works if you use the whole panglao database. If you use custom databases (such as combined lists of the marker repo), 
        you can ignore this parameter.
    inplace : boolean, default True
        Whether to add the annotations to the adata object in place.
    header : bool, default False
        Skip first line if header is True.
    min_hits : int, default 4
        Minimum number of hits required to consider a cell type for annotation.
    verbose : bool, default False
        Whether to print additional information.
    ignore_overwrite : bool, default False
        Whether to ignore the overwrite warning.

    Returns
    --------
    If inplace is True, the annotation is added to adata.obs in place.
    Else, a copy of the adata object is returned with the annotations added.
    """
    
    if db_path is None:
        raise ValueError("The parameter 'db_path' must be provided. This is the path to the cell type marker gene database file.")
    elif not os.path.exists(db_path):
        raise FileNotFoundError(f"The specified database file path '{db_path}' does not exist.")
    
    if cluster_column is None:
        raise ValueError("The parameter 'cluster_column' must be provided. This is the column in the .obs table of the AnnData object that contains the clustering information.")
    elif cluster_column not in genes_adata.obs.columns:
        raise KeyError(f"The specified cluster column '{cluster_column}' does not exist in the .obs table of the provided AnnData object.")
    
    if rank_genes_column is None:
        raise ValueError("The parameter 'rank_genes_column' must be provided. This is the column in the .uns table of the AnnData object that contains the rank genes scores.")
    elif rank_genes_column not in genes_adata.uns.keys():
        raise KeyError(f"The specified rank genes column '{rank_genes_column}' does not exist in the .uns table of the provided AnnData object.")

    go_on = True

    if not adata:
        adata = genes_adata

    if not inplace:
        adata = adata.copy()

    if output_path and db_path:
        cluster_path = f"{output_path}/ranked/clusters/{cluster_column}"
        ct_path = f"{output_path}/ranked/output/{cluster_column}"

        if not ignore_overwrite and os.path.exists(ct_path):
            print(f"Warning: The path {ct_path}/ already exists!\nAll annotation files will be overritten.")
            go_on = False

            if not go_on:
                go_on = input("Do you want to continue? (yes/no): ")
                go_on = True if go_on == "yes" else False

            if not go_on:
                print("Cell type annotation has been aborted.")

                return
            
        if verbose:
            print(f"Output folder: {ct_path}/", "\nDB file: " + db_path, f"\nCluster folder: {cluster_path}/",
                "\nTissue: " + tissue)
        if adata and genes_adata and cluster_column:
            # Create folders containing the annotation assignment table as well as the detailed scoring files per cluster
            if not os.path.exists(f'{ct_path}'):
                os.makedirs(f'{ct_path}')
                if verbose:
                    print(f'Created folder: {ct_path}')

            # Check if cluster_path exists
            if not ignore_overwrite and os.path.exists(cluster_path):
                user_input = input(f"The folder {cluster_path} already exists.\nDo you want to skip creating new ranked cluster files and keep the old ones? (yes/no): ")
                if user_input.lower() == 'yes':
                    print("Skipping the creation of new ranked cluster files.")
                else:
                    print("Creating new ranked cluster files.")
                    write_cluster_files(cluster_path, sample, adata, cluster_column, genes_adata, rank_genes_column)
            else:
                # Create folder if it doesn't exist and write files
                if not os.path.exists(cluster_path):
                    os.makedirs(cluster_path)
                    if verbose:
                        print(f'Created folder: {cluster_path}')

                write_cluster_files(cluster_path, sample, adata, cluster_column, genes_adata, rank_genes_column)

            # Perform the actual cell type annotation per clustering resolution
            if verbose:
                print("Starting cell type annotation.")
                print(output_path, ct_path, cluster_column)
            perform_cell_type_annotation(
                f"{ct_path}/", db_path, f"{cluster_path}/", tissue, species=species, header=header, min_hits=min_hits, verbose=verbose)

            # Add information to the adata object
            if verbose:
                print("Adding information to the adata object.")
            cta_dict = {}
            with open(f'{ct_path}/annotation.txt') as file:
                for line in file:
                    cluster, ct = line.split('\t')
                    cta_dict[cluster] = ct.rstrip()
            adata.obs[f'{ct_column}'] = adata.obs[f'{cluster_column}'].map(cta_dict)

            if verbose:
                print(f"Finished cell type annotation! The results are found in the .obs table {ct_column}.")

            if not inplace:
                return adata

        elif cluster_path:
            if verbose:
                print("Output folder: " + output_path, "\nDB file: " + db_path, "\nCluster folder: " + cluster_path,
                      "\nTissue: " + tissue)
            perform_cell_type_annotation(
                f"{output_path}/ranked/output/{cluster_column}/", db_path, cluster_path, tissue, header=header, verbose=verbose)
            if verbose:
                print(f"Cell type annotation of output path {ct_path}/ finished.")

        else:
            pass


def modify_ct(adata=None, annotation_dir=None, clustering_column="leiden_0.1", cell_type_column="cell_types_leiden_0.1", inplace=True):
    """
    This function can be used to make subsequent changes to cell types that were previously annotated with the annot_ct() function.
    For each annotated cluster, a choice of 10 possible alternative assignments is presented.

    Parameters
    ----------
    adata : anndata.AnnData, default None
        The anndata object containing cell type assignments from the annot_ct() function.
    annotation_dir : string, default None
        The path where the annotation files are being stored (should be the same path as the output_path parameter of the annot_ct function).
    clustering_column : string, default "leiden"
        The obs column containing the clustering information.
    cell_type_column : string, defaul "cell_types"
        The obs column containing the cell type annotation.
    inplace : boolean, default True
        Whether to add the new cell type assignments to the adata object in place.

    Returns
    --------
    If inplace is True, the modified annotation is added to adata.obs in place.
    Else, a copy of the adata object is returned with the annotations added.
    """

    if not inplace:
        adata = adata.copy()

    adata.obs[f'{cell_type_column}_mod'] = adata.obs[f'{cell_type_column}']

    modify = True
    while modify:
        cluster = int(input("Enter the number of the cluster you'd like to modify: "))
        df = pd.read_csv(f'{annotation_dir}/ranked/output/{clustering_column}/ranks/cluster_{cluster}', sep='\t', names=["Cell type", "Score", "Hits", "Number of marker genes", "Mean of UI"])
        display(df.head(10))
        new_ct = int(input("Please choose another cell type by picking a number of the corresponding index column: "))
        adata.obs[f'{cell_type_column}_mod'] = adata.obs[f'{cell_type_column}_mod'].cat.rename_categories({df.iat[0, 0]: df.iat[new_ct, 0]})
        print(f'Succesfully replaced {df.iat[0, 0]} with {df.iat[new_ct, 0]}.')
        umap = input("Would you like to see the updated UMAP? Enter yes or no: ")
        umap = True if umap == "yes" else False
        if umap:
            sc.pl.umap(adata, color=[f'{cell_type_column}_mod', f'{cell_type_column}'], wspace=0.5)
        modify = input("Would you like to modify another cluster? Enter yes or no: ")
        modify = True if modify == "yes" else False

    if not inplace:
        return adata


def calculate_normalized_diffs(df):
    """
    Calculate both normal and scaled normalized difference scores for adjacent cell type pairs in a DataFrame.

    This function first calculates the normal normalized differences between each cell type and 
    its immediate successor based on their scores. Then, it scales these differences within a range of 0 to 100,
    where the highest difference is scaled to 100 and the lowest to 0.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame sorted by 'Score', containing cell types and their corresponding scores.

    Returns
    -------
    tuple of dicts :
        The first dictionary contains the normal normalized differences, and the second dictionary 
        contains the scaled normalized differences, both keyed by cell type.
    """

    # Check if df has at least two rows
    if len(df) < 2:
        return {}, {}  

    # Calculate normal normalized diffs
    normal_diffs = {}
    for i in range(len(df) - 1):
        score_current = df.at[i, 'Score']
        score_next = df.at[i + 1, 'Score']
        normalized_diff = (score_current - score_next) / score_current if score_current != 0 else 0
        normal_diffs[df.at[i, 'Cell type']] = normalized_diff
    normal_diffs[df.at[len(df) - 1, 'Cell type']] = 0

    # Scale normalized diffs
    scaled_diffs = normal_diffs.copy()
    max_diff = max(scaled_diffs.values())
    min_diff = min(scaled_diffs.values())
    range_diff = max_diff - min_diff
    if range_diff != 0:
        for key in scaled_diffs:
            scaled_diffs[key] = round((scaled_diffs[key] - min_diff) / range_diff * 100)
    else:
        for key in scaled_diffs:
            scaled_diffs[key] = 0  # or 100 if all diffs are the same

    return normal_diffs, scaled_diffs


def show_tables(annotation_dir=None, n=5, clustering_column="leiden_0.1", show_diff=False):
    """
    Display dataframes for each cluster showing scores, hits, number of genes, mean of UI. Optionally, 
    it can also show both normal and scaled normalized differences of every potential cell type.

    For each cluster, this function reads the corresponding file, sorts the cell types based on their scores.
    If show_diff is True, it calculates both normal and scaled normalized differences between each cell type 
    and its immediate successor, and then displays the dataframe with these additional columns.

    Parameters
    ----------
    annotation_dir : string, default None
        The directory path where the annotation files are stored.
    n : int, default 5
        The maximum number of rows to display for each cluster.
    clustering_column : string, default "leiden_0.1"
        The clustering column used for cell type annotation.
    show_diff : bool, default False
        Whether to show the normalized differences in the output.
    """

    path = f'{annotation_dir}/ranked/output/{clustering_column}/ranks'

    files = os.listdir(path)
    for file in files:
        cluster = file.split("_")[1]
        ct_column = f"Cluster {cluster}"
        df = pd.read_csv(f'{path}/{file}', sep='\t', names=[ct_column, "Score", "Hits", "Number of marker genes", "Mean of UI"])

        if show_diff:
            # Calculate and add both normal and scaled diffs to the DataFrame if show_diff is True
            normal_diffs, scaled_diffs = calculate_normalized_diffs(df.rename(columns={ct_column: "Cell type"}))
            # df_sorted['Normalized Diff'] = df_sorted[ct_column].apply(lambda x: normal_diffs.get(x, 0))
            df['Adjacent Disparity'] = df[ct_column].apply(lambda x: scaled_diffs.get(x, 0))

        display(df.head(n))


def write_cluster_files(cluster_path, sample, adata, cluster_column, genes_adata, rank_genes_column, log=False):
    """
    Writes one file per cluster that contains gene IDs and their corresponding scores, sorted by ranking.

    Parameters
    ----------
    cluster_path : string, default None
        The path to the folder which contains the "cluster files": Tab-separated files containing the genes and
        the corresponding ranked scores. Use only if you already created your own cluster files.
    sample : string, default "sample"
        The name of the sample. E.g. "sample1" or "zebrafish". This will be used for naming the output files.
    adata : anndata.AnnData, default None
        The anndata object to add the annotations to.
    cluster_column : string, default None
        The column of the .obs table which contains the clustering information. E.g. "louvain" or "leiden".
    genes_adata : anndata.AnnData
        The anndata object which contains clustered data, gene ID as index as well as rank genes groups.
    rank_genes_column : string, default None
        The column of the .uns table which contains the rank genes scores. E.g. "rank_genes_groups".
    log : bool, default False
        Whether to log the progress.
    """

    clusters = adata.obs[f'{cluster_column}'].unique()
    total_clusters = len(clusters)

    for index, cluster in enumerate(clusters):
        with open(f'{cluster_path}/{sample}.cluster_{cluster}', 'w') as file:
            if log:
                print(f"Writing ranked cluster file {sample}.cluster_{cluster} ({index+1}/{total_clusters})")
            for i, gene in enumerate(genes_adata.uns[f'{rank_genes_column}']['names'][cluster]):
                score = genes_adata.uns[f'{rank_genes_column}']['scores'][cluster][i]
                file.write(f'{gene.split("_")[0]}\t{score}\n')


def parse_marker_database(file_path, tissue="all", species=None, header=False):
    """
    Parse a marker gene database file to extract marker genes for various cell types.

    The function can handle two types of general file formats:
    - Two-column: [Gene Symbol, Cell Type]
    - Three-column: [Gene Symbol, Cell Type, Score]
    Additionally, it can process a specialized six-column format specific to the PanglaoDB.

    Parameters
    ----------
    file_path : str
        The path to the marker gene file.
    tissue : str, default "all"
        The target tissue type. If set to "all", markers from all tissues will be included.
        Only applicable when processing the PanglaoDB format.
    species : str, default None
        The target species. If None, markers from all species will be included.
        Only applicable when processing the PanglaoDB format.
    skip_header : bool, default False
        Whether to skip the first line of the file as a header.

    Returns
    -------
    dict :
        A dictionary where each key is a cell type and the value is another dictionary. 
        The inner dictionary maps marker genes (keys) to their score (values), if available.
    """

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The specified file {file_path} was not found.")

    tissues = [tissue]
    panglao_dict = {}
    panglao_rank_dict = {}

    with open(file_path, "r") as panglao_file:
        if header:
            panglao_file.readline()
        for line in panglao_file.readlines():
            line_split = line.split("\t")
            if len(line_split) == 2:  # two-column file
                gene_symb, ct = line_split
                us = 1.0
                spec = ''
                if ct not in panglao_dict.keys():
                    panglao_dict[ct] = []
                panglao_dict[ct].append((us, gene_symb.strip()))
            elif len(line_split) == 3:  # three-column file
                gene_symb, ct, ub_i = line_split
                us = float(ub_i)
                if us != 0:
                    us = round(math.sqrt(1 / us))
                else:
                    us = 32
                spec = ''
                if ct not in panglao_dict.keys():
                    panglao_dict[ct] = []
                panglao_dict[ct].append((us, gene_symb.strip()))
            elif len(line_split) == 6:  # six-column file -> whole PanglaoDB
                spec, gene_symb, ct, n_genes, ub_i, organ = line_split
                us = float(ub_i)
                if us != 0:
                    us = round(math.sqrt(1 / us))
                else:
                    us = 32
                if species is None or species in spec:
                    if "all" in tissues:
                        if ct not in panglao_dict.keys():
                            panglao_dict[ct] = []
                        genes = [(us, gene_symb)]
                        if len(n_genes.split("|")) > 1:
                            for gene in n_genes.split("|"):
                                genes.append((us, gene.upper()))
                        elif n_genes != "NA":
                            genes.append((us, n_genes))
                        for gene in genes:
                            panglao_dict[ct].append(gene)
                    elif any(t in organ.lower() for t in tissues):
                        if ct not in panglao_dict.keys():
                            panglao_dict[ct] = []
                        genes = [(us, gene_symb)]
                        if len(n_genes.split("|")) > 1:
                            for gene in n_genes.split("|"):
                                genes.append((us, gene.upper()))
                        elif n_genes != "NA":
                            genes.append((us, n_genes))
                        for gene in genes:
                            panglao_dict[ct].append(gene)
            else:
                raise ValueError(f"The file format is not recognized. Please use a two-column, three-column, or PanglaoDB six-column format. Got {len(line_split)}.")

    for ct in panglao_dict.keys():
        rank_dict = {}
        for gene in panglao_dict[ct]:
            rank_dict[gene[1]] = gene[0]

        panglao_rank_dict[ct] = rank_dict

    return panglao_rank_dict


def calc_ranks(cell_marker_dict, cluster_gene_ranks, min_hits=4, verbose=True):
    """
    Calculate cell type annotation scores for each cluster based on differential gene scores and cell type marker genes.

    Parameters
    ----------
    cell_marker_dict : dict
        A dictionary representing the cell marker database. Keys are cell types, and values are dictionaries
        with gene names as keys and ubiquitousness scores as values.
    cluster_gene_ranks : dict
        A dictionary containing ranked gene scores for each gene in each cluster.
    min_hits : int, default 4
        The minimum number of marker genes required to consider a cell type for annotation in a cluster.
    verbose : bool, default True
        If True, prints additional information about the gene overlap between the database and input data.

    Returns
    -------
    dict :
        A dictionary with annotation scores, number of hits, total marker genes, and the average ubiquitousness
        score for each cell type in each cluster.
    """

    annotation_scores = {}
    input_genes = set()
    matched_genes = set()
    database_genes = set(cell_marker_dict.keys())

    # Initialize cluster keys in annotation scores dictionary
    for cluster in cluster_gene_ranks:
        annotation_scores[cluster] = {}

    # Calculate scores for each cell type and cluster
    for cell_type, markers in cell_marker_dict.items():
        num_markers = len(markers)

        for cluster, genes in cluster_gene_ranks.items():
            match_count = 0
            rank_scores = []
            ubiquity_scores = []

            for gene, rank_score in genes.items():
                input_genes.add(gene)

                if gene in markers:
                    matched_genes.add(gene)
                    ubiquity_score = markers[gene]
                    weighted_score = rank_score * ubiquity_score
                    rank_scores.append(weighted_score)
                    ubiquity_scores.append(ubiquity_score)
                    match_count += 1

            if match_count >= min_hits:
                avg_ubiquity_score = round(statistics.mean(ubiquity_scores))
                total_score = round(sum(rank_scores) / math.sqrt(num_markers))
                annotation_scores[cluster][cell_type.rstrip()] = [total_score, match_count, num_markers, avg_ubiquity_score]

    # Print summary if verbose is True
    if verbose:
        db_gene_count = len(database_genes)
        input_gene_count = len(input_genes)
        matched_gene_count = len(matched_genes)
        overlap_percentage = round(matched_gene_count / db_gene_count * 100)
        print(f"The database contains {db_gene_count} different genes. "
              f"The input data contains {input_gene_count} different genes. "
              f"The genes of the input data overlap with {matched_gene_count} genes in total, {overlap_percentage}%.")

    # Sort cell types by score and match count/total marker genes
    for cluster in annotation_scores:
        annotation_scores[cluster] = dict(sorted(annotation_scores[cluster].items(), 
                                                 key=lambda item: (-item[1][0], -item[1][1]/item[1][2])))


    return annotation_scores


def get_cell_types(cluster_path, db_path, tissue="all", species="Hs", header=False, min_hits=4, verbose=True):
    """
    Prepare database and clusters for upcoming ranking calculations.

    Parameters
    ----------
    cluster_path : string
        The path to the folder which contains the "cluster files": Tab-separated files containing the
        genes and the corresponding ranked scores.
    db_path : string
        The path to the cell type marker gene database file.
    tissue : string, default "all"
        If tissue is not "all", only marker genes found in the entered tissue will be taken into account.
    species : string, default "hs"
        The species of the data.
    header : bool, default False
        Skip first line if header is True.
    min_hits : int, default 4
        Minimum number of hits required to consider a cell type for annotation.
    verbose : bool, default True
        Whether to print additional information.

    Returns
    -------
    dictionary :
        The dictionary which contains the scores, the quantity of hits, the overall marker genes and
        the ubiquitousness index per cell type for each cluster.
    """

    db_dict = parse_marker_database(db_path, tissue=tissue, species=species, header=header)
    annotated_clusters = get_annotated_clusters(cluster_path=cluster_path)

    return calc_ranks(db_dict, annotated_clusters, min_hits=min_hits, verbose=verbose)


def get_annotated_clusters(cluster_path, show_duplicates=False):
    """
    Reads cluster files from a specified directory to aggregate positive ranked scores
    for each gene within each cluster.

    Parameters
    ----------
    cluster_path : str
        The path to the folder containing the cluster files.
    show_duplicates : bool, default False
        Whether to show duplicate genes with positive scores.

    Returns
    -------
    dict :
        A dictionary where keys are cluster names and values are dictionaries of genes with their 
        aggregated positive ranked scores.
    """

    annotated_clusters = {}
    files = os.listdir(cluster_path)
    for file in [x for x in files if not x.startswith(".")]:
        cname = file.split(".cluster_")[1]
        annotated_dict = {}        
        duplicate_genes = {}

        with open(os.path.join(cluster_path, file)) as cfile:
            annotated_dict[cname] = []
            lines = cfile.readlines()
            for line in lines:
                split = line.split("\t")
                if len(split) == 2:
                    gene = split[0].upper()
                    score = float(split[1].rstrip())
                    annotated_dict[cname].append([gene, score])

                    # Check for duplicate genes if show_duplicates is True
                    if show_duplicates and score > 0:
                        if gene in duplicate_genes:
                            duplicate_genes[gene].append(score)
                        else:
                            duplicate_genes[gene] = [score]

        # Output duplicate genes and their scores if show_duplicates is True
        if show_duplicates:
            for gene, scores in duplicate_genes.items():
                if len(scores) > 1:
                    print(f"Duplicate gene found in cluster {cname}: {gene} with scores {scores}")

        sum_dict = {}
        for gene in annotated_dict[cname]:
            if gene[1] > 0:  # Only consider positive values
                if gene[0] in sum_dict.keys():
                    sum_dict[gene[0]] += gene[1]
                else:
                    sum_dict[gene[0]] = gene[1]

        annotated_clusters[cname] = sum_dict

    return annotated_clusters


def perform_cell_type_annotation(output, db_path, cluster_path, tissue="all", species="Hs", header=False, min_hits=4, verbose=True):
    """
    Performs cell type identification, generate cell type assignment table
    and create ranks folder with files for further investigation (one per cluster).

    Parameters
    ----------
    output : string
        The path to the folder where the annotation file will be written and where the ranks folder will be created.
    db_path : string
        The path to the cell type marker gene database file.
    cluster_path : string
        The path to the folder which contains the "cluster files": Tab-separated files containing the genes and
        the corresponding ranked scores.
    tissue : string, default "all"
        If tissue is not "all", only marker genes found in the entered tissue will be taken into account.
    species : string, default "hs"
        The species of the data.
    header : bool, default False
        Skip first line if header is True.
    min_hits : int, default 4
        Minimum number of hits required to consider a cell type for annotation.
    verbose : bool, default True
        Whether to print additional information.
    """

    opath = output + "/ranks/"
    if not os.path.exists(opath):
        os.makedirs(opath)

    ct_dict = get_cell_types(cluster_path, db_path, tissue, species=species, header=header, min_hits=min_hits, verbose=verbose)
    write_annotation(ct_dict, output)


def write_annotation(ct_dict, output):
    """
    Writes a tab-separated file that contains exactly one cell type (the one with the highest score)
    for each cluster to which at least one cell type could be assigned ("cell type assignment table").

    Parameters
    ----------
    ct_dict : dictionary
        The dictionary which contains the scores, the quantity of hits, the overall marker genes and
        the ubiquitousness index per cell type for each cluster. This dictionary is being returned
        by the calc_ranks() method.
    output : string
        The path to the folder where the annotation file will be written.
    """

    with open(output + "/annotation.txt", "w") as c_file:
        for cluster in ct_dict.keys():
            if len(ct_dict[cluster]) > 0:
                highest_scored_cell_type = next(iter(ct_dict[cluster]))
                c_file.write(f"{cluster}\t{highest_scored_cell_type}\n")

            with open(f"{output}/ranks/cluster_{cluster}", "w") as d_file:
                for cell_type, values in ct_dict[cluster].items():
                    d_file.write(f"{cell_type}\t" + "\t".join(map(str, values)) + "\n")
    

def list_possible_settings(repo_path, lists_path=None, adata=None):
    """
    Lists all possible settings based on the repo and AnnData object.

    Parameters
    ----------
    repo_path : str
        Path to the marker repository.
    lists_path : str, default None
        The path of the folder which contains the marker lists (YAML files).
        If None, the lists folder of the repo_path is used.        
    adata : anndata.AnnData, default None
        The loaded AnnData object.
    """

    if not os.path.exists(repo_path):
        print(f"Repo path {repo_path} does not exist.")
        return

    print("Possible Settings:")
    print("-" * 40)

    combined_df_columns = list(combine_dfs(repo_path=repo_path, lists_path=lists_path).columns)
    print("1. Available columns to search in MarkerRepo:")
    for col in combined_df_columns:
        print(f"  - {col}")
    
    valid_organisms = read_whitelist("organism", repo_path=repo_path)['whitelist']
    print("\n2. Possible organisms or taxon IDs:")
    for org in valid_organisms:
        print(f"  - {org}")

    if adata:
        print("\n3. Available columns in adata.obs:")
        for col in adata.obs.columns:
            print(f"  - {col}")

        print("\n4. Available columns in adata.var:")
        for col in adata.var.columns:
            print(f"  - {col}")
    
    print("-" * 40)


def export_markers_from_anndata(adata, n=50, rank_genes_column='rank_genes_groups', file_name='ranked_markers'):
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

    Returns
    -------
    str :
        The path where the marker list is saved.
    """
    
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


def compare_cell_types(adata, column, obs_columns):
    """
    Group the observation DataFrame of an anndata object by a specific column and include relevant cell types from given obs columns.

    Parameters
    ----------
    adata : anndata.AnnData
        The anndata object containing the .obs DataFrame.
    column : str
        The column by which to group the .obs DataFrame.
    obs_columns : list of str
        List of column names in the .obs DataFrame to keep for the comparison.

    Returns
    -------
    DataFrame :
        A grouped DataFrame based on the specified column, incorporating relevant cell types from the obs columns.
    """
    
    if column in obs_columns:
        obs_columns.remove(column)

    obs_df = adata.obs
    columns_to_keep = [column] + obs_columns  # Include the grouping column and additional columns from obs_columns
    
    filtered_obs_df = obs_df[columns_to_keep]
    grouped_obs_df = filtered_obs_df.groupby(column).agg('first')
    
    return grouped_obs_df


def reformat_marker_list(input_path, suffix="_SCSA"):
    """
    Reformat a TSV marker list for compatibility with SCSA.

    Parameters
    ----------
    input_path : str
        Path to the input TSV marker list.
    suffix : str, "_SCSA"
        Suffix to append to output file name.

    Returns
    -------
    str
        Path where the reformatted marker list is saved.
    """

    df = pd.read_csv(input_path, header=None, sep='\t')
    df = df.iloc[:, :2]
    df.columns = ["marker", "cell_name"]
    output_path = f"{input_path}{suffix}"
    df.to_csv(output_path, index=False, sep='\t')

    return output_path


def update_adata_with_markers(adata, list_name, df, ignore_overwrite=False):
    """
    Updates the provided AnnData object (adata) with marker lists associated with a given file name.

    This function adds a new entry to the 'MarkerRepo' key in the .uns attribute of the AnnData object. 
    It creates a dictionary under 'MarkerRepo' if it does not exist and adds the marker lists from the 
    DataFrame 'df' under the specified 'list_name'. If 'list_name' already exists in 'MarkerRepo', 
    and ignore_overwrite is False, the user is repeatedly prompted to decide whether to overwrite, update, or leave unchanged 
    until a valid response is given. If ignore_overwrite is True, the existing entry is overwritten without prompt.

    Parameters
    ----------
    adata : AnnData
        The AnnData object to be updated.
    list_name : str
        The name under which the marker lists will be stored in the AnnData object.
    df : DataFrame
        The DataFrame containing the marker lists. The index of this DataFrame should be the list of markers 
        that will be stored in the AnnData object.
    ignore_overwrite : bool
        If set to True, any existing entry with the same list_name is overwritten without prompt.
    """

    if "MarkerRepo" not in adata.uns:
        adata.uns["MarkerRepo"] = {"marker_lists": {}}
    
    while list_name in adata.uns["MarkerRepo"]["marker_lists"] and not ignore_overwrite:
        response = input(f"'{list_name}' already exists. Overwrite (O), Update (U), or Leave unchanged (L)? [O/U/L]: ").strip().upper()
        if response == 'O':
            adata.uns["MarkerRepo"]["marker_lists"][list_name] = df.index.astype(str).tolist()
            break
        elif response == 'U':
            existing_list = set(adata.uns["MarkerRepo"]["marker_lists"][list_name])
            new_list = set(df.index.astype(str).tolist())
            combined_list = list(existing_list.union(new_list))
            adata.uns["MarkerRepo"]["marker_lists"][list_name] = combined_list
            break
        elif response == 'L':
            print("No changes made.")
            break
        else:
            print("Invalid response. Please enter O, U, or L.")

    if ignore_overwrite or list_name not in adata.uns["MarkerRepo"]["marker_lists"]:
        adata.uns["MarkerRepo"]["marker_lists"][list_name] = df.index.astype(str).tolist()
