import math
import os
import statistics
import pandas as pd
import scanpy as sc
from IPython.display import display
from .marker_repo import read_whitelist, get_whitelists, combine_dfs, export_marker_list


def annot_ct(genes_adata, adata=None, output_path=".", db_path=None, cluster_path=None, cluster_column=None, rank_genes_column=None, sample="sample", ct_column="cell_types", tissue="all", species="Hs", inplace=True, header=False):
    """
    If the script is called via a package (atactoolbox), please use this function.
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

    Returns
    --------
    If inplace is True, the annotation is added to adata.obs in place.
    Else, a copy of the adata object is returned with the annotations added.
    """
    
    go_on = True

    if not adata:
        adata = genes_adata

    if not inplace:
        adata = adata.copy()

    if output_path and db_path:
        cluster_path = f"{output_path}/ranked/clusters/{cluster_column}"
        ct_path = f"{output_path}/ranked/output/{cluster_column}"

        if os.path.exists(ct_path):
            print(f"Warning: The path {ct_path}/ already exists!\nAll annotation files will be overritten.")
            go_on = False

        if not go_on:
            go_on = input("Do you want to continue? (yes/no): ")
            go_on = True if go_on == "yes" else False

            if not go_on:
                print("Cell type annotation has been aborted.")

                return

        print(f"Output folder: {ct_path}/", "\nDB file: " + db_path, f"\nCluster folder: {cluster_path}/",
              "\nTissue: " + tissue)
        if adata and genes_adata and cluster_column:
            # Create folders containing the annotation assignment table as well as the detailed scoring files per cluster
            if not os.path.exists(f'{ct_path}'):
                os.makedirs(f'{ct_path}')
                print(f'Created folder: {ct_path}')

            # Check if cluster_path exists
            if os.path.exists(cluster_path):
                user_input = input(f"The folder {cluster_path} already exists.\nDo you want to skip creating new ranked cluster files and keep the old ones? (yes/no): ")
                if user_input.lower() == 'yes':
                    print("Skipping the creation of new ranked cluster files.")
                else:
                    print("Creating new ranked cluster files.")
                    write_cluster_files(cluster_path, sample, adata, cluster_column, genes_adata, rank_genes_column)
            else:
                # Create folder if it doesn't exist and write files
                os.makedirs(cluster_path)
                write_cluster_files(cluster_path, sample, adata, cluster_column, genes_adata, rank_genes_column)

            # Perform the actual cell type annotation per clustering resolution
            print("Starting cell type annotation.")
            print(output_path, ct_path, cluster_column)
            perform_cell_type_annotation(
                f"{ct_path}/", db_path, f"{cluster_path}/", tissue, species=species, header=header)

            # Add information to the adata object
            print("Adding information to the adata object.")
            cta_dict = {}
            with open(f'{ct_path}/annotation.txt') as file:
                for line in file:
                    cluster, ct = line.split('\t')
                    cta_dict[cluster] = ct.rstrip()
            adata.obs[f'{ct_column}'] = adata.obs[f'{cluster_column}'].map(cta_dict)

            print(f"Finished cell type annotation! The results are found in the .obs table {ct_column}.")

            if not inplace:
                return adata

        elif cluster_path:
            print("Output folder: " + output_path, "\nDB file: " + db_path, "\nCluster folder: " + cluster_path,
                  "\nTissue: " + tissue)
            perform_cell_type_annotation(
                f"{output_path}/ranked/output/{cluster_column}/", db_path, cluster_path, tissue, header=header)
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


def show_tables(annotation_dir=None, n=5, clustering_column="leiden_0.1"):
    """
    Show dataframes of each cluster which shows score, hits, number of genes and mean of the UI of every potential cell type.

    Parameters
    ----------
    annotation_dir : string, default None
        The path where the annotation files are being stored (should be the same path as the output_path parameter of the annot_ct function).
    n : int, default 5
        The maximum number of rows to show
    clustering_column : string, default "leiden"
        The clustering column of the obs table which has been used for cell type annotation.
    """

    path = f'{annotation_dir}/ranked/output/{clustering_column}/ranks'

    files = os.listdir(path)
    for file in files:
        cluster = file.split("_")[1]
        df = pd.read_csv(f'{path}/{file}', sep='\t', names=[f"Cluster {cluster}: Cell type", "Score", "Hits", "Number of marker genes", "Mean of UI"])
        display(df.head(n))


def write_cluster_files(cluster_path, sample, adata, cluster_column, genes_adata, rank_genes_column):
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
    """

    clusters = adata.obs[f'{cluster_column}'].unique()
    total_clusters = len(clusters)

    for index, cluster in enumerate(clusters):
        with open(f'{cluster_path}/{sample}.cluster_{cluster}', 'w') as file:
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
            else:  # six-column file -> whole PanglaoDB
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

    for ct in panglao_dict.keys():
        rank_dict = {}
        for gene in panglao_dict[ct]:
            rank_dict[gene[1]] = gene[0]

        panglao_rank_dict[ct] = rank_dict

    return panglao_rank_dict


def calc_ranks(cm_dict, annotated_clusters):
    """
    Identify cell types of each cluster by ranking each potential cell type using fitting genes, ranked scores,
    quantity of available marker genes per cell type aswell as using the panglao ubiquitousness index.

    Parameters
    ----------
    cm_dict : dictionary
        Dictionary which contains the cell marker database.
    annotated_clusters :
        Dictionary which contains the summed up ranked scores per gene for each cluster.

    Returns
    -------
    dictionary :
        The dictionary which contains the scores, the quantity of hits, the overall marker genes and
        the ubiquitousness index per cell type for each cluster.
    """

    ct_dict = {}
    data_genes = []
    data_hits = []
    db_genes = []

    for key in annotated_clusters.keys():
        ct_dict[key] = {}

    for celltype in cm_dict.keys():
        gene_count = len(cm_dict[celltype])
        for c in annotated_clusters.keys():
            count = 0
            ranks = []
            ub_scores = []

            for mgene in annotated_clusters[c].keys():
                data_genes.append(mgene)
                if mgene in cm_dict[celltype].keys():
                    data_hits.append(mgene)
                    gene_score, ub_score = annotated_clusters[c][mgene], cm_dict[celltype][mgene]
                    gene_score = gene_score * ub_score
                    ranks.append(gene_score)
                    ub_scores.append(ub_score)
                    count += 1

            # ranks = sorted(ranks, reverse=True)

            # if count >= 10:
            #     ranks = ranks[:10]

            if count > 4:
                ub_mean = round(statistics.mean(ub_scores))
                # TODO
                ct_dict[c][celltype.rstrip()] = [round(sum(ranks) / math.sqrt(gene_count)), count, gene_count,
                                        ub_mean]

    for ct in cm_dict.keys():
        for gene in cm_dict[ct].keys():
            db_genes.append(gene)

    data_genes = list(set(data_genes))
    data_hits = list(set(data_hits))
    db_genes = list(set(db_genes))

    print(f"The database contains {str(len(db_genes))} different genes.\
          \nThe input data contains {str(len(data_genes))} different genes.\
          \nThe genes of the input data overlap with {str(len(data_hits))} genes in total, {str(round(len(data_hits) / len(db_genes), 2) * 100)} percent.")

    return ct_dict


def get_cell_types(cluster_path, db_path, tissue="all", species="Hs", header=False):
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

    Returns
    -------
    dictionary :
        The dictionary which contains the scores, the quantity of hits, the overall marker genes and
        the ubiquitousness index per cell type for each cluster.
    """

    db_dict = parse_marker_database(db_path, tissue=tissue, species=species, header=header)
    annotated_clusters = get_annotated_clusters(cluster_path=cluster_path)

    return calc_ranks(db_dict, annotated_clusters)


def get_annotated_clusters(cluster_path):
    """
    Read cluster files and sum ranked scores if genes appear more than once per file.

    Parameters
    ----------
    cluster_path : string
        The path to the folder which contains the "cluster files": Tab-separated files containing the
        genes and the corresponding ranked scores.

    Returns
    -------
    dictionary :
        Dictionary which contains the summed up ranked scores per gene for each cluster.
    """

    annotated_clusters = {}
    files = os.listdir(cluster_path)
    for file in [x for x in files if not x.startswith(".")]:
        cname = file.split(".cluster_")[1]
        annotated_dict = {}
        with open(cluster_path + file) as cfile:
            annotated_dict[cname] = []
            lines = cfile.readlines()
            for line in lines:
                split = line.split("\t")
                if len(split) == 2:
                    annotated_dict[cname].append(
                        [split[0].upper(), float(split[1].rstrip())])

        sum_dict = {}
        for gene in annotated_dict[cname]:
            if gene[1] > 0:  # Only consider positive values
                if gene[0] in sum_dict.keys():
                    sum_dict[gene[0]] += gene[1]
                else:
                    sum_dict[gene[0]] = gene[1]

        annotated_clusters[cname] = sum_dict

    return annotated_clusters



def perform_cell_type_annotation(output, db_path, cluster_path, tissue="all", species="Hs", header=False):
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
    """

    opath = output + "/ranks/"
    if not os.path.exists(opath):
        os.makedirs(opath)

    ct_dict = get_cell_types(cluster_path, db_path, tissue, species=species, header=header)
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
        for dic in ct_dict.keys():
            sorted_dict = dict(
                sorted(ct_dict[dic].items(), key=lambda r: (r[1][0], r[1][1]), reverse=True))
            with open(output + "/ranks/" + "cluster_" + dic, "w") as d_file:
                for key in sorted_dict.keys():
                    d_file.write(key)
                    for value in sorted_dict[key]:
                        d_file.write("\t" + str(value))
                    d_file.write("\n")
            if len(sorted_dict.keys()) > 0:
                c_file.write(dic + "\t" + str(next(iter(sorted_dict))) + "\n")


def validate_settings(repo_path, adata, organism, rank_genes_column, genes_column, column, ensembl, col_to_search, search_terms):
    """
    Validates user settings including file paths, anndata object columns, and specified organism.

    Parameters
    ----------
    repo_path : str
        Path to the marker repository.
    adata : anndata.AnnData
        The loaded AnnData object.
    organism : str or int
        Organism name, taxon ID, or both. E.g., "mouse", 10090, or "mouse 10090".
    rank_genes_column : str, default None
        Column in .obs table where ranked genes are stored. None if no ranking performed yet.
    genes_column : str, default None
        Column in .var table where gene symbols or IDs are stored. None if index column already suits the need.
    column : str
        The column in .obs table of the clustering you want to annotate. E.g., "leiden" or "louvain".
    ensembl : bool
        True if the index of .var tables are Ensembl IDs, False otherwise.
    col_to_search : str
        Column to search in for marker list selection. None to search all columns.
    search_terms : list of str
        Search terms for marker list selection. "-" to exclude, "+" must contain.

    Returns
    -------
    bool :
        Returns True if all settings are valid, otherwise prints the errors and returns False.
    """ 

    errors = []

    # Validate if repo_path exists
    if not os.path.exists(repo_path):
        errors.append(f"Repo path {repo_path} does not exist.")

    if adata is None:
        errors.append("No AnnData object provided.")
        
    # List of valid organisms and their tax IDs
    get_whitelists(repo_path=repo_path)
    valid_organisms = read_whitelist("organism", repo_path=repo_path)['whitelist']

    # Validate the organism
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

    if column and column not in adata.obs.columns:
        formatted_obs_columns = "\n  - " + "\n  - ".join(adata.obs.columns)
        errors.append(f"Invalid column {column}.\nAvailable columns in adata.obs:{formatted_obs_columns}\n")

    # Validate col_to_search using the columns from the combined DataFrames in the repo
    combined_df_columns = list(combine_dfs(repo_path=repo_path).columns)
    if col_to_search and col_to_search not in combined_df_columns:
        formatted_combined_df_columns = "\n  - " + "\n  - ".join(combined_df_columns)
        errors.append(f"Invalid col_to_search {col_to_search}.\nAvailable columns:{formatted_combined_df_columns}\n")


    if errors:
        print("Validation failed due to the following errors:")
        print("-" * 40)
        for i, error in enumerate(errors, 1):
            print(f"{i}. {error}")
        print("-" * 40)

        return False
    else:
        print("All settings are valid.")
        print(f"Summary of settings:")
        print("-" * 40)
        print(f"  Repo path: {repo_path}")
        print(f"  Organism: {organism}")
        print(f"  Rank genes column: {rank_genes_column}")
        print(f"  Genes column: {genes_column}")
        print(f"  Clustering column: {column}")
        print(f"  Ensembl IDs: {ensembl}")
        print(f"  Column to search: {col_to_search}")
        print(f"  Search terms: {search_terms}")
        print("-" * 40)

        return True
    

def list_possible_settings(repo_path, adata):
    """
    Lists all possible settings based on the repo and AnnData object.

    Parameters
    ----------
    repo_path : str
        Path to the marker repository.
    adata : anndata.AnnData
        The loaded AnnData object.
    """

    if not os.path.exists(repo_path):
        print(f"Repo path {repo_path} does not exist.")
        return

    print("Possible Settings:")
    print("-" * 40)
    
    valid_organisms = read_whitelist("organism", repo_path=repo_path)['whitelist']
    print("1. Possible Organisms or Taxon IDs:")
    for org in valid_organisms:
        print(f"  - {org}")

    print("\n2. Available Columns in adata.obs:")
    for col in adata.obs.columns:
        print(f"  - {col}")

    print("\n3. Available Columns in adata.var:")
    for col in adata.var.columns:
        print(f"  - {col}")

    combined_df_columns = list(combine_dfs(repo_path=repo_path).columns)
    print("\n4. Available Columns to Search in Marker Repository:")
    for col in combined_df_columns:
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


def compare_cell_types(adata, column, marker_lists):
    """
    Group the observation DataFrame of an anndata object by a specific column and include relevant cell types from given marker lists.

    Parameters
    ----------
    adata : anndata.AnnData
        The anndata object containing the .obs DataFrame.
    column : str
        The column by which to group the .obs DataFrame.
    marker_lists : list of str
        List of paths to marker list files.

    Returns
    -------
    DataFrame :
        A grouped DataFrame based on the specified column, incorporating relevant cell types from the marker lists.
    """
    
    obs_df = adata.obs
    columns_to_keep = [column]
    
    for marker_list in marker_lists:
        name = marker_list.split("/")[-1]
        columns_to_keep.append(f"cell_types_{name}")
        
    filtered_obs_df = obs_df[columns_to_keep]
    grouped_obs_df = filtered_obs_df.groupby(column).agg('first')
    
    return grouped_obs_df