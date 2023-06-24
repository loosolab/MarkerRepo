import seaborn as sns
import matplotlib.pyplot as plt

def plot_gene_counts(gene_counts_df, file_name="gene_counts_barplot.png", dpi=300):
    """
    Plots a bar plot of the gene counts.

    Parameters
    ----------
    gene_counts_df : pd.DataFrame
        DataFrame containing the counts of target genes for each source gene.
    file_name : str, default "gene_counts_barplot.png"
        File name of plot.
    dpi : int, default 300
        The resolution of the exported plot.
    """

    sns.set_theme(style="darkgrid")
    gene_counts_df = gene_counts_df[gene_counts_df['Count target genes'] >= 2]

    # Group by the 'Count target genes' column and calculate the size of each group
    gene_counts = gene_counts_df.groupby('Count target genes').size().reset_index(name='Number of source genes')

    fig = sns.barplot(data=gene_counts, x='Count target genes', y='Number of source genes')
    plt.title("Number of source genes by count of target genes")
    
    # Save the figure with a dpi of 300
    fig.get_figure().savefig(file_name, dpi=dpi)

    plt.show()
