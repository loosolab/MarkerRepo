# Quickstart

The easiest way to get started with this project is by running the Jupyter notebooks available in the `notebooks` folder. Before you can use these notebooks, make sure to follow the steps outlined below to set up your Conda environment and install necessary packages.

## Pre-requisites: Set up your Conda environment

If you haven't already, make sure to set up your Conda environment by following the steps in the [Install the Conda environment and add it to Jupyter as a Kernel](#install-the-conda-environment-and-add-it-to-jupyter-as-a-kernel) section.

## Notebook Descriptions

This section provides an overview of the Jupyter notebooks available in the `notebooks` folder and what each notebook is designed to do.

### annotation.ipynb: Cell Type Annotation

The `annotation.ipynb` notebook serves as a comprehensive tool for annotating clustered h5ad files using the Marker Repo package. 

The notebook guides you through the entire annotation workflow. It starts by loading essential packages and setting up your Anndata object. Next, it offers options for gene ranking if not already done, followed by marker list creation and cell type annotation.

**Key Features:**

- Validating your initial settings and Anndata object.
- Ranking genes within clusters if not pre-ranked.
- Creating customized marker lists based on organism and other filtering criteria.
- Annotating cell types using the created marker lists.
- Visualizing annotation results on UMAP plots.

**When to Use:**

- When you have clustered h5ad files that require cell type annotation.
- When you want to utilize customized marker lists for annotation.

### guided_annotation.ipynb: Step-by-Step Cell Type Annotation

The `guided_annotation.ipynb` notebook is an interactive variant of the `annotation.ipynb` notebook. Unlike its counterpart, which is designed for a seamless, uninterrupted workflow, this notebook guides you through the annotation process step-by-step. It allows you to manually enter parameters at different stages, making it ideal for experimental setups where you may be unsure of the optimal parameters to use.

**When to Use:**
you're in the exploratory phase and wish to experiment with different parameters.
- If y
- If ou're not yet certain about the ideal settings for annotation and want a guided, interactive process.

**Note:** If you already know all the parameters you want to use and prefer a notebook that can be executed in one go, it's recommended to use the `annotation.ipynb` notebook instead.

### homology.ipynb: Transfer Marker Genes Across Organisms Using Homology

The `homology.ipynb` notebook enables the transfer of marker genes from a source organism to a target organism using homology-based methods. It offers two distinct approaches: one utilizing BioMart and Ensembl, and the other using the HomoloGene database. The notebook is structured to guide you through each step of the process, from selecting the appropriate method and organisms to the actual gene transfer.

**Key Features:**
- Interactive selection of source markers, target organisms, and other parameters.
- Wrapper function to streamline the gene transfer process.
- Offers the ability to select from a range of supported organisms, each backed by either Ensembl, or HomoloGene

**When to Use:**

- When you need to transfer marker genes from one organism to another, especially 
- when no marker lists are available for the organism you are analyzing, and you need a reliable method for generating them.

### scoring.ipynb: Scoring

The `scoring.ipynb` notebook allows for the weighting of individual markers using various functions. It enables the scoring of markers based on ubiquitousness index. 

**Key Features:**
- Enables comparison and scoring of markers from selected lists.
- Export functionality for scored marker lists.

**When to Use:**
- When you need to prioritize or weight marker genes for your analyses.
- When you are working with marker lists from multiple sources and need to generate a unified, weighted list.

### search_and_combine_lists.ipynb: Search and Combine Lists

The `search_and_combine_lists.ipynb` notebook provides a guide for searching and combining marker lists from the marker repository. It covers the process from selection to export, with wrapper functions streamlining the steps for various output formats and needs.

**Key Features:**
- Guided Search to compile a selection of marker lists based on metadata or other criteria.
- Combining and formatting selected lists into a unified DataFrame.
- Export options for different formats such as "two_column" or "score".
- Wrapper functions to automate the entire process for custom requirements.

**When to Use:**
- To familiarize yourself with the core functionalities of the Marker Repo.
- When you need a tailored marker list, assembled from multiple sources or based on specific criteria.
- When you require custom formatting styles for your marker lists.

### submit_lists.ipynb: Create and Upload Your Own Marker Lists

The `submit_lists.ipynb` notebook is designed to allow you to upload your own marker lists to the repository. This notebook takes you through every step, from entering metadata about your list to actually submitting it. It even supports validation and curation using whitelists, ensuring that your list is both accurate and useful.

**Key Features:**
- Step-by-step guide to entering metadata for your marker list.
- Utilizes whitelists for marker validation and curation.
- Allows for the automated entry of markers, which can be in the form of genes or genomic regions.
- Final validation step before publishing, ensuring that your list adheres to the repository's format and standards.

**When to Use:**
- When you have developed a unique marker list and wish to share it with the broader scientific community.
- When you want to ensure the quality and reliability of your marker list through guided validation steps.
- For a seamless process of contributing to the repository, making your work easily accessible and usable by others.

# Install the Conda environment and add it to Jupyter as a Kernel

Follow these steps to install the Conda environment and add it to Jupyter as a Kernel.

## Step 1: Install the Conda environment

1. Open your terminal
2. Navigate to the `environment.yaml` file in the repository
3. Run one of the following commands:

`conda env create -f environment.yaml`


`mamba env create -f environment.yaml`

This command creates a new Conda environment named `marker-repo` as described in the `environment.yaml` file.

## Step 2: Activate the Conda environment

You need to activate the environment before using it or adding it as a Kernel to Jupyter. To do this, run:

`conda activate marker-repo`


## Step 3: Add the environment to Jupyter as a Kernel

Add `marker-repo` as a Kernel in Jupyter by running:

`python -m ipykernel install --user --name=marker-repo`


After this step, `marker-repo` should appear as an available option when choosing a Kernel in Jupyter.

## Step 4: Deactivate the Conda environment

Once you are done, you can deactivate the `marker-repo` environment by running:

`conda deactivate`

That's it! You have successfully created a Conda environment, activated it, and added it as a Kernel in Jupyter. Select the `marker-repo` Kernel in the provided notebooks and begin working on your project.

# Installing the MarkerRepo Package

The MarkerRepo package allows you to use the functionalities of MarkerRepo in external environments. Follow these steps to install the package:

## Step 1: Activate the Conda environment

Before installing the package, make sure to activate the Conda environment where you want the package installed. Run:

`conda activate your-env`

## Step 2: Install the MarkerRepo package

With the Conda environment activated, navigate to the root directory of the MarkerRepo package and run:

`pip install .`

This command installs the MarkerRepo package into your `your-env` Conda environment.

Now you can import and use the MarkerRepo package in any Python script or notebook running in the `your-env` Conda environment.
