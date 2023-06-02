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
