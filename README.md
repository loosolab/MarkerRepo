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
