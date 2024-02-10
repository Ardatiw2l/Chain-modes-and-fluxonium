# Fluxonium Chain Mode Analysis

This repository contains the implementation of the ChainMode class, which models the modes in the inductive chain of a fluxonium qubit. The implementation is based on the theoretical framework outlined in the paper:

- Giovanni Viola and Gianluigi Catelani, "Collective modes in the fluxonium qubit." Physical Review B, vol. 92, no. 224511, 2015. DOI: 10.1103/PhysRevB.92.224511.

## Installation

To set up and use the ChainMode class for modeling fluxonium qubit modes, please follow these detailed steps:

### Step 1: Clone the Repository

First, you'll want to clone this repository to your local machine. To do so, open your terminal (or Command Prompt/PowerShell if you're on Windows) and run the following command:

```bash
git clone https://github.com/yourusername/yourrepository.git
```
Make sure to replace yourusername/yourrepository with the actual username and repository name where the code is hosted.
Step 2: Navigate to the Repository Directory

Once the repository has been cloned, navigate to the newly created directory by running:

```bash

cd yourrepository
```
You should now be in the root directory of the project.
Step 3: Create the Conda Environment

This project uses an environment.yml file to manage dependencies, ensuring that you have the right packages and versions to run the code. To create a Conda environment with these dependencies, execute:

```bash

conda env create -f environment.yml
```
This command will create a new environment and install all the necessary packages as specified in the environment.yml file.
Step 4: Activate the Conda Environment

After creating the environment, you need to activate it. You can do this with the following command:

```bash

conda activate environment-name
```
Replace environment-name with the name of the environment that was created from the environment.yml file. This name is typically found at the top of the environment.yml file under the name field.
Step 5: Verify Installation

To verify that the environment has been set up correctly and all dependencies are installed, you can run:

```bash

conda list
```
This command lists all the packages installed in the active Conda environment, allowing you to check that everything is in order.

You are now ready to use the ChainMode class and other modules provided in this repository.

## Contributing

Contributions to the project are welcome. Please ensure that your code adheres to the project's coding standards and includes appropriate tests. For major changes, please open an issue first to discuss what you would like to change.


## Citation

If you use this code in your research, please cite the original paper:
```bibtex
@article{ViolaCatelani2015,
  title={Collective modes in the fluxonium qubit},
  author={Giovanni Viola and Gianluigi Catelani},
  journal={Physical Review B},
  volume={92},
  number={224511},
  year={2015},
  publisher={American Physical Society},
  doi={10.1103/PhysRevB.92.224511}
}
```

For the YAML file tutorial, I've included instructions on how to create a Conda environment using the provided `environment.yml` file. This assumes you have such a file at the root of your repository to manage dependencies. If your setup process is different, you'll need to adjust the instructions accordingly.