### Scaffold hopping by holistic molecular descriptors in drug design
Francesca Grisoni,* Gisbert Schneider

ETH Zurich, Department of Chemistry and Applied Biosciences, RETHINK, Vladimir-Prelog-Weg 4, 8093, Zurich, Switzerland.

*Corresponding author: francesca.grisoni@pharma.ethz.ch 

## Preliminary steps
In order to download and run the code, ensure you have the following software on your machine:
•	Anaconda (for Python 3.7): www.anaconda.com
•	Git: https://www.atlassian.com/git/tutorials/install-git. 

## Getting the code
The open source code we will use for the current example is available as a GitHub repository at the following link: https://github.com/grisoniFr/scaffold_hopping_whales. GitHub is a version control platform (66) based on Git, which is useful to host and review code, manage projects, and build and share software. Users can clone the repository (i.e., get a local copy of the repository contents) with the following command on a Linux/Mac terminal or Windows command line:
git clone https://github.com/grisoniFr/scaffold_hopping_whales
A copy of the repository will be generated on your local machine, in the dedicated GitHub folder. Move to the repository with the following commands:
cd <path/to/folder> (from a Mac or Linux terminal)
cd <path\to\folder> (from Windows command line)
The folder, named “scaffold_hopping”, will contain the following elements:
•	README.md file, containing instructions on the repository on the installation. 
•	scaffold_hopping.yml, needed to install the necessary dependencies. 
•	code (folder), which contains all the necessary scripts to execute the calculations, along with the Jupyter notebook.
•	data (folder), containing the virtual screening dataset provided for this example.

## Setting up the virtual environment
Performing all the calculations within a virtual environment is recommended, to ensure that all the dependencies are properly installed, with no conflicts with previously installed Python packages. The set of necessary packages can be installed by creating the virtual environment with the provided “scaffold_hopping.yml” file:
conda env create -f scaffold_hopping.yml
To use the installed packages, activate the environment:
conda activate scaffold_hopping.yml

## Use the provided Jupyter notebook
Move to the “code” folder, where the Jupyter notebook file is contained, and launch Jupter Notebook, as follows:
``
cd code
``


