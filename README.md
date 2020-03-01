# Scaffold hopping by holistic molecular descriptors in drug design
Francesca Grisoni,* Gisbert Schneider

ETH Zurich, Department of Chemistry and Applied Biosciences, RETHINK, Vladimir-Prelog-Weg 4, 8093, Zurich, Switzerland.

*Corresponding author: francesca.grisoni@pharma.ethz.ch 

## Preliminary steps
In order to download and run the code, ensure you have the following software on your machine: <div>
*	Anaconda (for Python 3.7): www.anaconda.com <div>
*	Git: https://www.atlassian.com/git/tutorials/install-git. 

## Getting the code
Clone the repository, as follows:
``
git clone https://github.com/grisoniFr/scaffold_hopping_whales
``
A copy of the repository will be generated on your local machine, in the dedicated GitHub folder. Move to the donwloaded repository to start using it. 

## Setting up the virtual environment
Performing all the calculations within a virtual environment is recommended. You can import the environment information (as provided in the “scaffold_hopping.yml” file) as follows:

``
conda env create -f scaffold_hopping.yml
``
To use the installed packages, activate the environment:
``
conda activate scaffold_hopping
``

## Use the provided Jupyter notebook
Move to the [code](/code) folder, where the Jupyter notebook file is contained, and launch Jupter Notebook, as follows:
``
jupyter notebook
``
Click on the notebook file "virtual_screening_pipeline.jpynb". There, you will find additional information on the required calculation steps.

