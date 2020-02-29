# ======================================================================================================================
# * Weighted Holistic Atom Localization and Entity Shape (WHALES) descriptors *
#   v. 1, May 2018
# ----------------------------------------------------------------------------------------------------------------------
# This file contains all the necessary files to handle molecular properties and coordinates.
#
# Francesca Grisoni, May 2018, ETH Zurich & University of Milano-Bicocca, francesca.grisoni@unimib.it
# please cite as xxxx
# ======================================================================================================================

import numpy as np


def get_coordinates_and_prop(mol, property_name='partial_charges', do_charge=True):
    """
    Extracts all of the useful chemical information, i.e., the partial charge and the coordinates and formats it
    for atomic centred covariance matrix calculation.
    ====================================================================================================================
    :param
    mol: rdkit molecule
    do_charge: if True, the charges are computed
    do_geom: if True, it calculates MMF 3D coordinates
    :returns
    coords (n_atoms x 3): geometrical matrix (x-y-z coords)
    w (n_atoms x 1): partial charge array
    ====================================================================================================================
    Francesca Grisoni, 05/2018, v. beta
    ETH Zurich
    """

    # molecule preparation
    mol, property_name, err = prepare_mol(mol, property_name, do_charge)

    if err == 0:
        # pre-allocation
        n_at = mol.GetNumAtoms()  # num atoms
        coords = np.zeros((n_at, 3))  # init coords
        w = np.zeros((n_at, 1))  # init weights

        # coordinates and property
        for atom in range(n_at):  # loops over atoms, gets 3D coordinate matrix

            # gets atomic positions
            pos = mol.GetConformer().GetAtomPosition(atom)
            coords[atom, ] = [pos.x, pos.y, pos.z]

            # gets atomic properties
            w[atom] = mol.GetAtomWithIdx(atom).GetProp(property_name)
            
        # checks the weight values computed and throws and error if they are all 0
        if all(v == 0 for v in w):
            err = 1
    else:
        coords = []
        w = []

    return coords, w, err

# ----------------------------------------------------------------------------------------------------------------------


def prepare_mol(mol, property_name, do_charge):
    """
    Sets atomic properties if they are specified in the sdf, otherwise computes them. If specified, computes 3D coordinates
    using MMF.  The default number of iterations is 200, but it is progressively increased to 5000 (with a step of 500)
    in case convergence is not reached.
    ====================================================================================================================
    :param
    mol: molecule to be analyzed (from rdkit supplier)
    property_name: name of the property to be used
    do_charge: if True, partial charge is computed
    do_geom: if True, molecular geometry is optimized
    :return:
    mol: molecule with property and 3D coordinates (H depleted)
    property_name: updated on the basis of the settings
    ====================================================================================================================
    Francesca Grisoni, 12/2016, v. alpha
    ETH Zurich
    """

    from rdkit.Chem import AllChem as Chem
    err = 0

    # partial charges
    if do_charge is False:
        if property_name is not '':
            err = check_mol(mol, property_name, do_charge)
            if err == 0:
                # prepares molecule
                # mol = Chem.AddHs(mol)
                mol = Chem.RemoveHs(mol)
                n_at = mol.GetNumAtoms()
                # takes properties
                list_prop = mol.GetPropsAsDict()
                string_values = list_prop[property_name]  # extracts the property according to the set name
                string_values = string_values.split("\n")
                w = np.asarray(map(float, string_values))
        else:
            mol = Chem.AddHs(mol)
            n_at = mol.GetNumAtoms()
            w = np.ones((n_at, 1))/n_at
            w = np.asarray(map(float, w))
            property_name = 'equal_w'
            err = 0
        # extract properties
        for atom in range(n_at):
            mol.GetAtomWithIdx(atom).SetDoubleProp(property_name, w[atom])

        mol = Chem.RemoveHs(mol)

    # Gasteiger-Marsili Charges
    elif (do_charge is True) and (err is 0):
        Chem.ComputeGasteigerCharges(mol)
        property_name = '_GasteigerCharge'
        err = check_mol(mol, property_name, do_charge)

    return mol, property_name, err


# ----------------------------------------------------------------------------------------------------------------------
def check_mol(mol, property_name, do_charge):
    """
    checks if the property is annotated and gives 0 if it is
    """
    n_at = mol.GetNumAtoms()
    if do_charge is False:
        list_prop = mol.GetPropsAsDict()
        string_values = list_prop[property_name]  # extracts the property according to the set name
        if string_values == '' or string_values == ['']:
            err = 1
        else:
            err = 0
    else:
        from rdkit.Chem import AllChem as Chem
        err = 0
        atom = 0
        while atom < n_at:
            value = mol.GetAtomWithIdx(atom).GetProp(property_name)
            # checks for error (-nan, inf, nan)
            if value == '-nan' or value == 'nan' or value == 'inf':
                err = 1
                break

            atom += 1

    # checks for the number of atoms
    if n_at < 4:
        err = 1

    return err

