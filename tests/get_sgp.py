import numpy as np
from flare_pp._C_flare import SparseGP, NormalizedDotProduct, B2, Structure
from flare_pp.sparse_gp import SGP_Wrapper
from flare_pp.sparse_gp_calculator import SGP_Calculator
from flare.ase.atoms import FLARE_Atoms
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.build import make_supercell

# Define kernel.
sigma = 2.0
power = 1.0
kernel = NormalizedDotProduct(sigma, power)

# Define remaining parameters for the SGP wrapper.
sigma_e = 0.001
sigma_f = 0.05
sigma_s = 0.006
species_map = {6: 0, 8: 1}
single_atom_energies = {0: -5, 1: -6}
variance_type = "local"
max_iterations = 20
opt_method = "L-BFGS-B"
bounds = [(None, None), (sigma_e, None), (None, None), (None, None)]


def get_random_atoms(a=2.0, sc_size=2, numbers=[6, 8],
                     set_seed: int = None):

    """Create a random structure."""

    if set_seed:
        np.random.seed(set_seed)

    cell = np.eye(3) * a
    positions = np.array([[0, 0, 0], [a/2, a/2, a/2]])
    unit_cell = Atoms(cell=cell, positions=positions, numbers=numbers,
                      pbc=True)
    multiplier = np.identity(3) * sc_size
    atoms = make_supercell(unit_cell, multiplier)
    atoms.positions += (2 * np.random.rand(len(atoms), 3) - 1) * 0.1
    flare_atoms = FLARE_Atoms.from_ase_atoms(atoms)

    return flare_atoms


def get_isolated_atoms(numbers=[6, 8]):

    """Create a random structure."""

    a = 30.0
    cell = np.eye(3) * a
    positions = np.array([[0, 0, 0], [a/2, a/2, a/2]])
    unit_cell = Atoms(cell=cell, positions=positions, numbers=numbers,
                      pbc=True)
    atoms = unit_cell
    flare_atoms = FLARE_Atoms.from_ase_atoms(atoms)

    return flare_atoms


def get_empty_sgp(n_types=2, power=2, multiple_cutoff=False):
    kernel.power = power

    # Define B2 calculator.
    cutoff = 5.0
    cutoff_function = "quadratic"
    radial_basis = "chebyshev"
    radial_hyps = [0.0, cutoff]
    cutoff_hyps = []
    cutoff_matrix = cutoff * np.ones((n_types, n_types))
    if multiple_cutoff:
        cutoff_matrix += np.eye(n_types) - 1

    descriptor_settings = [n_types, 3, 2]
    b2_calc = B2(radial_basis, cutoff_function, radial_hyps, cutoff_hyps,
                 descriptor_settings, cutoff_matrix)

    empty_sgp = SGP_Wrapper(
        [kernel], [b2_calc], cutoff, sigma_e, sigma_f, sigma_s, species_map,
        single_atom_energies=single_atom_energies, variance_type=variance_type,
        opt_method=opt_method, bounds=bounds, max_iterations=max_iterations
    )

    return empty_sgp


def get_updated_sgp(n_types=2, power=2, multiple_cutoff=False):
    if n_types == 1:
        numbers = [6, 6]
    elif n_types == 2:
        numbers = [6, 8]

    sgp = get_empty_sgp(n_types, power, multiple_cutoff)

    # add a random structure to the training set
    training_structure = get_random_atoms(numbers=numbers)
    training_structure.calc = LennardJones()

    forces = training_structure.get_forces()
    energy = training_structure.get_potential_energy()
    stress = training_structure.get_stress()

    sgp.update_db(training_structure, forces, custom_range=(1, 2, 3, 4, 5),
                  energy=energy, stress=stress, mode="specific")

    # add an isolated atom to the training data
    training_structure = get_isolated_atoms(numbers=numbers)
    training_structure.calc = LennardJones()

    forces = training_structure.get_forces()
    energy = training_structure.get_potential_energy()
    stress = training_structure.get_stress()

    sgp.update_db(training_structure, forces, custom_range=(0,),
                  energy=energy, stress=stress, mode="specific")

    print("sparse_indices", sgp.sparse_gp.sparse_indices)

    return sgp


def get_sgp_calc(n_types=2, power=2, multiple_cutoff=False):
    sgp = get_updated_sgp(n_types, power, multiple_cutoff)
    sgp_calc = SGP_Calculator(sgp)

    return sgp_calc
