import time, os, shutil, glob, subprocess, sys
from copy import deepcopy
import pytest
import pkgutil, pyclbr
import importlib
import inspect
import numpy as np

from flare.learners.otf import OTF

from ase import units
import ase.calculators as ase_calculators
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation,
)
from ase import io
from ase.symbols import symbols2numbers

import yaml

with open(sys.argv[1], "r") as f:
    config = yaml.safe_load(f)

################################################################################
#                                                                              #
#                        Set up Supercell from ASE Atoms                       #
#                                                                              #
################################################################################

# parse parameters
atoms_config = config["supercell"]
atoms_file = atoms_config.get("file")
atoms_format = atoms_config.get("format")
replicate = ase_config.get("replicate", [1, 1, 1])
jitter = atoms_config.gt("jitter", 0)

super_cell = io.read(atoms_file, format=atoms_format)
super_cell *= replicate

# jitter positions to give nonzero force on first frame
for atom_pos in super_cell.positions:
    for coord in range(3):
        atom_pos[coord] += (2 * np.random.random() - 1) * jitter

################################################################################
#                                                                              #
#                           Set up ASE DFT calculator                          #
#                                                                              #
################################################################################

# find the module including the ASE DFT calculator class by name
dft_calc_name = config["dft_calc"]["name"]
dft_calc_kwargs = config["dft_calc"]["kwargs"]
dft_calc_params = config["dft_calc"]["params"]

dft_module_name = ""
for importer, modname, ispkg in pkgutil.iter_modules(ase_calculators.__path__):
    module_info = pyclbr.readmodule("ase.calculators." + modname)
    if dft_calc_name in module_info:
        dft_module_name = modname
        break

# import ASE DFT calculator module, and build a DFT calculator class object
dft_calc = None
dft_module = importlib.import_module("ase.calculators." + dft_module_name)
for name, obj in inspect.getmembers(dft_module, inspect.isclass):
    if name == dft_calc_name:
        dft_calc = obj(**dft_calc_kwargs)
dft_calc.set(**dft_calc_params)

################################################################################
#                                                                              #
#                           Set up ASE flare calculator                        #
#                                                                              #
################################################################################

flare_config = config["flare_calc"]
kernels = flare_config.get("kernels")
random_init_hyps = flare_config.get("random_init_hyps", True)
opt_algorithm = flare_config.get("opt_algorithm", "BFGS")
max_iterations = flare_config.get("max_iterations", 20)
bounds = flare_config.get("bounds", None)

if flare_config["gp"] == "GaussianProcess":
    from flare.bffs.gp import GaussianProcess
    from flare.bffs.mgp import MappedGaussianProcess
    from flare.bffs.gp.calculator import FLARE_Calculator
    from flare.utils.parameter_helper import ParameterHelper

    # create gaussian process model
    gp_parameters = flare_config.get("gp_parameters")
    n_cpus = flare_config.get("n_cpus", 1)
    use_mapping = flare_config.get("use_mapping", False)
    if use_mapping:
        grid_params = flare_config.get("grid_params")
        vap_map = flare_config.get("var_map", "pca")
    
    # set up GP hyperparameters
    pm = ParameterHelper(
        kernels=kernels,
        random=random_init_hyps,
        parameters=gp_parameters,
    )
    hm = pm.as_dict()
    
    gp_model = GaussianProcess(
        kernels=kernels,
        component="mc",
        hyps=hm["hyps"],
        cutoffs=hm["cutoffs"],
        hyp_labels=hm["hyp_labels"],
        opt_algorithm=opt_algorithm,
        n_cpus=n_cpus,
    )
    
    # create mapped gaussian process
    if use_mapping:
        unique_species = []
        for s in super_cell.symbols:
            if s not in unique_species:
                unique_species.append(s)
        coded_unique_species = symbols2numbers(unique_species)
        mgp_model = MappedGaussianProcess(
            grid_params=grid_params,
            unique_species=coded_unique_species,
            n_cpus=n_cpus,
            var_map=var_map,
        )
    else:
        mgp_model = None
    
    flare_calc = FLARE_Calculator(
        gp_model=gp_model,
        mgp_model=mgp_model,
        par=n_cpus > 1,
        use_mapping=use_mapping,
    )

elif flare_config["gp"] == "SGP_Wrapper":
    from flare_pp._C_flare import NormalizedDotProduct, B2, SquaredExponential
    from flare_pp.sparse_gp import SGP_Wrapper
    from flare_pp.sparse_gp_calculator import SGP_Calculator
    
    # Define kernels.
    kernels = []
    for k in flare_config["kernels"]:
        if k["name"] == "NormalizedDotProduct":
            kernels.append(NormalizedDotProduct(k["sigma"], k["power"]))
        elif k["name"] == "SquaredExponential":
            kernels.append(SquaredExponential(k["sigma"], k["ls"]))
        else:
            raise NotImplementedError(f"{k['name']} kernel is not implemented")
    
    # Define descriptor calculators.
    descriptors = []
    for d in flare_config["descriptors"]:
        if d["name"] == "B2":
            radial_hyps = [0.0, flare_config["cutoff"]]
            cutoff_hyps = []
            descriptor_settings = [len(flare_config["species"]), d["nmax"], d["lmax"]]
            b2_calc = B2(
                d["radial_basis"], 
                d["cutoff_function"], 
                radial_hyps, 
                cutoff_hyps,
                descriptor_settings
            )
            descriptors.append(b2_calc)
   
    # Define remaining parameters for the SGP wrapper.
    species_map = {flare_config["species"][i]: i for i in range(len(flare_config["species"]))}
    single_atom_energies = {i: flare_config["single_atom_energies"][i] for i in range(len(flare_config["single_atom_energies"]))}
    assert len(species_map) == len(single_atom_energies)

    sgp = SGP_Wrapper(
        kernels, 
        descriptors, 
        flare_config["cutoff"], 
        flare_config["energy_noise"], 
        flare_config["forces_noise"], 
        flare_config["stress_noise"], 
        species_map,
        single_atom_energies=single_atom_energies, 
        variance_type=flare_config["variance_type"],
        opt_method=opt_algorithm, 
        bounds=bounds,
        max_iterations=max_iterations,
    )

    flare_calc = SGP_Calculator(sgp)

else:
    raise NotImplementedError(f"{flare_config['gp']} is not implemented")

################################################################################
#                                                                              #
#                           Set up OTF training engine                         #
#                                                                              #
################################################################################

# intialize velocity
temperature = config["otf"].get("temperature", 0)
MaxwellBoltzmannDistribution(super_cell, temperature * units.kB)
Stationary(super_cell)  # zero linear momentum
ZeroRotation(super_cell)  # zero angular momentum

test_otf = OTF(
    super_cell,
    flare_calc=flare_calc,
    dft_calc=dft_calc,
    **config.get("otf"),
)

# run on-the-fly training
test_otf.run()
