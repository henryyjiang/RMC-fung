import numpy as np
import copy
from pymatgen.core.composition import Composition

def step(structure, latticeprob, σr, σl, σθ):
    new_struct = copy.deepcopy(structure)
    if np.random.rand() < latticeprob:
        lattice_step(new_struct, σl, σθ)
    else:
        positions_step(new_struct, σr)
    return new_struct

def lattice_step(structure, σl, σθ):
    structure.lattice = structure.lattice.from_parameters(
        np.maximum(0.0, structure.lattice.a + σl * np.random.randn()),
        np.maximum(0.0, structure.lattice.b + σl * np.random.randn()), 
        np.maximum(0.0, structure.lattice.c + σl * np.random.randn()), 
        structure.lattice.alpha + σθ * np.random.randn(), 
        structure.lattice.beta + σθ * np.random.randn(), 
        structure.lattice.gamma + σθ * np.random.randn()
    )

def positions_step(structure, σr):
    atom_i = np.random.choice(range(len(structure)))
    structure.sites[atom_i].a = (structure.sites[atom_i].a + 
        σr * np.random.rand() / structure.lattice.a) % 1
    structure.sites[atom_i].b = (structure.sites[atom_i].b + 
        σr * np.random.rand() / structure.lattice.b) % 1
    structure.sites[atom_i].c = (structure.sites[atom_i].c + 
        σr * np.random.rand() / structure.lattice.c) % 1

    if np.random.rand() < 0.05:
        comp = Composition("O")
        structure.append(comp.formula, [( np.random.rand() / structure.lattice.a) % 1, (np.random.rand() / structure.lattice.b) % 1, (np.random.rand() / structure.lattice.c) % 1])
    if np.random.rand() < 0.05:
        atom_b = np.random.choice(range(len(structure)))
        structure.remove_sites([atom_b])
def 𝛘2(exp, th, σ):
    #return np.mean((exp - th) ** 2) / (σ ** 2)
    return np.sum((exp - th) ** 2) / (σ ** 2)

def reject(structure):
    dists = structure.distance_matrix.flatten()
    return np.min(dists[dists > 0]) < 1

def rmc(optfunc, args, exp, σ, structure, N, latticeprob = 0.1, σr = 0.5, σl = 0.1, σθ = 1.0):
    structures = []
    𝛘2s = []
    accepts = []
    old_structure = structure
    old_𝛘2 = 𝛘2(exp, optfunc(old_structure, **(args)), σ)
    best = structure

    for _ in range(N):
        new_structure = step(structure, latticeprob, σr, σl, σθ)

        if 𝛘2(exp, optfunc(new_structure, **(args)), σ) < 𝛘2(exp, optfunc(best, **(args)), σ):
            best = new_structure

        new_𝛘2 = 𝛘2(exp, optfunc(new_structure, **(args)), σ)
        Δχ2 = new_𝛘2 - old_𝛘2
        #accept = np.random.rand() < np.exp(-Δχ2/2) and not reject(new_structure)
        #accept = Δχ2 < 0 or np.random.default_rng().random() > Δχ2*50
        accept = Δχ2 < 0

        structures.append(new_structure)
        𝛘2s.append(new_𝛘2)
        accepts.append(accept)

        if accept:
            #old_structure = copy.deepcopy(new_structure)
            old_structure = new_structure
            old_𝛘2 = new_𝛘2
        if _ % 1000 == 0:
            print(_)

    return structures, 𝛘2s, accepts, best
