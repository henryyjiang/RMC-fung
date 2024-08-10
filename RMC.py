import activestructopt.dataset.materialsproject
import activestructopt.dataset.rdf
import activestructopt.optimization.rmc.rmc
import activestructopt.optimization.mcmc.mcmc
from pymatgen.core.composition import Composition
import numpy as np
import matplotlib.pyplot as plt

target_structure = activestructopt.dataset.materialsproject.get_structure('mp-22526', 'yvC2UCUjpLKJgwZ4Vhx5xFKHvVfOiF7k')
#target_structure = activestructopt.dataset.materialsproject.get_structure('mp-13703', 'yvC2UCUjpLKJgwZ4Vhx5xFKHvVfOiF7k')
#target_structure = activestructopt.dataset.materialsproject.get_structure('mp-860963', 'yvC2UCUjpLKJgwZ4Vhx5xFKHvVfOiF7k')
starting_structure = target_structure.copy()
starting_structure.perturb(0.5)
print(target_structure)
rs = np.arange(0.5, 12.0, 0.01)
exp = activestructopt.dataset.rdf.get_rdf(target_structure, σ = 0.1)
exp2 = activestructopt.dataset.rdf.get_rdf(starting_structure, σ = 0.1)

plt.plot(rs, exp2, label="starting")
plt.plot(rs, exp, label="target")
plt.legend()
plt.show()

structures, 𝛘2s, accepts, best = activestructopt.optimization.rmc.rmc.rmc(
  activestructopt.dataset.rdf.get_rdf,
  {'σ': 0.1},#0.05
  exp,
  1,#0.5
  starting_structure,
  5000,
  latticeprob = 0.,
  σr = 0.05,#0.1
)

print(sum(accepts))

plt.plot(rs, exp, label = "target")
plt.plot(rs, activestructopt.dataset.rdf.get_rdf(structures[0], σ = 0.1), label = "starting random positions")
#plt.plot(rs, activestructopt.dataset.rdf.get_rdf(structures[np.nonzero(accepts)[-1][-1]], σ = 0.1), label = "last accepted positions")
plt.plot(rs, activestructopt.dataset.rdf.get_rdf(best, σ = 0.1), label = "best positions")
plt.legend()
plt.xlabel('r (Å)')
plt.ylabel('Radial Distribution Function')
plt.show()

plt.scatter(range(len(𝛘2s)), 𝛘2s, color = ['g' if a else 'r' for a in accepts])
plt.xlabel('steps')
plt.ylabel('X2')
plt.show()