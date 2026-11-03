"""Bimetallic Pt/Au chemical-ordering optimization showcase.

A two-element twin of ``examples/multimet_go.ipynb``. It compares how reliably and
how quickly two optimizers find the lowest-energy Pt/Au ordering, over several
independent restarts:

  - canonical Metropolis Monte Carlo (random swaps), and
  - guided TOP basin hopping (optimal-exchange swaps, ``model="TOP"``).

Pipeline:
  1. Build a set of Pt/Au truncated-octahedron nanoparticles and label them with
     single-point EMT energies.
  2. Fit two linear surrogate energy models on that *same* EMT data:
       - an ETOP model (Bayesian ridge regression over the global ETOP
         descriptor) that drives the Monte Carlo run -- only the ETOP descriptor
         supports the incremental feature updates ``run_monte_carlo`` needs; and
       - a TOP model (Bayesian ridge regression over the global TOP descriptor),
         converted to per-atom (TEC) environment energies that drive guided
         basin hopping.
  3. Optimize a fresh particle of the same size with both methods, repeated over
     several restarts.
  4. Because the two optimizers minimize *different* surrogates, their raw
     energies are not comparable. Each restart's best ordering is therefore
     re-scored with single-point EMT (the ground truth both surrogates
     approximate). The global minimum is the lowest EMT energy found by either
     method, and a restart "succeeds" if it reaches within ``--tol`` of it.
  5. Emit a single overview figure (lowest-energy structure + cumulative success
     rate vs. steps on a log axis, drawn as a staircase) and save the raw run
     data for reproducibility.

    python examples/ptau_go.py            # full, large-particle run (slow)
    python examples/ptau_go.py --quick     # fast smoke test
"""
import argparse
import os

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from ase.io import write as ase_write
from ase.visualize.plot import plot_atoms

from npl.core import Nanoparticle
from npl.calculators import (
    EMTCalculator,
    BayesianRRCalculator,
    compute_coefficients_for_linear_topological_model,
)
from npl.descriptors import ExtendedTopologicalFeaturesClassifier
from npl.descriptors.global_feature_classifier import testTopologicalFeatureClassifier
from npl.monte_carlo.monte_carlo_etop import run_monte_carlo
from npl.optimization import run_basin_hopping

matplotlib.use("Agg")  # force a headless backend before any figure is created

STOICHIOMETRY = {"Pt": 0.5, "Au": 0.5}
SYMBOLS = list(STOICHIOMETRY.keys())

FULL = dict(height=6, trunc=1, n_train=150, n_runs=10,
            mc_temperature=300, mc_max_steps=100000,
            bh_hopping_attempts=10000, bh_hops=10)
QUICK = dict(height=5, trunc=1, n_train=20, n_runs=3,
             mc_temperature=300, mc_max_steps=1000,
             bh_hopping_attempts=10, bh_hops=10)


def new_particle(height, trunc):
    p = Nanoparticle()
    p.truncated_octahedron(height, trunc, STOICHIOMETRY)
    return p


def build_training_set(n_particles, height, trunc, emt, etop, top):
    """Label n_particles random orderings with single-point EMT energy and both
    the ETOP and TOP global feature vectors, so both surrogates fit identical
    data."""
    training_set = []
    for _ in range(n_particles):
        p = new_particle(height, trunc)
        emt.compute_energy(p)
        etop.compute_feature_vector(p)
        top.compute_feature_vector(p)
        training_set.append(p)
    return training_set


def fit_mc_model(training_set, etop):
    """ETOP Bayesian ridge model that powers the Monte Carlo run."""
    calculator = BayesianRRCalculator(etop.get_feature_key())
    calculator.fit(training_set, "EMT", validation_set=0.1)
    return calculator


def fit_bh_model(training_set, top, n_atoms):
    """TOP Bayesian ridge model, converted to the per-atom (TEC) environment
    energies that drive guided basin hopping (the basin_hopping_top recipe)."""
    calculator = BayesianRRCalculator(top.get_feature_key())
    calculator.fit(training_set, "EMT")
    coefficients, environment_energies = compute_coefficients_for_linear_topological_model(
        calculator.get_coefficients(), SYMBOLS, n_atoms=n_atoms)
    calculator.set_coefficients(coefficients)
    calculator.set_feature_key("TEC")
    return calculator, environment_energies


def emt_score(particle, emt):
    """Single-point EMT energy of an ordering -- the common yardstick."""
    emt.compute_energy(particle)
    return particle.get_energy("EMT")


def optimize_with_mc(cfg, calculator, etop, emt):
    """Run independent MC restarts. Returns one (emt_energy, step, best_particle)
    record per restart, where `step` is when the run found its best ordering."""
    records = []
    for _ in range(cfg["n_runs"]):
        start = new_particle(cfg["height"], cfg["trunc"])
        best, accepted = run_monte_carlo(cfg["mc_temperature"], cfg["mc_max_steps"],
                                         start, calculator, etop)
        _, step = min(accepted, key=lambda es: es[0])
        records.append((emt_score(best, emt), step, best))
    return records


def optimize_with_bh(cfg, calculator, environment_energies, emt):
    """Run independent TOP basin-hopping restarts. Returns one
    (emt_energy, step, best_particle) record per restart."""
    records = []
    for _ in range(cfg["n_runs"]):
        start = new_particle(cfg["height"], cfg["trunc"])
        best, lowest_energies, _ = run_basin_hopping(
            start, calculator, environment_energies,
            cfg["bh_hopping_attempts"], cfg["bh_hops"], model="TOP")
        _, step = lowest_energies[-2]
        records.append((emt_score(best, emt), step, best))
    return records


def cumulative_success(records, global_best, tol, n_runs, x_max):
    """Staircase of the cumulative percentage of restarts that reached the global
    minimum (within `tol`), indexed by the step at which each such restart found
    its best ordering. Plotted with `step(where="post")` and a log x-axis."""
    success_steps = sorted(max(step, 1) for energy, step, _ in records
                           if energy <= global_best + tol)
    x = [1] + success_steps + [x_max]
    y = [0.0] + [(i + 1) * 100.0 / n_runs for i in range(len(success_steps))]
    y += [len(success_steps) * 100.0 / n_runs]
    return x, y


def split_hemispheres(atoms, gap=3.0):
    """Pull the two z-halves of the particle apart to expose the core ordering."""
    atoms = atoms.copy()
    atoms.center()
    z_mid = atoms.get_cell()[2][2] / 2 + 1.0
    lower = atoms[[a.index for a in atoms if a.position[2] < z_mid]]
    upper = atoms[[a.index for a in atoms if a.position[2] >= z_mid]]
    lower.translate((0.0, 0.0, -gap))
    upper.translate((0.0, 0.0, gap))
    return lower + upper


def make_figure(gm_particle, bh_curve, mc_curve, n_runs, x_max, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))

    plot_atoms(split_hemispheres(gm_particle.get_ase_atoms()),
               ax=axes[0], rotation="0x,75y,0z")
    axes[0].set_title("(a) Optimized Pt/Au ordering")
    axes[0].axis("off")

    ax = axes[1]
    ax.step(*bh_curve, where="post", label="BH")
    ax.step(*mc_curve, where="post", label="MC")
    ax.set_xscale("log")
    ax.set_xlim(1, x_max)
    ax.set_ylim(0, 105)
    ax.set_title("(b) Cumulative success rate")
    ax.set_xlabel("Steps (log scale)")
    ax.set_ylabel("Runs reaching global minimum (%)")
    ax.grid(True, which="both")
    ax.legend()

    fig.suptitle("Pt/Au chemical-ordering optimization over {} restarts".format(n_runs))
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true",
                        help="fast smoke-test config (small particle, few steps)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tol", type=float, default=1e-2,
                        help="EMT energy tolerance (eV) for counting a restart as having "
                             "reached the global minimum")
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__),
                                                         "ptau_go_results"))
    args = parser.parse_args()

    cfg = QUICK if args.quick else FULL
    np.random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    emt = EMTCalculator(relax_atoms=False)
    etop = ExtendedTopologicalFeaturesClassifier(SYMBOLS)
    top = testTopologicalFeatureClassifier(SYMBOLS)

    print("Building training set ({} particles, height={}, trunc={}) ...".format(
        cfg["n_train"], cfg["height"], cfg["trunc"]))
    training_set = build_training_set(cfg["n_train"], cfg["height"], cfg["trunc"], emt, etop, top)
    n_atoms = training_set[0].get_n_atoms()
    print("Particle size: {} atoms. Fitting ETOP (MC) and TOP (BH) models on EMT data ...".format(
        n_atoms))

    mc_calculator = fit_mc_model(training_set, etop)
    bh_calculator, environment_energies = fit_bh_model(training_set, top, n_atoms)

    print("Optimizing chemical ordering: Monte Carlo ({} restarts) ...".format(cfg["n_runs"]))
    mc_records = optimize_with_mc(cfg, mc_calculator, etop, emt)
    print("Optimizing chemical ordering: TOP basin hopping ({} restarts) ...".format(
        cfg["n_runs"]))
    bh_records = optimize_with_bh(cfg, bh_calculator, environment_energies, emt)

    all_records = mc_records + bh_records
    global_best = min(energy for energy, _, _ in all_records)
    gm_particle = min(all_records, key=lambda r: r[0])[2]
    be_mc = min(energy for energy, _, _ in mc_records)
    be_bh = min(energy for energy, _, _ in bh_records)
    print("Best EMT energy  MC: {:.3f}  BH: {:.3f}  -> global minimum {:.3f} eV".format(
        be_mc, be_bh, global_best))

    x_max = max(max(step for _, step, _ in all_records), 1)
    mc_curve = cumulative_success(mc_records, global_best, args.tol, cfg["n_runs"], x_max)
    bh_curve = cumulative_success(bh_records, global_best, args.tol, cfg["n_runs"], x_max)

    fig_path = os.path.join(args.outdir, "ptau_go_overview.png")
    make_figure(gm_particle, bh_curve, mc_curve, cfg["n_runs"], x_max, fig_path)

    data_path = os.path.join(args.outdir, "ptau_go_data.npz")
    np.savez(data_path,
             emt_energies_mc=[r[0] for r in mc_records], steps_mc=[r[1] for r in mc_records],
             emt_energies_bh=[r[0] for r in bh_records], steps_bh=[r[1] for r in bh_records],
             global_best=global_best, n_atoms=n_atoms, **cfg)
    gm_path = os.path.join(args.outdir, "global_minimum.xyz")
    ase_write(gm_path, gm_particle.get_ase_atoms())

    print("\nWrote:")
    print("  figure : {}".format(fig_path))
    print("  data   : {}".format(data_path))
    print("  GM xyz : {}".format(gm_path))


if __name__ == "__main__":
    main()
