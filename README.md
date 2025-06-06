# Relax-Fix-and-Exclude (RFE) for MINLP problems with Multilinear Interpolations

This repository contains the code necessary to reproduce the results reported in
> Pacheco, B. M., Antunes, P. M., Camponogara, E., Seman, L. O., Rosa, V. R., Vieira, B. F., & Longhi, C. (2025). A relax-fix-and-exclude algorithm for an MINLP problem with multilinear interpolations. arXiv. https://doi.org/10.48550/ARXIV.2502.21249

## Requirements

All package requirements are contained in the Julia environment in `Project.toml`. Note that we use our custom package [Oil.jl](github.com/gos-ufsc/Oil.jl/) to manipulate the oil production components. Because Oil.jl is not a listed package, you have to instantiate the environment manually.

Naturally, licenses to Gurobi and BARON are required to use the solvers.

## Problem formulations

The construction of the MINLP problem and the MILP relaxation is implemented in `models.jl`. Note that we use Gurobi as our default solver for those, such that a license is required. Nevertheless, it shoudl be easy to change from Gurobi to any combination of NLP solver and MILP solver both in `rfe.jl` and in `batch_rfe.jl`.

## Data

Our instances rely on a VLP curve... **TODO**

## Experiments

The instances used in our paper are available under `scenarios/`. A direct reproduction of our experiments can be achieved by running the `batch_*.jl` files. To run your own experiments with RFE, refer to `rfe.jl`, in which we have a detailed implementation of the algorithm for a sample problem.
