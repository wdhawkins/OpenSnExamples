#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HEU_MET_FAST_003 benchmark
"""

import sys
import os
import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    # Append parent directory to locate the pyopensn modules
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, PETScGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased

if __name__ == "__main__":

    meshgen = FromFileMeshGenerator(
        # filename="./HEU_MET_FAST_003.msh",
        filename="./HEU_MET_FAST_003_coarse.msh",
        partitioner=PETScGraphPartitioner(type='parmetis'),
    )
    grid = meshgen.Execute()
    # grid.ExportToPVTU("HEU_MET_FAST_003_mesh_only")

    # get "measured" volumes
    # Block ID: 0 Volume: 1286.10113382
    # Block ID: 1 Volume: 5636.96006452
    volumes_per_block = grid.ComputeVolumePerBlockID()

    # define radii per block-ID
    radii = {1: 6.7820, 2: 11.8620}

    # sort block IDs by increasing radius
    block_ids = sorted(volumes_per_block.keys())

    # check:
    missing = [blk for blk in block_ids if blk not in radii]
    if missing:
        raise KeyError(f"No radius provided for block IDs: {missing}")

    # compute exact volumes
    exact_volumes_per_block = {}
    prev_R = 0.0
    for blk in block_ids:
        R = radii[blk]
        # shell from prev_R to R
        exact_volumes_per_block[blk] = (4.0 / 3.0) * np.pi * (R**3 - prev_R**3)
        prev_R = R

    # build the ratios array
    ratios = np.array([
        exact_volumes_per_block[blk] / volumes_per_block[blk]
        for blk in block_ids
    ])

    # only rank 0 prints
    if rank == 0:
        for idx, blk in enumerate(block_ids):
            print(f"Block {blk}: measured = {volumes_per_block[blk]:.11e}, "
                  f" exact = {exact_volumes_per_block[blk]:.11e}, "
                  f" ratio = {ratios[idx]:.6f}")

    # load and modify XS
    xs_oralloy = MultiGroupXS()
    xs_oralloy.LoadFromOpenMC("Oralloy_LANL30g.h5", "Oralloy", 294.0)
    xs_oralloy.SetScalingFactor(ratios[0])

    xs_tuballoy = MultiGroupXS()
    xs_tuballoy.LoadFromOpenMC("Tuballoy_LANL30g.h5", "Tuballoy", 294.0)
    xs_tuballoy.SetScalingFactor(ratios[1])

    num_groups = xs_tuballoy.num_groups
    if rank == 0:
        print("num groups =", num_groups)

    # Solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=8,
                    n_azimuthal=16,
                    scattering_order=3
                ),
                "inner_linear_method": "petsc_gmres",
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "l_max_its": 500,
                "l_abs_tol": 1.0e-6,
            },
        ],
        xs_map=[
            {"block_ids": [1], "xs": xs_oralloy},
            {"block_ids": [2], "xs": xs_tuballoy},
        ],
        options={
            "scattering_order": 3,
            "use_precursors": False,
            "verbose_inner_iterations": True,
            "verbose_outer_iterations": True,
        },
    )

    # k_solver = PowerIterationKEigenSolver(
    #     lbs_problem=phys,
    #     k_tol=1e-6,
    # )
    k_solver = NonLinearKEigenSolver(
        lbs_problem=phys,
        nl_max_its=500,
        nl_abs_tol=1.0e-8,
    )
    k_solver.Initialize()
    k_solver.Execute()
    k = k_solver.GetEigenvalue()
    # only rank 0 prints
    if rank == 0:
        print(f"Computed k-eigenvalue: {k}")

    # export
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)
    vtk_basename = "HEU_MET_FAST_003"
    # export only the flux of group g (first []), moment 0 (second [])
    FieldFunctionGridBased.ExportMultipleToVTK(
        [fflist[g][0] for g in range(num_groups)],
        vtk_basename
    )
