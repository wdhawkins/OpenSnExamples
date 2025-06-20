#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenSn text in block example problem
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    # Append parent directory to locate the pyopensn modules
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import MeshGenerator, FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS, CreateSimpleOneGroup
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased

if __name__ == "__main__":

    # Mesh
    meshgen = FromFileMeshGenerator(filename="opensn.msh")
    grid = meshgen.Execute()

    # Cross sections
    Ng = 1
    xs_source = MultiGroupXS()
    xs_source.CreateSimpleOneGroup(0.1, 0.1)
    xs_block = MultiGroupXS()
    xs_block.CreateSimpleOneGroup(1.0e-6, 0.001)

    # Source
    src_strength = [0.0 for _ in range(Ng)]
    src_strength[0] = 100.0
    src = VolumetricSource(block_ids=[1], group_strength=src_strength)

    # Quadrature
    pquad = GLCProductQuadrature3DXYZ(n_polar=8, n_azimuthal=16, scattering_order=0)

    # Solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=Ng,
        groupsets=[
            {
                "groups_from_to": (0, Ng - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
            }
        ],
        xs_map=[
            {"block_ids": [1], "xs": xs_source},
            {"block_ids": [2], "xs": xs_block},
        ],
        options={
            "scattering_order": 0,
            "volumetric_sources": [src],
        }
    )
    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Export results to VTK
    fflist = phys.GetScalarFieldFunctionList()
    FieldFunctionGridBased.ExportMultipleToVTK(fflist, "flux")
