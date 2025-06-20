#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Urban source example problemS
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
    meshgen = FromFileMeshGenerator(filename="urban_source.msh")
    grid = meshgen.Execute()

    # Cross sections
    Ng = 1
    xs_concrete = MultiGroupXS()
    xs_concrete.CreateSimpleOneGroup(0.05, 0.98)
    xs_source = MultiGroupXS()
    xs_source.CreateSimpleOneGroup(1.0, 0.5)
    xs_air = MultiGroupXS()
    xs_air.CreateSimpleOneGroup(0.0001, 0.999)

    # Source
    src_strength = [0.0 for _ in range(Ng)]
    src_strength[0] = 1000000.0
    src = VolumetricSource(block_ids=[2], group_strength=src_strength)

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
                "l_max_its": 500,
            },
        ],
        xs_map=[
            {"block_ids": [1], "xs": xs_concrete},
            {"block_ids": [2], "xs": xs_source},
            {"block_ids": [3], "xs": xs_air},
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
