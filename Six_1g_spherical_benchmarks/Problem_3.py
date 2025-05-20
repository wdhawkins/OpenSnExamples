#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Problem 3
Chang, B. (2021a, February 19).
Six 1D polar transport test problems for the deterministic and monte-carlo method.
LLNL-TR-819680
https://www.osti.gov/servlets/purl/1769096
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
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased, GetFieldFunctions
    from pyopensn.logvol import SphereLogicalVolume, BooleanLogicalVolume

if __name__ == "__main__":

    meshgen = FromFileMeshGenerator(
        filename="./meshes/three_region_sphere.msh",
        partitioner=PETScGraphPartitioner(type='parmetis'),
    )
    grid = meshgen.Execute()

    # create XS
    num_groups = 1
    xs_void = MultiGroupXS()
    xs_void.CreateSimpleOneGroup(sigma_t=0.0,c=0.0)
    # xs_mat = MultiGroupXS()
    # xs_mat.CreateSimpleOneGroup(sigma_t=1.0,c=0.0)

    # volumetric src
    mg_src = VolumetricSource(block_ids=[2], group_strength=[1.])

    # angular quadrature
    nazimu = 4
    npolar = 2
    pquad = GLCProductQuadrature3DXYZ(npolar, nazimu)

    # Solver
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "l_max_its": 500,
                "l_abs_tol": 1.0e-6,
                "gmres_restart_interval": 30
            },
        ],
        xs_map=[
            {"block_ids": [1,2,3], "xs": xs_void},
        ],
        options={
            "scattering_order": 0,
            "volumetric_sources": [mg_src],
            "boundary_conditions": [
                {"name":"xmin", "type":"vacuum"},
                {"name":"xmax", "type":"vacuum"},
                {"name":"ymin", "type":"vacuum"},
                {"name":"ymax", "type":"vacuum"},
                {"name":"zmin", "type":"vacuum"},
                {"name":"zmax", "type":"vacuum"}
            ]
        },

    )

    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # export
    fflist = phys.GetFieldFunctions()
    """
        vtk_basename = "Problem_3_"
        # export only the scalar flux of group g
        FieldFunctionGridBased.ExportMultipleToVTK(
            [fflist[g] for g in range(num_groups)],
            vtk_basename
        )
    """
    # compute the average flux over a logical volume
    def average_flx_logvol(logvol):
        ffvol = FieldFunctionInterpolationVolume()
        ffvol.SetOperationType("avg")
        ffvol.SetLogicalVolume(logvol)
        ffvol.AddFieldFunction(fflist[0])
        ffvol.Initialize()
        ffvol.Execute()
        avgval = ffvol.GetValue()
        return avgval

    # compute the average flux over spherical shells
    def compute_average_flux(N_logvols, r_max):
        r_vals = np.linspace(0, r_max, N_logvols + 1)
        avg_flx = np.zeros(N_logvols)
        for i in range(N_logvols):
            if i != 0:
                inner_vol = SphereLogicalVolume(r=r_vals[i])
                outer_vol = SphereLogicalVolume(r=r_vals[i + 1])
                logvol = BooleanLogicalVolume(parts=[{"op":True,"lv":outer_vol},{"op":False,"lv":inner_vol}])
            else:
                logvol = SphereLogicalVolume(r=r_vals[i + 1])
            avg_flx[i] = average_flx_logvol(logvol)
        return r_vals, avg_flx

    # perform the calculation
    N_logvols = 15
    r_vals, avg_flx = compute_average_flux(N_logvols, 1)
    if rank == 0:
        for r,v in zip(r_vals[1:],avg_flx):
            print("end-radius: {:.2f}, avg-value: {:.6f}".format(r, v))

    # compute the analytical answer
    def get_phi(r, q, a, b):
        a2 = a**2
        b2 = b**2
        r2 = r**2
        phi = q*(b+(b2-r2)/(2*r)*np.log((b+r)/(np.abs(r-b))))-q*(a+(a2-r2)/(2*r)*np.log((a+r)/(np.abs(r-a))))
        return phi

    if rank==0:

        q = 0.5
        a = 1./3.
        b = 2./3.
        eps = 1e-3
        r_fine = np.linspace(0.+eps, 1.0-eps, 100)
        exact_phi = np.zeros_like(r_fine)
        for i,r in enumerate(r_fine):
            exact_phi[i] = get_phi(r, q, a, b)
            print(r_fine[i], exact_phi[i])


        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(r_fine, exact_phi, label='Exact')
        plt.step(r_vals, np.append(avg_flx, avg_flx[-1]), where='post', label='Radial Average Flux in Spherical Shells')
        plt.xlabel("Radius")
        plt.ylabel("Average Flux")
        plt.title("Flux as a function of radius")
        plt.grid(True)
        plt.legend()
        # Save the figure as a PNG file
        plt.savefig("Problem_3.png", dpi=300, bbox_inches='tight')
        plt.show()