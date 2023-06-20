# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

import code_aster
from code_aster.Commands import *

import os
import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True

except ImportError:
    HAS_MATPLOTLIB = False

from math import pi, sqrt, log

code_aster.init("--test")

test = code_aster.TestCase()

###################################################################################
#
#   Analytical solution
#   -laplacian(u) = f
#   u = sin(Pi*x)*sin(Pi*y)
#   f = 2.0*Pi*Pi*sin(Pi*x)*sin(Pi*y)
#
#   Weak form: (grad u, grad v) = (f,v)
#   HHO unknowns : huT = (uT, udT)
#
####################################################################################

# number of refinement
nb_reff = 6

# define analytical solution
u_ana = FORMULE(VALE="sin(pi*X)*sin(pi*Y)", NOM_PARA=("X", "Y"))

# define load function
f_load = FORMULE(VALE="2.0*pi*pi*sin(pi*X)*sin(pi*Y)", NOM_PARA=("X", "Y"))

# error save
error = {}
conv_order = {}

# initial_mesh - domain [0,1]^2 - 4 quads
mesh0_quad = code_aster.Mesh.buildSquare(refine=1)
# split in triangular mesh
mesh0_tri = CREA_MAILLAGE(
    MAILLAGE=mesh0_quad, MODI_MAILLE=_F(TOUT="OUI", OPTION="QUAD_TRIA3")
)
# convert for hho-cells
mesh0_hho = CREA_MAILLAGE(MAILLAGE=mesh0_tri, MODI_HHO=_F(TOUT="OUI"))

# size of a triangle
h0 = sqrt(2) / 2

for order in ("LINEAIRE", "QUADRATIQUE"):
    error[order] = {"h": [], "L2": [], "H1": []}
    mesh = mesh0_hho
    h = h0
    for i_reff in range(nb_reff):
        ## DEFINE PROBLEM
        # create mesh - refine previous mesh
        mesh = mesh.refine(1)
        h = h / 2

        # size of a cell
        error[order]["h"].append(sqrt(2) / ((i_reff + 1) ** 2))

        # define material
        coeff = DEFI_MATERIAU(THER=_F(LAMBDA=1, RHO_CP=1), HHO=_F(COEF_STAB=1))

        # apply material on mesh
        mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

        # define finite element model
        model = AFFE_MODELE(
            MAILLAGE=mesh,
            AFFE=_F(
                TOUT="OUI",
                MODELISATION="PLAN_HHO",
                FORMULATION=order,
                PHENOMENE="THERMIQUE",
            ),
        )

        # define Dirichlet BC
        bc = AFFE_CHAR_CINE(
            MODELE=model,
            THER_IMPO=_F(GROUP_MA=("RIGHT", "LEFT", "TOP", "BOTTOM"), TEMP=0.0),
        )

        # define external load
        load = AFFE_CHAR_THER_F(
            MODELE=model, SOURCE=_F(GROUP_MA=("SURFACE"), SOUR=f_load)
        )

        ## COMPUTE DISCRETE SOLUTION
        # define physical problem
        phys_pb = code_aster.PhysicalProblem(model, mater)
        phys_pb.addLoad(load)
        phys_pb.addDirichletBC(bc)

        # compute DOF numbering
        phys_pb.computeDOFNumbering()

        # create discrete computation
        disc_comp = code_aster.DiscreteComputation(phys_pb)

        # compute rigidity matrix: (lambda * GkT(huT), GkT(hvT))_T + lambda * stab(huT, hvT)
        rigidity = disc_comp.getLinearStiffnessMatrix(assembly=True)

        # compute external load: (f, vT)_T
        rhs = disc_comp.getNeumannForces()

        # compute Dirichlet BC to apply
        diriBCs = disc_comp.getDirichletBC()

        # define linear solver - MUMPS
        mySolver = code_aster.MumpsSolver()

        # factorize and solve
        mySolver.factorize(rigidity)
        u_hho = mySolver.solve(rhs, diriBCs)

        ## COMPUTE ERROR
        # create hho handler
        hho = code_aster.HHO(phys_pb)

        # Project analytical solution on HHO space
        u_proj = hho.projectOnHHOSpace(u_ana)

        # compute difference
        u_diff = u_hho - u_proj

        # Compute mass matrix M = (vT, wT) - RHO_CP == 1 in DEFI_MATERIAU
        mass = disc_comp.getMassMatrix(assembly=True)

        # compute L2 and H1-errors
        error[order]["L2"].append(sqrt((mass * u_diff).dot(u_diff)))
        error[order]["H1"].append(sqrt((rigidity * u_diff).dot(u_diff)))

# compute convergence order
for order in ("LINEAIRE", "QUADRATIQUE"):
    conv_order[order] = {}
    for norm in ("L2", "H1"):
        xlog = [log(h) for h in error[order]["h"]]
        ylog = [log(y) for y in error[order][norm]]
        A = np.vstack([xlog, np.ones(len(xlog))]).T
        m, c = np.linalg.lstsq(A, ylog, rcond=None)[0]
        conv_order[order][norm] = [m, c]

    print(
        "Convergence order %s -> L2-norm: %f and H1-norm: %f"
        % (order, conv_order[order]["L2"][0], conv_order[order]["H1"][0])
    )

    # test convergence order
    test.assertAlmostEqual(
        conv_order[order]["L2"][0],
        {"LINEAIRE": 2.839235, "QUADRATIQUE": 3.784263}[order],
        delta=1e-4,
    )
    test.assertAlmostEqual(
        conv_order[order]["H1"][0],
        {"LINEAIRE": 1.889858, "QUADRATIQUE": 2.835326}[order],
        delta=1e-4,
    )

# plot figure with matplotlib
if HAS_MATPLOTLIB and os.getenv("DISPLAY"):
    import aster_core

    # disable floating point exceptions from matplotlib
    aster_core.matfpe(-1)

    ylim = {"L2": [1e-10, 0.1], "H1": [1e-7, 1.0]}
    the_conv = {
        "LINEAIRE": {"L2": {"o": 3, "c": 2e-2}, "H1": {"o": 2, "c": 0.2}},
        "QUADRATIQUE": {"L2": {"o": 4, "c": 1e-3}, "H1": {"o": 3, "c": 0.02}},
    }
    # plot the data
    for norm in ("L2", "H1"):
        plt.plot(1, 1, 1)

        plt.xscale("log")
        plt.yscale("log")

        # set the limits
        plt.xlim([0.02, 2])
        plt.ylim(ylim[norm])
        plt.xlabel("mesh-size")
        plt.ylabel("%s-error" % norm)

        plt.plot(
            error["LINEAIRE"]["h"],
            error["LINEAIRE"][norm],
            marker="o",
            color="tab:blue",
            label="k=1, computed",
        )
        m, c = conv_order["LINEAIRE"][norm]
        plt.plot(
            error["LINEAIRE"]["h"],
            [
                the_conv["LINEAIRE"][norm]["c"] * h ** the_conv["LINEAIRE"][norm]["o"]
                for h in error["LINEAIRE"]["h"]
            ],
            "k--",
            label="k=1, theorical",
            color="tab:blue",
        )

        plt.plot(
            error["QUADRATIQUE"]["h"],
            error["QUADRATIQUE"][norm],
            marker="o",
            color="tab:orange",
            label="k=2, computed",
        )
        m, c = conv_order["QUADRATIQUE"][norm]
        plt.plot(
            error["QUADRATIQUE"]["h"],
            [
                the_conv["QUADRATIQUE"][norm]["c"]
                * h ** the_conv["QUADRATIQUE"][norm]["o"]
                for h in error["LINEAIRE"]["h"]
            ],
            "k--",
            label="k=2, theorical",
            color="tab:orange",
        )

        plt.legend()
        plt.title("%s convergence error for HHO" % norm)

        # save plot
        savedir = "/tmp/" or os.getcwd()
        plt.savefig(os.path.join(savedir, "%s_error.png" % norm))
        plt.clf()

# close
code_aster.close()
