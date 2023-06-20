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

from math import pi, sqrt, log

code_aster.init()

# initial_mesh - domain [0,1]^2 - 4 quads
mesh0_quad = code_aster.Mesh.buildSquare(refine=1)
# split in triangular mesh
mesh0_tri = CREA_MAILLAGE(
    MAILLAGE=mesh0_quad, MODI_MAILLE=_F(TOUT="OUI", OPTION="QUAD_TRIA3")
)
# convert for hho-cells
mesh0_hho = CREA_MAILLAGE(MAILLAGE=mesh0_tri, MODI_HHO=_F(TOUT="OUI"))

# define finite element model
model = AFFE_MODELE(
    MAILLAGE=mesh0_hho,
    AFFE=_F(
        TOUT="OUI",
        MODELISATION="PLAN_HHO",
        FORMULATION="LINEAIRE",
        PHENOMENE="THERMIQUE",
    ),
)

# define Dirichlet BC
bc = AFFE_CHAR_CINE(
    MODELE=model, THER_IMPO=_F(GROUP_MA=("RIGHT", "LEFT", "TOP", "BOTTOM"), TEMP=0.0)
)

# define load function
f_load = FORMULE(VALE="2.0*pi*pi*sin(pi*X)*sin(pi*Y)", NOM_PARA=("X", "Y"))

# define external load
load = AFFE_CHAR_THER_F(MODELE=model, SOURCE=_F(GROUP_MA=("SURFACE"), SOUR=f_load))

# define material
coeff = DEFI_MATERIAU(THER=_F(LAMBDA=1, RHO_CP=1), HHO=_F(COEF_STAB=1))

# apply material on mesh
mater = AFFE_MATERIAU(MODELE=model, AFFE=_F(TOUT="OUI", MATER=coeff))

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

# create hho handler
hho = code_aster.HHO(phys_pb)

# define analytical solution
u_ana = FORMULE(VALE="sin(pi*X)*sin(pi*Y)", NOM_PARA=("X", "Y"))

# Project analytical solution on HHO space
u_proj = hho.projectOnHHOSpace(u_ana)

# compute difference
u_diff = u_hho - u_proj

# Compute mass matrix M = (vT, wT) - RHO_CP == 1 in DEFI_MATERIAU
mass = disc_comp.getMassMatrix(assembly=True)

# compute L2 and H1-errors
l2_error = sqrt((mass * u_diff).dot(u_diff))
h1_error = sqrt((rigidity * u_diff).dot(u_diff))

print("L2-norm: %f and H1-norm %f" % (l2_error, h1_error))

# Project analytical solution on HHO space
u_proj = hho.projectOnLagrangeSpace(u_hho)

# Save solution in MED-format
u_proj.printMedFile("/tmp/u_proj.med")

code_aster.close()
