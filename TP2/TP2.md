# TP1: Solve laplacian problem with Hybrid High-Order methods

## Problem to solve

The main goal is to solve the following weak fork of the Poisson problem:
Find $u \in H^1_0(\Omega)$ such that
$$ (\nabla u, \nabla v)_\Omega = (f,v)_\Omega, \: \forall  v \in H^1_0(\Omega) $$
with $f \in L^2(\Omega)$. This problem is well-posed.

## HHO discretization

The continuous weak form is discretized with HHO methods. Let $\mathcal{T}_h$ be a mesh of $\Omega$, with $T$ a generic cell and $\partial T$ its boundary. For $k \geq 0$, we introduce the local discrete space $U^k_T$ such that

$$ U^k_T := P^k_d(T) \times P^k_{d-1}(\partial T) $$
where $P^k_d(T)$ is a polynomial space of degree at most $k$.

The local HHO unknowns is $\hat{u}_T := (u_T, u_{\partial T}) \in U^k_T$. The discrete weak form to solve: Find $\hat{u}_h \in U^k_{h,0}$ such that

$$
 \sum_{T  \in  \mathcal{T}_h}(G_T(\hat{u}_T), G_T(\hat{v}_T))_T + \beta (S_{\partial T}(\hat{u}_T), S_{\partial T}(\hat{v}_T))_{\partial T} = \sum_{T  \in  \mathcal{T}_h} (f, v_T)_T, \: \forall  \hat{v}_T \in U^k_{h,0}
$$
where $G_T$ is the HHO reconstructed gradient, $S_T$ the stabilization operator and $\beta$ the stabilization coefficient. This discrete problem is well-posed for $\beta > 0$.

For more theoretical details, see the following book "Hybrid high-order methods. A primer with application to solid mechanics" (<https://arxiv.org/abs/2106.09348>)

## A first simple example

To launch a test, you need at least two files:

- the .export file with runtime parameters and external data
- the .comm or .py file: main file with code_aster or python syntax

There is a basic example below:

- tp2.export

```none
P time_limit 30
P memory_limit 512
P ncpus 1
P mpi_nbcpu 1
P mpi_nbnoeud 1

F comm tp2.py D 1
F mess mess.txt R 6
```

- tp2.py

```python
import code_aster
from code_aster.Commands import *

from math import pi, sqrt, log

code_aster.init()

assert True

code_aster.close()
```

To run your test inside the container (we suppose that all files are in current working folder)

```bash
/opt/public/code_aster/bin/run_aster tp2.export
```

The run is a success if you have the following output

```none
------------------------------------------------------------
--- DIAGNOSTIC JOB : OK
------------------------------------------------------------
```

Moreover, there a complete output file named mess.txt. The online documentation is here <https://code-aster.org/V2/spip.php?rubrique19>

## A test with a known solution

The domain is the unit square $\Omega=[0,1]^2$ with analytical solution $u(x,y)=sin(\pi x)*sin(\pi y)$ and source load $f(x,y)=2\pi^2sin(\pi x)*sin(\pi y)$. Hence, we have homogeneous Dirichlet boundary conditions.

In a code_aster command file, we can mix python and DSL (Domain Specific Language). Firstly, we define a mesh.

### Mesh definition

```python
# initial_mesh - domain [0,1]^2 - 4 quads
mesh0_quad = code_aster.Mesh.buildSquare(refine=1)

# split in triangular mesh
mesh0_tri = CREA_MAILLAGE(MAILLAGE=mesh0_quad, MODI_MAILLE=_F(TOUT="OUI", OPTION="QUAD_TRIA3"))

# convert for hho-cells
mesh0_hho = CREA_MAILLAGE(MAILLAGE=mesh0_tri, MODI_HHO=_F(TOUT="OUI"))
```

Because of the architecture of code_aster, a special treatment is needed on the mesh to support HHO finite element.

### Finite element definition

So, we define the finite element space chosen

```python
# define finite element model
model = AFFE_MODELE(
            MAILLAGE=mesh0_hho,
            AFFE=_F(TOUT="OUI", MODELISATION="PLAN_HHO", FORMULATION="LINEAIRE", PHENOMENE="THERMIQUE"),
        )
```

There is two choices for HHO formulation in code_aster: "LINEAIRE" (k=1) and "QUADRATIQUE" (k=2).

### Boundary conditions and load

We use homogeneous boundary conditions

```python
# define Dirichlet BC
bc = AFFE_CHAR_CINE(
            MODELE=model, THER_IMPO=_F(GROUP_MA=("RIGHT", "LEFT", "TOP", "BOTTOM"), TEMP=0.0)
        )
```

and the source term

```python
# define load function
f_load = FORMULE(VALE="2.0*pi*pi*sin(pi*X)*sin(pi*Y)", NOM_PARA=("X", "Y"))

# define external load
load = AFFE_CHAR_THER_F(MODELE=model, SOURCE=_F(GROUP_MA=("SURFACE"), SOUR=f_load))
```

### Define material

```python
# define material
coeff = DEFI_MATERIAU(THER=_F(LAMBDA=1, RHO_CP=1), HHO=_F(COEF_STAB=1))

# apply material on mesh
mater = AFFE_MATERIAU(MODELE=model, AFFE=_F(TOUT="OUI", MATER=coeff))
```

### Elementary computations

Firstly, we define the physical problem with model, material and loads then we compute DoFs numbering

```python
# define physical problem
phys_pb = code_aster.PhysicalProblem(model, mater)
phys_pb.addLoad(load)
phys_pb.addDirichletBC(bc)

# compute DOF numbering
phys_pb.computeDOFNumbering()
```

Now, we define the DiscreteComputation object. Then, we compute the rigidity matrix (with stabilization term included) and external load

```python
# create discrete computation
disc_comp = code_aster.DiscreteComputation(phys_pb)

# compute rigidity matrix: (lambda * GkT(huT), GkT(hvT))_T + lambda * stab(huT, hvT)
rigidity = disc_comp.getLinearStiffnessMatrix(assembly=True)

# compute external load: (f, vT)_T
rhs = disc_comp.getNeumannForces()
```

### Resolution

Finally, the linear system is solved

```python
# compute Dirichlet BC to apply
diriBCs = disc_comp.getDirichletBC()

# define linear solver - MUMPS
mySolver = code_aster.MumpsSolver()

# factorize and solve
mySolver.factorize(rigidity)
u_hho = mySolver.solve(rhs, diriBCs)
```

### Compute $L^2$ and $H^1$-norm errors

Now, we can compute the $L^2$ and $H^1$-norm errors compare to the analytical solution

```python
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

print("L2-norm: %f and H1-norm %f"%(l2_error, h1_error))
```

### Post-processing

Since HHO solution is not directly viewable by Paraview, we project the HHO solution on Lagrange finite element space

```python
# Project analytical solution on HHO space
u_proj = hho.projectOnLagrangeSpace(u_hho)

# Save solution in MED-format
u_proj.printMedFile("/tmp/u_proj.med")
```

The med file can be opened with Paraview inside salome_meca.

### Other possibilities

You can change some parameters and see the differences

- refine the mesh
- change the stabilization parameter "COEF_STAB"
- change approximation order from "LINEAIRE" to "QUADRATIQUE"

The correction is given in tp2_1.py

## Convergence order

The goal of this section is to perform a mesh refinement in order to compute converge order. We recall that for HHO, we have a convergence in $h^{k+1}$ in $H^1$-norm and $h^{k+2}$ in $L^2$-norm for a solution smooth enougth.

For that, we define the number of refinement that we want with one python on mesh refinement and one on approximation order. The only modification is that a different mesh is used at every step

```python
# number of refinement
nb_reff = 6

# initial_mesh - domain [0,1]^2 - 4 quads
mesh0_quad = code_aster.Mesh.buildSquare(refine=1)
# split in triangular mesh
mesh0_tri = CREA_MAILLAGE(
    MAILLAGE=mesh0_quad, MODI_MAILLE=_F(TOUT="OUI", OPTION="QUAD_TRIA3")
)
# convert for hho-cells
mesh0_hho = CREA_MAILLAGE(MAILLAGE=mesh0_tri, MODI_HHO=_F(TOUT="OUI"))

# size of a triangle on the original mesh
h0 = sqrt(2) / 2

for order in ("LINEAIRE", "QUADRATIQUE"):
    mesh = mesh0_hho
    h = h0
    for i_reff in range(nb_reff):
        ## DEFINE PROBLEM
        # create mesh - refine previous mesh
        mesh = mesh.refine(1)
        h = h / 2

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
            MODELE=model, THER_IMPO=_F(GROUP_MA=("RIGHT", "LEFT", "TOP", "BOTTOM"), TEMP=0.0))

        # rest of the resolution is identical
        # ...
```

You can saved error like that after declaration of variable:

```python
# compute L2 and H1-errors
error[order]["h"].append(sqrt(2) / ((i_reff + 1) ** 2))
error[order]["L2"] = l2_error
error[order]["H1"] = h1_error
```

### Linear regression

The asymptotic convergence order are computed using a linear regression in $\log-\log$ scale.

```python
import numpy as np

# compute convergence order
conv_order = {}
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
```

Have you the correct convergence order ?

### Visualization

We use matplotlib to plot the different curves with the following lines

```python
import os
import matplotlib as plt
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
            the_conv["QUADRATIQUE"][norm]["c"] * h ** the_conv["QUADRATIQUE"][norm]["o"]
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
```

The correction is given in tp2_2.py. As before, you can change the stabilization coefficient.
