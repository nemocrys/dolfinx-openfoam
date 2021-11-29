"""
Solid plate participant in flow-over-plate tutorial using FEniCS
"""

# import ufl
# import dolfinx
# import dolfinx.io
# from mpi4py import MPI
# from petsc4py import PETSc


from mpi4py import MPI
from dolfinx import Function, FunctionSpace, VectorFunctionSpace, Constant, DirichletBC
from dolfinx.fem import LinearProblem, locate_dofs_geometrical
from dolfinx.generation import RectangleMesh
# from dolfinx.mesh import CellType  # this should work according to https://jorgensd.github.io/dolfinx-tutorial/chapter2/diffusion_code.html?highlight=rectangular
from dolfinx.io import XDMFFile

from ufl import FiniteElement, VectorElement, TrialFunction, TestFunction, grad, inner, dot, dx
# from dolfinx import interpolate, near, MeshFunction, rhs, lhs

# from fenicsprecice import Adapter
import numpy as np


# class ComplementaryBoundary(SubDomain):
#     """Determines if a point is at the complementary boundary with tolerance of
#     1E-14.
#     :func inside(): returns True if point belongs to the boundary, otherwise
#                     returns False
#     """

#     def __init__(self, subdomain):
#         self.complement = subdomain
#         SubDomain.__init__(self)

#     def inside(self, x, on_boundary):
#         tol = 1E-14
#         if on_boundary and not self.complement.inside(x, on_boundary):
#             return True
#         else:
#             return False

# class TopBoundary(SubDomain):
#     """Determines if the point is at the top boundary with tolerance of 1E-14.
#     :func inside(): returns True if point belongs to the boundary, otherwise
#                     returns False
#     """

#     def inside(self, x, on_boundary):
#         tol = 1E-14
#         if on_boundary and near(x[1], y_top, tol):
#             return True
#         else:
#             return False
def top_boundary(x):
    tol = 1E-14
    return np.isclose(x[1], y_top, tol)


# class BottomBoundary(SubDomain):
#     """Determines if the point is at the bottom boundary with tolerance of
#     1E-14.

#     :func inside(): returns True if point belongs to the boundary, otherwise
#                     returns False
#     """

#     def inside(self, x, on_boundary):
#         tol = 1E-14
#         if on_boundary and near(x[1], y_bottom, tol):
#             return True
#         else:
#             return False
def bottom_boundary(x):
    tol = 1E-14
    return np.isclose(x[1], y_bottom, tol)


def determine_heat_flux(V_g, u, k, flux):
    """
    compute flux following http://hplgit.github.io/INF5620/doc/pub/fenics_tutorial1.1/tu2.html#tut-poisson-gradu
    :param V_g: Vector function space
    :param u: solution where gradient is to be determined
    :param k: thermal conductivity
    :param flux: returns calculated flux into this value
    """

    w = TrialFunction(V_g)
    v = TestFunction(V_g)

    a = inner(w, v) * dx
    L = inner(-k * grad(u), v) * dx
    # solve(a == L, flux)
    problem = LinearProblem(a, L) #, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    flux = problem.solve()


# Create mesh and define function space
nx = 100
ny = 25
nz = 1

fenics_dt = 0.01  # time step size
dt_out = 0.2  # interval for writing VTK files
y_top = 0
y_bottom = y_top - .25
x_left = 0
x_right = x_left + 1


mesh = RectangleMesh(MPI.COMM_WORLD, np.array([[x_left, y_bottom, 0],[x_right, y_top, 1]]), [nx,ny])#, cell_type=CellType.triangle)
scalar_element = FiniteElement("P", mesh.ufl_cell(), 1)
vector_element = VectorElement("P", mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, scalar_element)
V_g = FunctionSpace(mesh, vector_element)
# V = FunctionSpace(mesh, 'P', 1)
# V_g = VectorFunctionSpace(mesh, 'P', 1)

alpha = 1  # m^2/s, https://en.wikipedia.org/wiki/Thermal_diffusivity
k = 100  # kg * m / s^3 / K, https://en.wikipedia.org/wiki/Thermal_conductivity

# Define boundary condition
u_D = Function(V)
with u_D.vector.localForm() as loc:  # according to https://jorgensd.github.io/dolfinx-tutorial/chapter2/ns_code2.html#boundary-conditions
    loc.set(310)

# We will only exchange flux in y direction on coupling interface. No initialization necessary.
V_flux_y = V_g.sub(1)

# coupling_boundary = TopBoundary()
coupling_boundary = top_boundary
# bottom_boundary = BottomBoundary()

# Define initial value
u_n = Function(V, name="T")
u_n.interpolate(u_D)

# Adapter definition and initialization
# precice = Adapter(adapter_config_filename="precice-adapter-config.json")  # TODO

# precice_dt = precice.initialize(coupling_boundary, read_function_space=V, write_object=V_flux_y)  # TODO

# Create a FEniCS Expression to define and control the coupling boundary values
# coupling_expression = precice.create_coupling_expression()  # TODO

# Assigning appropriate dt
dt = Constant(mesh, 1)
# dt.assign(np.min([fenics_dt, precice_dt])) TODO

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
# F = u * v / dt * dx + alpha * dot(grad(u), grad(v)) * dx - u_n * v / dt * dx
# a, L = lhs(F), rhs(F)
a = alpha * dot(grad(u), grad(v)) * dx + u * v / dt * dx
L = u_n * v / dt * dx

# apply constant Dirichlet boundary condition at bottom edge
# apply Dirichlet boundary condition on coupling interface
# bcs = [DirichletBC(V, coupling_expression, coupling_boundary), DirichletBC(V, u_D, bottom_boundary)]  # TODO
bcs = []
dofs_bottom = locate_dofs_geometrical(V, bottom_boundary)
bcs.append(DirichletBC(u_D, dofs_bottom))

dofs_coupling = locate_dofs_geometrical(V, coupling_boundary)
value_coupling = Function(V)
with value_coupling.vector.localForm() as loc:
    loc.set(200)
bcs.append(DirichletBC(value_coupling, dofs_bottom))



# Time-stepping
u_np1 = Function(V)
t = 0
u_D.t = t + dt


with XDMFFile(MPI.COMM_WORLD, "result.xdmf", "w") as xdmf:
    xdmf.write_mesh(mesh)
    xdmf.write_function(u_n, 0)

    # trying to solve the problem
    problem = LinearProblem(a, L, bcs=bcs) #, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    u_np1 = problem.solve()
    u_n.interpolate(u_np1)
    xdmf.write_function(u_n, 1)
    exit()
    
    print("output vtk for time = {}".format(float(t)))
    n = 0

    fluxes = Function(V_g, name="Fluxes")

    while precice.is_coupling_ongoing():

        if precice.is_action_required(precice.action_write_iteration_checkpoint()):  # write checkpoint
            precice.store_checkpoint(u_n, t, n)

        read_data = precice.read_data()

        # Update the coupling expression with the new read data
        precice.update_coupling_expression(coupling_expression, read_data)

        dt.assign(np.min([fenics_dt, precice_dt]))

        # Compute solution
        # solve(a == L, u_np1, bcs)
        problem = LinearProblem(a, L, bcs=bcs) #, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        u_np1 = problem.solve()

        # Dirichlet problem obtains flux from solution and sends flux on boundary to Neumann problem
        determine_heat_flux(V_g, u_np1, k, fluxes)
        fluxes_y = fluxes.sub(1)  # only exchange y component of flux.
        precice.write_data(fluxes_y)

        precice_dt = precice.advance(dt(0))

        if precice.is_action_required(precice.action_read_iteration_checkpoint()):  # roll back to checkpoint
            u_cp, t_cp, n_cp = precice.retrieve_checkpoint()
            u_n.assign(u_cp)
            t = t_cp
            n = n_cp
        else:  # update solution
            u_n.assign(u_np1)
            t += float(dt)
            n += 1

        if precice.is_time_window_complete():
            tol = 10e-5  # we need some tolerance, since otherwise output might be skipped.
            if abs((t + tol) % dt_out) < 2 * tol:  # output if t is a multiple of dt_out
                print("output vtk for time = {}".format(float(t)))
                xdmf.write_function(u_n, t)

        # Update dirichlet BC
        u_D.t = t + float(dt)

# Hold plot
precice.finalize()
