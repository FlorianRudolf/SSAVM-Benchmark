#!/usr/bin/python

import sys
from dolfin import *
from math import pi



def make_boundary_facet_function(mesh, index):
  facetfunction = FacetFunction('size_t', mesh)
  facetfunction.set_all(0)

  for f in facets(mesh):
    if f.exterior():
      facetfunction[f] = index

  return facetfunction



def mark_facet_function(mesh, vertexfunc, ff, lines, index):
  for f in facets(mesh):
    if not f.exterior():
      continue

    found = True

    for v in vertices(f):
      val = vertexfunc[v]
      if not val in lines:
        found = False
        break

    if found:
      ff[f] = index



def make_facet_function(mesh, vertexfunc, u0_lines, u0_index, u1_lines, u1_index, u0_lines_dual, u0_index_dual, u1_lines_dual, u1_index_dual):
  facetfunction = FacetFunction('size_t', mesh)
  facetfunction.set_all(0)

  mark_facet_function(mesh, vertexfunc, facetfunction, u0_lines, u0_index)
  mark_facet_function(mesh, vertexfunc, facetfunction, u1_lines, u1_index)
  mark_facet_function(mesh, vertexfunc, facetfunction, u0_lines_dual, u0_index_dual)
  mark_facet_function(mesh, vertexfunc, facetfunction, u1_lines_dual, u1_index_dual)

  return facetfunction


def do_fem(filename, index):
  # Mesh and function space
  #mesh = UnitSquareMesh(2, 2)

  ofn = filename[filename.find("/")+1:]

  mesh = Mesh(filename + ".xml")
  V = FunctionSpace(mesh, 'Lagrange', 1)

  tdim = mesh.topology().dim()
  mesh.init(tdim-1, tdim)

  vertexfunc = MeshFunction('size_t', mesh, filename + "_bnd.xml")




  if index == 0:
    facetfunction = make_facet_function(mesh, vertexfunc, [9,0,1,19,10], 1, [4,5,6,7,14,15,16], 2, [7,8,9,17,18], 3, [1,2,3,4,11,12,13], 4)
  elif index == 1:
    facetfunction = make_facet_function(mesh, vertexfunc, [1,2,3,11,12], 1, [6,7,8,9,16,17,18], 2, [9,0,1,19,10], 3, [3,4,5,6,13,14,15], 4)
  elif index == 2:
    facetfunction = make_facet_function(mesh, vertexfunc, [3,4,5,13,14], 1, [8,9,0,1,18,19,10], 2, [1,2,3,11,12], 3, [5,6,7,8,15,16,17], 4)
  elif index == 3:
    facetfunction = make_facet_function(mesh, vertexfunc, [5,6,7,15,16], 1, [0,1,2,3,10,11,12], 2, [3,4,5,13,14], 3, [7,8,9,0,17,18,19], 4)
  elif index == 4:
    facetfunction = make_facet_function(mesh, vertexfunc, [7,8,9,17,18], 1, [2,3,4,5,12,13,14], 2, [5,6,7,15,16], 3, [9,0,1,2,19,10,11], 4)





  bc1 = DirichletBC(V, Constant(0), facetfunction, 1)
  bc2 = DirichletBC(V, Constant(1), facetfunction, 2)

  bcs = [bc1, bc2]

  # Variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  f = Constant(0.0)
  a = inner(nabla_grad(u), nabla_grad(v))*dx
  L = f*v*dx

  # compute solution
  u = Function(V)
  solve(a == L, u, bcs)

  flux = project(-grad(u), VectorFunctionSpace(mesh, 'Lagrange', 1))
  flux_normal = dot(flux, FacetNormal(mesh))

  file = File("results/solution_" + str(index) + "_" + ofn + ".pvd")
  file << u





  bc0_dual = DirichletBC(V, Constant(0), facetfunction, 3)
  bc1_dual = DirichletBC(V, Constant(1), facetfunction, 4)

  bcs_dual = [bc0_dual, bc1_dual]

  ### Variational problem
  u_dual = TrialFunction(V)
  a_dual = inner(nabla_grad(u_dual), nabla_grad(v))*dx
  L_dual = f*v*dx

  # compute solution
  u_dual = Function(V)
  solve(a_dual == L_dual, u_dual, bcs_dual)

  file = File("results/solution_dual_" + str(index) + "_" + ofn +  ".pvd")
  file << u_dual



set_log_active(False)

info(parameters, True)
prm = parameters['krylov_solver'] # short form
prm['absolute_tolerance'] = 1E-16
prm['relative_tolerance'] = 1E-12
prm['maximum_iterations'] = 10000
set_log_level(PROGRESS)

#filenames = ["gs_c0.262144", "gs_c0.131072", "gs_c0.065536", "gs_c0.032768", "gs_c0.016384", "gs_c0.008192", "gs_c0.004096", "gs_c0.002048", "gs_c0.001024", "gs_c0.000512", "gs_c0.000256", "gs_c0.000128", "gs_c6.4e-05", "gs_c3.2e-05", "gs_c1.6e-05", "gs_c8e-06", "gs_c4e-06", "gs_c2e-06", "gs_c1e-06"]
#filenames = ["gs_c0.262144", "gs_c0.131072", "gs_c0.065536", "gs_c0.032768", "gs_c0.016384", "gs_c0.008192", "gs_c0.004096"]
#filenames = ["gs_c0.000128"]
#filenames = ["gs_c0.008192"]


filenames = ["gs_c0.016384", "gs_c0.001024", "gs_c0.000128", "gs_c1.6e-05"]
#filenames = ["gs_c0.016384"]

def ftos(x):
  s = "%.1E" % x
  [a,b] = s.split("E")
  if b[1] == '0':
    b = "-" + b[2] + "\phantom{0}"
    #b[1] = b[0]
    #b[0] = '0'
    #b = b.replace("0", "\phantom{0}")
  return a + " \\times 10^{" + b + "}"


for filename in filenames:

  print(filename)

  mesh_nonsym = Mesh("data/" + filename + "_nonsym.xml")
  mesh_rotsym = Mesh("data/" + filename + "_rotsym.xml")
  mesh_reflrotsym = Mesh("data/" + filename + "_reflrotsym.xml")

  for i in range(0,5):
    do_fem("data/" + filename + "_nonsym", i)
    do_fem("data/" + filename + "_rotsym", i)
    do_fem("data/" + filename + "_reflrotsym", i)




