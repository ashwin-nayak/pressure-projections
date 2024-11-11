# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh

# Generate the mesh for a fluid domain: a quarter of an annulus of interior radius 1.0 and exterior radius 2.0
function generate_mesh(N)
    # Initialize Gmsh
    gmsh.initialize()
    gmsh.model.add("fluid")
    # Define the geometry
    # Origin point
    p0 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    # Corners of the quarter of the annular domain
    p1 = gmsh.model.geo.addPoint(1.0, 0., 0.)
    p2 = gmsh.model.geo.addPoint(2.0, 0., 0.)
    p3 = gmsh.model.geo.addPoint(0.0, 2.0, 0.)
    p4 = gmsh.model.geo.addPoint(0.0, 1.0, 0.)
    # Lines of the quarter of the annular domain
    gmsh.model.geo.addLine(p1, p2, 1)
    gmsh.model.geo.addLine(p3, p4, 2)
    # Circular curve between p2 and p3
    gmsh.model.geo.addCircleArc(p2, p0, p3, 3)
    #gmsh.model.geo.addLine(p2, p3, 3)
    # Circular curve between p3 and p4
    gmsh.model.geo.addCircleArc(p4, p0, p1, 4)
    #gmsh.model.geo.addCircleArc(p4, p1, 4)
    # Curve loop and surface
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    # Get transfinite mesh
    gmsh.model.geo.mesh.setTransfiniteCurve(1, N)
    gmsh.model.geo.mesh.setTransfiniteCurve(2, N)
    gmsh.model.geo.mesh.setTransfiniteCurve(3, ceil(π*N))
    gmsh.model.geo.mesh.setTransfiniteCurve(4, ceil(π*N))
    gmsh.model.geo.mesh.setTransfiniteSurface(1, "Alternate")
    # Synchronize the model to update the geometry
    gmsh.model.geo.synchronize()
    # Physical groups: 1 for the absorbing boundary and 2 for the rest of the boundary sides
    gmsh.model.addPhysicalGroup(1, [3], 1)
    gmsh.model.setPhysicalName(1, 1, "AbsorbingBoundary")
    gmsh.model.addPhysicalGroup(1, [1, 2], 2)
    gmsh.model.setPhysicalName(1, 2, "RigidBoundary")
    gmsh.model.addPhysicalGroup(1, [4], 3)
    #gmsh.model.addPhysicalGroup(1, [3, 4], 3)
    gmsh.model.setPhysicalName(1, 3, "WallPressure")
    # Physical group for the fluid domain
    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.setPhysicalName(2, 1, "FluidDomain")
    # Physical group for the rigid vertices
    gmsh.model.addPhysicalGroup(0, [p1, p4], 1)
    gmsh.model.setPhysicalName(0, 1, "RigidVertices")
    # Generate 2D mesh
    gmsh.model.mesh.generate(2)
    # Write mesh to file
    gmsh.write("./meshes/mesh_$(N).msh")
    gmsh.finalize()
    # This steps are to visualize the mesh and his labels using Gridap's API for Gmsh
    model = GmshDiscreteModel("./meshes/mesh_$(N).msh")
    # writevtk(model,"./meshes//mesh_$(N)")
    return model
end
