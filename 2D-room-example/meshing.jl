# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

using Gmsh
using Gridap
using Gridap.Io
using GridapGmsh

# Generate the mesh for a fluid domain:
# a vertical cross-section of a auditorium-like geometry with staircase bottom
function generate_mesh(N)
    # Room parameters
    # room  width and height: room dimensions
    room_width = 8.0
    room_height = 4.5
    # dias : dias dimensions
    dias_anchor = [0.0, 0.0]
    dias_width = 1.5
    dias_height = 0.25
    # steps : distance to the first step, steps dimensions
    first_step_dist = 3.0
    steps_width = 1.0
    steps_height = 0.25
    steps_number = 5
    # tol : tolerance
    tol=1e-10

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.model.add("fluid")

    # Corners of the domain
    p = zeros(Int, 6+2*steps_number)
    p[1] = gmsh.model.geo.addPoint(0., dias_height, 0.)
    p[2] = gmsh.model.geo.addPoint(dias_width, dias_height, 0.)
    p[3] = gmsh.model.geo.addPoint(dias_width, 0., 0.)
    for j=1:steps_number
        p[4+2*(j-1)] = gmsh.model.geo.addPoint(first_step_dist+(j-1)*steps_width, (j-1)*steps_height, 0.)
        p[5+2*(j-1)] = gmsh.model.geo.addPoint(first_step_dist+(j-1)*steps_width, j*steps_height, 0.)
    end
    p[4+2*steps_number] = gmsh.model.geo.addPoint(room_width, steps_number*steps_height, 0.)
    p[5+2*steps_number] = gmsh.model.geo.addPoint(room_width, room_height, 0.)
    p[6+2*steps_number] = gmsh.model.geo.addPoint(0., room_height, 0.)

    # Lines of the domain
    for j=1:6+2*steps_number-1
        gmsh.model.geo.addLine(p[j], p[j+1], j)
    end
    gmsh.model.geo.addLine(p[6+2*steps_number], p[1], 6+2*steps_number)
    # Curve loop and surface
    gmsh.model.geo.addCurveLoop(Vector(1:6+2*steps_number), 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    # Synchronize the model to update the geometry
    gmsh.model.geo.synchronize()
    # Physical groups: 1 for the absorbing boundary and 2 for the rest of the boundary sides
    gmsh.model.addPhysicalGroup(1, Vector(1:3+steps_number*2), 1)
    gmsh.model.setPhysicalName(1, 1, "AbsorbingBoundary")
    gmsh.model.addPhysicalGroup(1, Vector(5+steps_number*2:6+2*steps_number), 2)
    gmsh.model.setPhysicalName(1, 2, "RigidBoundary")
    gmsh.model.addPhysicalGroup(1, Vector([4+steps_number*2]), 3)
    gmsh.model.setPhysicalName(1, 3, "WallPressure")
    # Physical group for the fluid domain
    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.setPhysicalName(2, 1, "FluidDomain")
    # Physical group for the rigid vertices
    gmsh.model.addPhysicalGroup(0, p[4+2*steps_number:5+2*steps_number], 1)
    gmsh.model.setPhysicalName(0, 1, "RigidVertices")
    # Generate 2D mesh
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", room_width/N)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Netgen")
    # Write mesh to file
    gmsh.write("./meshes/mesh_$(N).msh")
    gmsh.finalize()
    # This steps are to visualize the mesh and his labels using Gridap's API for Gmsh
    model = GmshDiscreteModel("./meshes/mesh_$(N).msh")
    writevtk(model,"./meshes//mesh_$(N)")
    return model
end
