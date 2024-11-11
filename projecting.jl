# Copyright (c) 2024 Ashwin. S. Nayak, Andrés Prieto, Daniel Fernández Comesaña
#
# This file is part of pressure-projections
#
# SPDX-License-Identifier:  MIT

"""
Computation of the different projections of the solution of the vector-valued Helmholtz equation
attending to the different finite element spaces (continuous and discontinuous Galerkin approaches) used in the variational formulation
for the  and L^2 and H^1-projections. Additionally, the scalar Helmholtz equation is solved using the continuous Galerkin approach
to obtain the pressure field in H^1. All the projections are computed using the Gridap.jl package.
"""
# Load the packages
using Gridap
using Gridap.Fields
using Gridap.Geometry

# Computation of the displacement field using a Raviart-Thomas finite element discretization
function RT_solve(k, orderFE, model, boundary_data)
    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain
    xp = get_physical_coordinate(Ω)

    # Define the finite element space: Raviart-Thomas of orderFE-1 (the lowest order is ZERO!!)
    reffe = ReferenceFE(raviart_thomas, Float64, orderFE-1)
    V = TestFESpace(model, reffe, conformity=:Hdiv, dirichlet_tags=["RigidBoundary"], vector_type=Vector{ComplexF64})

    # Define the trial function with null Dirichlet boundary conditions
    uD = boundary_data["RigidBoundary"]
    U = TrialFESpace(V, uD)

    # Define the measure for the fluid and porous domains
    degree = 2*orderFE # Quadrature degree
    dΩ = Measure(Ω, degree)

    # Define the boundary terms in the absorbing boundary
    absorbing_tags = ["AbsorbingBoundary"]
    Za = boundary_data["AbsorbingBoundary"]
    Γa = BoundaryTriangulation(model, tags=absorbing_tags)
    dΓa = Measure(Γa, degree)
    na = get_normal_vector(Γa) # Normal vector to the boundary

    # Define the boundary term in the wall pressure boundary
    wallpressure_tags = ["WallPressure"]
    pD = boundary_data["WallPressure"]
    Γp = BoundaryTriangulation(model, tags=wallpressure_tags)
    dΓp = Measure(Γp, degree)
    np = get_normal_vector(Γp) # Normal vector to the boundary

    # Define the variational problem: first factor in the cdot product is the conjugated test function
    a(u, v) = ∫(divergence(v) ⋅ divergence(u))dΩ - k^2*∫(v ⋅ u)dΩ - 1im*k*Za*∫((na ⋅ v) ⋅ (na ⋅ u))dΓa
    b(v) = ∫(-(np ⋅ v) ⋅ (pD ∘ xp))dΓp

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    uh = solve(op)

    return uh
end

# Computation of the pressure field using a standard continuous Galerkin finite element discretization
function CG_solve(k, orderFE, model, boundary_data)
    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain
    xp = get_physical_coordinate(Ω)

    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"DirichletPressure",["RigidVertices", "WallPressure"])

    # Define the finite element space: scalar Discontinuous polynomials of a given order
    reffe = ReferenceFE(lagrangian, Float64, orderFE)
    V = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="DirichletPressure", vector_type=Vector{ComplexF64})

    # Define the trial function with null Dirichlet boundary conditions
    pD = boundary_data["WallPressure"]
    U = TrialFESpace(V, pD)

    # Define the measure for the fluid and porous domains
    degree = 2*orderFE # Quadrature degree
    dΩ = Measure(Ω, degree)

    # Define the boundary terms in the absorbing boundary
    absorbing_tags = ["AbsorbingBoundary"]
    Za = boundary_data["AbsorbingBoundary"]
    Γa = BoundaryTriangulation(model, tags=absorbing_tags)
    dΓa = Measure(Γa, degree)
    na = get_normal_vector(Γa) # Normal vector to the boundary

    # Define the boundary term in the rigid wall boundary
    rigidwall_tags = ["RigidBoundary"]
    uD = boundary_data["RigidBoundary"]
    Γp = BoundaryTriangulation(model, tags=rigidwall_tags)
    dΓp = Measure(Γp, degree)
    np = get_normal_vector(Γp) # Normal vector to the boundary

    # Define the variational problem: first factor in the cdot product is the conjugated test function
    a(u, v) =  ∫(∇(v) ⋅ ∇(u))dΩ - k^2*∫(v ⋅ u)dΩ - 1im*k/Za*∫(v ⋅ u)dΓa
    b(v) = k^2*∫(v ⋅ (np ⋅ (uD ∘ xp)))dΓp

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    ph = solve(op)

    return ph
end

# Computation of the L^2-projection of the pressure field using a discontinuous Galerkin finite element discretization
function L2_DG_projection(uh, orderFE, model)
    # Define the finite element space: scalar Discontinuous polynomials of a given order
    reffe = ReferenceFE(lagrangian, Float64, orderFE-1)
    V = TestFESpace(model, reffe, conformity=:L2, vector_type=Vector{ComplexF64})
    U = TrialFESpace(V)

    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain

    # Define the measure for the fluid and porous domains
    degree = 2*orderFE # Quadrature degree
    dΩ = Measure(Ω, degree)

    # Define the variational problem
    ph = -divergence(uh)
    a(u, v) = ∫(v ⋅ u)dΩ
    b(v) = ∫(v ⋅ ph)dΩ

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    ph = solve(op)

    return ph
end

# Computation of the L^2-projection of the pressure field using a continuous Galerkin finite element discretization
function L2_CG_projection(uh, orderFE, model)
    # Define the finite element space: scalar Discontinuous polynomials of a given order
    reffe = ReferenceFE(lagrangian, Float64, orderFE)
    V = TestFESpace(model, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
    U = TrialFESpace(V)

    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain

    # Define the measure for the fluid and porous domains
    degree = 2*orderFE # Quadrature degree
    dΩ = Measure(Ω, degree)

    # Define the variational problem
    ph = -divergence(uh)
    a(u, v) = ∫(v ⋅ u)dΩ
    b(v) = ∫(v ⋅ ph)dΩ

    # Define alternative variational problem using the divergence theorem on uh
    # Γ = BoundaryTriangulation(model)
    # dΓ = Measure(Γ, degree)
    # nb = get_normal_vector(Γ) # Normal vector to the transducer boundary
    # a(u, v) = ∫(v ⋅ u)dΩ
    # b(v) = ∫(∇(v) ⋅ uh)dΩ - ∫(v ⋅ (nb ⋅ uh))dΓ

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    ph = solve(op)

    return ph
end

# Computation of the H^1-projection of the pressure field using a continuous Galerkin finite element discretization
function H1_projection(k, uh, orderFE, model)
    # Define the finite element space: scalar Lagrange polynomials of a given order
    reffe = ReferenceFE(lagrangian, Float64, orderFE)
    V = TestFESpace(model, reffe, conformity=:H1, vector_type=Vector{ComplexF64})
    U = TrialFESpace(V)

    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain

    # Define the measure for the fluid and porous domains
    degree = 2*orderFE # Quadrature degree
    dΩ = Measure(Ω, degree)

    # Define the variational problem
    ph = -divergence(uh)
    a(u, v) = k^2*∫(v ⋅ u)dΩ + ∫(∇(v)⋅∇(u))dΩ
    b(v) = k^2*∫(v ⋅ ph)dΩ + k^2*∫(∇(v) ⋅ uh)dΩ

    # Define alternative variational problem using the divergence theorem on uh using H^1-norm
    # Γ = BoundaryTriangulation(model)
    # dΓ = Measure(Γ, degree)
    # nb = get_normal_vector(Γ) # Normal
    # ph = -divergence(uh)
    # a(u, v) = (1+k^2)*∫(v ⋅ u)dΩ
    # b(v) = ∫(v ⋅ ph)dΩ + k^2*∫(∇(v) ⋅ uh)dΩ - k^2*∫(v ⋅ (nb ⋅ uh))dΓ

    # # Define alternative variational problem using the divergence theorem on uh using the energy norm
    # Γ = BoundaryTriangulation(model)
    # dΓ = Measure(Γ, degree)
    # nb = get_normal_vector(Γ) # Normal
    # ph = -divergence(uh)
    # a(u, v) = 2*∫(v ⋅ u)dΩ
    # b(v) = ∫(v ⋅ ph)dΩ + ∫(∇(v) ⋅ uh)dΩ - ∫(v ⋅ (nb ⋅ uh))dΓ

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    ph = solve(op)

    return ph
end

# Computation of the H^1-projection of the pressure field using a discontinuous Galerkin finite element discretization
function H2_projection_DG(k, uh, orderFE, model, N)
    # Define the finite element space: P1+Raviart-Thomas of a given order
    reffeₚ = ReferenceFE(lagrangian, Float64, orderFE)
    V = TestFESpace(model, reffeₚ, conformity=:L2, vector_type=Vector{ComplexF64})
    U = TrialFESpace(V)

    # Define the fluid domain
    Ω = Triangulation(model) # Computational domain

    # Define the measure for the fluid and porous domains
    degree = 2*orderFE # Quadrature degree
    dΩ = Measure(Ω, degree)

    # Define boundary of the triangulation
    Γ = BoundaryTriangulation(model)
    dΓ = Measure(Γ, degree)
    n_Γ = get_normal_vector(Γ) # Normal

    # Define the inner boundary of the triangulation
    Λ = SkeletonTriangulation(model)
    dΛ = Measure(Λ,degree)
    n_Λ = get_normal_vector(Λ) # Normal

    # Compute the pressure from the divergence
    ph = -divergence(uh)
    # Using that ∇(p) = k^2*u, and Δ(p)=k^2*div(u)

    # Define the variational problem
    ph = -divergence(uh)
    γ = orderFE*(orderFE+1.0) # Penalty parameter =8 (originally)
    # Compute mesh size
    h = 1.0/N
    # Variational terms in the interior of the elements and on their edges
    a_Ω(u, v) = k^2*∫(v ⋅ u)dΩ + ∫(∇(v) ⋅ ∇(u))dΩ #+ ∫(Δ(v) ⋅ Δ(u))dΩ
    a_Λ(u, v) = ∫( (γ/h)*jump(v) ⋅ jump(u) )dΛ

    # Variational formulation
    a(u, v) = a_Ω(u, v) + a_Λ(u, v)
    b(v) = k^2*∫(v ⋅ ph)dΩ + k^2*∫(∇(v) ⋅ uh)dΩ #- k^2*∫(Δ(v) ⋅ ph)dΩ

    # Assembly the system
    op = AffineFEOperator(a, b, U, V)
    # Solve the system
    ph = solve(op)

    return ph
end
