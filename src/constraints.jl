abstract type BoundaryCondition end

struct Pinned <: BoundaryCondition end
struct Fixed <: BoundaryCondition end
struct Free <: BoundaryCondition end
struct Roller <: BoundaryCondition end
struct Hinge <: BoundaryCondition end
struct Joint <: BoundaryCondition end

function boundary_condition(::Pinned, AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    irow += 1
    _displ!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    irow += 1
    _bending!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    return irow
end

function boundary_condition(::Fixed, AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    irow += 1
    _displ!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    irow += 1
    _rotation!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    return irow
end

function boundary_condition(::Free, AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    irow += 1
    _bending!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    irow += 1
    _shear!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    return irow
end

function boundary_condition(::Roller, AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    if xe == 0.0
        sgn = -1.0
        irow += 1
        _displ!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn)
        _rotation!(AA, fii, T, irow - 2, jcol, Λ, L, Γ, xe, sgn)
        _bending!(AA, fii, T, irow - 1, jcol, Λ, L, Γ, xe, sgn)
    else
        irow += 1
        _displ!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
        irow += 1
        _rotation!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
        irow += 1
        _bending!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    end
    return irow
end

function boundary_condition(::Hinge, AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    if xe == 0.0
        sgn = -1.0
        _shear!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn)
        _displ!(AA, fii, T, irow - 2, jcol, Λ, L, Γ, xe, sgn)
        irow += 1
        _bending!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn)
    else
        irow += 1
        _displ!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
        irow += 1
        _bending!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
        irow += 1
        _shear!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    end
    return irow
end

function boundary_condition(::Joint, AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    if xe == 0.0
        sgn = -1.0
        _displ!(AA, fii, T, irow - 3, jcol, Λ, L, Γ, xe, sgn)
        _rotation!(AA, fii, T, irow - 2, jcol, Λ, L, Γ, xe, sgn)
        _bending!(AA, fii, T, irow - 1, jcol, Λ, L, Γ, xe, sgn)
        _shear!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn)
    else
        irow += 1
        _displ!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
        irow += 1
        _rotation!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
        irow += 1
        _bending!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
        irow += 1
        _shear!(AA, fii, T, irow, jcol, Λ, L, Γ, xe)
    end
    return irow
end

const boundary_conditions = Dict(
    "pinned" => Pinned(),
    "fixed" => Fixed(),
    "free" => Free(),
    "roller" => Roller(),
    "hinge" => Hinge(),
    "joint" => Joint()
)

function _displ!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn=1.0)
    if xe == 0.0
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+2] += sgn
        T[irow, jcol+4] += sgn
        push!(AA, copy(T))
        push!(fii, λ -> one(λ))
    else
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+1] += sgn
        push!(AA, copy(T))
        push!(fii, λ -> sin(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+2] += sgn
        push!(AA, copy(T))
        push!(fii, λ -> cos(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+3] += sgn
        push!(AA, copy(T))
        push!(fii, λ -> sinh(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+4] += sgn
        push!(AA, copy(T))
        push!(fii, λ -> cosh(λ * Λ))
    end
end

function _rotation!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn=1.0)
    if xe == 0.0
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+1] += sgn * (Λ / L)
        T[irow, jcol+3] += sgn * (Λ / L)
        push!(AA, copy(T))
        push!(fii, λ -> one(λ))
    else
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+1] += sgn * (Λ / L)
        push!(AA, copy(T))
        push!(fii, λ -> cos(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+2] -= sgn * (Λ / L)
        push!(AA, copy(T))
        push!(fii, λ -> sin(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+3] += sgn * (Λ / L)
        push!(AA, copy(T))
        push!(fii, λ -> cosh(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+4] += sgn * (Λ / L)
        push!(AA, copy(T))
        push!(fii, λ -> sinh(λ * Λ))
    end
end

function _bending!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn=1.0)
    if xe == 0.0
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+2] -= sgn * (Λ / L)^2 * Γ
        T[irow, jcol+4] += sgn * (Λ / L)^2 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> one(λ))
    else
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+1] -= sgn * (Λ / L)^2 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> sin(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+2] -= sgn * (Λ / L)^2 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> cos(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+3] += sgn * (Λ / L)^2 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> sinh(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+4] += sgn * (Λ / L)^2 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> cosh(λ * Λ))
    end
end

function _shear!(AA, fii, T, irow, jcol, Λ, L, Γ, xe, sgn=1.0)
    if xe == 0.0
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+1] -= sgn * (Λ / L)^3 * Γ
        T[irow, jcol+3] += sgn * (Λ / L)^3 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> one(λ))
    else
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+1] -= sgn * (Λ / L)^3 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> cos(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+2] += sgn * (Λ / L)^3 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> sin(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+3] += sgn * (Λ / L)^3 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> cosh(λ * Λ))
        # add new matrix and function
        T .= 0.0
        T[irow, jcol+4] += sgn * (Λ / L)^3 * Γ
        push!(AA, copy(T))
        push!(fii, λ -> sinh(λ * Λ))
    end
end