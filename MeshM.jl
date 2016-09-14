module MeshM

export Mesh1d

" A simple 1d mesh class, assumes constant spacing."
type Mesh1d
    "Number of cells."
    nCells::Int64

    "The step size."
    h::Float64

    "The time-step."
    dt::Float64

    "Location of the cell centres."
    xCells::Vector{Float64}

    " Location of the faces."
    xFaces::Vector{Float64}

    """
    Constructor.

    Parameters:
    nCells::Int64
        Number of cells
    xStart::Real
        Beginning of the interval
    xEnd::Real
        End of the interval
    """
    function Mesh1d(nCells::Int64, xStart::Real, xEnd::Real, dt::Real)
        if xEnd < xStart
            error("Ill-defined interval!")
        else
            xFaces = linspace(xStart, xEnd, nCells+1)
            h = step(xFaces)
            xCells = (xFaces + 0.5*h)[1:end-1]
            return new(nCells, h, dt, xCells, xFaces)
        end
    end
end
end
