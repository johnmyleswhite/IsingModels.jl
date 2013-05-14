using Distributions

module IsingModels

using Distributions

export Ising

type Ising <: Distribution
    state::Array{Int8}
    isferromagnetic::Bool
    isperiodic::Bool
    temperature::Float64
    applied_field::Array{Int8}
    iterations::Int
end

function Ising(n::Integer;
               m::Integer = -1,
               isferromagnetic::Bool = true,
               isperiodic::Bool = true,
               temperature::Real = 1.0,
               iterations::Integer = -1,
               applied_field::Matrix = Array(Int8, 0, 0))

    if m == -1
        m = n
    end

    if iterations == -1
        iterations = n * n * 1_000
    end

    state = Array(Int8, n, m)
    fill!(state, -1)

    Ising(state,
          isferromagnetic,
          isperiodic,
          temperature,
          applied_field,
          iterations)
end

function determine_j(d::Ising)
    if d.isferromagnetic
        return int8(1)
    else
        return int8(-1)
    end
end

function determine_local_field(d::Ising, r::Integer, c::Integer, n::Integer, m::Integer, j::Int8)
    # Sum neighbors
    s = int8(0)

    # Upper neighbor
    if r > 1
        s += d.state[r - 1, c] * j
    else
        if d.isperiodic
            s += d.state[n, c] * j
        end
    end

    # Lower neighbor
    if r + 1 <= n
        s += d.state[r + 1, c] * j
    else
        if d.isperiodic
            s += d.state[1, c] * j
        end
    end
  
    # Left neighbor
    if c > 1
        s += d.state[r, c - 1] * j
    else
        if d.isperiodic
            s += d.state[r, m] * j
        end
    end

    # Right neighbor
    if c + 1 <= m
        s += d.state[r, c + 1] * j
    else
        if d.isperiodic
            s += d.state[r, c] * j
        end
    end

    # Apply external field
    if !isempty(d.applied_field)
        s += d.applied_field[r, c]
    end

    return s
end

function flipbit!(d::Ising, n::Integer, m::Integer, j::Int8)
    r, c = rand(1:n), rand(1:m)

    b = determine_local_field(d, r, c, n, m, j)

    p = 1 / (1 + exp(-2 * b / d.temperature))

    if !isfinite(p)
        d.state[r, c] = int8(1)
    else
        if rand() < p
            d.state[r, c] = int8(1)
        else
            d.state[r, c] = int8(-1)
        end
    end

    return
end

function Distributions.rand!(d::Ising)
    # Determine constants
    j = determine_j(d)
    n, m = size(d.state)

    # Randomly initialize state
    for r in 1:n
        for c in 1:m
            if rand() < 0.5
                d.state[r, c] = int8(-1)
            else
                d.state[r, c] = int8(1)
            end
        end
    end

    # Flip bits for a fixed number of iterations
    for iteration in 1:d.iterations
        flipbit!(d, n, m, j)
    end

    return
end

function Distributions.rand(d::Ising)
    rand!(d)
    return copy(d.state)
end

Distributions.mean(d::Ising) = NaN
Distributions.var(d::Ising) = NaN

end # module
