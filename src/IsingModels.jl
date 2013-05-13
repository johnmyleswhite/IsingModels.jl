using Distributions

module IsingModels

using Distributions

export Ising

type Ising <: Distribution
    state::Array{Int}
    isferromagnetic::Bool
    isperiodic::Bool
    temperature::Float64
    applied_field::Array{Int}
    iterations::Int
end

function Ising(n::Integer;
               isferromagnetic::Bool = true,
               isperiodic::Bool = true,
               temperature::Real = 1.0,
               iterations::Integer = -1)
    state = Array(Int, n, n)
    fill!(state, -1)
    if iterations == -1
        iterations = n * n * 1_000
    end
    Ising(state,
          isferromagnetic,
          isperiodic,
          temperature,
          Array(Int, 0, 0),
          iterations)
end

function determine_j(d::Ising)
    if d.isferromagnetic
        return 1
    else
        return -1
    end
end

function determine_local_field(d::Ising,
                               row_index::Integer,
                               col_index::Integer)

    # Determine constants
    j = determine_j(d)
    n, m = size(d.state)

    # Sum neighbors
    s = 0

    # Upper neighbor
    if row_index > 1
        s += d.state[row_index - 1, col_index] * j
    else
        if d.isperiodic
            s += d.state[n, col_index] * j
        end
    end

    # Lower neighbor
    if row_index + 1 <= n
        s += d.state[row_index + 1, col_index] * j
    else
        if d.isperiodic
            s += d.state[1, col_index] * j
        end
    end
  
    # Left neighbor
    if col_index > 1
        s += d.state[row_index, col_index - 1] * j
    else
        if d.isperiodic
            s += d.state[row_index, m] * j
        end
    end

    # Right neighbor
    if col_index + 1 <= m
        s += d.state[row_index, col_index + 1] * j
    else
        if d.isperiodic
            s += d.state[row_index, col_index] * j
        end
    end

    if !isempty(d.applied_field)
        s += d.applied_field[row_index, col_index]
    end

    return s
end

function flipbit!(d::Ising)
    n, m = size(d.state)

    row_index, col_index = rand(1:n), rand(1:m)

    b = determine_local_field(d, row_index, col_index)

    p = 1 / (1 + exp(-2 * b / d.temperature))

    if !isfinite(p)
        d.state[row_index, col_index] = 1
    else
        if rand() < p
            d.state[row_index, col_index] = 1
        else
            d.state[row_index, col_index] = -1    
        end
    end

    return
end

function Distributions.rand!(d::Ising)
    # Randomly initialize state
    n, m = size(d.state)
    for r in 1:n
        for c in 1:m
            if rand() < 0.5
                d.state[r, c] = -1
            else
                d.state[r, c] = 1
            end
        end
    end

    # Flip bits for a fixed number of iterations
    for iteration in 1:d.iterations
        flipbit!(d)
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
