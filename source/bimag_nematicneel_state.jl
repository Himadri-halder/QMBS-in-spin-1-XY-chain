# Constructs scar state from given basis list
function scar_state_1D(N::Int64,n::Int64,basis_list::Vector{Vector{Int64}})
    
    scar_state::Vector{Int64} = Vector{Int64}()
    norm_constant::Float64 = sqrt((factorial(big(N-n))*factorial(big(n)))/factorial(big(N)))
    
    for basis::Vector{Int64} in basis_list
        if all(x::Int64 -> x != 0, basis) && count(x::Int64 -> x == 1, basis) == n  # Only bimagnon states
            phase::Int64 = 1
            for (r::Int64,i::Int64) in enumerate(basis)
                if i == 1
                    phase *= (-1)^r
                end
            end
            push!(scar_state,phase)
        else
            push!(scar_state,0::Int64)
        end
    end
    
    return (norm_constant*scar_state)
end

# Constructs scar state from given basis list
function scar_state_2D(Nx::Int64,Ny::Int64,n::Int64,site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},
                        basis_list::Vector{Vector{Int64}})
    
    N = Nx*Ny
    scar_state::Vector{Int64} = Vector{Int64}()
    norm_constant::Float64 = sqrt((factorial(big(N-n))*factorial(big(n)))/factorial(big(N)))
    
    for basis::Vector{Int64} in basis_list
        if all(x::Int64 -> x != 0, basis) && count(x::Int64 -> x == 1, basis) == n  # Only bimagnon states
            phase::Int64 = 1
            for (r::Int64,i::Int64) in enumerate(basis)
                if i == 1
                    coords = site_coords[r-1]
                    x,y = coords[1],coords[2]
                    phase *= (-1)^(x+y)
                end
            end
            push!(scar_state,phase)
        else
            push!(scar_state,0::Int64)
        end
    end
    
    return (norm_constant*scar_state)
end

# Constructs nematic neel state from given basis list
function nematic_neel_1D(N::Int64,basis_list::Vector{Vector{Int64}})
    
    neel_state::Vector{Int64} = Vector{Int64}()
    norm_constant::Float64 = 1/2^(N/2)
    
    for basis::Vector{Int64} in basis_list
        # Only bimagnon states are present in nematic neel state
        if all(x::Int64 -> x != 0, basis) 
            phase::Int64 = 1
            for (r::Int64,i::Int64) in enumerate(basis)
                if i == -1
                    phase *= (-1)^(r+1)
                end
            end
            push!(neel_state,phase)
        else
            push!(neel_state,0::Int64)
        end
    end
    
    return(norm_constant*neel_state)
end

# Constructs nematic neel state from given basis list
function nematic_neel_2D(Nx::Int64,Ny::Int64,site_coords::OrderedDict{Int64,Tuple{Int64,Int64}},
                          basis_list::Vector{Vector{Int64}})
    
    N = Nx*Ny
    neel_state::Vector{Int64} = Vector{Int64}()
    norm_constant::Float64 = 1/2^(N/2)
    
    for basis::Vector{Int64} in basis_list
        # Only bimagnon states are present in nematic neel state
        if all(x::Int64 -> x != 0, basis)
            phase::Int64 = 1
            for (r::Int64,i::Int64) in enumerate(basis)
                if i == -1
                    coords = site_coords[r-1]
                    x,y = coords[1],coords[2]
                    phase *= (-1)^(x+y)
                end
            end
            push!(neel_state,phase)
        else
            push!(neel_state,0::Int64)
        end
    end
    
    return(norm_constant*neel_state)
end

# Constructs nematic ferro state from given basis list
function nematic_ferro_1D(N::Int64,basis_list::Vector{Vector{Int64}})
    
    ferro_state::Vector{Int64} = Vector{Int64}()
    norm_constant::Float64 = 1/2^(N/2)
    
    for basis::Vector{Int64} in basis_list
        # Only bimagnon states are present in nematic ferro state
        if all(x::Int64 -> x != 0, basis) 
            phase::Int64 = 1
            push!(ferro_state,phase)
        else
            push!(ferro_state,0::Int64)
        end
    end
    
    return(norm_constant*ferro_state)
end