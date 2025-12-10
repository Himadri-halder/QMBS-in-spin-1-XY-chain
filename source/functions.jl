# Finds index of a basis in a basis list
function find_state_index(basis::Vector{Int64},basis_list::Vector{Vector{Int64}})
    
    for (i::Int64,s::Vector{Int64}) in enumerate(basis_list)
        if s == basis
            return i
        end
    end
    return -1
end

# Performs spin-flip
function flip(l0::Int64,l1::Int64)
    if l0 == l1
        if l0 == 0
            return 1,[[1,-1],[-1,1]]
        else
            return 0,[[l0,l1]]
        end
    else
        if l0 == -1 && l1 == 0
            return 1,[[0,-1]]
        elseif l0 == -1 && l1 == 1
            return 1,[[0, 0]]
        elseif l0 == 0 && l1 == -1
            return 1,[[-1, 0]]
        elseif l0 == 0 && l1 == 1
            return 1,[[1, 0]]
        elseif l0 == 1 && l1 == -1
            return 1,[[0, 0]]
        elseif l0 == 1 && l1 == 0
            return 1,[[0, 1]]
        end
    end
end

# Generates time evolved state of a given state
function state_evolution(exp_H::Matrix{ComplexF64},energy::Vector{Float64},eigenvector::Matrix{Float64},
                         eigenvector_inv::Matrix{Float64},time::Float64,state::Vector{T},
                         evolved_state::Vector{ComplexF64}) where {T<:Union{Float64,ComplexF64}}
    
    exp_H .= eigenvector*Diagonal(exp.(-1im*time*energy))*eigenvector_inv
    mul!(evolved_state,exp_H,state)

    return evolved_state
end

# Computes von neumann bipartite entanglement entropy of a system (full system)
function bipartite_entang_entropy_full(state::Vector{T},s::Int64,N_A::Int64,
                                       N_B::Int64) where {T<:Union{Float64,ComplexF64}}   
    
    dim_A::Int64 = convert(Int64,(2*s+1)^N_A)
    dim_B::Int64 = convert(Int64,(2*s+1)^N_B)
    psi_matrix::Matrix{T} = reshape(state,(dim_A,dim_B))
    singular_values::Vector{Float64} = svdvals(psi_matrix)
    p::Vector{Float64} = singular_values.^2
    p_non_zero::Vector{Float64} = filter(x -> x>1e-15, p)
    entropy::Float64 = -sum(p_non_zero.*log.(p_non_zero))
    
    return entropy
end

# Splits a list at k point into two lists (Useful at splitting full basis into subsystems' basis)
function split_list(list::Vector{T},k::Int64) where T<:Union{Int64,Float64,ComplexF64}
    return list[1:k],list[k+1:end]
end

# Finds index of a basis of a subsystem in full basis
function find_subsys_basis_index(k::Int64,basis::Vector{Int64},
                                 total_basis_list::Vector{Tuple{Vector{Int64},Vector{Int64}}})
        
    index_list::Vector{Int64} = Vector{Int64}(undef,0)
    for (i::Int64,s::Tuple{Vector{Int64},Vector{Int64}}) in enumerate(total_basis_list)
        if s[k] == basis
            push!(index_list,i)
        end
    end
    return index_list
end

# Computes von neumann bipartite entanglement entropy of a symmetry resolved state
function bipartite_entang_entropy_sym(state::Vector{T},basis_list::Vector{Vector{Int64}},N_A::Int64,
                                       N_B::Int64) where {T<:Union{Float64,ComplexF64}} 
   
    """Feed basis_list and state corresponding to that particular symmetry sector"""
    """N_A and N_B mean no. of sites in A and B subsytems"""
    A_B_basis_list::Vector{Tuple{Vector{Int64},Vector{Int64}}} =
                     Vector{Tuple{Vector{Int64},Vector{Int64}}}(undef,length(basis_list))
    A_basis_list::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,0)
    for (i,basis) in enumerate(basis_list)
        A_basis,B_basis = split_list(basis,N_A)
        push!(A_basis_list,A_basis)
        A_B_basis_list[i] = (A_basis,B_basis)
    end
    unique!(A_basis_list)
    
    H_A::Int64 = length(A_basis_list)
    rho_mat_A::Matrix{T} = zeros(T,H_A,H_A)
    A_basis_to_index = Dict{Vector{Int64},Int64}()
    for (k,A_basis) in enumerate(A_basis_list)
        A_basis_to_index[A_basis] = k
    end
    for (ind_c,A_basis) in enumerate(A_basis_list)
        A_index_list = find_subsys_basis_index(1,A_basis,A_B_basis_list)
        for i in A_index_list
            B_basis = A_B_basis_list[i][2]
            B_index_list = find_subsys_basis_index(2,B_basis,A_B_basis_list)
            for j in B_index_list
                A_basis_1 = A_B_basis_list[j][1]
                ind_r = A_basis_to_index[A_basis_1]
                rho_mat_A[ind_c,ind_r] += state[i]*conj(state[j])
            end
        end
    end    
    p::Vector{Float64} = eigvals(Hermitian(rho_mat_A))
    p_non_zero::Vector{Float64} = filter(x -> x>1e-15, p)
    entropy::Float64 = -sum(p_non_zero.*log.(p_non_zero)) 
    
    return entropy
end

# Moving Window Variance
function moving_variance(F::Vector{Float64},window_size::Int64,step_size::Int64)
    
    var_list::Vector{Float64} = Vector{Float64}(undef,0)
    var_time_indices::Vector{Float64} = Vector{Float64}(undef,0)
    for i in 1:step_size:length(F)-window_size
        segment = F[i:i+window_size]
        push!(var_list,var(segment))  # Compute variance
        push!(var_time_indices,i)
    end
    
    return var_time_indices,var_list
end

# Finds first peak in fidelity
function fidelity_first_peak(time::Vector{Float64},fidelity::Vector{Float64})
    
    peak_time::Float64 = -1.0
    peak_fidelity::Float64 = -1.0
    for i in 2:length(fidelity)
        if (fidelity[i-1] < fidelity[i]) && (fidelity[i] > fidelity[i+1])
            peak_time = time[i]
            peak_fidelity = fidelity[i]
            break
        end
    end
    
    return (peak_time,peak_fidelity)
end    

# Detect Noisy Region: Variance Falls Below a Threshold
function detect_noisy_region(var_list::Vector{Float64},var_time_indices::Vector{Float64},
                             time_list::Vector{Float64},threshold::Float64,h::Float64)
    
    var_noise_start_idx::Int64 = findfirst(x -> x < threshold*maximum(var_list),var_list)
    noise_start_idx::Int64 = var_time_indices[var_noise_start_idx]
    noise_start_time::Float64 = h*time_list[noise_start_idx]
    
    return noise_start_idx,noise_start_time
end

# Calculates centered moving average
function moving_average(F::Vector{Float64},noise_start_idx::Int64,window::Int64)
    
    new_F = similar(F)
    n::Int64 = length(F)
    for i in 1:n
        if i < noise_start_idx+window
            new_F[i] = F[i]
        elseif i >= noise_start_idx+window
            left = max(1,i-window)  #To handle edge cases
            right = min(n,i+window)
            new_F[i] = mean(F[left:right])
        end
    end
    
    return new_F
end

# Detects the time when fidelity first hits a value below the tolerance
function decay_time(time_list::Vector{Float64},decay_parameter_list::Vector{Float64},tol::Float64,
                    h::Float64)
    
    ht = h*time_list
    min_p = Inf
    pos::Int64 = 1
    a = Inf
    
    for (i::Int64,p::Float64) in enumerate(decay_parameter_list)
        if p < tol
            if p < min_p
                min_p = p
                pos = i
                a = p
            else
                break
            end
        end
    end
    
    return a,ht[min(pos,length(ht))]::Float64
end