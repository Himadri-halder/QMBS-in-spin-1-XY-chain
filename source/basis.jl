## Generates basis list for a particular mz block (1D Lattice)
function gen_basis_mz_1D(N::Int64,mz::Int64)
    
    basis = [collect(b) for b in product(fill(-1:1,N)...) if sum(b) == mz]
    
    return basis::Vector{Vector{Int64}}
end

## Generates basis list for a particular mz block (2D Lattice)
function gen_basis_mz_2D(Nx::Int64,Ny::Int64,mz::Int64)
    
    basis = [collect(b) for b in product(fill(-1:1,Nx*Ny)...) if sum(b) == mz]
    
    return basis::Vector{Vector{Int64}}
end

## Generates basis list for a particular mz block (1D Lattice) [optimized] 
function gen_basis_mz_1D_op(N::Int64,mz::Int64)
    
    spin_vals = [1,0,-1]
    basis_states::Vector{Vector{Int64}} = Vector{Vector{Int64}}()

    # Stack-based iterative version to avoid recursive function overhead
    initial_state::Vector{Int64} = zeros(Int64,N)
    stack = [(initial_state,1,0)]  # (state,position,current_mz)

    while !isempty(stack)
        (state,pos,current_mz) = pop!(stack)

        # Check if we've assigned all spins
        if pos > N
            if current_mz == mz
                push!(basis_states,copy(state))  # Copy to avoid mutation issues
            end
            continue
        end

        # Early pruning: calculate remaining possible spins to see if reaching mz is feasible
        remaining_spins = N-pos+1
        max_possible_mz = current_mz + remaining_spins  # All +1 spins
        min_possible_mz = current_mz - remaining_spins  # All -1 spins

        # If reaching mz is impossible, prune this branch
        if mz < min_possible_mz || mz > max_possible_mz
            continue
        end

        # Otherwise, push possible spin assignments onto the stack
        for s in spin_vals
            new_state = copy(state)
            new_state[pos] = s
            push!(stack,(new_state,pos+1,current_mz+s))
        end
    end

    return basis_states
end

## Generates basis list for a particular mz block (2D Lattice) [optimized]
function gen_basis_mz_2D_op(Nx::Int64,Ny::Int64,mz::Int64)
    
    N = Nx*Ny
    spin_vals = [1,0,-1]
    basis_states::Vector{Vector{Int64}} = Vector{Vector{Int64}}()

    # Stack-based iterative version to avoid recursive function overhead
    initial_state::Vector{Int64} = zeros(Int64,N)
    stack = [(initial_state,1,0)]  # (state,position,current_mz)

    while !isempty(stack)
        (state,pos,current_mz) = pop!(stack)

        # Check if we've assigned all spins
        if pos > N
            if current_mz == mz
                push!(basis_states,copy(state))  # Copy to avoid mutation issues
            end
            continue
        end

        # Early pruning: calculate remaining possible spins to see if reaching mz is feasible
        remaining_spins = N-pos+1
        max_possible_mz = current_mz + remaining_spins  # All +1 spins
        min_possible_mz = current_mz - remaining_spins  # All -1 spins

        # If reaching mz is impossible, prune this branch
        if mz < min_possible_mz || mz > max_possible_mz
            continue
        end

        # Otherwise, push possible spin assignments onto the stack
        for s in spin_vals
            new_state = copy(state)
            new_state[pos] = s
            push!(stack,(new_state,pos+1,current_mz+s))
        end
    end

    return basis_states
end

## Generates full basis list for 1D Chain
function gen_basis_full_1D(N::Int64)
    
    # Create the basis as an array of arrays
    basis = [collect(b) for b in product(fill(-1:1,N)...)]
    # Ensure the Basis is a 1D array of tuples
    basis_list = collect(Iterators.flatten([basis]));
    
    return basis_list::Vector{Vector{Int64}}
end

## Generates full basis list for 2D Lattice
function gen_basis_full_2D(Nx::Int64,Ny::Int64)

    # Create the basis as an array of arrays
    basis = [collect(b) for b in product(fill(-1:1,Nx*Ny)...)]
    # Ensure the Basis is a 1D array of tuples
    basis_list = collect(Iterators.flatten([basis]));
    
    return basis_list::Vector{Vector{Int64}}
end