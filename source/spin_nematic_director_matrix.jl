## (S_r^+)^2 Matrix
function S_plus_sq(r::Int64,basis_list::Vector{Vector{Int64}},basis_new::Vector{Int64})
    
    H_size = length(basis_list)
    data,rows,cols = Vector{Int64}(),Vector{Int64}(),Vector{Int64}()
    
    for basis in basis_list
        ind1 = find_state_index(basis,basis_list)        
        sr = basis[r]
        basis_new .= copy(basis)
        
        if sr == -1
            basis_new[r] = 1
            ind2 = find_state_index(basis_new,basis_list)
            if ind2 >= 0
                push!(data,2)
                push!(rows,ind1)
                push!(cols,ind2) 
            end
        end
    end
        
    H_sparse = sparse(rows,cols,data,H_size,H_size)
    
    return (transpose(H_sparse),H_size)
end

## (S_r^-)^2 Matrix
function S_minus_sq(r::Int64,basis_list::Vector{Vector{Int64}},basis_new::Vector{Int64})
    
    H_size = length(basis_list)
    data,rows,cols = Vector{Int64}(),Vector{Int64}(),Vector{Int64}()
    
    for basis in basis_list
        ind1 = find_state_index(basis,basis_list)        
        sr = basis[r]
        basis_new .= copy(basis)
        
        if sr == 1
            basis_new[r] = -1
            ind2 = find_state_index(basis_new,basis_list)
            if ind2 >= 0
                push!(data,2)
                push!(rows,ind1)
                push!(cols,ind2) 
            end
        end
    end
        
    H_sparse = sparse(rows,cols,data,H_size,H_size)
    
    return (transpose(H_sparse),H_size)
end

## J^+ Operator Matrix
function J_plus_matrix(N::Int64,basis_list::Vector{Vector{Int64}},basis_new::Vector{Int64})
    
    H_size = length(basis_list)
    data,rows,cols = Vector{Int64}(),Vector{Int64}(),Vector{Int64}()
    
    for basis in basis_list
        ind1 = find_state_index(basis,basis_list)        
        
        for r in 1:N
            basis_new .= copy(basis)
            sr = basis[r]
            if sr == -1
                basis_new[r] = 1
#                 println(basis," ",r," ",basis_new)
                ind2 = find_state_index(basis_new,basis_list)
                if ind2 >= 0
                    push!(data,(-1)^r)
                    push!(rows,ind1)
                    push!(cols,ind2) 
                end
            end
        end
    end
    
    H_sparse = sparse(rows,cols,data,H_size,H_size)
    
    return (transpose(H_sparse),H_size)
end

## J^- Operator Matrix
function J_minus_matrix(N::Int64,basis_list::Vector{Vector{Int64}},basis_new::Vector{Int64})
    
    H_size = length(basis_list)
    data,rows,cols = Vector{Int64}(),Vector{Int64}(),Vector{Int64}()
    
    for basis in basis_list
        ind1 = find_state_index(basis,basis_list)        
        
        for r in 1:N
            basis_new .= copy(basis)
            sr = basis[r]
            if sr == 1
                basis_new[r] = -1
#                 println(basis," ",r," ",basis_new)
                ind2 = find_state_index(basis_new,basis_list)
                if ind2 >= 0
                    push!(data,(-1)^r)
                    push!(rows,ind1)
                    push!(cols,ind2) 
                end
            end
        end
    end
    
    H_sparse = sparse(rows,cols,data,H_size,H_size)
    
    return (transpose(H_sparse),H_size)
end

## J^-J^+ Operator Matrix
function odlro_matrix(N::Int64,basis_list::Vector{Vector{Int64}},basis_new1::Vector{Int64},
                         basis_new2::Vector{Int64})
    
    H_size = length(basis_list)
    data,rows,cols = Vector{Int64}(),Vector{Int64}(),Vector{Int64}()
    
    for basis in basis_list
        ind1 = find_state_index(basis,basis_list)  
        
        for r1 in 1:N
            sr1 = basis[r1]
            if sr1 == -1
                basis_new1 .= copy(basis)
                basis_new1[r1] = 1
                
                for r2 in 1:N
                    sr2 = basis_new1[r2]
                    if sr2 == 1
                        basis_new2 .= copy(basis_new1)
                        basis_new2[r2] = -1
                        ind2 = find_state_index(basis_new2,basis_list)
                        if ind2 >= 0
                            push!(data,(-1)^(r1+r2))
                            push!(rows,ind1)
                            push!(cols,ind2)
                        end
                    end
                end
            end
        end
    end
    
    H_sparse = sparse(rows,cols,data,H_size,H_size)
    
    return (transpose(H_sparse),H_size)
end