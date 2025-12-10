## Define Spin-1 XY 1D Hamiltonian
function Ham_XY_1D(J1::Float64,J3::Float64,D::Float64,h::Float64,N::Int64,basis_list::Vector{Vector{Int64}},
                 basis_new1::Vector{Int64},basis_new2::Vector{Int64},nearest_bond_list::Vector{Vector{Int64}},
                 third_bond_list::Vector{Vector{Int64}},use_sparse::Bool)
    
    H_size::Int64 = length(basis_list)
    
    if use_sparse == true
        data::Vector{Float64} = Vector{Float64}()        # Stores the matrix element values
        rows::Vector{Int64} = Vector{Int64}()            # Stores the row indices for sparse matrix
        cols::Vector{Int64} = Vector{Int64}()            # Stores the column indices for sparse matrix
    else
        H_dense::Matrix{Float64} = zeros(Float64,H_size,H_size) # Matrix for non-sparse case
    end
         
    for basis::Vector{Int64} in basis_list
        diag_energy::Float64 = h*sum(basis)+D*(sum(basis.^2))
        ind1::Int64 = find_state_index(basis,basis_list)
        
        for (bond1,bond2) in zip(nearest_bond_list,third_bond_list)
            
            #For NN
            s1_1 = basis[bond1[1]+1]     
            s2_1 = basis[bond1[2]+1]
            flag1,kk1 = flip(s1_1,s2_1)
            basis_new1 .= copy(basis)
            
            if flag1 != 0                 
                for j1 in kk1
                    basis_new1[bond1[1]+1] = j1[1]   
                    basis_new1[bond1[2]+1] = j1[2]
                    ind2 = find_state_index(basis_new1,basis_list)
                    if ind2 >= 0
                        if use_sparse != true
                            H_dense[ind1,ind2] += J1
                        else
                            push!(data,J1)
                            push!(rows,ind1)
                            push!(cols,ind2)
                        end
                    end
                end
            end            
                        
            #For second NN
            s1_2 = basis[bond2[1]+1]     
            s2_2 = basis[bond2[2]+1]
            flag2,kk2 = flip(s1_2,s2_2)   
            basis_new2 .= copy(basis)
            
            if flag2 != 0                
                for j2 in kk2
                    basis_new2[bond2[1]+1] = j2[1]   
                    basis_new2[bond2[2]+1] = j2[2]
                    ind3 = find_state_index(basis_new2,basis_list)
                    if ind3 >= 0
                        if use_sparse != true
                            H_dense[ind1,ind3] += J3
                        else
                            push!(data,J3)
                            push!(rows,ind1)
                            push!(cols,ind3)
                        end
                    end
                end
            end
        end
        
        if use_sparse != true
            H_dense[ind1,ind1] = diag_energy
        else
            push!(data,diag_energy)
            push!(rows,ind1)
            push!(cols,ind1)
        end
    end  
    
    if use_sparse == true
        H_sparse = sparse(rows,cols,data,H_size,H_size)
    end
    
    if use_sparse == true
        return (H_sparse,H_size)
    else
        return (H_dense,H_size)
    end
end 

## Define Spin-1 XY 1D Hamiltonian with next nearest neighbour interactions
function Ham_2nd_nbr_1D(N::Int64,basis_list::Vector{Vector{Int64}},basis_new::Vector{Int64},
                    sec_nbr_bond_list::Vector{Vector{Int64}},use_sparse::Bool)
    
    H_size::Int64 = length(basis_list)
    
    if use_sparse == true
        data::Vector{Float64} = Vector{Float64}()  # Stores the matrix element values
        rows::Vector{Int64} = Vector{Int64}()            # Stores the row indices for sparse matrix
        cols::Vector{Int64} = Vector{Int64}()            # Stores the column indices for sparse matrix
    else
        H_dense::Matrix{Float64} = zeros(Float64,H_size,H_size) # Matrix for non-sparse case
    end
        
    for basis::Vector{Int64} in basis_list
        ind1::Int64 = find_state_index(basis,basis_list)
        
        for bond in sec_nbr_bond_list
            s1 = basis[bond[1]+1]
            s2 = basis[bond[2]+1]
            
            flag,kk = flip(s1,s2)
            basis_new .= copy(basis)
            if flag != 0
                for j in kk
                    basis_new[bond[1]+1] = j[1]
                    basis_new[bond[2]+1] = j[2]
                    ind2 = find_state_index(basis_new,basis_list)
                    if ind2 >= 0
                        if use_sparse != true
                            H_dense[ind1,ind2] += 1.0
                        else
                            push!(data,1.0)
                            push!(rows,ind1)
                            push!(cols,ind2)
                        end
                    end
                end
            end
        end
    end
        
    if use_sparse == true
        H_sparse = sparse(rows,cols,data,H_size,H_size)
    end
    
    if use_sparse == true
        return (H_sparse,H_size)
    else
        return (H_dense,H_size)
    end
end 

## Define Spin-1 XY 2D Hamiltonian
function Ham_XY_2D(J1::Float64,h::Float64,Nx::Int64,Ny::Int64,basis_list::Vector{Vector{Int64}},
                    basis_new::Vector{Int64},nearest_bond_list::Vector{Vector{Int64}},use_sparse::Bool)
    
    H_size::Int64 = length(basis_list)
    
    if use_sparse == true
        data::Vector{Float64} = Vector{Float64}()  # Stores the matrix element values
        rows::Vector{Int64} = Vector{Int64}()            # Stores the row indices for sparse matrix
        cols::Vector{Int64} = Vector{Int64}()            # Stores the column indices for sparse matrix
    else
        H_dense::Matrix{Float64} = zeros(Float64,H_size,H_size) # Matrix for non-sparse case
    end
        
    for basis::Vector{Int64} in basis_list
        diag_energy::Float64 = h*sum(basis)
        ind1::Int64 = find_state_index(basis,basis_list)
        
        for bond in nearest_bond_list
            s1 = basis[bond[1]+1]
            s2 = basis[bond[2]+1]
            
            flag,kk = flip(s1,s2)
            basis_new .= copy(basis)
            if flag != 0
                for j in kk
                    basis_new[bond[1]+1] = j[1]
                    basis_new[bond[2]+1] = j[2]
                    ind2 = find_state_index(basis_new,basis_list)
                    if ind2 >= 0
                        if use_sparse != true
                            H_dense[ind1,ind2] += J1
                        else
                            push!(data,J1)
                            push!(rows,ind1)
                            push!(cols,ind2)
                        end
                    end
                end
            end
        end
        
        if use_sparse != true
            H_dense[ind1,ind1] = diag_energy
        else
            push!(data,diag_energy)
            push!(rows,ind1)
            push!(cols,ind1)
        end
    end  
    
    if use_sparse == true
        H_sparse = sparse(rows,cols,data,H_size,H_size)
    end
    
    if use_sparse == true
        return (H_sparse,H_size)
    else
        return (H_dense,H_size)
    end
end 

## Define Spin-1 XY 2D Hamiltonian with diagonal interactions
function Ham_diag_interac_2D(Nx::Int64,Ny::Int64,basis_list::Vector{Vector{Int64}},
                      basis_new::Vector{Int64},diagonal_bond_list::Vector{Vector{Int64}},use_sparse::Bool)
    
    H_size::Int64 = length(basis_list)
    
    if use_sparse == true
        data::Vector{Float64} = Vector{Float64}()  # Stores the matrix element values
        rows::Vector{Int64} = Vector{Int64}()            # Stores the row indices for sparse matrix
        cols::Vector{Int64} = Vector{Int64}()            # Stores the column indices for sparse matrix
    else
        H_dense::Matrix{Float64} = zeros(Float64,H_size,H_size) # Matrix for non-sparse case
    end
        
    for basis::Vector{Int64} in basis_list
        ind1::Int64 = find_state_index(basis,basis_list)
                
        for bond in diagonal_bond_list
            s1 = basis[bond[1]+1]
            s2 = basis[bond[2]+1]
            
            flag,kk = flip(s1,s2)
            basis_new .= copy(basis)
            if flag != 0
                for j in kk
                    basis_new[bond[1]+1] = j[1]
                    basis_new[bond[2]+1] = j[2]
                    ind2 = find_state_index(basis_new,basis_list)
                    if ind2 >= 0
                        if use_sparse != true
                            H_dense[ind1,ind2] += 1.0
                        else
                            push!(data,1.0)
                            push!(rows,ind1)
                            push!(cols,ind2)
                        end
                    end
                end
            end
        end
    end  
    
    if use_sparse == true
        H_sparse = sparse(rows,cols,data,H_size,H_size)
    end
    
    if use_sparse == true
        return (H_sparse,H_size)
    else
        return (H_dense,H_size)
    end
end

## Hp = sum_r S^z_r.S^z_(r+1)
function Hp_diag(basis_list::Vector{Vector{Int64}},nearest_bond_list::Vector{Vector{Int64}},use_sparse::Bool)
                     
    H_size::Int64 = length(basis_list)
    
    if use_sparse == true
        data::Vector{Float64} = Vector{Float64}()  # Stores the matrix element values
        rows::Vector{Int64} = Vector{Int64}()            # Stores the row indices for sparse matrix
        cols::Vector{Int64} = Vector{Int64}()            # Stores the column indices for sparse matrix
    else
        H_dense::Matrix{Float64} = zeros(Float64,H_size,H_size) # Matrix for non-sparse case
    end
        
    for (ind1::Int64,basis::Vector{Int64}) in enumerate(basis_list)
        mat_elem::Float64 = 0.0
        
        for bond1 in nearest_bond_list
            s1_1 = basis[bond1[1]+1]     
            s2_1 = basis[bond1[2]+1]
            mat_elem += s1_1*s2_1
        end
        if use_sparse != true
            H_dense[ind1,ind1] += mat_elem
        else
            push!(data,mat_elem)
            push!(rows,ind1)
            push!(cols,ind1)
        end
    end
        
    if use_sparse == true
        H_sparse = sparse(rows,cols,data,H_size,H_size)
    end
    
    if use_sparse == true
        return (H_sparse,H_size)
    else
        return (H_dense,H_size)
    end
end