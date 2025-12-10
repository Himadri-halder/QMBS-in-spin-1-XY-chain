## Calculates non-degenerate perturbative energy upto given 'order'
function nondeg_pert_energy(H_p::Matrix{Float64},energy::Vector{Float64},eigenvector::Matrix{Float64},
                    state_index::Int64,state::Vector{Float64},state_energy::Float64,outvect_n::Vector{Float64},
                    outvect_m::Vector{Float64},outvect_k::Vector{Float64},order::String)
    
    #Index of the 'state' of which perturbed version is wanted, is denoted as 'n'
    order_int::Int64 = parse(Int64,order)
    E_n::Float64 = state_energy
    n::Int64 = state_index
    
    mul!(outvect_n,H_p,state)
    H_nn::Float64 = dot(state,outvect_n)
    E_corr_1::Float64 = H_nn
    
    E_corr_2::Float64 = E_corr_3::Float64 = 0.0
    
    if order_int > 1
        for m::Int64 in 1:size(eigenvector,2)
            if m != n
                E_m::Float64 = energy[m]
                psi_m::Vector{Float64} = eigenvector[:,m]
                mul!(outvect_m,H_p,psi_m)
                H_nm::Float64 = dot(state,outvect_m)                            
                H_mn::Float64 = conj(H_nm)
                E_corr_2 += abs2(H_mn)/(E_n-E_m)

                if order_int > 2
                    for k::Int64 in 1:size(eigenvector,2)
                        if k != n
                            E_k::Float64 = energy[k]
                            psi_k::Vector{Float64} = eigenvector[:,k]
                            mul!(outvect_k,H_p,psi_k)
                            H_mk::Float64 = dot(psi_m,outvect_k)
                            H_kn::Float64 = dot(psi_k,outvect_n)
                            E_corr_3 += (H_kn*H_mk*H_nm)/((E_n-E_m)*(E_n-E_k))
                        end
                    end
                
                    E_corr_3 -= H_nn*abs2(H_mn)/(E_n-E_m)^2
                end                   
            end
        end
    end
    
    E_corr_list::Vector{Float64} =  Vector{Float64}()
    if order_int == 1
        E_corr_list = [E_corr_1]
    elseif order_int == 2
        E_corr_list = [E_corr_1,E_corr_2]
    elseif order_int == 3
        E_corr_list= [E_corr_1,E_corr_2,E_corr_3]
    end
    
    return E_corr_list
end


function nondeg_4th_pert_energy(H_p::Matrix{Float64},energy::Vector{Float64},eigenvector::Matrix{Float64},
                    state_index::Int64,state::Vector{Float64},state_energy::Float64,outvect_n::Vector{Float64},
                    outvect_m::Vector{Float64},outvect_k::Vector{Float64},outvect_l::Vector{Float64})
    
    #Index of the 'state' of which perturbed version is wanted, is denoted as 'n'
    E_n::Float64 = state_energy
    n::Int64 = state_index
    
    mul!(outvect_n,H_p,state)
    H_nn::Float64 = dot(state,outvect_n)                    # H_nn

    E_corr_4::Float64 = 0.0
    for m::Int64 in 1:size(eigenvector,2)
        if m != n
            E_m::Float64 = energy[m]
            psi_m::Vector{Float64} = eigenvector[:,m]
            mul!(outvect_m,H_p,psi_m)
            H_nm::Float64 = dot(state,outvect_m)            # H_nm & H_mn                         
            H_mn::Float64 = conj(H_nm)
            E_corr_4 += ((H_nn^2)*abs2(H_mn))/((E_n-E_m)^3)
            
            for k::Int64 in 1:size(eigenvector,2)
                if k != n
                    E_k::Float64 = energy[k]
                    psi_k::Vector{Float64} = eigenvector[:,k]    
                    mul!(outvect_k,H_p,psi_k)
                    H_mk::Float64 = dot(psi_m,outvect_k)         # H_mk
                    H_kn::Float64 = dot(psi_k,outvect_n)         # H_kn
                    E_corr_4 -= (H_nn*H_kn*H_mk*H_nm)/((E_n-E_m)*(E_n-E_k)^2)
                    
                    for l::Int64 in 1:size(eigenvector,2)
                        if l != n
                            E_l::Float64 = energy[l]
                            psi_l::Vector{Float64} = eigenvector[:,l]
                            mul!(outvect_l,H_p,psi_l)
                            H_kl::Float64 = dot(psi_k,outvect_l)    # H_kl
                            H_ml::Float64 = dot(psi_m,outvect_l)    # H_ml
                            H_ln::Float64 = dot(psi_l,outvect_n)    # H_ln
                            E_corr_4 +=  (H_ln*H_kl*H_mk*H_nm)/((E_n-E_m)*(E_n-E_k)*(E_n-E_l))       
                            E_corr_4 -= (H_nn*H_ln*H_ml*H_nm)/((E_n-E_l)*(E_n-E_m)^2)
                        end
                    end
                end
            end           
        end
    end
    
    return E_corr_4
end


## Second order non-degenerate perturbative wave-function        
function nondeg_pert_wf_overlap(lambda::Float64,H_p::Matrix{Float64},energy::Vector{Float64},
                eigenvector::Matrix{Float64},state_index::Int64,state::Vector{Float64},state_energy::Float64,
                outvect_n::Vector{Float64},outvect_m::Vector{Float64},outvect_k::Vector{Float64},order::String)
    
    #Index of the 'state' of which perturbed version is wanted, is denoted as 'n'
    order_int::Int64 = parse(Int64,order)
    E_n::Float64 = state_energy
    n::Int64 = state_index
    
    psi_11::Float64 = psi_12::Float64 = psi_21::Float64 = psi_22::Float64 = 0.0
    mul!(outvect_n,H_p,state)
    H_nn::Float64 = dot(state,outvect_n)
    
    for m::Int64 in 1:size(eigenvector,2)
        if m != n
            psi_m::Vector{Float64} = eigenvector[:,m]
            mul!(outvect_n,H_p,state)
            H_mn::Float64 = dot(psi_m,outvect_n)
            psi_11 += abs2(H_mn)/(E_n-energy[m])^2
            
            psi_12_k::Float64 = psi_21_k::Float64 = 0.0
            if order_int > 1
                for k::Int64 in 1:size(eigenvector,2)
                    if k != n
                        psi_k::Vector{Float64} = eigenvector[:,k]
                        mul!(outvect_k,H_p,psi_k)
                        H_mk::Float64 = dot(psi_m,outvect_k)
                        H_kn::Float64 = dot(psi_k,outvect_n)
                        H_km::Float64 = conj(H_mk)
                        H_nk::Float64 = conj(H_kn)
                        psi_12_k += (H_kn*H_mk)/((E_n-energy[m])*(E_n-energy[k]))
                        psi_21_k += (H_nk*H_km)/((E_n-energy[m])*(E_n-energy[k]))
                    end
                end
            
                H_nm::Float64 = conj(H_mn)
                psi_12_tail::Float64 = ((H_nn*H_mn)/(E_n-energy[m])^2)
                psi_21_tail::Float64 = ((H_nn*H_nm)/(E_n-energy[m])^2)

                psi_12 += (H_nm/(E_n-energy[m]))*(psi_12_k-psi_12_tail)
                psi_21 += (psi_21_k-psi_21_tail)*(H_mn/(E_n-energy[m]))

                psi_22 += (psi_12_k-psi_12_tail)*(psi_21_k-psi_21_tail)
            end
        end
    end
    
    #|<psi|psi_pert>|^2/<psi_pert|psi_pert>
    fidelity_list::Vector{Float64} =  Vector{Float64}()
    fidelity_1::Float64 = 1/(1+((lambda^2)*psi_11))
    if order_int == 1
        fidelity_list = [fidelity_1]
    elseif order_int == 2
        fidelity_2::Float64 = 1/(1+((lambda^2)*psi_11)+((lambda^3)*psi_12)+((lambda^3)*psi_21)+      ((lambda^4)*psi_22))
        fidelity_list = [fidelity_1,fidelity_2]
    end
    
    return fidelity_list
end  


## Creates the degenerate perturbation matrix from eigenvectors
function degn_1st_pert_matrix(H_p::Matrix{Float64},eigenvector::Matrix{Float64},row_V::Vector{Float64},
                              col_V::Vector{Float64},outvect::Vector{Float64},H_size::Int64)
    
    H_deg::Matrix{Float64} = zeros(Float64,H_size,H_size)
    
    for row in 1:H_size
        row_V .= eigenvector[:,row]
        for col in 1:H_size
            col_V .= eigenvector[:,col]
            mul!(outvect,H_p,col_V)
            H_deg[row,col] = dot(row_V,outvect)
        end
    end
    
    return H_deg::Matrix{Float64}
end