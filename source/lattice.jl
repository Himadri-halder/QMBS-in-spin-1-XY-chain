# Generates all nearest-neighbour bonds for 1D chain
function gen_bonds_1D(N::Int64,boundary::String)
    
    if N == 2 
        Jbondlist::Vector{Vector{Int64}} = [[0,1]]         # Special Case for N=2
    else
        Jbondlist = Vector{Vector{Int64}}(undef,0)  # Initialzing as empty array
        if boundary == "PBC"
            for i in 0:(N-1)
                xbond = [i,mod(i+1,N)]
                push!(Jbondlist,xbond)
            end
        elseif boundary == "OBC"
            for i in 0:(N-2)
                xbond = [i,i+1]
                push!(Jbondlist,xbond)
            end
        end
    end
    
    return Jbondlist
end

# Generates all 2nd nearest-neighbour bonds for 1D chain
function gen_2nd_nbr_bonds_1D(N::Int64,boundary::String)
    
    if N == 2
        error("Please put N>2 to generate 2nd nearest neighbours")
    else 
        bondlist::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,0)
        if boundary == "OBC"
            for site in 0:N-3
                xbond = [site,site+2]
                push!(bondlist,xbond)
            end
        elseif boundary == "PBC"
            for site in 0:N-1
                xbond = [site,mod(site+2,N)]
                push!(bondlist,xbond)
            end
        end
    end
    
    return bondlist
end

# Generates all 2nd nearest-neighbour bonds for 1D chain
function gen_3rd_nbr_bonds_1D(N::Int64,boundary::String)
    
    if N <= 3
        error("Please put N>3 to generate 3rd nearest neighbours")
    else 
        bondlist::Vector{Vector{Int64}} = Vector{Vector{Int64}}(undef,0)
        if boundary == "OBC"
            for site in 0:N-4
                xbond = [site,site+3]
                push!(bondlist,xbond)
            end
        elseif boundary == "PBC"
            for site in 0:N-1
                xbond = [site,mod(site+3,N)]
                push!(bondlist,xbond)
            end
        end
    end
    
    return bondlist
end

# Generates lattice sites and their (x,y) co-ordinates
function lattice_coords_2D(Nx::Int64,Ny::Int64)
    
    coord::OrderedDict{Int64,Tuple{Int64,Int64}} = OrderedDict{Int64,Tuple{Int64,Int64}}()
    for i::Int64 in 0:(Nx*Ny-1)
        x = mod(i,Nx)
        y = div(i,Nx)
        coord[i] = (x,y)
    end
    return coord
end

# Generates all nearest-neighbour bonds for 2D lattice
function gen_bonds_2D(Nx::Int64,Ny::Int64,boundary::String)
    
    N = Nx*Ny
    if Nx == 2 && Ny == 2
        Jbondlist::Vector{Vector{Int64}} = [[0,1],[0,2],[1,3],[2,3]]         # Special Case for N=4
    else
        Jbondlist = Vector{Vector{Int64}}(undef,0)  # Initialzing as empty array
        if boundary == "PBC"
            for i in 0:(N-1)
                x = div(i,Nx)            
                xbond = [i,Nx*x+mod(i+1,Nx)]
                ybond = [i,mod(i+Nx,N)]
                push!(Jbondlist,xbond)
                push!(Jbondlist,ybond)
            end
        elseif boundary == "OBC"
            for i in 0:(N-1)
                if mod(i+1,Nx) != 0
                    xbond = [i,i+1]
                    push!(Jbondlist,xbond)
                end
                if div(((i+Nx)%N),Nx) != 0
                    ybond = [i,i+Nx]
                    push!(Jbondlist,ybond)
                end
            end
        end
    end
    
    return Jbondlist
end

# Generates diagonal interaction for 2D lattice
function gen_diagonal_bonds_2D(Nx::Int64,Ny::Int64,lattice_coords::OrderedDict{Int64,Tuple{Int64,Int64}},
                            boundary::String)
    
    N = Nx*Ny
    
    if Nx == 2 && Ny == 2
        bondlist::Vector{Vector{Int64}} = [[0,3],[1,2]]
    else
        bondlist = Vector{Vector{Int64}}(undef,0)
        if boundary == "OBC"
            for site in 0:(N-1)
                coord = lattice_coords[site]
                x,y = coord[1],coord[2]
                x1,y1 = x+1,y+1              
                if 0 <= x1 <= Nx-1 && 0 <= y1 <= Ny-1
                    new_site1 = (y1*Ny)+x1
                    push!(bondlist,[site,new_site1])
                end
                x2,y2 = x+1,y-1
                if 0 <= x2 <= Nx-1 && 0 <= y2 <= Ny-1
                    new_site2 = (y2*Ny)+x2
                    push!(bondlist,[site,new_site2])
                end
            end

        elseif boundary == "PBC"
            for site in 0:(N-1)
                coord = lattice_coords[site]
                x,y = coord[1],coord[2]
                x1,y1 = mod((x+1),Nx),mod((y+1),Ny)
                x2,y2 = mod((x+1),Nx),mod((y-1),Ny)
                new_site1 = (y1*Ny)+x1
                new_site2 = (y2*Ny)+x2

                append!(bondlist,[[site,new_site1],[site,new_site2]])
            end
        end
    end

    return bondlist
end