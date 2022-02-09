using Base.Cartesian
using BitBasis
using Combinatorics
using SparseArrays
using Arpack
using LinearAlgebra
using Plots



function occupation_ed(specs,Nₙs,temp::Float64)

    N=length(Nₙs)
    occ=zeros(Float64,length(Nₙs[1]))
    z=0.0
    for i in 1:N
        Nₙ=Nₙs[i]
        bolzman=sum(exp.(-(1/temp).*specs[i]))
        occ.+=bolzman.*Nₙ
        z+=bolzman
    end
    return occ./z
end

function energy_ed(specs,Nₙs,temp::Float64)

    N=length(Nₙs) 
    enr=0.0
    z=0.0
    for i in 1:N   
        enr+=sum(exp.(-(1/temp).*specs[i]) .*specs[i] )
        z+=sum(exp.(-(1/temp).*specs[i])  )
    end
    return enr./z
end

function energy_obs(Dₙ::Vector{Float64} ,β::Float64)

    D = Diagonal(Dₙ)
    expD = Diagonal(exp.(-β*Dₙ))
    enr = tr(D*expD)    
    return enr
end

function number_obs(Dₙ::Vector{Float64},Nₙ ,β::Float64)

    D = Diagonal(Dₙ)
    expD = Diagonal(exp.(-β*Dₙ))
    N=sum(Nₙ)
    enr = tr(N*expD)    
    return enr
end

function partition_obs(Dₙ::Vector{Float64} ,β::Float64)

    expD = Diagonal(exp.(-β*Dₙ))
    par = tr(expD)    
    return par

end

function SC_corr_obs(spacial_dims::NTuple{DIMS2,Int64},index_i::CartesianIndex{DIMS},index_j::CartesianIndex{DIMS} , Nₙ , Dₙ::Vector{Float64} , Uₙ::Matrix{Float64} , β::Float64) where {DIMS} where {DIMS2}
    P = SC_corr_sub(spacial_dims,Nₙ, index_i , index_j) 
    Pd = Uₙ'*P*Uₙ
    expD = Diagonal(exp.(-β*Dₙ))
    par = tr(Pd*expD)    
    return par
end 

function SC_corr_sw_tau_obs(spacial_dims::NTuple{DIMS2,Int64},num_species::Int64,index_i::CartesianIndex{DIMS},index_j::CartesianIndex{DIMS} 
    , Nₙ , Dₙ::Vector{Float64} , Uₙ::Matrix{Float64} , Nₘ , Dₘ::Vector{Float64} , Uₘ::Matrix{Float64}
    , β::Float64,taus::Vector{Float64} ) where {DIMS} where {DIMS2}
    
    corr = zeros(Float64,length(taus))
    for c in 1:num_species
        for g in 1:2
            if g==1
                index_i_sp=CartesianIndex(index_i,c) #(index_i...,c)#
                index_j_sp=CartesianIndex(index_j,c) #(index_j...,c)#
            else
                index_i_sp=CartesianIndex(index_j,c) #(index_j...,c)#
                index_j_sp=CartesianIndex(index_i,c)
            end

            Δᵢ⁺ = SC_create_sub(spacial_dims,Nₙ,Nₘ, index_i_sp) 
            Δⱼ  = SC_annihilate_sub(spacial_dims,Nₘ,Nₙ, index_j_sp) 
            Δᵢ⁺d = Uₘ' *  Δᵢ⁺ * Uₙ
            Δⱼd =  Uₙ' *  Δⱼ  * Uₘ   
           
            for i in 1:length(taus)
                τ=taus[i]
                e⁻⁽ᵝ⁻ᵀ⁾ᴰⁿ = Diagonal(exp.(-(β-τ).*Dₙ))
                e⁻ᵀᴰᵐ =     Diagonal(exp.(-τ*Dₘ))
                corr[i] += tr(e⁻⁽ᵝ⁻ᵀ⁾ᴰⁿ * Δⱼd * e⁻ᵀᴰᵐ * Δᵢ⁺d)                     
            end

        end
    end
    return corr
end

function SC_corr_dw_tau_obs(spacial_dims::NTuple{DIMS2,Int64},num_species::Int64,rᵢ::CartesianIndex{DIMS},rⱼ::CartesianIndex{DIMS} 
    , Nₙ , Dₙ::Vector{Float64} , Uₙ::Matrix{Float64} , Nₘ , Dₘ::Vector{Float64} , Uₘ::Matrix{Float64}
    , β::Float64,taus::Vector{Float64} ) where {DIMS} where {DIMS2}
    
    dims = (spacial_dims... ,num_species)
    corr = zeros(Float64,length(taus))
    for c in 1:num_species
        for g in 1:2
            if g==1
                index_i_sp=CartesianIndex(rᵢ,c) #(index_i...,c)#
                index_j_sp=CartesianIndex(rⱼ,c) #(index_j...,c)#
            else
                index_i_sp=CartesianIndex(rⱼ,c) #for the hemition conjugate of ΔᵢΔⱼ^† 
                index_j_sp=CartesianIndex(rᵢ,c)
            end
            Δᵢ⁺=[]
            Δⱼ=[]
           # Δᵢ⁺d=[]
           # Δⱼd=[]
            η=0.0
            for dir in 1:length(spacial_dims)
                for lr in 1:2
                    index_i = index_i_sp 
                    index_iₓ = hop(index_i,dir,lr,dims)                    
                    index_j = index_j_sp
                    index_jₓ = hop(index_j,dir,lr,dims)
                    lr==1 ? η=0.25 : η=-0.25   
                    if  (dir==1 && lr==1) #!isempty(Δᵢ⁺)
                        Δᵢ⁺  = η.*SC_create_sub(spacial_dims,Nₙ,Nₘ, index_i,index_iₓ)# c^†_{index_i,↑} c^†_{index_iₓ,↓} 
                        Δᵢ⁺ += η.* SC_create_sub(spacial_dims,Nₙ,Nₘ, index_iₓ,index_i)# c^†_{index_iₓ,↑} c^†_{index_i,↓} 
                        Δⱼ   = η.*SC_annihilate_sub(spacial_dims,Nₘ,Nₙ, index_j,index_jₓ) # c_{index_j,↓} c^_{index_jₓ,↑}
                        Δⱼ  += η.*SC_annihilate_sub(spacial_dims,Nₘ,Nₙ, index_jₓ,index_j) # c_{index_jₓ,↓} c^_{index_j,↑}
                    else
                        Δᵢ⁺ += η.*SC_create_sub(spacial_dims,Nₙ,Nₘ, index_i,index_iₓ)# c^†_{index_i,↑} c^†_{index_iₓ,↓} 
                        Δᵢ⁺ += η.* SC_create_sub(spacial_dims,Nₙ,Nₘ, index_iₓ,index_i)# c^†_{index_iₓ,↑} c^†_{index_i,↓} 
                        Δⱼ  += η.*SC_annihilate_sub(spacial_dims,Nₘ,Nₙ, index_j,index_jₓ) # c_{index_j,↓} c^_{index_jₓ,↑}
                        Δⱼ  += η.*SC_annihilate_sub(spacial_dims,Nₘ,Nₙ, index_jₓ,index_j) # c_{index_jₓ,↓} c^_{index_j,↑}
                    end
                end
            end

            Δᵢ⁺d = Uₘ' *  Δᵢ⁺ * Uₙ
            Δⱼd =  Uₙ' *  Δⱼ  * Uₘ   
           
            for i in 1:length(taus)
                τ=taus[i]
                e⁻⁽ᵝ⁻ᵀ⁾ᴰⁿ = Diagonal(exp.(-(β-τ).*Dₙ))
                e⁻ᵀᴰᵐ =     Diagonal(exp.(-τ*Dₘ))
                corr[i] += tr(e⁻⁽ᵝ⁻ᵀ⁾ᴰⁿ * Δⱼd * e⁻ᵀᴰᵐ * Δᵢ⁺d)                     
            end

        end
    end
    return corr
end

function observables(spacial_dims::NTuple{DIMS,Int64},num_species::Int64, t::Float64, U::Float64,μ::Float64
    ,βs::Vector{Float64},taus::Vector{Float64}) where {DIMS}
   
    num_spin=2
    N_max=prod(spacial_dims)
    nt=0:N_max
    n=[]
    for i in 1:num_spin*num_species
        push!(n,nt)
    end
    
    Nₙs=collect(Iterators.product(n...))[:]
    b=length(βs)
    tt=length(taus)
    partitions = zeros(Float64,b)
    energys = zeros(Float64,b)
    numbers = zeros(Float64,b)
    SC_corrs=zeros(Float64,spacial_dims...,b)
    sw_tau=zeros(Float64,spacial_dims...,tt,b)
    dw_tau=zeros(Float64,spacial_dims...,tt,b)
    
    
    for Nₙ in Nₙs 
        for i in 1:b #runs over βs
            β=βs[i] 
            if !isempty(taus)
                Nₘ =  (Nₙ[1]+1,Nₙ[2],Nₙ[3]+1,Nₙ[4])
                Dₘ =[]
                Uₘ=[]
                if Nₘ in Nₙs
                    Hₘ=hamiltonian_sub(spacial_dims,Nₘ,t,U,μ)
                    Dₘ , Uₘ = eigen(Hₘ)
                end
            end
            Hₙ=hamiltonian_sub(spacial_dims,Nₙ,t,U,μ)
            Dₙ , Uₙ = eigen(Hₙ)
            partitions[i] += partition_obs(Dₙ ,β)
            energys[i]    += energy_obs(Dₙ ,β)
            numbers[i]    += number_obs(Dₙ , Nₙ, β)
            for   xt in 1:spacial_dims[1]
                for yt  in 1:spacial_dims[2]
                    x=1
                    y=1
                    SC_corrs[xt,yt,i]+=SC_corr_obs(spacial_dims,CartesianIndex(x,y) , CartesianIndex(xt,yt) , Nₙ , Dₙ , Uₙ, β)
                   # if !isempty(taus) && !isempty(Dₘ)
                    #    sw_tau[xt,yt,:,i]+=SC_corr_sw_tau_obs(spacial_dims,num_species,CartesianIndex(x,y),CartesianIndex(xt,yt) 
                     #   , Nₙ , Dₙ  , Uₙ , Nₘ , Dₘ  , Uₘ , β,taus ) 
                   #     dw_tau[xt,yt,:,i]+=SC_corr_dw_tau_obs(spacial_dims,num_species,CartesianIndex(x,y),CartesianIndex(xt,yt) 
                      #  , Nₙ , Dₙ  , Uₙ , Nₘ , Dₘ  , Uₘ , β,taus ) 
                   # end
                end
            end           
        end
    end

    for i in 1:b
        SC_corrs[:,:,i] = SC_corrs[:,:,i]./partitions[i]
        sw_tau[:,:,:,i] = sw_tau[:,:,:,i]./partitions[i]
        dw_tau[:,:,:,i] = dw_tau[:,:,:,i]./partitions[i]
    end
    return partitions, energys./partitions, numbers./partitions,SC_corrs, sw_tau, dw_tau
end


function spectrum(spacial_dims::NTuple{DIMS,Int64},num_species::Int64, t::Float64, U::Float64,μ::Float64) where {DIMS}
    num_spin=2
    N_max=prod(spacial_dims)
    nt=0:N_max
    n=[]
    for i in 1:num_spin*num_species
        push!(n,nt)
    end
    Dₙ=[]
    #Uₙ=[]
    Nₙs=collect(Iterators.product(n...))[:]
    for Nₙ in Nₙs
        ham=hamiltonian_sub(spacial_dims,Nₙ,t,U,μ)
        E=eigen(ham)
        push!(Dₙ,E.values)
        #push!(Uₙ,E.vectors)
    end
    return Dₙ  , Nₙs
end

function hamiltonian_sub(spacial_dims::NTuple{DIMS,Int64},Nₙ, t::Float64, U::Float64,μ::Float64) where {DIMS}

   states , state_num =  build_states(spacial_dims,Nₙ)
   N=length(states)
   H_sub=zeros(Float64,N,N)
   for ket in 1:N  
        state_ket=states[ket]
        state_sp , co = hamiltonian_operator(state_ket,t,U,μ)
        for j in 1:length(state_sp)
            state_bra=state_sp[j]
            bra = findall(x->x==packbits(state_bra[:]),state_num)
            H_sub[ket,bra...]+=co[j]
        end
    end
    return H_sub
end

function SC_corr_sub(spacial_dims::NTuple{DIMS2,Int64},Nₙ, index_i::CartesianIndex{DIMS} , index_j::CartesianIndex{DIMS}) where {DIMS} where {DIMS2} 

    states ,state_num =  build_states(spacial_dims,Nₙ)
    N=length(states)
   
    SC_sub=zeros(Float64,N,N)
    for ket in 1:N  
     state_ket=states[ket]
         state_sp , co = SC_corr_operator(state_ket,index_i,index_j)
         for j in 1:length(state_sp)
             state_bra=state_sp[j]
             bra = findall(x->x==packbits(state_bra[:]),state_num)
             SC_sub[ket,bra].+=co[j]
         end
     end
     return SC_sub
 end
 

function SC_corr_operator(state::Array{Int8,DIMS2} , index_i::CartesianIndex{DIMS}, index_j::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    # index_i contains only the spacial index 
    dims = size(state)
    spacial_dims = length(index_i)
    num_species = dims[spacial_dims+1]
    state_sp = Array{Int8,length(dims)}[]#state#
    co=Float64[]#0.0#

    for c in 1:num_species 
        for g in 1:2
            if g==1
                index_i_sp=CartesianIndex(index_i,c) #(index_i...,c)#
                index_j_sp=CartesianIndex(index_j,c) #(index_j...,c)#
            else
                index_i_sp=CartesianIndex(index_j,c) #(index_j...,c)#
                index_j_sp=CartesianIndex(index_i,c)
            end
            state_temp , ex1 =  SC_create_operatr(state,index_i_sp)
            ex2=0
            ex1 == 0 ? state_f = state_temp : (state_f , ex2) = SC_annihilate_operatr(state_temp, index_j_sp)
            if ex1*ex2 != 0.0
                push!(state_sp,state_f)
                push!(co,ex1*ex2) 
            end  
        end
    end
    return state_sp , co
end


function hamiltonian_operator(state::Array{Int8,DIMS2} ,  t::Float64,U::Float64,μ::Float64) where {DIMS2}
    
    dims = size(state)
    spacial_dims = length(dims)-2
    state_sp = Array{Int8,length(dims)}[]
    co=Float64[]

    diag_term = U*sum( sum(state[ntuple(k->:,spacial_dims+1)...,1].-state[ntuple(k->:,spacial_dims+1)...,2] ,dims=spacial_dims+1 ).^2 ) #interaction term
    diag_term = -diag_term - μ*sum(state) #chemical potential term
    push!(state_sp,state)
    push!(co,diag_term)

    #hopping term 
    occupations = findall(x->x==1, state) # indices which are  occupied 
    for occupation in occupations
        for dir in 1:spacial_dims
            for lr in 1:2
                index_i = occupation
                index_j = hop(index_i,dir,lr,dims)
                state_f , co_temp = hopping_operatr(state,index_i,index_j)
                if co_temp != 0
                    push!(state_sp,state_f)
                    push!(co,-t*co_temp)
                end
            end
        end
    end
    return state_sp , co
end

function create(state::Array{Int8,DIMS2} ,index::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    # creation operator at (x,y,species,spin)
    new_state = copy(state)
    ex = 0  # denotes if the state was annihilated by the creation operator. 
    if state[index] == 0
        new_state[index] = 1 
        ex = check_fermion_comm(state ,index)#1
    end
        return new_state , ex
end 

function annihilate(state::Array{Int8,DIMS2} ,index::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    # annihilation operator at (x,y,spin,species)
    new_state = copy(state)
    ex = 0  # denotes if the state was annihilated by the annihilation operator. 
    if state[index] == 1
        new_state[index] = 0
        ex = check_fermion_comm(state ,index)#1
    end
        return new_state , ex
end

function check_fermion_comm(state::Array{Int8,DIMS2} ,index::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    state_l=LinearIndices(state)
    index_l=state_l[index]
    sign = sum(state[:][index_l+1:end])
    ex = (-1)^(sign-1)#(-1)^(sign-1)
    return ex 
end

function hopping_operatr(state::Array{Int8,DIMS2} ,index_i::CartesianIndex{DIMS},index_j::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    #hopping from  (x_i,y_i,species_i,spin_i) to  (x_j,y_j,species_j,spin_j)
    state_temp , ex1 = annihilate(state,index_i)
    ex1 == 0 ? state_f = state_temp : (state_f , ex2) = create(state_temp, index_j)
    return state_f , ex1*ex2
       
end

function SC_create_operatr(state::Array{Int8,DIMS2} ,index_i::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    #sc operator c^†_{index_i,↑} c^†_{index_i,↓} where   index_i=(x_i,y_i,species_i) 
    up_s   = CartesianIndex(index_i,1) 
    down_s = CartesianIndex(index_i,2)

     state_temp , ex1 = create(state,down_s)
    ex2=0
    ex1 == 0 ? state_f = state_temp : (state_f , ex2) = create(state_temp, up_s)
    return state_f , ex1*ex2   
end

function SC_create_operatr(state::Array{Int8,DIMS2} ,index_i::CartesianIndex{DIMS},index_j::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    #sc operator c^†_{index_i,↑} c^†_{index_j,↓} where   index_i=(x_i,y_i,species_i) 
    up_s   = CartesianIndex(index_i,1) 
    down_s = CartesianIndex(index_j,2)

     state_temp , ex1 = create(state,down_s)
    ex2=0
    ex1 == 0 ? state_f = state_temp : (state_f , ex2) = create(state_temp, up_s)
    return state_f , ex1*ex2   
end

function SC_create_sub(spacial_dims::NTuple{DIMS2,Int64},Nₙ, Nₘ, index_i::CartesianIndex{DIMS} ) where {DIMS} where {DIMS2} 

    statesₙ ,state_numₙ =  build_states(spacial_dims,Nₙ)
    statesₘ ,state_numₘ =  build_states(spacial_dims,Nₘ)
    N=length(statesₙ)
    M=length(statesₘ)
    mat=zeros(Float64,M,N)#spzeros(Float64,N,N)
    for ket in 1:N  
        state_ket=statesₙ[ket]
        state_bra , co = SC_create_operatr(state_ket,index_i)
        bra = findall(x->x==packbits(state_bra[:]),state_numₘ)
        mat[bra,ket].+=co
     end
     return mat
 end

 function SC_create_sub(spacial_dims::NTuple{DIMS2,Int64},Nₙ, Nₘ, index_i::CartesianIndex{DIMS}, index_j::CartesianIndex{DIMS} ) where {DIMS} where {DIMS2} 
    # matrix representation of the operator c^†_{index_i,↑} c^†_{index_j,↓} where   index_i=(x_i,y_i,species_i)
    statesₙ ,state_numₙ =  build_states(spacial_dims,Nₙ)
    statesₘ ,state_numₘ =  build_states(spacial_dims,Nₘ)
    N=length(statesₙ)
    M=length(statesₘ)
    mat=zeros(Float64,M,N)
    for ket in 1:N  
        state_ket=statesₙ[ket]
        state_bra , co = SC_create_operatr(state_ket,index_i,index_j)
        bra = findall(x->x==packbits(state_bra[:]),state_numₘ)
        mat[bra,ket].+=co
     end
     return mat
 end


function SC_annihilate_operatr(state::Array{Int8,DIMS2} ,index_i::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    #sc operator c^†_{index_i,↓} c^_{index_i,↑} where   index_i=(x_i,y_i,species_i) 
    up_s   = CartesianIndex(index_i,1) 
    down_s = CartesianIndex(index_i,2)
    ex2=0
    state_temp , ex1 = annihilate(state,up_s)
    ex1 == 0 ? state_f = state_temp : (state_f , ex2) = annihilate(state_temp, down_s)
    return state_f , ex1*ex2   
end

function SC_annihilate_operatr(state::Array{Int8,DIMS2} ,index_i::CartesianIndex{DIMS},index_j::CartesianIndex{DIMS}) where {DIMS} where {DIMS2}
    #sc operator c_{index_i,↓} c^_{index_j,↑} where   index_i=(x_i,y_i,species_i) 
    up_s   = CartesianIndex(index_j,1) 
    down_s = CartesianIndex(index_i,2)
    ex2=0
    state_temp , ex1 = annihilate(state,up_s)
    ex1 == 0 ? state_f = state_temp : (state_f , ex2) = annihilate(state_temp, down_s)
    return state_f , ex1*ex2   
end

function SC_annihilate_sub(spacial_dims::NTuple{DIMS2,Int64},Nₙ, Nₘ, index_i::CartesianIndex{DIMS} ) where {DIMS} where {DIMS2} 

    statesₙ ,state_numₙ =  build_states(spacial_dims,Nₙ)
    statesₘ ,state_numₘ =  build_states(spacial_dims,Nₘ)
    N=length(statesₙ)
    M=length(statesₘ)
    mat=zeros(Float64,M,N)#spzeros(Float64,N,N)
    for ket in 1:N  
        state_ket=statesₙ[ket]
        state_bra , co = SC_annihilate_operatr(state_ket,index_i)
        bra = findall(x->x==packbits(state_bra[:]),state_numₘ)
        mat[bra,ket].+=co
     end
     return mat
 end

 function SC_annihilate_sub(spacial_dims::NTuple{DIMS2,Int64},Nₙ, Nₘ, index_i::CartesianIndex{DIMS}, index_j::CartesianIndex{DIMS} ) where {DIMS} where {DIMS2} 

    statesₙ ,state_numₙ =  build_states(spacial_dims,Nₙ)
    statesₘ ,state_numₘ =  build_states(spacial_dims,Nₘ)
    N=length(statesₙ)
    M=length(statesₘ)
    mat=zeros(Float64,M,N)#spzeros(Float64,N,N)
    for ket in 1:N  
        state_ket=statesₙ[ket]
        state_bra , co = SC_annihilate_operatr(state_ket,index_i,index_j)
        bra = findall(x->x==packbits(state_bra[:]),state_numₘ)
        mat[bra,ket].+=co
     end
     return mat
 end

function hop(index::CartesianIndex{DIMS},dir::Int64,lr::Int64,dims::NTuple{DIMS,Int64}) where {DIMS}
    # update index
    if (lr==1)
      hop_index= index[dir]==dims[dir] ? 1 : index[dir]+1
    else
      hop_index= index[dir]==1 ? dims[dir] : index[dir]-1
    end
    # generate a new CartesianIndex with updated index
    CartesianIndex(Base.setindex(Tuple(index), hop_index, dir))

end

function build_states(dims::NTuple{DIMS,Int64}) where {DIMS}#(N,L_x,L_y,num_species,num_spin)
    # returen all states in the Fock space,  each state is reshaped to a an array of dims=(L_x,L_y,num_spin,num_species)
    Nₙ_n=prod(dims)
    N=2^Nₙ_n
    states=[];
    for c in 0:N-1
        state = dec2bin(c,Nₙ_n)
        push!(states,reshape(state , dims...))
    end
    return states
end

function build_states(dims::NTuple{DIMS,Int64},Nₙ::NTuple{DIMS2,Int64}) where {DIMS} where {DIMS2}
    #Nₙ is a vector of  Nₙber operators for the different types of electrons (classified by their specie and spin)
    #dims is a vector of the spacial dimensions.
    # returen all states in the sub-space of the Fock space 
    Nₙ_sites=prod(dims) 
    states=[] # TO DO -set type to integer
    for Nsp in Nₙ
        st = [ones(Int8,Nsp);zeros(Int8,Nₙ_sites-Nsp)]
        bit_combs=unique(permutations(st) )
        Nₙ_comb = Int64[] 
        for bit_comb in bit_combs 
            push!(Nₙ_comb,packbits(bit_comb)  )
        end
        push!(states,Nₙ_comb)
    end
    states_comb = collect(Iterators.product(states... ))[:]
    states_bin=[]
    for sc in states_comb
        push!(states_bin,bin_states(sc,dims,Int8(length(Nₙ)/2)))
    end
    state_num=Int64[]
    N=length(states_bin)
    for c in 1:N#state in states
     state_c = states_bin[c]
         push!(state_num,packbits(state_c[:]))
    end
       return states_bin , state_num
end

function bin_states(state::NTuple{DIMS2,Int64},dims::NTuple{DIMS,Int64},Nspecies::Int8) where {DIMS} where {DIMS2}
    Nspin=2
     bs=zeros(Int8,dims...,Nspecies,Nspin)
     for ns in 1:Nspecies
         for np in 1:Nspin
            bs[ntuple(k->:,length(dims))...,ns,np]=reshape(bitarray(state[ns+(np-1)*Nspecies],prod(dims)),dims)
         end
    end
    return bs
end


