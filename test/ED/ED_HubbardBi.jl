using Base.Cartesian
using BitBasis

spacial_dims=Int8(2)
num_spin=2;
num_species=2;
L_x=2;
L_y=1;
num_sites = L_x*L_y
num_n=num_spin*num_species*num_sites;
N=2^num_n;
t=1.0; 
U=1.0;
mu=1.0;
## TO DO - assign types to all variables and functions

function hamiltonian_sub(spacial_dims::NTuple{DIMS,Int64},Num::Int64, t::Float64, U::Float64,μ::Float64) where {DIMS}

   states =  build_states(spacial_dims,Num)
   N=length(states)
   state_num=Int[]
   for state in states
        push(state_num,packbits(state[:]))
   end
   H_sub=zeros(Float64,N,N)
   for state in states
        H_temp=[]


   end
end
function build_states(dims::NTuple{DIMS,Int64}) where {DIMS}#(N,L_x,L_y,num_species,num_spin)
    # returen all states in the Fock space,  each state is reshaped to a an array of dims=(L_x,L_y,num_spin,num_species)
    num_n=prod(dims)
    N=2^num_n
    states=[];
    for c in 0:N-1
        state = dec2bin(c,num_n)
        push!(states,reshape(state , dims...))
    end
    return states
end

function build_states(dims::NTuple{DIMS,Int64},Num) where {DIMS}
    #Num is a vector of  number operators for the different types of electrons (classified by their specie and spin)
    #dims is a vector of the spacial dimensions.
    # returen all states in the sub-space of the Fock space 
    num_sites=prod(dims) 
    states=[] # TO DO -set type to integer
    for Nsp in Num
        st = [ones(Int8,Nsp);zeros(Int8,num_sites-Nsp)]
        bit_combs=unique(permutations(st) )
        num_comb = [] # TO DO -set type to integer
        for bit_comb in bit_combs 
            push!(num_comb,packbits(bit_comb)  )
        end
        push!(states,num_comb)
    end
    states_comb = collect(Iterators.product(states... ))[:]
    states_bin=[]
    for sc in states_comb
        push!(states_bin,bin_states(sc,dims,Int64(length(Num)/2)))
    end
    return states_bin
end

function bin_states(state,dims::NTuple{DIMS,Int64},Nspecies) where {DIMS}
    Nspin=2
     bs=zeros(Int8,dims...,Nspecies,Nspin)
     for ns in 1:Nspecies
         for np in 1:Nspin
            bs[ntuple(k->:,length(dims))...,ns,np]=reshape(bitarray(state[ns+(np-1)*Nspecies],prod(dims)),dims)
         end
    end
    return bs
end



function T_operator(state,  t::Float64)
    #hopping part of hamiltonian
    dims = size(state)
    spacial_dims=length(dims)-2
    state_sp=Array{Int8,length(dims)}[]#state
    co=Float64[]#0.0
    occupations = findall(x->x==1, state) # indices which are  occupied 

    for occupation in occupations
        for dir in 1:spacial_dims
            for lr in 1:2
                index_i = occupation
                index_j = hop(index_i,dir,lr,dims)
                state_f , co_temp = hopping_operatr(state,index_i,index_j)
                if co_temp != 0
                    push!(state_sp,state_f)
                    push!(co,-t)
                end
            end
        end
    end
    return state_sp , co
end

function U_operator(state, U::Float64)
    spacial_dims = ndims(state)-2
    num = U*sum( sum(state[ntuple(k->:,spacial_dims+1),1].-state[ntuple(k->:,spacial_dims+1),2] ,dims=spacial_dims+1 ).^2 )
    return num
end

function create(state ,index::CartesianIndex{DIMS}) where {DIMS}
    # creation operator at (x,y,spin,species)
    new_state = state
    ex = 0.0  # denotes if the state was annihilated by the creation operator. 
    if state[index] == 0
        new_state[index] = 1 
        ex = 1
    end
        return new_state , ex
end

function annihilate(state ,index::CartesianIndex{DIMS}) where {DIMS}
    # annihilation operator at (x,y,spin,species)
    new_state = state
    ex = 0  # denotes if the state was annihilated by the annihilation operator. 
    if state[index] == 1
        new_state[index] = 0
        ex = 1
    end
        return new_state , ex
end

function hopping_operatr(state_i,index_i::CartesianIndex{DIMS},index_j::CartesianIndex{DIMS}) where {DIMS}
    #hopping from  (x_i,y_i,spin_i,species_i) to  (x_j,y_j,spin_j,species_f)
    state_temp , ex = annihilate(state_i,index_i)
    ex == 0 ? state_f = state_temp : (state_f , ex) = create(state_temp, index_j)
    return state_f , ex   
end

function number(state,index::CartesianIndex{DIMS}) where {DIMS}
    #number operator at  (x,y,spin,species)
    ex = state[index]
    return state , ex # passing and returning state is unnecessary - perhaps remove
end

function dec2bin(x::Int64,pad::Int64)

    bin = zeros(Int8,pad)
    for i in 1:pad
        bin[i]=x%2
        x = x÷2
    end
    return bin
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


