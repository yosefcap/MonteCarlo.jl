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
function Hamiltonian_operator()

end

function build_states(N,L_x,L_y,num_spin,num_species)
    # returen all states in the Fock space,  each state is reshaped to a an array of (L_x,L_y,num_spin,num_species)
    states=zeros(N,L_x,L_y,num_spin,num_species);
    for c=0:N-1
        state = dec2bin(c,num_n)
        states[c+1,:,:,:] = reshape(state , L_x,L_y,num_spin,num_species)   
    end
    return states
end

function create(state ,x,y,spin,species)
    # creation operator at (x,y,spin,species)
    new_state = state
    ex = 0.0  # denotes if the state was annihilated by the creation operator. 
    if state[x,y,spin,species] == 0
        new_state[x,y,spin,species] = 1
        ex = 1
    end
        return new_state , ex
end

function annihilate(state ,x,y,spin,species)
    # annihilation operator at (x,y,spin,species)
    new_state = state
    ex = 0  # denotes if the state was annihilated by the annihilation operator. 
    if state[x,y,spin,species] == 1
        new_state[x,y,spin,species] = 0
        ex = 1
    end
        return new_state , ex
end

function hop(state_i,x_i,y_i,spin_i,species_i,x_j,y_j,spin_j,species_j)
    #hopping from  (x_i,y_i,spin_i,species_i) to  (x_j,y_j,spin_j,species_f)
    state_temp , ex = annihilate(state_i, x_i,y_i,spin_i,species_i)
    ex == 0 ? state_f = state_temp : (state_f , ex) = create(state_temp, x_j,y_j,spin_j,species_j)
    return state_f , ex   
end

function number(state,x,y,spin,species)
    #number operator at  (x,y,spin,species)
    ex = state[x,y,spin,species]
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

