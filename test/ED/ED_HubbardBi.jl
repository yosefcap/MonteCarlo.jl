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

state_p=zeros(N,L_x,L_y,num_spin,num_species);
for c=0:N-1
    state = dec2bin(c,num_n)
    state_p[c+1,:,:,:] = reshape(state , L_x,L_y,num_spin,num_species)   
end

function dec2bin(x::Int64,pad::Int64)

    bin = zeros(pad)
    for i in 1:pad
        bin[i]=x%2
        x = x÷2
    end
    return bin
end

function hop(state_i , x_i,y_i,x_j,y_j  , spin,species)
    #hopping from site x_i,y_i to site x_j,y_j
    if (state_i[x_i,y_i,spin,species]==1 && state_i[x_j,y_j,spin,species]==0 )
        new_state  = 1


end