"""
    HubbardModelBi(lattice; params...)
    HubbardModelBi(L, dims; params...)
    HubbardModelBi(params::Dict)
    HubbardModelBi(params::NamedTuple)
    HubbardModelBi(; params...)

Defines a repulsive (positive `U`) Hubbard model on a given (or derived) 
`lattice`. If a linear system size `L` and dimensionality `dims` is given, the
`lattice` will be a Cubic lattice of fitting size.

Additional parameters (keyword arguments) include:
* `l::AbstractLattice = lattice`: The lattice the model uses. The keyword 
argument takes precedence over the argument `lattice`.
* `U::Float64 = 1.0 > 0.0` is the absolute value of the Hubbard Interaction.
* `t::Float64 = 1.0` is the hopping strength.

Internally, a discrete Hubbard Stratonovich transformation (Hirsch 
transformation) is used in the spin/magnetic channel to enable DQMC. The 
resulting Hubbard Stratonovich fiels is real.
To reduce computational cost a specialized `BlockDiagonal` representation of the
greens matrix is used. 
"""



@with_kw_noshow struct HubbardModelBi{LT<:AbstractLattice} <: HubbardModel
    # user optional
    U::Float64 = 0.6
    @assert U >= 0. "U must be positive."
    t::Float64 = 1.0
    μ::Float64 = -0.4
    a::Float64 = 0.8 # coeficient of the \lambda_1 term.
    # mandatory (this or (L, dims))
    l::LT

    # non-user fields
    flv::Int = 2
    # to avoid allocations (TODO always real?)
    IG::Matrix{Float64}  = Matrix{Float64}(undef, 2length(l), 2)
    IGR::Matrix{Float64} = Matrix{Float64}(undef, 2length(l), 2)
    R::Matrix{Float64}   = Matrix{Float64}(undef, 2, 2)
    Δ::Diagonal{Float64, Vector{Float64}} = Diagonal(Vector{Float64}(undef, 2))
    RΔ::Matrix{Float64}  = Matrix{Float64}(undef, 2, 2)
end


HubbardModelBi(params::Dict{Symbol, T}) where T = HubbardModelBi(; params...)
HubbardModelBi(params::NamedTuple) = HubbardModelBi(; params...)
function HubbardModelBi(lattice::AbstractLattice; kwargs...)
    HubbardModelBi(l = lattice; kwargs...)
end
function HubbardModelBi(L, dims; kwargs...)
    l = choose_lattice(HubbardModelBi, dims, L)
    HubbardModelBi(l = l; kwargs...)
end

# cosmetics
import Base.summary
import Base.show
Base.summary(model::HubbardModelBi) = "repulsive Hubbard model - Bilayer" 
function Base.show(io::IO, model::HubbardModelBi)
    print(io, "repulsive Hubbard model bilayer, $(length(model.l)) sites")
end
Base.show(io::IO, m::MIME"text/plain", model::HubbardModelBi) = print(io, model)


# Convenience
@inline parameters(m::HubbardModelBi) = (N = length(m.l), t = m.t, U = m.U)

# optional optimization
hopping_matrix_type(::Type{DQMC}, ::HubbardModelBi) = BlockDiagonal{Float64, 2, Matrix{Float64}}
greens_matrix_type( ::Type{DQMC}, ::HubbardModelBi) = BlockDiagonal{Float64, 2, Matrix{Float64}}

const HubbardBiDistribution = (Int8(-1), Int8(1),Int8(-2), Int8(2))

@inline function Base.rand(::Type{DQMC}, m::HubbardModelBi, nslices::Int)
    rand(HubbardBiDistribution, length(m.l), nslices)
end


"""
    hopping_matrix(mc::DQMC, m::HubbardModelBi)

Calculates the hopping matrix \$T_{i, j}\$ where \$i, j\$ are
site indices.

# TODO
# Note that since we have a time reversal symmetry relating spin-up
# to spin-down we only consider one spin sector (one flavor) for the repulsive
# Hubbard model in the DQMC simulation.

This isn't a performance critical method as it is only used once before the
actual simulation.
"""
function hopping_matrix(mc::DQMC, model::HubbardModelBi)
    N = length(model.l)
    T = zeros(N, N)

    # Nearest neighbor hoppings
    @inbounds @views begin
        for (src, trg) in neighbors(model.l, Val(true))
            trg == -1 && continue
            T[trg, src] += -model.t
        end
    end

    return BlockDiagonal(T, copy(T))
    end 


"""
    interaction_matrix_exp!(
        mc::DQMC, model::HubbardModelBi, result, conf, slice, power
    )

Calculate the interaction matrix exponential `expV = exp(- power * delta_tau * V(slice))`
and store it in `result::Matrix`.

This is a performance critical method.
"""
@inline @bm function interaction_matrix_exp!(mc::DQMC, model::HubbardModelBi,
            result::Diagonal, conf::HubbardConf, slice::Int, power::Float64=1.)
    Δτ = mc.parameters.delta_tau
    a = [model.a 1.0-model.a] 
    b  = sqrt(model.a/(1-model.a))
    x  = exp(0.5Δτ * model.U)
    α = [acosh(x-(1/sqrt(2)/b)*(x^2-1)), acosh(x+(b/sqrt(2))*(x^2-1))/2] #devide by 2 because  conf[i, slice] of \lambda_2 was labeled with a 2.
    #lambda = acosh(exp(0.5 * model.U * dtau))
    N = length(lattice(model))
    
    @inbounds for i in 1:N
        result.diag[i] = a[abs(conf[i, slice])]*exp(sign(power) * α[abs(conf[i, slice])] * conf[i, slice])
    end
    @inbounds for i in 1:N
        result.diag[i+N] = a[abs(conf[i, slice])]*exp(sign(power) * α[abs(conf[i, slice])] * conf[i, slice])
    end
    nothing
    end 

# TODO
# both of these could be a little bit faster if we used the zeros imposed by
# BlockDiagonal (i.e. G[i, i+N] = 0)

@inline @inbounds @bm function propose_local(
        mc::DQMC, model::HubbardModelBi, i::Int, slice::Int, conf::HubbardConf
    )
    N = length(model.l)
    G = mc.stack.greens
    Δτ = mc.parameters.delta_tau
    Δ = model.Δ 
    R = model.R

    #α = acosh(exp(0.5Δτ * model.U))
    # not sure it is the write place to do this
    a  = [model.a , 1.0-model.a]
    b  = sqrt(model.a/(1-model.a))
    x  = exp(0.5Δτ * model.U)
    α = [acosh(x-(1/sqrt(2)/b)*(x^2-1)) , acosh(x+(b/sqrt(2))*(x^2-1))]
    
    old_conf = conf[i, slice]
    new_conf = old_conf==-1 ? rand((Int8(1),Int8(-2),Int8(2))) :
               old_conf==1  ? rand((Int8(-1),Int8(-2),Int8(2))) :
               old_conf==-2 ? rand((Int8(-1),Int8(1),Int8(2))) :
                                    rand((Int8(-1),Int8(1),Int8(-2)))
            
    ΔE_Boson = α[abs(new_conf)]*sign(new_conf) - α[abs(old_conf)]*sign(old_conf)
    co_ratio = a[abs(new_conf)] / a[abs(old_conf)]
    Δ[1, 1] = co_ratio*exp(ΔE_Boson) - 1.0
    Δ[2, 2] = co_ratio*exp(-ΔE_Boson) - 1.0

    # Unrolled R = I + Δ * (I - G)
    R[1, 1] = 1.0 + Δ[1, 1] * (1.0 - G[i, i])
    R[1, 2] = - Δ[1, 1] * G[i, i+N]
    R[2, 1] = - Δ[2, 2] * G[i+N, i]
    R[2, 2] = 1.0 + Δ[2, 2] * (1.0 - G[i+N, i+N])

    # Calculate det of 2x2 Matrix
    # det() vs unrolled: 206ns -> 2.28ns
    detratio =  (R[1, 1] * R[2, 2] - R[1, 2] * R[2, 1])^2

    # There is no bosonic part (exp(-ΔE_Boson)) to the partition function.
    # Therefore pass 0.0
    return detratio,  0.0, new_conf #,nothing
end

@inline @inbounds @bm function accept_local!(
        mc::DQMC, model::HubbardModelBi, i::Int, slice::Int, 
        conf::HubbardConf, detratio, args...
    )

    @bm "accept_local (init)" begin
        N = length(model.l)
        G = mc.stack.greens
        IG = model.IG
        IGR = model.IGR
        Δ = model.Δ
        R = model.R 
        RΔ = model.RΔ
    end

    # inverting R in-place, using that R is 2x2, i.e.:
    # M^-1 = [a b; c d]^-1 = 1/det(M) [d -b; -c a]
    @bm "accept_local (inversion)" begin
        detratio2 =  R[1, 1] * R[2, 2] - R[1, 2] * R[2, 1]
        inv_div = 1.0 / detratio2
        R[1, 2] = -R[1, 2] * inv_div
        R[2, 1] = -R[2, 1] * inv_div
        x = R[1, 1]
        R[1, 1] = R[2, 2] * inv_div
        R[2, 2] = x * inv_div
    end

    # Compute (I - G) R^-1 Δ
    @bm "accept_local (IG, R)" begin
        @inbounds for m in axes(IG, 1)
            IG[m, 1] = -G[m, i]
            IG[m, 2] = -G[m, i+N]
        end
        IG[i, 1] += 1.0
        IG[i+N, 2] += 1.0      
        vmul!(RΔ, R, Δ)
        vmul!(IGR, IG, RΔ)
    end

    # Slowest part, don't think there's much more we can do
    @bm "accept_local (finalize computation)" begin
        # G = G - IG * (R * Δ) * G[i:N:end, :]

        # @turbo for m in axes(G, 1), n in axes(G, 2)
        #     mc.s.greens_temp[m, n] = IGR[m, 1] * G[i, n] + IGR[m, 2] * G[i+N, n]
        # end

        # @turbo for m in axes(G, 1), n in axes(G, 2) 
        #     G[m, n] = G[m, n] - mc.s.greens_temp[m, n]
        # end

        # BlockDiagonal version
        G1 = G.blocks[1]
        G2 = G.blocks[2]
        temp1 = mc.stack.greens_temp.blocks[1]
        temp2 = mc.stack.greens_temp.blocks[2]

        @turbo for m in axes(G1, 1), n in axes(G1, 2)
            temp1[m, n] = IGR[m, 1] * G1[i, n]
        end
        @turbo for m in axes(G2, 1), n in axes(G2, 2)
            temp2[m, n] = IGR[m+N, 2] * G2[i, n]
        end

        @turbo for m in axes(G1, 1), n in axes(G1, 2)
            G1[m, n] = G1[m, n] - temp1[m, n]
        end
        @turbo for m in axes(G2, 1), n in axes(G2, 2)
            G2[m, n] = G2[m, n] - temp2[m, n]
        end

        # Always
        conf[i, slice] = args[2]
    end

    nothing
end


@inline function energy_boson(mc::DQMC, m::HubbardModelBi, hsfield = conf(mc))
    # There is no purely bosonic part in the partition function
    0.0
end

greens(mc::DQMC, model::HubbardModelBi) = greens(mc)

function save_model(
        file::JLDFile,
        m::HubbardModelBi,
        entryname::String="Model"
    )
    write(file, entryname * "/VERSION", 1)
    write(file, entryname * "/tag", "HubbardModelBi")

    write(file, entryname * "/U", m.U)
    write(file, entryname * "/t", m.t)
    save_lattice(file, m.l, entryname * "/l")
    write(file, entryname * "/flv", m.flv)

    nothing
end

#     load_parameters(data, ::Type{<: DQMCParameters})
#
# Loads a DQMCParameters object from a given `data` dictionary produced by
# `JLD.load(filename)`.
function _load(data, ::Val{:HubbardModelBi})
    if !(data["VERSION"] == 1)
        throw(ErrorException("Failed to load HubbardModelBi version $(data["VERSION"])"))
    end

    l = _load(data["l"], to_tag(data["l"]))
    HubbardModelBi(
        U = data["U"],
        t = data["t"],
        l = l,
        flv = data["flv"]
    )
end
to_tag(::Type{<: HubbardModelBi}) = Val(:HubbardModelBi)



################################################################################
### Measurement kernels
################################################################################



function intE_kernel(mc, model::HubbardModelBi, G::GreensMatrix)
    # up-down zero
    model.U * sum((diag(G.val.blocks[1]) .- 0.5) .* (diag(G.val.blocks[2]) .- 0.5))
end
