
using NonteCarlo

betas = (1.0)#(1.0, 2.0, 5.0, 6.0, 7.0)
#lattice = SquareLattice(4)
mus = (-0.4 )#vcat(-2.0:0.5:-0.5, -0.1:0.1:1.1, 1.25, 1.5, 2.0)

dqmcs = []
    counter = 0
    N = length(mus) * length(betas)
    @time for beta in betas, mu in mus
        counter += 1
        print("\r[", lpad("$counter", 2), "/$N]")
        m = HubbardModelBi(2,1)#(l = lattice, t = 1.0, U = 4.0, μ = mu)
        dqmc = DQMC(
            m, beta = beta, delta_tau = 0.125, safe_mult = 8, 
            thermalization = 1000, sweeps = 1000, measure_rate = 1,
            recorder = Discarder()
        )
        dqmc[:occ] = occupation(dqmc, m)
        #dqmc[:PC] = pairing_correlation(dqmc, m, kernel = MonteCarlo.pc_kernel)
        run!(dqmc, verbose = false)

        # for simplicity we just keep the whole simulation around
        push!(dqmcs, dqmc)
    end
mean(mean(dqmcs[1][:occ]))