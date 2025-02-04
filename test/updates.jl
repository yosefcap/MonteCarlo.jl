@testset "Update Type Tree" begin
    @test MonteCarlo.AbstractLocalUpdate <: MonteCarlo.AbstractUpdate
    @test MonteCarlo.AbstractGlobalUpdate <: MonteCarlo.AbstractUpdate
    @test MonteCarlo.Adaptive <: MonteCarlo.AbstractUpdate

    @test LocalSweep <: MonteCarlo.AbstractLocalUpdate

    @test MonteCarlo.AbstractParallelUpdate <: MonteCarlo.AbstractGlobalUpdate
    @test GlobalShuffle <: MonteCarlo.AbstractGlobalUpdate
    @test GlobalFlip <: MonteCarlo.AbstractGlobalUpdate
    @test SpatialShuffle <: MonteCarlo.AbstractGlobalUpdate
    @test TemporalShuffle <: MonteCarlo.AbstractGlobalUpdate
    @test Denoise <: MonteCarlo.AbstractGlobalUpdate
    @test DenoiseFlip <: MonteCarlo.AbstractGlobalUpdate
    @test StaggeredDenoise <: MonteCarlo.AbstractGlobalUpdate

    @test ReplicaExchange <: MonteCarlo.AbstractParallelUpdate
    @test ReplicaPull <: MonteCarlo.AbstractParallelUpdate
end

struct BadUpdate <: MonteCarlo.AbstractLocalUpdate end
MonteCarlo.name(::BadUpdate) = "BadUpdate"
MonteCarlo.update(::BadUpdate, args...) = 0.0

struct GoodUpdate <: MonteCarlo.AbstractLocalUpdate end
MonteCarlo.name(::GoodUpdate) = "GoodUpdate"
MonteCarlo.update(::GoodUpdate, args...) = 1.0

@testset "Scheduler" begin
    model = HubbardModelAttractive(2,2)

    # These do not advance sweeps and should therefore error
    @test_throws ErrorException SimpleScheduler()
    @test_throws ErrorException SimpleScheduler(Adaptive())
    @test_throws ErrorException SimpleScheduler(GlobalFlip())

    @test_throws MethodError AdaptiveScheduler()
    @test_throws ErrorException AdaptiveScheduler(tuple(), tuple())
    @test_throws ErrorException AdaptiveScheduler((Adaptive(),), tuple())
    @test_throws ErrorException AdaptiveScheduler((GlobalFlip(),), tuple())

    # Check constructor
    schedulers = (
        SimpleScheduler(LocalSweep(), GlobalFlip(), LocalSweep(2)),
        AdaptiveScheduler((LocalSweep(), GlobalFlip(), LocalSweep(2)), tuple()),
    )
    for scheduler in schedulers
        Random.seed!(123)
        @test scheduler.sequence == MonteCarlo.AcceptanceStatistics.(
            (LocalSweep(), GlobalFlip(), LocalSweep(), LocalSweep())
        )
        @test scheduler.idx == 0

        mc = DQMC(model, beta=1.0, scheduler = scheduler)
        MonteCarlo.reverse_build_stack(mc, mc.stack)
        MonteCarlo.propagate(mc)

        @test mc.scheduler.idx == 0
        @test mc.last_sweep == 0

        MonteCarlo.update(mc.scheduler, mc, mc.model)
        @test mc.scheduler.sequence[1].total == 1
        @test mc.scheduler.idx == 1
        @test mc.last_sweep == 1
        
        MonteCarlo.update(mc.scheduler, mc, mc.model)
        @test mc.scheduler.sequence[2].total == 1
        @test mc.scheduler.idx == 2
        @test mc.last_sweep == 2

        MonteCarlo.update(mc.scheduler, mc, mc.model)
        @test mc.scheduler.sequence[3].total == 1
        @test mc.scheduler.idx == 3
        @test mc.last_sweep == 3
    end

    schedulers = (
        SimpleScheduler(GoodUpdate(), GoodUpdate(), BadUpdate()),
        AdaptiveScheduler((GoodUpdate(), GoodUpdate(), BadUpdate()), tuple())
    )
    for scheduler in schedulers
        mc = DQMC(model, beta=1.0, scheduler = scheduler)
        MonteCarlo.reverse_build_stack(mc, mc.stack)
        MonteCarlo.propagate(mc)

        for _ in 1:300
            MonteCarlo.update(mc.scheduler, mc, mc.model)
        end

        io = IOBuffer()
        MonteCarlo.show_statistics(io, mc.scheduler)
        @test String(take!(io)) == "Update statistics (since start):\n\tBadUpdate              0.0% accepted   (  0 / 100)\n\tGoodUpdate           100.0% accepted   (200 / 200)\n\t--------------------------------------------------\n\tTotal                 66.7% accepted   (200 / 300)\n"
    end
end

@testset "AdaptiveScheduler" begin
    scheduler = AdaptiveScheduler(
        (BadUpdate(), Adaptive()), 
        (GoodUpdate(), BadUpdate())
    )

    # check defaults
    @test scheduler.grace_period == 99
    @test scheduler.minimum_sampling_rate == 0.01
    @test scheduler.adaptive_rate == 9.0

    model = HubbardModelAttractive(2,2)
    mc = DQMC(model, beta=1.0, scheduler = scheduler)

    # Checks without adaptive corrections
    # should be save to run, 1 / 2^50 chance of already adapting
    for _ in 1:200
        MonteCarlo.update(scheduler, mc, model)
    end

    @test scheduler.sampling_rates == [0.5, 0.5, 1e-10]
    @test scheduler.sequence[1].accepted == 0.0
    @test scheduler.sequence[1].total == 100
    # no adaptive stuff yet, both updates should be picked at 50%
    @test scheduler.adaptive_pool[1].accepted == scheduler.adaptive_pool[1].total
    @test 30 <= scheduler.adaptive_pool[1].total <= 70
    @test scheduler.adaptive_pool[2].accepted == 0
    @test 30 <= scheduler.adaptive_pool[2].total <= 70

    # Checks with adaptive corrections
    for _ in 1:800
        MonteCarlo.update(scheduler, mc, model)
    end
    @test scheduler.sampling_rates[1] > 0.925
    @test scheduler.sampling_rates[2] < 0.075
    @test scheduler.sampling_rates[3] == 1e-10
    @test scheduler.sequence[1].accepted == 0.0
    @test scheduler.sequence[1].total == 500

    # With 0 acceptance rate we should have the hard limit for total:
    i_min = ceil(Int, log(
        scheduler.adaptive_rate / (scheduler.adaptive_rate+1), 
        scheduler.minimum_sampling_rate / 0.5
    )) + 99 + 1 # +1 to be save
    @test scheduler.adaptive_pool[1].accepted == scheduler.adaptive_pool[1].total
    @test scheduler.adaptive_pool[1].total > 500 - i_min
    @test scheduler.adaptive_pool[2].accepted == 0
    @test scheduler.adaptive_pool[2].total < i_min
end


using MonteCarlo: conf, current_slice, nslices

@testset "Global update" begin
    # Verify that split greens calculation matches the normal one
    pivot = Vector{Int}(undef, 16)
    tempv = Vector{Float64}(undef, 16)

    Ul = Matrix{Float64}(undef, 16, 16)
    Dl = Vector{Float64}(undef, 16)
    Tl = Matrix{Float64}(undef, 16, 16)
    
    Ur = Matrix{Float64}(undef, 16, 16)
    Dr = Vector{Float64}(undef, 16)
    Tr = Matrix{Float64}(undef, 16, 16)

    for i in 1:10
        Bl = rand(16, 16)
        copyto!(Tl, Bl)
        MonteCarlo.udt_AVX_pivot!(Ul, Dl, Tl, pivot, tempv)
        Br = rand(16, 16)
        copyto!(Tr, Br)
        MonteCarlo.udt_AVX_pivot!(Ur, Dr, Tr, pivot, tempv)
        G1 = Matrix{Float64}(undef, 16, 16)
        MonteCarlo.calculate_greens_AVX!(Ul, Dl, Tl, Ur, Dr, Tr, G1, pivot, tempv)

        copyto!(Tl, Bl)
        MonteCarlo.udt_AVX_pivot!(Ul, Dl, Tl)
        copyto!(Tr, Br)
        MonteCarlo.udt_AVX_pivot!(Ur, Dr, Tr)
        G2 = Matrix{Float64}(undef, 16, 16)
        MonteCarlo.calculate_inv_greens_udt(Ul, Dl, Tl, Ur, Dr, Tr, G2, pivot, tempv)
        MonteCarlo.finish_calculate_greens(Ul, Dl, Tl, Ur, Dr, Tr, G2, pivot, tempv)

        @test G1 ≈ G2
    end

    models = (HubbardModelAttractive(2,2,mu=0.5), HubbardModelRepulsive(2,2))
    for model in models
        @testset "$(typeof(model))" begin
            mc1 = DQMC(model, beta=2.0)
            mc2 = DQMC(model, beta=2.0)
            
            # Verify probabilities and greens in global updates
            # This is the backbone check for all global and parallel updates
            for _ in 1:10
                # re-initialize everything with a random conf
                mc1.conf .= rand(DQMC, model, nslices(mc1))
                mc2.conf .= mc1.conf
                for mc in (mc1, mc2)
                    MonteCarlo.init!(mc)
                    MonteCarlo.reverse_build_stack(mc, mc.stack)
                    MonteCarlo.propagate(mc)
                end

                # global update
                temp_conf = shuffle(deepcopy(conf(mc1)))
                detratio, ΔE_boson, passthrough = MonteCarlo.propose_global_from_conf(mc1, model, temp_conf)
                global_p = exp(- ΔE_boson) * detratio
                MonteCarlo.accept_global!(mc1, model, temp_conf, passthrough)

                # global update through successive local updates
                local_p = 1.0
                for t in 1:nslices(mc2)
                    for i in 1:length(lattice(mc2))
                        if mc2.conf[i, current_slice(mc2)] != temp_conf[i, current_slice(mc2)]
                            detratio, ΔE_boson, passthrough = MonteCarlo.propose_local(
                                mc2, model, i, current_slice(mc2), conf(mc2)
                            )
                            local_p *= real(exp(- ΔE_boson) * detratio)
                            MonteCarlo.accept_local!(
                                mc2, model, i, current_slice(mc2), conf(mc2), detratio, ΔE_boson, passthrough
                            )
                        end
                    end
                    MonteCarlo.propagate(mc2)
                end

                # Move to correct time slice
                for t in 1:nslices(mc2)
                    MonteCarlo.propagate(mc2)
                end

                # Verify
                @test local_p ≈ global_p
                @test mc1.conf == mc2.conf
                @test current_slice(mc1) == current_slice(mc2)
                @test mc1.stack.greens ≈ mc2.stack.greens
            end
        end
    end

    function setup()
        model = HubbardModelAttractive(8,2)
        mc = DQMC(model, beta=1.0)
        MonteCarlo.init!(mc)
        MonteCarlo.reverse_build_stack(mc, mc.stack)
        MonteCarlo.propagate(mc)
        c = deepcopy(conf(mc))
        return mc, model, c
    end

    # Check config adjustments for global updates
    mc, model, c = setup()
    MonteCarlo.update(GlobalFlip(), mc, model)
    @test mc.temp_conf == -c

    # make this unlikely to fail randomly via big conf size
    mc, model, c = setup()
    MonteCarlo.update(GlobalShuffle(), mc, model)
    @test sum(mc.temp_conf) == sum(c)
    @test mc.temp_conf != c

    mc, model, c = setup()
    u = SpatialShuffle()
    MonteCarlo.init!(mc, u)
    MonteCarlo.update(u, mc, model)
    @test mc.temp_conf == c[u.indices, :]

    mc, model, c = setup()
    u = TemporalShuffle()
    MonteCarlo.init!(mc, u)
    MonteCarlo.update(u, mc, model)
    @test mc.temp_conf == c[:, u.indices]

    #=
    Lattice with conf:

        7 8 9                    + - -
        -----                    -----
    3 | 1 2 3 | 1            + | - + + | -
    6 | 4 5 6 | 4            - | + + - | +
    9 | 7 8 9 | 7            - | + - - | +
        -----                    -----
        1 2 3                    - + +


    Denoise                   DenoiseFlip             StaggeredDenoise
    1 = - -> +                1 = - -> -              1 = - -> -
    2 = + -> +                2 = + -> -              2 = + -> +
    3 = + -> -                3 = + -> +              3 = + -> +
    4 = + -> +                4 = + -> -              4 = + -> +
    5 = + -> +                5 = + -> -              5 = + -> -
    6 = - -> +                6 = - -> -              6 = - -> +
    7 = + -> -                7 = + -> +              7 = + -> +
    8 = - -> +                8 = - -> -              8 = - -> +
    9 = - -> -                9 = - -> +              9 = - -> +
    to surrounding            to -surrounding         mult by +1 for even, -1 for odd sites
    =#
    spatial = [-1, +1, +1, +1, +1, -1, +1, -1, -1]
    _conf = [spatial[i] for i in 1:9, slice in 1:10]

    model = HubbardModelAttractive(3, 2)
    mc = DQMC(model, beta=1.0)
    MonteCarlo.init!(mc)
    mc.conf .= _conf
    MonteCarlo.reverse_build_stack(mc, mc.stack)
    MonteCarlo.propagate(mc)
    MonteCarlo.update(Denoise(), mc, model)
    @test mc.temp_conf == [[+1, +1, -1, +1, +1, +1, -1, +1, -1][i] for i in 1:9, slice in 1:10]

    model = HubbardModelAttractive(3, 2)
    mc = DQMC(model, beta=1.0)
    MonteCarlo.init!(mc)
    mc.conf .= _conf
    MonteCarlo.reverse_build_stack(mc, mc.stack)
    MonteCarlo.propagate(mc)
    MonteCarlo.update(DenoiseFlip(), mc, model)
    @test mc.temp_conf == [[-1, -1, +1, -1, -1, -1, +1, -1, +1][i] for i in 1:9, slice in 1:10]

    model = HubbardModelAttractive(3, 2)
    mc = DQMC(model, beta=1.0)
    MonteCarlo.init!(mc)
    mc.conf .= _conf
    MonteCarlo.reverse_build_stack(mc, mc.stack)
    MonteCarlo.propagate(mc)
    MonteCarlo.update(StaggeredDenoise(), mc, model)
    @test mc.temp_conf == [[-1, +1, +1, +1, -1, +1, +1, +1, +1][i] for i in 1:9, slice in 1:10]
end

using Distributed
@testset "Parallel" begin
    @test isempty(MonteCarlo.connected_ids)
    
    # backend
    MonteCarlo.add_worker!(myid())
    @test MonteCarlo.connected_ids == [myid()]
    MonteCarlo.add_worker!(myid())
    @test MonteCarlo.connected_ids == [myid()]
    MonteCarlo.add_worker!(-1)
    @test MonteCarlo.connected_ids == [myid(), -1]
    MonteCarlo.remove_worker!(-1)
    @test MonteCarlo.connected_ids == [myid()]
    MonteCarlo.remove_worker!(myid())
    @test MonteCarlo.connected_ids == Int[]

    connect(myid()) 
    yield()
    @test MonteCarlo.connected_ids == [myid()]
    disconnect(myid())
    yield()
    @test MonteCarlo.connected_ids == Int[]
    connect([myid()]) 
    yield()
    @test MonteCarlo.connected_ids == [myid()]
    disconnect([myid()])
    yield()
    @test MonteCarlo.connected_ids == Int[]

    @test !isready(MonteCarlo.weight_probability)
    MonteCarlo.put_weight_prob!(myid(), 0.9, 0.6)
    yield()
    @test isready(MonteCarlo.weight_probability)
    @test take!(MonteCarlo.weight_probability) == (myid(), 0.9, 0.6)
end

@testset "Utilities/Printing" begin
    @test MonteCarlo.name(NoUpdate()) == "NoUpdate"
    @test MonteCarlo.name(LocalSweep()) == "LocalSweep"
    @test MonteCarlo.name(GlobalFlip()) == "GlobalFlip"
    @test MonteCarlo.name(GlobalShuffle()) == "GlobalShuffle"
    @test MonteCarlo.name(ReplicaExchange(1)) == "ReplicaExchange"
    @test MonteCarlo.name(ReplicaPull()) == "ReplicaPull"
    @test MonteCarlo.name(Adaptive()) == "Adaptive"

    # Idk if I want to keep these methods...
    model = HubbardModelAttractive(2,2) 
    mc = DQMC(model, beta=1.0)
    for T in (LocalSweep, GlobalFlip, GlobalShuffle, NoUpdate, ReplicaPull)
        @test T(mc, model) == T()
    end
    @test ReplicaExchange(mc, model, 1) == ReplicaExchange(1)
end