# Saving and loading happens in a nested fashion:
# save(filename, mc) calls
#   save_mc(filename, mc) calls
#       save_model(filename, mc.model) calls
#           save_lattice(filename, mc.lattice) (TODO, not used yet)
#       save_measurements(filename, mc) calls
#           save_measurement(filename, measurement)
#
# loading follows the same structure
# > Each level (beyond the outermost save) also has an argument `entryname::String`.
# > Each save_something() should write a "VERSION" and a "type" to the file. The
#   former simplifies updates down the line, the latter allows for dispatch.
# > Each load_something should be dispatched on a type, e.g.
#   `load_model(data, ::Type(HubbardModelAttractive))`


"""
    save(filename, mc; force_overwrite=false, allow_rename=true)

Saves the given MonteCarlo simulation `mc` to a JLD-file `filename`.

If `allow_rename = true` the filename will be adjusted if it already exists. If
`force_overwrite = true` it will be overwritten. In this case a temporary backup
will be created. If neither are true an error will be thrown.
"""
function save(filename, mc::MonteCarloFlavor; force_overwrite=false, allow_rename=true)
    endswith(filename, ".jld") || (filename *= ".jld")

    # handle ranming and overwriting
    isfile(filename) && !force_overwrite && !allow_rename && throw(ErrorException(
        "Cannot save because \"$filename\" already exists. Consider setting " *
        "`allow_rename = true` to adjust the filename or `force_overwrite = true`" *
        " to overwrite the file."
    ))
    if isfile(filename) && !force_overwrite && allow_rename
        filename = _generate_unqiue_JLD_filename(filename)
    end

    if force_overwrite
        parts = splitpath(filename)
        parts[end] = "." * parts[end]
        temp_filename = _generate_unqiue_JLD_filename(joinpath(parts...))
        mv(filename, temp_filename)
    end

    mode = isfile(filename) ? "r+" : "w"
    file = jldopen(filename, mode)
    write(file, "VERSION", 1)
    save_mc(file, mc, "MC")
    close(file)

    if force_overwrite
        rm(temp_filename)
    end

    return filename
end

# Something like
# existing_file.jld -> existing_file_aJ3c.jld
function _generate_unqiue_JLD_filename(filename)
    isfile(filename) || return filename
    # those map to 0-9,    A-Z,        a-z
    x = rand("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
    s = "_$(Char(x))"
    filename = filename[1:end-4] * s * ".jld"
    while isfile(filename)
        x = rand("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz")
        s = string(Char(x))
        filename = filename[1:end-4] * s * ".jld"
    end
    filename
end

"""
    load(filename)

Loads a MonteCarlo simulation from the given JLD-file `filename`.
"""
function load(filename)
    data = JLD.load(filename)
    if !(data["VERSION"] == 1)
        throw(ErrorException("Failed to load $filename version $(data["VERSION"])"))
    end
    load_mc(data)
end


"""
    resume!(filename[; kwargs...])

Resumes a Monte Carlo simulation from a savefile generated by `run!(mc)`. Takes
the same keyword arguments as `run!`. Returns the simulation and the state
returned by `run!`.

See also: [`run!`](@ref)
"""
function resume!(filename; kwargs...)
    data = JLD.load(filename)

    if !(data["VERSION"] == 1)
        throw(ErrorException("Failed to load $filename version $(data["VERSION"])"))
    end
    @assert haskey(data, "last_sweep")

    mc = load_mc(data)
    last_sweep = data["last_sweep"]
    resume_init!(mc)
    load_rng!(data)

    state = run!(mc, start = last_sweep + 1; kwargs...)
    mc, state
end


function save_mc(filename::String, mc::MonteCarloFlavor, entryname::String="MC")
    mode = isfile(filename) ? "r+" : "w"
    file = jldopen(filename, mode)
    save_mc(file, mc, entryname)
    close(file)
    nothing
end
load_mc(data) = load_mc(data["MC"], data["MC"]["type"])



#     save_model(filename, model, entryname)
#
# Save (minimal) information necessary to reconstruct the given `model` in a
# jld-file `filename` under group `entryname`.
#
# By default the full model object is saved. When saving a simulation, the
# entryname defaults to `MC/Model`.
function save_model(filename::String, model, entryname::String)
    mode = isfile(filename) ? "r+" : "w"
    file = jldopen(filename, mode)
    save_model(file, model, entryname)
    close(file)
    nothing
end
function save_model(file::JLD.JldFile, model, entryname::String)
    write(file, entryname * "/VERSION", 0)
    write(file, entryname * "/type", typeof(model))
    write(file, entryname * "/data", model)
    nothing
end

#     load_model(data, ::Type{Model})
#
# Loads a model from a given `data` dictionary produced by `JLD.load(filename)`.
# The second argument can be used for dispatch between different models.
function load_model(data, ::DataType)
    @assert data["VERSION"] == 0
    data["data"]
end


# TODO: Not used currently
#     save_lattice(filename, lattice, entryname)
#
# Save (minimal) information necessary to reconstruct the given `lattice` in a
# jld-file `filename` under group `entryname`.
#
# By default the full lattice object is saved. When saving a simulation, the
# entryname defaults to `MC/Model/Lattice`.
function save_lattice(filename::String, lattice::AbstractLattice, entryname::String)
    mode = isfile(filename) ? "r+" : "w"
    file = jldopen(filename, mode)
    save_lattice(file, lattice, entryname)
    close(file)
    nothing
end
function save_lattice(file::JLD.JldFile, lattice::AbstractLattice, entryname::String)
    write(file, entryname * "/VERSION", 0)
    write(file, entryname * "/type", typeof(lattice))
    write(file, entryname * "/data", lattice)
    nothing
end

#     load_lattice(data, ::Type{Lattice})
#
# Loads a lattice from a given `data` dictionary produced by `JLD.load(filename)`.
function load_lattice(data, ::DataType)
    @assert data["VERSION"] == 0
    data["data"]
end


const _GLOBAL_RNG = VERSION < v"1.3.0" ? Random.GLOBAL_RNG : Random.default_rng()

"""
    save_rng(filename [; rng = _GLOBAL_RNG, entryname = "RNG"])

Saves the current state of Julia's random generator (`Random.GLOBAL_RNG`) to the
given `filename`.
"""
function save_rng(
        filename::String;
        rng::MersenneTwister = _GLOBAL_RNG,
        entryname::String="RNG"
    )
    mode = isfile(filename) ? "r+" : "w"
    file = jldopen(filename, mode)
    save_rng(file, rng=rng, entryname=entryname)
    close(file)
end
function save_rng(
        file::JLD.JldFile;
        rng::MersenneTwister = _GLOBAL_RNG,
        entryname::String="RNG"
    )
    try
        write(file, entryname, rng)
    catch e
        error("Error while saving RNG state: ", e)
    end
end

"""
  load_rng!(data[; rng = _GLOBAL_RNG, entryname = "RNG"])

Loads an RNG from a given `data` dictinary generated by `JLD.load()` to `rng`.
"""
function load_rng!(
        data::Dict;
        rng = _GLOBAL_RNG,
        entryname::String="RNG"
    )
    try
        copy!(rng, data[entryname])
    catch e
        error("Error while restoring RNG state: ", e)
    end
end
