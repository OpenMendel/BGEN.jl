function length(x::Samples)
    return length(x.samples)
end
function getindex(x::Samples, i)
    return x.samples[i]
end
function setindex(x::Samples, v, i)
    x.samples[i] = v
end
function firstindex(x::Samples)
    return firstindex(x.samples)
end
function lastindex(x::Samples)
    return lastindex(x.samples)
end

function Samples(io::IOStream, n_samples::Integer)
    sample_header_length = read(io, UInt32)
    sample_n_check = read(io, UInt32)
    if sample_n_check != 0
        @assert n_samples == sample_n_check "Inconsistent number of samples"
    else
        @warn "Sample names unavailable. Do you have a separate '.sample' file?"
    end
    samples = String[]
    for i in 1:n_samples
        id_length = read(io, UInt16)
        id = read(io, id_length)
        push!(samples, String(id))
    end
    Samples(samples)
end

function Samples(path::String, n_samples::Integer)
    @assert endswith(path, ".sample") "Extension of the file should be .sample"
    io = open(path)
    readline(io) # header
    readline(io) # types
    samples = readlines(io)
    @assert length(samples) == n_samples "Inconsistent number of samples"
    Samples(samples)
end

function Samples(n_samples::Integer)
    samples = [string(i) for i in 1:n_samples]
    Samples(samples)
end
