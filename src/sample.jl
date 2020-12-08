struct Samples
    samples :: Vector{String}
end

function Samples(io::IOStream, n_samples::Integer)
    sample_header_length = read(io, UInt32)
    sample_n_check = read(io, UInt32)
    @assert n_samples == sample_n_check "Inconsistent number of samples"
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
