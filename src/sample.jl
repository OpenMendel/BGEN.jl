function get_samples(io::IOStream, n_samples::Integer)
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
    samples
end

function get_samples(path::String, n_samples::Integer)
    @assert endswith(path, ".sample") "Extension of the file should be .sample"
    io = open(path)
    keys = split(readline(io)) # header
    key_idx = ("ID_1" in keys && "ID_2" in keys) ? 2 : 1
    readline(io) # types
    samples = readlines(io)
    samples = map(x -> split(x, " ")[key_idx], samples)
    @assert length(samples) == n_samples "Inconsistent number of samples"
    close(io)
    samples
end

function get_samples(n_samples::Integer)
    samples = [string(i) for i in 1:n_samples]
    samples
end
