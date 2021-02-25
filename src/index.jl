function Index(path::AbstractString)
    db = SQLite.DB(path)
    Index(path, db, [], [], [], [])
end

@inline function _check_idx(b::Bgen)
    @assert b.idx !== nothing "bgen index (.bgi) is needed"
end

function select_region(idx::Index, chrom::AbstractString;
        start=nothing, stop=nothing)
    if start === nothing && stop === nothing
        q = "SELECT file_start_position FROM Variant WHERE chromosome=?"
        params = (chrom,)
    elseif stop === nothing
        q = "SELECT file_start_position FROM Variant" *
                    " WHERE chromosome=? AND position>=?"
        params = (chrom, start)
    else
        q = "SELECT file_start_position FROM Variant" *
                    " WHERE chromosome=? AND position>=? AND position<=?"
        params = (chrom, start, stop)
    end
    r = (DBInterface.execute(idx.db, q, params) |> columntable)[1]
end

"""
    select_region(bgen, chrom; start=nothing, stop=nothing)
Select variants from a region. Returns variant start offsets on the file.
Returns a `VariantIteratorFromOffsets` object.
"""
function select_region(b::Bgen, chrom::AbstractString;
        start=nothing, stop=nothing)
    _check_idx(b)
    offsets = select_region(b.idx, chrom; start=start, stop=stop)
    VariantIteratorFromOffsets(b, offsets)
end

function variant_by_rsid(idx::Index, rsid::AbstractString)
    q = "SELECT file_start_position FROM Variant WHERE rsid= ?"
    params = (rsid,)
    r = (DBInterface.execute(idx.db, q, params) |> columntable)[1]
    if length(r) == 0
        @error "variant with rsid $rsid not found"
    elseif length(r) > 1
        @error "multiple variant matches with $rsid"
    end
    return r[1]
end

"""
    variant_by_rsid(bgen, rsid)
Find a variant by rsid
"""
function variant_by_rsid(b::Bgen, rsid::AbstractString)
    _check_idx(b)
    offset = variant_by_rsid(b.idx, rsid)
    return Variant(b, offset)
end

function variant_by_pos(idx::Index, pos::Integer)
    q = "SELECT file_start_position FROM Variant WHERE position= ?"
    params = (pos,)
    r = (DBInterface.execute(idx.db, q, params) |> columntable)[1]
    if length(r) == 0
        @error "variant match at $pos not found"
    elseif length(r) > 1
        @error "multiple variant matches at $pos"
    end
    return r[1]
end

"""
    variant_by_pos(bgen, pos)
Get the variant of bgen variant given `pos` in the index file
"""
function variant_by_pos(b::Bgen, pos::Integer)
    _check_idx(b)
    offset = variant_by_pos(b.idx, pos)
    return Variant(b, offset)
end

function variant_by_index(idx::Index, index::Integer)
    q = "SELECT file_start_position FROM Variant LIMIT 1 OFFSET ?"
    params = (index - 1,)
    r = (DBInterface.execute(idx.db, q, params) |> columntable)[1]
    return r[1]
end

"""
    variant_by_index(bgen, n)
get the `n`-th variant (1-based).
"""
function variant_by_index(b::Bgen, index::Integer)
    _check_idx(b)
    offset = variant_by_index(b.idx, index)
    return Variant(b, offset)
end

function offsets(idx::Index)
    if length(idx.offsets) != 0
        return idx.offsets
    end
    q = "SELECT file_start_position FROM Variant"
    r = (DBInterface.execute(idx.db, q) |> columntable)[1]
    resize!(idx.offsets, length(r))
    idx.offsets .= r
    return r
end

function rsids(idx::Index)
    if length(idx.rsids) != 0
        return idx.rsids
    end
    q = "SELECT rsid FROM Variant"
    r = (DBInterface.execute(idx.db, q) |> columntable)[1]
    resize!(idx.rsids, length(r))
    idx.rsids .= r
    return r
end

function chroms(idx::Index)
    if length(idx.chroms) != 0
        return idx.chroms
    end
    q = "SELECT chromosome FROM Variant"
    r = (DBInterface.execute(idx.db, q) |> columntable)[1]
    resize!(idx.chroms, length(r))
    idx.chroms .= r
    return r
end

function positions(idx::Index)
    if length(idx.positions) != 0
        return idx.positions
    end
    q = "SELECT position FROM Variant"
    r = (DBInterface.execute(idx.db, q) |> columntable)[1]
    resize!(idx.positions, length(r))
    idx.positions .= r
    return r
end
