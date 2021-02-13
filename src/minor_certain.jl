"""
    minor_certain(freq, n_checked, z)
Check if minor allele is certain.

- `freq`: frequency of minor or major allele
- `n_checked`: number of individuals checked so far
- `z`: cutoff of "z" value, defaults to `5.0`
"""
function minor_certain(freq::Float64, n_checked::Integer, z::Float64=5.0)
    delta = z * sqrt(freq * (1-freq) / n_checked)
    return !(freq - delta < 0.5 && freq + delta > 0.5)
end
