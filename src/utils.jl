# Distances with periodic boundaries, random displacements and the Metropolis criterion
using StaticArrays
using Random
function euclidean_distance(ri::MVector{3,Float64},rj::MVector{3,Float64})::Float64
    # computes the euclidean distance between ri and rj.
    # outputs r = √(Δx^2 + Δy^2 + Δz^2)
    Δx = ri[1]-rj[1]
    Δy = ri[2]-rj[2]
    Δz = ri[3]-rj[3]
    rij  = sqrt(Δx*Δx + Δy*Δy+ Δz*Δz)
    return(rij)
end # euclidean distance

function min_config_distance(r::Vector{MVector{3,Float64}})::Tuple{Float64,Int64,Int64}
    # checks each pairwise distance in r and returns the minimum distance and the indices of the particles that that distance is between
    # so the user can identify if there is overlap or not. Needs to be called after renormalization to box units
    # and periodic boundary assurances are applied

    N = length(r)
    m = length(r[1])
    @assert m==3 "position matrix doesn't seem to be 3xN"

    min_distance = typemax(Float64)
    i_min = 0
    j_min = 0
    @inbounds for i in 1:(N-1)
        for j in (i+1):N
            distance = euclidean_distance(r[i],r[j])
            if distance < min_distance
                min_distance = distance
                i_min = i 
                j_min = j 
            end
        end
    end # ij pair loop

    return(min_distance,i_min,j_min)
end #min_config_distance(r)

@inline function euclidean_distance_squared_pbc(ri_box::MVector{3,Float64},rj_box::MVector{3,Float64})::Float64
    # computes the euclidean squared distance between ri and rj, written out to avoid allocations
    # does so assuming ri and rj are in box units for periodic boundary conditions
    # outputs r2 = Δx^2 + Δy^2 + Δz^2
    # (squared, so no square root is needed for the LJ potential)
    Δx = ri_box[1]-rj_box[1]
    Δy = ri_box[2]-rj_box[2]
    Δz = ri_box[3]-rj_box[3]

    # now doing PBC in box=1 units
    Δx = Δx - round(Δx)
    Δy = Δy - round(Δy)
    Δz = Δz - round(Δz)

    rij_squared  = Δx*Δx + Δy*Δy+ Δz*Δz
    return(rij_squared)
end # euclidean distance squared


function translate_by_random_vector!(r::MVector{3,Float64},δr_max::Float64,rng::MersenneTwister,c::SimCache)
    # displaces r in place by a uniform random amount in [-δr_max, δr_max] along each axis; PBC are not applied
    @inbounds for i in 1:3
        r[i] = r[i] + δr_max* (2* rand(rng) - 1.0)
    end
    return nothing
end # translate_by_random_vector

function metropolis(ΔE::Float64,T_σ::Float64,rng::MersenneTwister=MersenneTwister())::Bool
    # in the case of a translational move, the acceptance criteria reduces to the standard metropolis one. Assumes ΔE and T_σ have LJ units (both have energy units)
    exponent = -1*ΔE/T_σ 

    if exponent < -75 # exp(-75) ≈ 3e-33: reject without evaluating
        return(false)
    elseif exponent > 0 # downhill, accept without evaluating due to min(1,e^exponent) and e^exponent > 1 if exponent > 0 
        return(true)
    else 
        ζ = rand(rng) # random number between zero and 1
        accept = (exp(exponent) > ζ)   #boolean
        return(accept)
    end 

end # metropolis